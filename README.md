# LESwHT (Large Eddy Simulation with Heat Transfer)

LESwHT is a large eddy simulation (LES) workflow for incompressible turbulent flow with immersed boundaries (IB) and heat transfer (HT). The repository includes:

- Python-based grid generation and IBM pre-processing
- Fortran LES solver (OpenMP shared-memory parallelism; see [Parallel runs](#parallel-runs-openmp))
- Post-processors for instantaneous and averaged fields

## Current Repository Layout

- `01_pre_processor/`: grid generation (`grid.py`) and IBM pre-processing (`preprocessing.py`)
- `02_solver/`: main LES/IBM/heat-transfer solver (`solver_exec` built from `f90/*.f90`)
- `03_post_processor/`: Python post-processing package
	- `inst/`: instantaneous field post-processing
	- `avg/`: averaged field post-processing
- `geometry/`: body geometry function source (`geometry/funcbody.f90`)
- `output/`: simulation outputs (`field`, `field_avg`, `post_inst`, `post_avg`, `grid`, `ibmpre`, `ftr`)
- `run.sh`: full pipeline runner (reset → pre-process → solver)
- `reset.sh`: cleanup and re-initialize output directories

## Requirements

- Linux (recommended)
- GNU Fortran compiler (`gfortran`)
- GNU C compiler (`gcc`) for the f2py extension build path in pre-processing
- Python 3 (tested: 3.10.12)
- NumPy (tested: 2.2.6; with `numpy.f2py`)
- OpenMP runtime (`libgomp`)

> Notes
>
> - `01_pre_processor/Makefile` builds `lib_ibm_body` using f2py with `--fcompiler=gnu95`.
> - The preprocessor and solver run scripts export `CC/CXX/FC/F77/F90` as `gcc/g++/gfortran` to prevent accidental compiler switching from environment defaults.
> - Scripts set `OMP_NUM_THREADS` from `SLURM_CPUS_PER_TASK` when available; otherwise they default to `4`.

## Quick Start (Full Workflow)

From the project root:

```bash
bash run.sh
```

This executes:

1. `bash reset.sh`
2. `01_pre_processor/run_grid.sh`
3. `01_pre_processor/run_preprocessing.sh`
4. `02_solver/run_solver.sh`

## Stage-by-Stage Run

### 1) Reset

```bash
bash reset.sh
```

Actions:

- cleans build artifacts in `01_pre_processor` and `02_solver`
- recreates `output/` subdirectories

### 2) Grid Generation

```bash
cd 01_pre_processor
bash run_grid.sh
```

Inputs:

- `01_pre_processor/grid.input`

Outputs:

- `output/grid/grid.dat`
- optional debug files in `output/grid/` (depending on `grid.input` debug options)

### 3) IBM Pre-processing

```bash
cd 01_pre_processor
bash run_preprocessing.sh
```

Inputs:

- `01_pre_processor/preprocessing.input`
- `output/grid/grid.dat`

Outputs:

- IBM preprocessed binaries in `output/ibmpre/`

### 4) Solver

```bash
cd 02_solver
bash run_solver.sh
```

Inputs:

- `02_solver/settings.input`
- `02_solver/boundary.input`
- preprocessed data in `output/ibmpre/`

Outputs:

- instantaneous fields in `output/field/`
- averaged fields in `output/field_avg/` (if enabled in solver settings)

## Parallel runs (OpenMP)

The solver runs on one shared-memory node. `02_solver/run_solver.sh` sets the thread count from
`SLURM_CPUS_PER_TASK` (4 otherwise) and the runtime settings below:

```bash
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close OMP_STACKSIZE=1G
export MALLOC_MMAP_THRESHOLD_=33554432 MALLOC_TRIM_THRESHOLD_=68719476736 MALLOC_MMAP_MAX_=0
numactl --interleave=all ./solver_exec    # used by run_solver.sh when numactl is available
```

- The `MALLOC_` settings let glibc reuse the work arrays that the ADI sweeps and the dynamic SGS model
  allocate at every call, instead of mapping and unmapping them each time.
- `numactl --interleave=all` spreads the arrays over all memory channels of a node with several NUMA
  domains (two-socket nodes, for example).
- Loops run over (k, j) pairs or line blocks, so 2D runs (one z-plane) and thin 3D domains use all
  threads, not only `n3m` of them.
- Sums are added in a fixed order, so the fields and the history files (`output/ftr`) are the same for
  any thread count.

### Poisson solver

Choose the Poisson solver at build time with `POISS` (`make new POISS=poiss_z_ft`, or
`POISS=poiss_z_ft bash run_solver.sh`):

| `POISS` | Periodic directions | Method |
|---|---|---|
| `poiss_ftft` (default) | x and z | FFT in x and z, TDMA in y |
| `poiss_z_ft` | z only (inlet/outlet in x) | FFT in z, multigrid in x-y |
| `poiss_x_ft` | x only (walls in z) | FFT in x, multigrid in y-z |
| `poiss_mgrd` | any | 3D multigrid with BiCGSTAB |
| `poiss_fd` | any, y non-periodic | direct solver (no iterations); a periodic direction must be uniform |

`poiss_fd` diagonalizes x and z with the eigenvectors of their 1D operators and solves one tridiagonal
system in y for each (x-mode, z-mode) combination, so every line is independent. Its cost rises as
`n1m^2` and `n3m^2`, so it suits moderate x and z sizes.

## Immersed bodies at the domain boundaries

- A body may cross the inlet or the outlet (for example channel walls that run through both ends, or a
  body cut by the outlet). At the convective outlet the solid part takes the body velocity (zero), the
  convective speed is the mean over the fluid part, and the fluid part is scaled so that the outflow
  equals the inflow.
- A body may touch or come within two cells of a non-periodic end (a wall, the inlet or the outlet). An
  interpolation stencil that would leave the domain skips that direction. IBM files written by an
  older preprocessor can hold such stencils; the solver then stops and asks for a new preprocessing run.
- A body may cross a periodic boundary in x or z. `FUNCBODY` must then return the periodic image itself,
  for example with the distance `xd = x - xc - lx*anint((x - xc)/lx)` for a body centred at `xc`.
- `BC_Y*`/`BC_Z*` = 1 (far field) sets the normal velocity to zero and the tangential velocity to zero
  gradient at that boundary.

After updating from an older version, rerun `01_pre_processor` (`run_grid.sh` and `run_preprocessing.sh`)
before the solver. Older IBM files lack the forcing points on periodic boundaries and list the
zero-eddy-viscosity cells one index off; the solver only detects the out-of-range stencils.

## Known limits

- Moving bodies (`IMOVING = 1`) are not available; the solver stops at start-up.
- `ICH = 1` (constant flow rate) needs a periodic x; with an inlet use `ICH = 0`.
- A y-periodic domain is not available.
- LES with one z-plane: `FILTER = 3` runs as `FILTER = 1`, and `FILTER = 2` stops.

## Post-processing

### Instantaneous (`03_post_processor/inst`)

```bash
cd 03_post_processor/inst
bash postprocess.sh
```

Uses:

- `03_post_processor/inst/post_inst.input`

### Averaged (`03_post_processor/avg`)

From `03_post_processor/avg/`:

```bash
bash postavgprocess.sh
```

Uses:

- `03_post_processor/avg/post_avg.input`

## Common Input Files to Edit

- `01_pre_processor/grid.input`: grid resolution/domain and spacing options
- `01_pre_processor/preprocessing.input`: IBM/heat-transfer preprocessing options
- `02_solver/settings.input`: solver physics and run controls
- `02_solver/boundary.input`: boundary conditions
- `03_post_processor/inst/post_inst.input`: instantaneous post-processing setup
- `03_post_processor/avg/post_avg.input`: averaged post-processing setup

## Cleaning and Re-running

- Full cleanup + output reset:

```bash
bash reset.sh
```

- Clean solver build only:

```bash
cd 02_solver && make clean
```

- Clean preprocessor extension build only:

```bash
cd 01_pre_processor && make clean
```