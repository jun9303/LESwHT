# LESwHT (Large Eddy Simulation with Heat Transfer)

LESwHT is a CFD workflow for incompressible turbulent flow with immersed boundaries and optional heat transfer. The repository includes:

- Python-based grid generation and IBM pre-processing
- Fortran solver (OpenMP)
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

- Linux
- Intel Fortran compiler (`ifort`)
- Intel C compiler (`icc`) for the f2py extension build path in pre-processing
- Python 3
- NumPy (with `numpy.f2py`)
- OpenMP runtime (`libiomp5`)

> Notes
>
> - `01_pre_processor/Makefile` builds `lib_ibm_body` using f2py with `--fcompiler=intelem`.
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