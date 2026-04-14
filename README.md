# LESwHT (Large Eddy Simulation with Heat Transfer)

LESwHT is a large eddy simulation (LES) workflow for incompressible turbulent flow with immersed boundaries (IB) and heat transfer (HT). The repository includes:

- Python-based grid generation and IBM pre-processing
- Fortran LES solver (OpenMP shared-memory parallelism)
- Post-processors for instantaneous and averaged fields

## Current Repository Layout

- `01_pre_processor/`: grid generation (`grid.py`) and IBM pre-processing (`preprocessing.py`)
- `02_solver/`: main LES/IBM/heat-transfer solver (`solver_exec` built from `f90/*.f90`)
- `03_post_processor/`: Python post-processing package
	- `inst/`: instantaneous field post-processing
	- `avg/`: averaged field post-processing
- `geometry/`: body geometry function source (`geometry/funcbody.f90`)
- `output/`: per-case simulation outputs and active case pointer file
- `run.sh`: primary orchestration script for the full workflow (geometry → pre-process → solver → post-process)
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

This executes the orchestrated flow in `run.sh`:

1. Geometry generation via `00_dimple_generator/main.py` (`-t -d -s -x`)
2. Grid generation via `01_pre_processor/run_grid.sh`
3. IBM pre-processing via `01_pre_processor/run_preprocessing.sh`
4. Automatic CFR target update in `02_solver/settings.input` (`UDRV_I`)
5. Solver execution via `02_solver/run_solver.sh`
6. Full post-processing via `03_post_processor/run_postprocess_all.sh`

`run.sh` should be treated as the default entrypoint for end-to-end runs.

## Active Running Case File

The active case is tracked in `output/running_case.env`. This file is intentionally simple so you can edit it manually to point post-generation tools to any existing case folder.

Current format:

```bash
LESWHT_CASE_NAME=5_0.0_1.0_1.0_20260410_121706
LESWHT_OUTPUT_ROOT=/anvil/projects/x-phy250071/LESwHT/output/5_0.0_1.0_1.0_20260410_121706
LESWHT_GEOMETRY_STL=/anvil/projects/x-phy250071/LESwHT/output/5_0.0_1.0_1.0_20260410_121706/geometry/dimple_slab_ascii.stl
```

How it is used:

1. `00_dimple_generator/main.py` writes/updates `output/running_case.env` after geometry generation.
2. `01_pre_processor/run_grid.sh` sources `output/running_case.env` and writes grid outputs into `${LESWHT_OUTPUT_ROOT}/grid`.
3. `01_pre_processor/run_preprocessing.sh` sources `output/running_case.env` and writes preprocessing outputs into `${LESWHT_OUTPUT_ROOT}/ibmpre`.

How to switch active case manually:

1. Edit `output/running_case.env`.
2. Set `LESWHT_OUTPUT_ROOT` and `LESWHT_GEOMETRY_STL` to an existing case directory.
3. Re-run `run_grid.sh` and/or `run_preprocessing.sh`.

## Stage-by-Stage Run

If you are not using the orchestrated `run.sh` flow, you can run stages manually as below.

### 1) Geometry Generation

From `00_dimple_generator/`:

```bash
python3 main.py -t 0 -d 45.0 -s 1.0 -x 1.0 --compact
```

Outputs:

- case directory under `output/<topology>_<depth>_<scale>_<stretch>_<timestamp>/`
- geometry STL at `.../geometry/dimple_slab_ascii.stl`
- active case file `output/running_case.env`

### 2) Reset (optional)

```bash
bash reset.sh
```

Actions:

- cleans build artifacts in `01_pre_processor` and `02_solver`
- recreates only the `output/` root directory (case subdirectories are created by geometry generation)

### 3) Grid Generation

```bash
cd 01_pre_processor
bash run_grid.sh
```

Inputs:

- `output/running_case.env` (active case selection)
- case parameters from `LESWHT_CASE_NAME` when available (fallback: `01_pre_processor/grid.input`)

Outputs:

- `output/grid/grid.dat`
- optional debug files in `output/grid/` (depending on `grid.input` debug options)

### 4) IBM Pre-processing

```bash
cd 01_pre_processor
bash run_preprocessing.sh
```

Inputs:

- `output/running_case.env` (active case selection)
- `01_pre_processor/preprocessing.input`
- `${LESWHT_OUTPUT_ROOT}/grid/grid.dat`

Outputs:

- IBM preprocessed binaries in `output/ibmpre/`

### 5) Solver

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