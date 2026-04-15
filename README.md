# LESwHT (Large Eddy Simulation with Heat Transfer)

LESwHT is a large eddy simulation (LES) workflow for incompressible turbulent flow with immersed boundaries (IB) and heat transfer (HT).

## dimple_surf Branch Status

This branch documents and supports the experimental dimple-surface pipeline. The default end-to-end entrypoint is `run.sh`, which now includes:

1. Dimple geometry generation
2. Grid generation + IBM pre-processing
3. CFR update step that rewrites `UDRV_I` in `02_solver/settings.input` from STL-derived effective height
4. Solver run
5. Automatic post-processing for all discovered instantaneous and averaged fields

## Repository Layout (Current)

- `00_dimple_generator/`: dimple slab STL generator and per-case output bootstrap
- `01_pre_processor/`: grid generation (`grid.py`) and IBM pre-processing (`preprocessing.py`)
- `02_solver/`: LES/IBM/HT solver (`solver_exec`, Fortran)
- `03_post_processor/`: post-processing orchestration (`run_postprocess_all.sh`) and per-mode runners
- `geometry/`: body geometry source (`funcbody.f90`)
- `output/`: generated case directories and active case pointer (`running_case.env`)
- `run.sh`: full experimental pipeline for dimple cases
- `reset.sh`: clean build artifacts and reset root output directory

## Requirements

- Linux
- GNU toolchain: `gcc`, `g++`, `gfortran`
- Python 3.10+
- Python packages from `requirements.txt` (notably `numpy`, `trimesh`, `manifold3d`)
- OpenMP runtime (`libgomp`)

Recommended setup from project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick Start (Experimental dimple_surf Workflow)

From project root:

```bash
bash run.sh
```

Optional parameters:

```bash
bash run.sh -t 0 -d 45.0 -s 1.0 -x 1.0
```

Arguments in `run.sh`:

- `-t`: topology ID (`0` to `6`)
- `-d`: dimple depth (wall units)
- `-s`: isotropic scaling
- `-x`: streamwise stretch/scaling

`run.sh` stage flow:

1. `[1/4]` Generate dimple slab geometry by running `00_dimple_generator/main.py`.
2. `[2/4]` Run `01_pre_processor/run_grid.sh` and `01_pre_processor/run_preprocessing.sh`.
3. `[2.5/4]` Load STL from `output/running_case.env`, compute effective dimple volume/height, and update `UDRV_I` in `02_solver/settings.input`.
4. `[3/4]` Build and run solver with `02_solver/run_solver.sh`.
5. `[4/4]` Run `03_post_processor/run_postprocess_all.sh`.

Treat `run.sh` as the default command for experimental dimple runs.

## Dimple Generator Description

The generator in `00_dimple_generator/` builds a watertight periodic dimple slab STL from baseline templates and writes the run context for downstream scripts.

Primary entry:

```bash
cd 00_dimple_generator
python3 main.py -t 0 -d 30.0 -s 1.0 -x 1.0 --compact
```

Supported topologies:

- `0`: Cylinder
- `1`: Spherical
- `2`: Teardrop Down
- `3`: Teardrop Up
- `4`: Diamond
- `5`: Triangle Down
- `6`: Triangle Up

Important behavior:

- Uses baseline STL templates under `00_dimple_generator/baselineSTL/`.
- Topology `0` uses an analytic cylindrical baseline path.
- Boolean carving is performed with trimesh/manifold backend when available.
- `--compact` (default) creates a compact periodic cell (`z in [-5, 5]`, bounded streamwise extent).

Generator outputs:

- Case directory: `output/<topology>_<depth>_<scale>_<stretch>_<timestamp>/`
- Geometry STL: `.../geometry/dimple_slab_ascii.stl`
- Active case pointer: `output/running_case.env`

## Active Case File (running_case.env)

`output/running_case.env` is the interface between geometry generation and all downstream stages.

Example:

```bash
LESWHT_CASE_NAME=0_45.0_1.0_1.0_20260415_101500
LESWHT_OUTPUT_ROOT=/anvil/projects/x-phy250071/LESwHT/output/0_45.0_1.0_1.0_20260415_101500
LESWHT_GEOMETRY_STL=/anvil/projects/x-phy250071/LESwHT/output/0_45.0_1.0_1.0_20260415_101500/geometry/dimple_slab_ascii.stl
```

How it is used:

1. `00_dimple_generator/main.py` writes/updates this file.
2. `01_pre_processor/run_grid.sh` and `01_pre_processor/run_preprocessing.sh` source it.
3. `run.sh` CFR update stage reads `LESWHT_GEOMETRY_STL` from it.
4. Solver and post-process scripts use `LESWHT_OUTPUT_ROOT` for case-local I/O.

## Manual Stage-by-Stage Run

If you are not using full orchestration:

1. Generate geometry:

```bash
cd 00_dimple_generator
python3 main.py -t 0 -d 45.0 -s 1.0 -x 1.0 --compact
```

2. Run pre-processing:

```bash
cd ../01_pre_processor
bash run_grid.sh
bash run_preprocessing.sh
```

3. Run solver:

```bash
cd ../02_solver
bash run_solver.sh
```

4. Run all post-processing:

```bash
cd ../03_post_processor
bash run_postprocess_all.sh
```

Direct postprocessor runners:

- Instantaneous: `03_post_processor/inst/run_postprocess.sh`
- Averaged: `03_post_processor/avg/run_postprocess.sh`

## Common Inputs to Edit

- `01_pre_processor/grid.input`
- `01_pre_processor/preprocessing.input`
- `02_solver/settings.input`
- `02_solver/boundary.input`
- `03_post_processor/inst/post_inst.input`
- `03_post_processor/avg/post_avg.input`

## Cleaning and Re-running

Full reset:

```bash
bash reset.sh
```

Targeted clean:

```bash
cd 01_pre_processor && make clean
cd 02_solver && make clean
```