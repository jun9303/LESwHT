# LESwHT: dimple_surf branch

This branch runs LESwHT (large eddy simulation with immersed boundaries and heat transfer) for turbulent
channel flow over dimpled walls. Its pipeline starts with a **dimple generator**: a periodic dimpled slab is
built as an STL surface, and the solver's `geometry/funcbody.f90` reads that STL to tell solid from fluid.
Everything after the generator (grid, IBM pre-processing, solver, post-processing) runs on the generated
case.

The solver, its OpenMP settings and its bug fixes are the same as on `master`; see the `master` README for
the general solver notes (Poisson solver choice, bodies at domain boundaries, known limits).

## What differs from master

- `00_dimple_generator/`: builds the dimpled slab STL and a case directory, and writes
  `output/running_case.env`, which the later stages read.
- `geometry/funcbody.f90`: inside/outside test against the generated STL (not an analytic body).
- Case-local output: every stage writes under `output/<case>/` (set by `LESWHT_OUTPUT_ROOT`).
- `run.sh`: runs every stage in order, including a step that sets the constant-flow-rate target `UDRV_I`
  from the dimple volume.
- `03_post_processor/run_postprocess_all.sh`: post-processes every instantaneous and averaged field.
- Solver options used here: an optional eighth value `IDELTA_SGS` on the SGS line of `settings.input`
  (filter width of the thermal SGS model; 1: cube root of the cell volume, 2: largest cell size), a CFL
  step controller, a drag history split into IBM, inertia and wall-shear parts, and ICH = 1/2 runs that
  start from rest.

## Requirements

- Linux, GNU compilers (`gcc`, `g++`, `gfortran`) with OpenMP (`libgomp`)
- Python 3.10 or newer with the packages in `requirements.txt` (`numpy`, `scipy`, `trimesh`, `manifold3d`)

Set up once, from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Activate the environment (`source .venv/bin/activate`) in each new shell before `run.sh` or a stage
script. `run.sh` itself uses `.venv/bin/python` for the generator and the `UDRV_I` step, but the grid, IBM
and post-processing scripts call `python3` from `PATH`.

## Run the full pipeline

From the repository root:

```bash
bash run.sh -t 0 -d 30.0 -s 1.0 -x 1.0
```

| Option | Meaning | Default |
|---|---|---|
| `-t` | dimple topology (see below) | 0 |
| `-d` | dimple depth in wall units | 30.0 |
| `-s` | isotropic scale factor | 1.0 |
| `-x` | streamwise stretch factor | 1.0 |

Topologies: 0 cylinder, 1 spherical, 2 teardrop down, 3 teardrop up, 4 diamond, 5 triangle down,
6 triangle up.

## Run step by step

1. **Generate the dimpled slab.**

   ```bash
   cd 00_dimple_generator
   python3 main.py -t 0 -d 30.0 -s 1.0 -x 1.0 --compact
   ```

   This creates `output/<topology>_<depth>_<scale>_<stretch>_<timestamp>/` with
   `geometry/dimple_slab_ascii.stl`, and writes `output/running_case.env`, for example:

   ```bash
   LESWHT_CASE_NAME=0_30.0_1.0_1.0_20260415_101500
   LESWHT_OUTPUT_ROOT=/path/to/LESwHT/output/0_30.0_1.0_1.0_20260415_101500
   LESWHT_GEOMETRY_STL=/path/to/LESwHT/output/0_30.0_1.0_1.0_20260415_101500/geometry/dimple_slab_ascii.stl
   ```

2. **Build the grid and the IBM data.** Both scripts read `output/running_case.env`.

   ```bash
   cd ../01_pre_processor
   bash run_grid.sh
   bash run_preprocessing.sh
   ```

   Inputs: `grid.input` and `preprocessing.input`.

3. **Set the flow-rate target.** `run.sh` does this between steps 2 and 4: it computes the effective dimple
   height `h_eff` (dimple volume divided by the slab footprint) from the STL and writes
   `UDRV_I = 2 / (2 + h_eff)` into `02_solver/settings.input`. When running by hand, set `UDRV_I` yourself.

4. **Build and run the solver.**

   ```bash
   cd ../02_solver
   bash run_solver.sh
   ```

   Inputs: `settings.input` and `boundary.input`. The script reads `output/running_case.env` as well, and the
   solver stops with a message if it cannot read the STL named there. The thread count is `SLURM_CPUS_PER_TASK` under Slurm and 64
   otherwise (set `OMP_NUM_THREADS` in the script to change it). `POISS=poiss_fd bash run_solver.sh` builds
   with another Poisson solver; the default `poiss_ftft` fits the x- and z-periodic slab.

5. **Post-process.**

   ```bash
   cd ../03_post_processor
   bash run_postprocess_all.sh
   ```

   Single modes: `inst/run_postprocess.sh` (instantaneous) and `avg/run_postprocess.sh` (averaged), with
   `inst/post_inst.input` and `avg/post_avg.input`.

Results are under `output/<case>/`: `grid/`, `ibmpre/`, `field/`, `field_avg/`, `ftr/` (histories),
`post_inst/` and `post_avg/`.

## Clean up

```bash
bash reset.sh                      # build files and the output directory
cd 01_pre_processor && make clean  # pre-processor extension only
cd 02_solver && make clean         # solver only
```
