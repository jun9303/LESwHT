# LESwHT (Large Eddy Simulation with Heat Transfer)

**LESwHT** is a computational fluid dynamics (CFD) solver designed for unsteady, incompressible, turbulent flows. It extends the LICA code to include complex heat transfer capabilities, making it suitable for simulating multi-physics problems involving fluid flow, heat transfer, and solid-fluid interactions.

## Key Features
* **Large Eddy Simulation (LES):** Resolves large-scale turbulent structures.
* **Immersed Boundary Method (IBM):** Handles complex geometries and solid boundaries on a Cartesian grid.
* **Conjugate Heat Transfer (CHT):** Solves coupled heat conduction (solid) and convection (fluid) equations simultaneously.
* **Subgrid-Scale Models:**
    * **Turbulence:** Vreman model.
    * **Heat Flux:** Linear Eddy viscosity model.
* **Buoyancy:** Includes Boussinesq approximation (Grashof number) for natural convection.

## Directory Structure
* `01_pre_processor/`: Grid generation and geometry setup (IBM).
* `02_solver/`: Main CFD engine (Fortran 90).
* `03_post_inst_processor/`: Visualization of instantaneous flow fields.
* `04_post_avg_processor/`: Calculation of time-averaged statistics.
* `global_lib/`: Shared libraries and geometry definitions.

## Prerequisites
* **OS:** Linux (Ubuntu 22.04 LTS recommended)
* **Compilers:**
    * Intel Fortran (`ifort`)
    * Python (3.10.12 recommended) with NumPy (1.26.4)
* **Hardware:** Multi-core CPU support via OpenMP

## Installation & Setup

### 1. Environment Setup
Add the following alias to your shell configuration (`.bashrc`) to standardise the compilation command:
```bash
alias cf90='ifort -r8 -i8 -O3 -xT -mcmodel=medium -i-dynamic -vec-report0 -warn nounused'
```
### 2. Compilation
The repository includes Makefiles in respective directories. Alternatively, you can compile manually using the alias above. Ensure you include the `-openmp` flag for parallel support.

#### Usage Workflow
##### Step 1: Pre-processing
Navigate to `01_pre_processor` and run the generation scripts to create the grid and IBM data:
```bash
cd 01_pre_processor
./makegrid.sh
./makeibmpre.sh
```
##### Step 2: Running the Solver
1. Navigate to `02_solver`.
2. Configure physics in settings.in (Re, Pr, Time steps) and boundaries in boundary.in.
3. Execute the solver script:
```bash
./run_solver.sh
```

##### Step 3: Post-processing
Navigate to `03_post_inst_processor` to convert binary output to Tecplot format:

```bash
cd 03_post_inst_processor
./postprocess.sh
```