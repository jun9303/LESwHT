#!/bin/bash

set -e

# 1. Resource Limits
ulimit -s unlimited
# ulimit -v unlimited  # Let scheduler handle memory

# 2. Environment Setup (Uncomment/Edit for your cluster)
# module purge
# module load gcc/<version>            # Load GNU compilers
# module load python/3.x    # Load Python 3
# source /path/to/venv/bin/activate # Activate virtual env if needed

# 2-1. GNU Toolchain (hardened)
export FC=gfortran
export F77=gfortran
export F90=gfortran
export CC=gcc
export CXX=g++

# 3. OpenMP Settings
#    If running via Slurm, use $SLURM_CPUS_PER_TASK, otherwise default to 4.
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=4
fi

export OMP_STACKSIZE=1G
export OMP_SCHEDULE="dynamic"
export OMP_DYNAMIC=TRUE

# 4. Compilation
echo "Compiling Solver (make new)..."
make new

# 5. Execution
echo "Running Solver..."
rm -f *.o *.mod *.e core coredir.* fort.*
./solver_exec

echo "Solver finished."
