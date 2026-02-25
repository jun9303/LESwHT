#!/bin/bash

set -e

# 1. Resource Limits
ulimit -s unlimited
# ulimit -v unlimited  # Let scheduler handle memory

# 2. OpenMP Settings
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=4
fi

export OMP_STACKSIZE=1G
export OMP_SCHEDULE="dynamic"
export OMP_DYNAMIC=TRUE

# 3. GNU Toolchain (hardened)
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran
export F90=gfortran

# 4. Clean Previous Output
echo "Cleaning old IBM pre-processing files..."
make clean
rm -rf ../output/ibmpre/*

# 5. Compilation (f2py)
#    Ensure standard output/error is visible if compilation fails
echo "Compiling Fortran Extension (make new)..."
make flib

# 6. Execution
echo "Running IBM Pre-processor..."
python3 preprocessing.py

echo "Done."