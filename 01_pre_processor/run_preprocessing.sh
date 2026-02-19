#!/bin/bash

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

# 3. Clean Previous Output
echo "Cleaning old IBM pre-processing files..."
make clean
rm -rf ../output/ibmpre/*

# 4. Compilation (f2py)
#    Ensure standard output/error is visible if compilation fails
echo "Compiling Fortran Extension (make new)..."
make flib
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

# 5. Execution
echo "Running IBM Pre-processor..."
python3 preprocessing.py

echo "Done."