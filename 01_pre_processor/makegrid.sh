#!/bin/bash

# 1. Resource Limits
ulimit -s unlimited
# ulimit -v unlimited  # Let scheduler handle memory

# 2. Environment Setup (Uncomment/Edit for your cluster)
# module purge
# module load intel/19.1.3.304         # Load Intel Compilers
# module load python/3.x    # Load Python 3
# source /path/to/venv/bin/activate # Activate virtual env if needed

# 3. OpenMP Settings
#    If running via Slurm, use $SLURM_CPUS_PER_TASK, otherwise default to 4.
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=4
fi

export OMP_STACKSIZE=512M      # Increased for safety on large grids
export OMP_SCHEDULE="dynamic"
export OMP_DYNAMIC=TRUE

# 4. Clean Previous Output
#    (No root commands like 'drop_caches' allowed here)
echo "Cleaning old grid files..."
rm -rf ../output/grid/*

# 5. Execution
echo "Running Grid Generator..."
python3 01_grid.py

echo "Done."