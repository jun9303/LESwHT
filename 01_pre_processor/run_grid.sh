#!/bin/bash

set -e

if [ -f ../output/running_case.env ]; then
    # shellcheck source=/dev/null
    source ../output/running_case.env
fi

export LESWHT_CASE_NAME="${LESWHT_CASE_NAME:-}"
export LESWHT_OUTPUT_ROOT="${LESWHT_OUTPUT_ROOT:-../output}"
export LESWHT_GEOMETRY_STL="${LESWHT_GEOMETRY_STL:-}"

# 1. Resource Limits
ulimit -s unlimited
# ulimit -v unlimited  # Let scheduler handle memory

# 2. OpenMP Settings
#    If running via Slurm, use $SLURM_CPUS_PER_TASK, otherwise default to 64.
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=64
fi

export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_STACKSIZE=1G
export OMP_SCHEDULE="dynamic"
export OMP_DYNAMIC=FALSE
echo "OpenMP config: OMP_NUM_THREADS=${OMP_NUM_THREADS}, OMP_DYNAMIC=${OMP_DYNAMIC}, OMP_PROC_BIND=${OMP_PROC_BIND}, OMP_PLACES=${OMP_PLACES}"

# 3. Clean Previous Output
echo "Cleaning old grid files..."
mkdir -p "${LESWHT_OUTPUT_ROOT}/grid"
rm -rf "${LESWHT_OUTPUT_ROOT}/grid"/*

# 4. Execution
echo "Running Grid Generator..."
python3 grid.py

echo "Done."