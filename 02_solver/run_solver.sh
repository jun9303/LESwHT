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
export OMP_DYNAMIC=FALSE        # keep the requested number of threads
export OMP_PLACES=cores         # one thread per core, threads kept close to each other
export OMP_PROC_BIND=close

# glibc: serve the work arrays that the solver allocates and frees at every call (ADI planes,
# SGS arrays) from the heap instead of mmap/munmap, which page-faults on every call
export MALLOC_MMAP_THRESHOLD_=33554432
export MALLOC_TRIM_THRESHOLD_=68719476736
export MALLOC_MMAP_MAX_=0

# 4. Compilation
#    POISS selects the Poisson solver (see Makefile), e.g. POISS=poiss_z_ft bash run_solver.sh
POISS=${POISS:-poiss_ftft}
echo "Compiling Solver (make new POISS=$POISS)..."
make new POISS=$POISS

# 5. Execution
#    On a node with several NUMA domains, numactl --interleave=all spreads the arrays over all
#    memory channels (the arrays are first written by one thread). Used when numactl is available.
RUN="./solver_exec"
if command -v numactl > /dev/null 2>&1 && numactl --interleave=all true > /dev/null 2>&1; then
    RUN="numactl --interleave=all ./solver_exec"
fi
echo "Running Solver ($RUN, $OMP_NUM_THREADS threads)..."
rm -f *.o *.mod *.e core coredir.* fort.*
$RUN

echo "Solver finished."
