#!/bin/sh
ulimit -s unlimited
ulimit -v unlimited
export OMP_STACKSIZE=99999
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=10
export OMP_DYNAMIC=TRUE
make new
./solver_exec
