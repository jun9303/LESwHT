#!/bin/bash
ulimit -s unlimited
ulimit -v unlimited
export OMP_STACKSIZE=99999
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=TRUE
echo 3 > /proc/sys/vm/drop_caches
rm -rf ../output/grid/*
python3 01_grid.py #python ver. >= 3