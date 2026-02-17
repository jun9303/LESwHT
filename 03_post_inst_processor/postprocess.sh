#!/bin/bash
echo 3 > /proc/sys/vm/drop_caches
ifort -r8 -i8 -O3 -mcmodel=medium -warn nounused post_inst.f90 -o postproc
./postproc
rm -rf postproc *.mod
