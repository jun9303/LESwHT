#!/bin/bash
ifort -r8 -i8 -O3 -mcmodel=medium -warn nounused post_inst.f90 ../geometry/f90/funcbody.f90 -o postproc
./postproc
rm -rf postproc *.mod
