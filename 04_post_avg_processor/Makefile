#FC = mpiifort -r8 -i8 -O3 -pg -g -mcmodel=medium
#FC = mpiifort -r8 -i8 -O3 -openmp -pg -g -mcmodel=medium
#FC = mpif90 -O3 -pg -g -mcmodel=medium
#FC = mpif90 -O3 -openmp -pg -g -mcmodel=medium
#FC = ifort -r8 -i8 -O3 -xT -mcmodel=medium -i-dynamic -vec-report0 -warn nounused
FC = ifort -r8 -i8 -O3 -xT -mcmodel=medium -i-dynamic -vec-report0 -warn nounused
GEOM_SRC = ../geometry/f90/funcbody.f90

post_avg: post_avg.f90
	$(FC) $@.f90 $(GEOM_SRC) -o postavgproc $(LIB)

clean:
	rm -f *.o *.mod *.e core
	rm -rf coredir.*

post:
	make clean
	make
	./postavgproc
	rm -rf postavgproc
	make clean
