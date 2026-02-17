FC = ifort -r8 -i8 -O3 -qopenmp -mcmodel=medium -warn nounused

OBJS =  mod_common.o \
        mod_flowarray.o \
        mod_poiss.o \
        misc_init.o \
        misc_hist.o \
        sgs.o \
        slv_mmtm.o \
        slv_engy.o \
        slv_cont.o \
        slv_auxi.o \
        poiss_ftft.o \
        ibm_body.o \
        funcbody.o \
        
solver: $(OBJS)
	$(FC) ./f90/$@.f90 -o solver_exec $(OBJS) $(LIB)

mod_common.o: ./f90/mod_common.f90
	$(FC) ./f90/mod_common.f90 -c $(LIB)

mod_flowarray.o: ./f90/mod_flowarray.f90
	$(FC) ./f90/mod_flowarray.f90 -c $(LIB)

mod_poiss.o: ./f90/mod_poiss.f90
	$(FC) ./f90/mod_poiss.f90 -c $(LIB)

misc_init.o: ./f90/lib/misc_init.f90
	$(FC) ./f90/lib/misc_init.f90 -c $(LIB)

misc_hist.o: ./f90/lib/misc_hist.f90
	$(FC) ./f90/lib/misc_hist.f90 -c $(LIB)

sgs.o: ./f90/lib/sgs.f90
	$(FC) ./f90/lib/sgs.f90 -c $(LIB)

slv_mmtm.o: ./f90/lib/slv_mmtm.f90
	$(FC) ./f90/lib/slv_mmtm.f90 -c $(LIB)

slv_engy.o: ./f90/lib/slv_engy.f90
	$(FC) ./f90/lib/slv_engy.f90 -c $(LIB)

slv_cont.o: ./f90/lib/slv_cont.f90
	$(FC) ./f90/lib/slv_cont.f90 -c $(LIB)

slv_auxi.o: ./f90/lib/slv_auxi.f90
	$(FC) ./f90/lib/slv_auxi.f90 -c $(LIB)

poiss_mgrd.o: ./f90/lib/poiss/poiss_mgrd.f90
	$(FC) ./f90/lib/poiss/poiss_mgrd.f90 -c $(LIB)

poiss_x_ft.o: ./f90/lib/poiss/poiss_x_ft.f90
	$(FC) ./f90/lib/poiss/poiss_x_ft.f90 -c $(LIB)

poiss_z_ft.o: ./f90/lib/poiss/poiss_z_ft.f90
	$(FC) ./f90/lib/poiss/poiss_z_ft.f90 -c $(LIB)

poiss_ftft.o: ./f90/lib/poiss/poiss_ftft.f90
	$(FC) ./f90/lib/poiss/poiss_ftft.f90 -c $(LIB)

poiss.o: ./f90/lib/poiss/poiss.f90
	$(FC) ./f90/lib/poiss/poiss.f90 -c $(LIB)

ibm_body.o: ../global_lib/f90/ibm_body.f90
	$(FC) ../global_lib/f90/ibm_body.f90 -c $(LIB)

funcbody.o: ../global_lib/f90/funcbody.f90
	$(FC) ../global_lib/f90/funcbody.f90 -c $(LIB)

clean:
	rm -f *.o *.mod *.e core
	rm -rf coredir.*

new:
	make clean
	make
	make clean
