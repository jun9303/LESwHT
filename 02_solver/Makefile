# Compiler & Flags
FC = ifort -r8 -i8 -O3 -qopenmp -mcmodel=medium -warn nounused
LIB = 

# Directory for FFTE
FFTE_DIR = ./f90/lib/poiss/ffte-7.0

# FFTE Object List
FFTE_OBJS = $(FFTE_DIR)/zfft1d.o \
            $(FFTE_DIR)/zfft2d.o \
            $(FFTE_DIR)/zfft3d.o \
            $(FFTE_DIR)/fft235.o \
            $(FFTE_DIR)/kernel.o \
            $(FFTE_DIR)/mfft235.o \
            $(FFTE_DIR)/factor.o

# Solver Object List
# For the Poisson solver, choose one among poiss_[mgrd, x_ft, z_ft, ftft].o depending on simulation setups
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
        funcbody.o \
        ibm_body.o \
        ibm_forcing.o

# Main Solver Executable Linkage
solver: $(FFTE_OBJS) $(OBJS)
	$(FC) ./f90/solver.f90 -o solver_exec $(OBJS) $(FFTE_OBJS) $(LIB)

#-----------------------------------------------------------------------
#     FFTE Compilation Rules (Fortran 77)
#-----------------------------------------------------------------------
# The -I$(FFTE_DIR) flag ensures param.h is found during compilation
$(FFTE_DIR)/%.o: $(FFTE_DIR)/%.f
	$(FC) -fPIC -I$(FFTE_DIR) -c $< -o $@

#-----------------------------------------------------------------------
#     Module Dependencies
#-----------------------------------------------------------------------
mod_common.o: ./f90/mod_common.f90
	$(FC) ./f90/mod_common.f90 -c $(LIB)

mod_flowarray.o: ./f90/mod_flowarray.f90
	$(FC) ./f90/mod_flowarray.f90 -c $(LIB)

mod_poiss.o: ./f90/mod_poiss.f90
	$(FC) ./f90/mod_poiss.f90 -c $(LIB)

#-----------------------------------------------------------------------
#     Solver Libraries
#-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
#     Poisson Solver
#-----------------------------------------------------------------------
poiss_mgrd.o: ./f90/lib/poiss/poiss_mgrd.f90
	$(FC) ./f90/lib/poiss/poiss_mgrd.f90 -c $(LIB)

poiss_x_ft.o: ./f90/lib/poiss/poiss_x_ft.f90
	$(FC) ./f90/lib/poiss/poiss_x_ft.f90 -c $(LIB)

poiss_z_ft.o: ./f90/lib/poiss/poiss_z_ft.f90
	$(FC) ./f90/lib/poiss/poiss_z_ft.f90 -c $(LIB)

poiss_ftft.o: ./f90/lib/poiss/poiss_ftft.f90
	$(FC) ./f90/lib/poiss/poiss_ftft.f90 -c $(LIB)

#-----------------------------------------------------------------------
#     External & Geometry Linkage
#-----------------------------------------------------------------------
funcbody.o: ../geometry/funcbody.f90
	$(FC) ../geometry/funcbody.f90 -c $(LIB)

ibm_body.o: ./f90/lib/ibm_body.f90
	$(FC) ./f90/lib/ibm_body.f90 -c $(LIB)

ibm_forcing.o: ./f90/lib/ibm_forcing.f90 mod_common.o mod_flowarray.o
	$(FC) ./f90/lib/ibm_forcing.f90 -c $(LIB)

#-----------------------------------------------------------------------
#     Utilities
#-----------------------------------------------------------------------
clean:
	rm -f *.o *.mod *.e core coredir.* fort.*
	rm -f solver_exec
	rm -f $(FFTE_DIR)/*.o

new:
	make clean
	make solver
