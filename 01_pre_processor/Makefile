MKLROOT = /opt/intel/compilers_and_libraries_2018.1.163/linux
FC = f2py3 --fcompiler=intelem -L$(MKLROOT)/mkl/lib/intel64 -lmkl_rt -L$(MKLROOT)/lib/intel64 -liomp5
OUTPUT = lib_ibm_body

OBJS =  ../global_lib/f90/ibm_body.f90 \
        ../global_lib/f90/funcbody.f90 \

flib: $(OBJS)
	  $(FC) -c $(OBJS) -m $(OUTPUT)

clean:
	rm -f *.o *.mod *.e core
	rm -rf coredir.*

new:
	make clean
	make
	mv $(OUTPUT).* ../global_lib/$(OUTPUT).so
