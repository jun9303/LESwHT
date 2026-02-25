# Compiler Settings
PYTHON ?= python3
F2PY = $(PYTHON) -m numpy.f2py

# Intel Compiler Configuration
FCOMPILER ?= intelem
COMPILER ?= intel
INTEL_CC ?= icc
INTEL_CXX ?= icpc
INTEL_FC ?= ifort

# Optimization & OpenMP Flags
OPT_FLAGS = -O3 -qopenmp -fpp

# Linker Flags
LIBS = -liomp5

# Target Definitions
OUTPUT_NAME = lib_ibm_body

SRCS = ../02_solver/f90/lib/ibm_body.f90 \
	../geometry/funcbody.f90

.PHONY: all flib clean

all: flib

flib:
	@echo "----------------------------------------------------------"
	@echo "Building $(OUTPUT_NAME) using $(FCOMPILER) via f2py..."
	@echo "----------------------------------------------------------"
	CC=$(INTEL_CC) CXX=$(INTEL_CXX) FC=$(INTEL_FC) F77=$(INTEL_FC) F90=$(INTEL_FC) LDSHARED="$(INTEL_FC)" $(F2PY) -c $(SRCS) \
		-m $(OUTPUT_NAME) \
		--compiler=$(COMPILER) \
		--fcompiler=$(FCOMPILER) \
		--opt="$(OPT_FLAGS)" \
		$(LIBS)
	@echo "----------------------------------------------------------"
	@echo "Build Complete."
	@echo "----------------------------------------------------------"

clean:
	@echo "Cleaning up build artifacts..."
	rm -f *.o *.mod *.so
	rm -f $(OUTPUT_NAME).*.so