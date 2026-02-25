# Compiler Settings
PYTHON ?= python3
F2PY = $(PYTHON) -m numpy.f2py

# GNU Compiler Configuration (hardened)
override FCOMPILER := gnu95
override COMPILER := unix
override CC := gcc
override CXX := g++
override FC := gfortran

# Optimization & OpenMP Flags (legacy compatibility)
OPT_FLAGS ?= -O3 -fopenmp -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -ffree-line-length-none -fallow-argument-mismatch

# Linker Flags
LDFLAGS ?= -fopenmp
LIBS ?=

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
	CC=$(CC) CXX=$(CXX) FC=$(FC) F77=$(FC) F90=$(FC) LDFLAGS="$(LDFLAGS)" LDSHARED="$(FC) -shared" $(F2PY) -c $(SRCS) \
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