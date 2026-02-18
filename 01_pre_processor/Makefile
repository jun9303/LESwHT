# Compiler Settings
PYTHON ?= python3
F2PY = $(PYTHON) -m numpy.f2py

# Intel Compiler Configuration
FCOMPILER = intelem
CC = icc

# Optimization & OpenMP Flags
OPT_FLAGS = -O3 -qopenmp -fpp

# Linker Flags
LIBS = -liomp5

# Target Definitions
OUTPUT_NAME = lib_ibm_body

SRCS = ../02_solver/f90/lib/ibm_body.f90 \
       ../geometry/f90/funcbody.f90

.PHONY: all flib clean

all: flib

flib:
	@echo "----------------------------------------------------------"
	@echo "Building $(OUTPUT_NAME) using $(FCOMPILER) via f2py..."
	@echo "----------------------------------------------------------"
	CC=$(CC) $(F2PY) -c $(SRCS) \
		-m $(OUTPUT_NAME) \
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