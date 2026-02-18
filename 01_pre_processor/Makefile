# Compiler Settings
PYTHON ?= python3
F2PY = $(PYTHON) -m numpy.f2py

# Intel Compiler Configuration
FCOMPILER = intelem
# Force C Compiler to Intel C (icc)
CC = icc

# Optimization & OpenMP Flags
OPT_FLAGS = -O3 -qopenmp -fpp

# Linker Flags
LIBS = -liomp5

# Target Definitions
OUTPUT_NAME = lib_ibm_body
GLOBAL_LIB_DIR = ../global_lib
SRCS = ../global_lib/f90/ibm_body.f90 \
       ../global_lib/f90/funcbody.f90

.PHONY: all flib clean

all: flib

flib:
	@echo "----------------------------------------------------------"
	@echo "Building $(OUTPUT_NAME) using $(FCOMPILER) via f2py..."
	@echo "----------------------------------------------------------"
	# CC=$(CC) passes the 'icc' command to distutils/f2py
	CC=$(CC) $(F2PY) -c $(SRCS) \
		-m $(OUTPUT_NAME) \
		--fcompiler=$(FCOMPILER) \
		--opt="$(OPT_FLAGS)" \
		$(LIBS)
	@echo "----------------------------------------------------------"
	@echo "Installing $(OUTPUT_NAME).so to $(GLOBAL_LIB_DIR)..."
	@echo "----------------------------------------------------------"
	mv $(OUTPUT_NAME).*.so $(GLOBAL_LIB_DIR)/$(OUTPUT_NAME).so

clean:
	@echo "Cleaning up build artifacts..."
	rm -f *.o *.mod *.so
	rm -f $(GLOBAL_LIB_DIR)/$(OUTPUT_NAME).so