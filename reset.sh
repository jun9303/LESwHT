#!/bin/bash

# Define the target output directory
OUT_DIR="./output"

# CLEAN UP PREPROCESSING AND SOLVER BUILDS
echo "Cleaning up compilation files..."

if [ -d "01_pre_processor" ]; then
    echo "Cleaning 01_pre_processor..."
    make -C 01_pre_processor clean
else
    echo "Warning: 01_pre_processor directory not found."
fi

if [ -d "02_solver" ]; then
    echo "Cleaning 02_solver..."
    make -C 02_solver clean
else
    echo "Warning: 02_solver directory not found."
fi

# COMPLETE RESET: Check if the directory exists, and if so, wipe it.
if [ -d "$OUT_DIR" ]; then
    echo "Resetting: Removing existing '$OUT_DIR' directory and old data..."
    rm -rf "$OUT_DIR"
fi

echo "Initializing clean '$OUT_DIR' directory structure..."

# INITIALIZATION: Create the new directory hierarchy
# The -p flag ensures parent directories are created as needed without errors
mkdir -p "$OUT_DIR/ftr"
mkdir -p "$OUT_DIR/field"
mkdir -p "$OUT_DIR/field_avg"
mkdir -p "$OUT_DIR/post"
mkdir -p "$OUT_DIR/post_avg"
mkdir -p "$OUT_DIR/grid"
mkdir -p "$OUT_DIR/ibmpre"

echo "Initialization complete. Ready for new simulation."