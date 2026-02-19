#!/bin/bash

# Exit immediately if any command fails
set -e

# Reset before running the full workflow
echo ""
echo "[1/3] Cleaning up old builds and output data..."
bash reset.sh

echo "========================================"
echo "                 LESwHT                 "
echo "========================================"

# Run the Preprocessor
echo ""
echo "[2/3]Executing Preprocessing Scripts..."
cd 01_pre_processor
bash run_grid.sh
bash run_preprocessing.sh
cd ..

# Run the Solver
echo ""
echo "[3/3] Executing Solver..."
cd 02_solver
bash run_solver.sh
cd ..

echo ""
echo "========================================"
echo "               Completed!               "
echo "========================================"
