#!/bin/bash
set -euo pipefail

PYTHON=${PYTHON:-python3}

echo "Running instantaneous postprocessor..."
$PYTHON post_inst.py

echo "Done. Outputs are in ../../output/post_inst"
