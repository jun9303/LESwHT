#!/bin/bash
set -euo pipefail

PYTHON=${PYTHON:-python3}

echo "Running averaged postprocessor..."
$PYTHON post_avg.py

echo "Done. Outputs are in ../../output/post_avg"
