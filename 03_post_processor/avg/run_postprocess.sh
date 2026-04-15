#!/bin/bash
set -euo pipefail

if [ -f ../../output/running_case.env ]; then
	# shellcheck source=/dev/null
	source ../../output/running_case.env
fi

export LESWHT_OUTPUT_ROOT="${LESWHT_OUTPUT_ROOT:-../../output}"

PYTHON=${PYTHON:-python3}

echo "Running averaged postprocessor..."
$PYTHON post_avg.py

echo "Done. Outputs are in ${LESWHT_OUTPUT_ROOT}/post_avg"
