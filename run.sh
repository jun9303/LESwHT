#!/bin/bash

# Exit immediately if any command fails
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

TOPOLOGY=0
DEPTH=30.0
SCALE=1.0
STRETCH=1.0

while getopts "t:d:s:x:h" opt; do
	case "$opt" in
		t) TOPOLOGY="$OPTARG" ;;
		d) DEPTH="$OPTARG" ;;
		s) SCALE="$OPTARG" ;;
		x) STRETCH="$OPTARG" ;;
		h)
			echo "Usage: $0 [-t topology] [-d depth] [-s scale] [-x stretch]"
			echo "  -t  Topology ID (0-6)"
			echo "  -d  Dimple depth in wall units"
			echo "  -s  Isotropic scale factor"
			echo "  -x  Streamwise stretch factor"
			exit 0
			;;
		*)
			echo "Invalid option. Use -h for help."
			exit 1
			;;
	esac
done


echo "========================================"
echo "                 LESwHT                 "
echo "========================================"

# Run the Dimple Geometry Generator
echo ""
echo "[1/4] Generating dimple slab geometry..."
if [ -x "$PROJECT_ROOT/.venv/bin/python" ]; then
	PYTHON_CMD="$PROJECT_ROOT/.venv/bin/python"
else
	PYTHON_CMD="python3"
fi

cd "$PROJECT_ROOT/00_dimple_generator"
"$PYTHON_CMD" main.py -t "$TOPOLOGY" -d "$DEPTH" -s "$SCALE" -x "$STRETCH" --compact
cd "$PROJECT_ROOT"

# Run the Preprocessor
echo ""
echo "[2/4] Executing Preprocessing Scripts..."
cd "$PROJECT_ROOT/01_pre_processor"
bash run_grid.sh
bash run_preprocessing.sh
cd "$PROJECT_ROOT"

# Update UDRV_I for CFR using effective channel height from dimple volume.
echo ""
echo "[2.5/4] Updating UDRV_I for CFR from effective channel height..."
if [ -f "$PROJECT_ROOT/output/running_case.env" ]; then
	# shellcheck source=/dev/null
	source "$PROJECT_ROOT/output/running_case.env"
else
	echo "Error: running_case.env not found at $PROJECT_ROOT/output/running_case.env"
	exit 1
fi

export LESWHT_OUTPUT_ROOT="${LESWHT_OUTPUT_ROOT:-$PROJECT_ROOT/output}"
export LESWHT_GEOMETRY_STL="${LESWHT_GEOMETRY_STL:-}"

if [ -z "$LESWHT_GEOMETRY_STL" ] || [ ! -f "$LESWHT_GEOMETRY_STL" ]; then
	echo "Error: LESWHT_GEOMETRY_STL is not set to a valid STL file."
	exit 1
fi

export LESWHT_SETTINGS_FILE="$PROJECT_ROOT/02_solver/settings.input"
if [ ! -f "$LESWHT_SETTINGS_FILE" ]; then
	echo "Error: settings file not found at $LESWHT_SETTINGS_FILE"
	exit 1
fi

"$PYTHON_CMD" - <<'PY'
import os
import trimesh

stl_path = os.environ['LESWHT_GEOMETRY_STL']
settings_file = os.environ['LESWHT_SETTINGS_FILE']

mesh = trimesh.load_mesh(stl_path, process=False)
if isinstance(mesh, trimesh.Scene):
	mesh = trimesh.util.concatenate(tuple(mesh.geometry.values()))

# V_dimple definition:
#   total STL-domain volume (axis-aligned bounds) minus slab-solid volume.
# h_eff definition:
#   V_dimple / (X_stl * Z_stl)
slab_volume = abs(float(mesh.volume))
b = mesh.bounds
total_domain_volume = float((b[1, 0] - b[0, 0]) * (b[1, 1] - b[0, 1]) * (b[1, 2] - b[0, 2]))
v_dimple = max(0.0, total_domain_volume - slab_volume)
x_stl = float(b[1, 0] - b[0, 0])
z_stl = float(b[1, 2] - b[0, 2])

if x_stl <= 0.0 or z_stl <= 0.0:
	raise RuntimeError('Computed non-positive X_stl or Z_stl from STL bounds')

h_eff = v_dimple / (x_stl * z_stl)
if h_eff < 0.0:
	raise RuntimeError('Computed negative h_eff from STL-based dimple volume')

udrv_i = 2.0 / (2.0 + h_eff)
udrv_i = max(0.0, min(1.0, udrv_i))

with open(settings_file, 'r', encoding='utf-8') as f:
	lines = f.read().splitlines()

updated = False
for i, line in enumerate(lines):
	if line.strip().upper().startswith('IRESET'):
		if i + 1 >= len(lines):
			raise RuntimeError('Malformed settings.input around IRESET section')
		vals = lines[i + 1].split()
		if len(vals) < 6:
			raise RuntimeError('Expected 6 values on IRESET data line')
		vals[5] = f'{udrv_i:.10f}'
		lines[i + 1] = '       '.join(vals)
		updated = True
		break

if not updated:
	raise RuntimeError('Could not find IRESET header in settings.input')

with open(settings_file, 'w', encoding='utf-8') as f:
	f.write('\n'.join(lines) + '\n')

print('CFR UDRV_I update complete:')
print(f'  X_stl                     = {x_stl:.10f}')
print(f'  Z_stl                     = {z_stl:.10f}')
print(f'  V_total_domain(STL bbox)  = {total_domain_volume:.10e}')
print(f'  V_slab(STL solid)         = {slab_volume:.10e}')
print(f'  V_dimple                  = {v_dimple:.10e}')
print(f'  h_eff (=V_dimple/X_stl/Z_stl) = {h_eff:.10f}')
print(f'  UDRV_I (=2/(2+h_eff), <=1)    = {udrv_i:.10f}')
PY

# Run the Solver
echo ""
echo "[3/4] Executing Solver..."
cd "$PROJECT_ROOT/02_solver"
bash run_solver.sh
cd "$PROJECT_ROOT"

# Run all Postprocessors
echo ""
echo "[4/4] Executing Postprocessing (all fields)..."
cd "$PROJECT_ROOT/03_post_processor"
bash run_postprocess_all.sh
cd "$PROJECT_ROOT"

echo ""
echo "========================================"
echo "               Completed!               "
echo "========================================"
