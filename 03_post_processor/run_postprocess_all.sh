#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"

if [ -f "${PROJECT_ROOT}/output/running_case.env" ]; then
    # shellcheck source=/dev/null
    source "${PROJECT_ROOT}/output/running_case.env"
fi

export LESWHT_CASE_NAME="${LESWHT_CASE_NAME:-}"
export LESWHT_OUTPUT_ROOT="${LESWHT_OUTPUT_ROOT:-${PROJECT_ROOT}/output}"
export LESWHT_GEOMETRY_STL="${LESWHT_GEOMETRY_STL:-}"

INST_DIR="${LESWHT_OUTPUT_ROOT}/field"
AVG_DIR="${LESWHT_OUTPUT_ROOT}/field_avg"

if [ ! -d "${INST_DIR}" ]; then
    echo "No instantaneous field directory found: ${INST_DIR}"
    exit 0
fi

if [ ! -d "${AVG_DIR}" ]; then
    echo "No averaged field directory found: ${AVG_DIR}"
    exit 0
fi

mapfile -t INST_FIELDS < <(find "${INST_DIR}" -maxdepth 1 -type f -name 'fld*' -printf '%f\n' | sort)

AVG_MANIFEST="${AVG_DIR}/fav_manifest.dat"
AVG_FIELDS=()
DISCOVERED_AVG_FIELDS=()

if [ -s "${AVG_MANIFEST}" ]; then
    while IFS= read -r line; do
        line="${line%%#*}"
        line="${line%%!*}"
        line="$(echo "${line}" | xargs)"
        if [ -n "${line}" ]; then
            AVG_FIELDS+=("${line}")
        fi
    done < "${AVG_MANIFEST}"
fi

while IFS= read -r fpath; do
    rel="${fpath#${AVG_DIR}/}"
    DISCOVERED_AVG_FIELDS+=("${rel}")
done < <(find "${AVG_DIR}" -type f \( -name 'fav*' -o -name 'field_avg.bin' \) ! -name 'fav_manifest.dat' | sort)

if [ ${#AVG_FIELDS[@]} -eq 0 ]; then
    AVG_FIELDS=("${DISCOVERED_AVG_FIELDS[@]}")
elif [ ${#DISCOVERED_AVG_FIELDS[@]} -gt 0 ]; then
    # Keep manifest ordering first, then append any additional discovered files
    # (e.g., field_avg/combined/*) while removing duplicates.
    mapfile -t AVG_FIELDS < <(printf '%s\n' "${AVG_FIELDS[@]}" "${DISCOVERED_AVG_FIELDS[@]}" | awk 'NF { if (!seen[$0]++) print $0 }')
fi

if [ ${#INST_FIELDS[@]} -eq 0 ]; then
    echo "No instantaneous field files found under ${INST_DIR}."
else
    echo "Found ${#INST_FIELDS[@]} instantaneous field file(s)."
fi

if [ ${#AVG_FIELDS[@]} -eq 0 ]; then
    echo "No averaged field files found under ${AVG_DIR}."
else
    echo "Found ${#AVG_FIELDS[@]} averaged field file(s)."
fi

if [ ${#INST_FIELDS[@]} -eq 0 ] && [ ${#AVG_FIELDS[@]} -eq 0 ]; then
    echo "Nothing to postprocess."
    exit 0
fi

make_input_with_fields() {
    local template_path="$1"
    local output_path="$2"
    local ind_film="$3"
    shift 3

    python3 - "$template_path" "$output_path" "$ind_film" "$@" <<'PY'
import sys
from pathlib import Path

template = Path(sys.argv[1])
out = Path(sys.argv[2])
ind_film = int(sys.argv[3])
fields = sys.argv[4:]

lines = template.read_text(encoding='utf-8').splitlines()
anchor = None
for i, line in enumerate(lines):
    if line.strip().startswith('ANIFLD'):
        anchor = i
        break

if anchor is None:
    raise RuntimeError(f'ANIFLD section not found in {template}')

new_lines = lines[: anchor + 1]
new_lines.append(f"{len(fields)}       {ind_film}")
new_lines.extend(fields)
out.write_text("\n".join(new_lines) + "\n", encoding='utf-8')
PY
}

INST_INPUT="${ROOT_DIR}/inst/post_inst.input"
AVG_INPUT="${ROOT_DIR}/avg/post_avg.input"
INST_TMP="$(mktemp)"
AVG_TMP="$(mktemp)"
INST_BAK=""
AVG_BAK=""

cleanup() {
    rm -f "${INST_TMP}" "${AVG_TMP}"
    if [ -n "${INST_BAK}" ] && [ -f "${INST_BAK}" ]; then
        cp "${INST_BAK}" "${INST_INPUT}"
        rm -f "${INST_BAK}"
    fi
    if [ -n "${AVG_BAK}" ] && [ -f "${AVG_BAK}" ]; then
        cp "${AVG_BAK}" "${AVG_INPUT}"
        rm -f "${AVG_BAK}"
    fi
}
trap cleanup EXIT

if [ ${#INST_FIELDS[@]} -gt 0 ]; then
    make_input_with_fields "${INST_INPUT}" "${INST_TMP}" 1 "${INST_FIELDS[@]}"
    INST_BAK="$(mktemp)"
    cp "${INST_INPUT}" "${INST_BAK}"
    cp "${INST_TMP}" "${INST_INPUT}"

    echo "Running instantaneous postprocessing for all fields..."
    (
        cd "${ROOT_DIR}/inst"
        bash run_postprocess.sh
    )
fi

if [ ${#AVG_FIELDS[@]} -gt 0 ]; then
    make_input_with_fields "${AVG_INPUT}" "${AVG_TMP}" 1 "${AVG_FIELDS[@]}"
    AVG_BAK="$(mktemp)"
    cp "${AVG_INPUT}" "${AVG_BAK}"
    cp "${AVG_TMP}" "${AVG_INPUT}"

    echo "Running averaged postprocessing for all fields..."
    (
        cd "${ROOT_DIR}/avg"
        bash run_postprocess.sh
    )
fi

echo "Postprocessing orchestration complete for case root: ${LESWHT_OUTPUT_ROOT}"
