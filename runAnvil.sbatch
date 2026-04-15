#!/bin/bash
#SBATCH --job-name=dimple
#SBATCH --partition=shared
#SBATCH --account=phy250071
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=96:00:00
#SBATCH --output=./.sout/slurm-%x-%j.out
#SBATCH --error=./.sout/slurm-%x-%j.err

set -euo pipefail

# Under Slurm, the batch script can execute from a spool path (e.g. /var/spool/slurm/job*/).
# Use submission directory as project root so relative project files resolve correctly.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  SUBMIT_ROOT="$SLURM_SUBMIT_DIR"
else
  SUBMIT_ROOT="$(pwd)"
fi

if [ ! -f "$SUBMIT_ROOT/run.sh" ]; then
  echo "Error: run.sh not found in SUBMIT_ROOT=$SUBMIT_ROOT"
  echo "Hint: submit from repo root, e.g., cd /anvil/projects/x-phy250071/LESwHT && sbatch run_anvil_shared.sh ..."
  exit 1
fi

JOB_LOCAL_ROOT="${SLURM_TMPDIR:-/tmp/${USER}/slurm_${SLURM_JOB_ID:-manual}}"
PROJECT_BASENAME="$(basename "$SUBMIT_ROOT")"
WORK_ROOT="$JOB_LOCAL_ROOT/$PROJECT_BASENAME"

echo "========================================"
echo "       LESwHT on Anvil (Shared)         "
echo "========================================"
echo "Job ID            : ${SLURM_JOB_ID:-N/A}"
echo "Partition         : ${SLURM_JOB_PARTITION:-N/A}"
echo "Nodes             : ${SLURM_JOB_NUM_NODES:-1}"
echo "CPUs per task     : ${SLURM_CPUS_PER_TASK:-64}"
echo "Requested walltime: 96:00:00"
echo "Submit root       : $SUBMIT_ROOT"
echo "Job local root    : $WORK_ROOT"

# OpenMP policy for full shared-node usage.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=dynamic
export OMP_STACKSIZE=1G

ulimit -s unlimited

echo "Staging project to job-local storage..."
mkdir -p "$WORK_ROOT"
rsync -a --delete \
  --exclude '.git/' \
  --exclude '.venv/' \
  --exclude '.sout/' \
  --exclude 'output/' \
  "$SUBMIT_ROOT/" "$WORK_ROOT/"
mkdir -p "$WORK_ROOT/output"

sync_output_back() {
  echo "Syncing output back to submit directory..."
  mkdir -p "$SUBMIT_ROOT/output"
  rsync -a "$WORK_ROOT/output/" "$SUBMIT_ROOT/output/"
}

JOB_STATUS=0
cd "$WORK_ROOT"

# Always use virtual environment from submit directory (do not stage/copy .venv).
if [ -f "$SUBMIT_ROOT/.venv/bin/activate" ]; then
  # shellcheck source=/dev/null
  source "$SUBMIT_ROOT/.venv/bin/activate"
fi

bash "$WORK_ROOT/run.sh" "$@" || JOB_STATUS=$?
sync_output_back

exit "$JOB_STATUS"
