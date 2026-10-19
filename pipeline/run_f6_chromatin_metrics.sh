#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f6_chromatin_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f6_chromatin_%j.err
#SBATCH -J f6_chromatin

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f6 chromatin metrics across all stages ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
python3 pipeline/f6_chromatin_metrics.py "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
