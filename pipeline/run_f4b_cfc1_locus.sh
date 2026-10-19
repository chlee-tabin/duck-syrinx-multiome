#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4b_cfc1_locus_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4b_cfc1_locus_%j.err
#SBATCH -J f4b_cfc1_locus

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f4b CFC1 locus coverage ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
python3 pipeline/f4b_cfc1_locus_coverage.py "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
