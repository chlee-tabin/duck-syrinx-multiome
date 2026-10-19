#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_HH13_floor20_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_HH13_floor20_%j.err
#SBATCH -J g2c_HH13_f20

# Re-run g2c_aggregated_descriptive.py at HH13 INFERENCE_POOL using the
# floor=20 pseudobulk counts file (rebuilt by g1 rerun, 2026-05-12).  Goal:
# boost HH13 depth from 1L+2R (floor=50) to ≥3L+3R for the cross-stage UpSet.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== g2c_aggregated HH13 (ArchR, floor=20) ==="
date
STAGE="HH13" \
    DA_DIR="$SCRATCH/phase2/atac/da_archr" \
    COUNTS_PATH="$SCRATCH/phase2/atac/da_archr/pseudobulk_counts_floor20.tsv.gz" \
    python3 pipeline/g2c_aggregated_descriptive.py

echo ""
echo "=== DONE ==="
date
