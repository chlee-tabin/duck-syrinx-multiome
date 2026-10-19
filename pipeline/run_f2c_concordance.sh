#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2c_concordance_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2c_concordance_%j.err
#SBATCH -J f2c_concordance

# Phase 1.3 (per fig04_chromatin_exploratory_plan.md): SnapATAC2 (f2)
# vs ArchR union peak set cross-validation.  Runs after f2
# emits union_peaks.bed.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f2c_archr_snapatac_concordance ==="
date

python3 pipeline/f2c_archr_concordance.py

echo "=== DONE ==="
date
