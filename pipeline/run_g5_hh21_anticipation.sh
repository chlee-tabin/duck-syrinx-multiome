#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g5_hh21_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g5_hh21_%j.err
#SBATCH -J g5_hh21_anticipation

# Phase 5 (per fig04_chromatin_exploratory_plan.md): cross-modal
# intersection of HH21 Ventral spl. DEGs ∩ HH17/16 L>R ATAC peaks.
# Conditional on g2_lpm_hh1716 having produced significant peaks.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== g5_hh21_wave2_anticipation ==="
date

python3 pipeline/g5_hh21_wave2_anticipation.py

echo ""
echo "=== DONE ==="
date
