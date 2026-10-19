#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h8_side_shuffle_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h8_side_shuffle_%j.err
#SBATCH -J h8_side_shuffle

# Laterality (side-shuffle) null for the HH17/16 gkm-SVM: re-derive the L-biased
# peak set under 100 global side permutations (18L/10R preserved) and retrain,
# testing whether the real side assignment yields a more learnable signature
# than random splits. AUC_real should land ~0.72 (matches h2 random-neg) as a
# reconstruction check.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"

# h8 has no MODE knob — its cutoffs are fixed module-level constants (P_POS etc.).
N_SHUFFLE=100 python3 pipeline/h8_side_shuffle_null.py

date
echo "=== DONE ==="
