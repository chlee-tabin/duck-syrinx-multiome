#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH13_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH13_%j.err
#SBATCH -J h4_perm_HH13

# Permutation null at HH13 INFERENCE_POOL (Early LPM + Cardiopharyngeal LPM).
# Mirrors run_h4_HH25.sh / run_h4_permutation.sh (HH17/16).
#
# Per-peak empirical p is the primary deliverable: drives the permutation-
# corrected genome-wide L-vs-R UpSet summary, as the
# replacement for the binom_p-based g2c classification (which is mis-calibrated
# under the asymmetric pseudobulk depth at HH25).
#
# Note: HH13 has thin pseudobulk depth at floor=50 (only 1L+2R Early_LPM, fails
# the "need at least 2 per side" guard).  Use floor=30 instead, which gives
# 6L+5R Early_LPM = 11 pseudobulks (still no Cardiopharyngeal_LPM at HH13 in
# either floor, but Early_LPM is the dominant population at this stage).
#
# Expected wall: ~5-15 min.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# COUNTS_PATH override → use floor=30 counts (instead of default floor=50);
# h4_permutation_test.py honors COUNTS_PATH (see env block at top of script).
N_PERM=${N_PERM:-1000} STAGE=${STAGE:-HH13} \
COUNTS_PATH="$SCRATCH/phase2/atac/da_archr/pseudobulk_counts_floor30.tsv.gz" \
  python3 pipeline/h4_permutation_test.py

date
echo "=== DONE ==="
