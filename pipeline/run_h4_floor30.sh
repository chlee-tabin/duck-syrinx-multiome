#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_floor30_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_floor30_%j.err
#SBATCH -J h4_perm_floor30

# Re-run h4 permutation test against the floor=30 pseudobulk set to test
# whether the directional L-bias claim survives the more permissive
# cell-count cutoff.  Floor=50 gave 18L/10R imbalance + empirical p=0.35;
# floor=30 gives 25L/17R (more balanced).  If directional claim is real,
# empirical p should drop substantially; if it's a methodology artefact,
# p stays ≈0.5.
# (2026-08-11: p=0.35 was the pre-fix tail-biased estimator — see
# pipeline/_lr_binom.py; corrected HH17/16 values in
# figures/data/h4_summary_HH17_16.tsv: count-level emp_p L=0.91 / R=0.61.)

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

# h4_permutation_test.py picks up COUNTS_PATH from env; write outputs to
# a separate dir so we don't overwrite floor=50 results.
N_PERM=1000 \
    STAGE=HH17_16 \
    COUNTS_PATH="$SCRATCH/phase2/atac/da_archr/pseudobulk_counts_floor30.tsv.gz" \
    OUT_DIR="$SCRATCH/phase2/atac/perm_test_floor30" \
    python3 pipeline/h4_permutation_test.py

date
echo "=== DONE ==="
