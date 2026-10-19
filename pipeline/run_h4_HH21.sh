#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH21_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH21_%j.err
#SBATCH -J h4_perm_HH21

# Permutation null at HH21 INFERENCE_POOL (Ventral spl. + Syringeal + BARX1+).
# Mirrors run_h4_HH25.sh / run_h4_permutation.sh (HH17/16).
#
# Per-peak empirical p is the primary deliverable: drives the permutation-
# corrected genome-wide L-vs-R UpSet summary, as the
# replacement for the binom_p-based g2c classification.
#
# HH21 has the largest pseudobulk N at INFERENCE_POOL: ~19L + 18R (37 total)
# which is roughly balanced. Permutation should mostly confirm the binom-
# based 13K L-biased / 0 R-biased call.  (Pre-fix tail-biased estimator — see
# pipeline/_lr_binom.py; the 13K/0 monodirectionality was the artefact.
# Corrected re-run in figures/data/h4_summary_HH21.tsv: 12,135 L / 9,722 R at
# raw p<0.05, count-level emp_p L=0.37 / R=0.83.)
#
# Expected wall: ~30-60 min.

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

N_PERM=${N_PERM:-1000} STAGE=${STAGE:-HH21} python3 pipeline/h4_permutation_test.py

date
echo "=== DONE ==="
