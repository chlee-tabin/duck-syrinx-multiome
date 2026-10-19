#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_%j.err
#SBATCH -J h4_permutation

# Item (3) of overnight chromatin queue: empirical permutation null for
# the (pre-2026-08-11) 9,202 L-biased vs 0 R-biased monodirectional claim at HH17/16 —
# since shown to be an artifact of the tail-biased binomial (see
# h4_permutation_test.py header); the script now asks whether nominal L/R counts exceed
# the shuffle null in either direction
# INFERENCE_POOL (raw p<0.05, 283,818 ArchR peaks).
# Expected wall: 1-3h.

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

N_PERM=${N_PERM:-1000} STAGE=${STAGE:-HH17_16} python3 pipeline/h4_permutation_test.py

date
echo "=== DONE ==="
