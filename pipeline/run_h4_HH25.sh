#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH25_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h4_perm_HH25_%j.err
#SBATCH -J h4_perm_HH25

# HH25 sensitivity arm:
# Re-evaluate the HH25 R-bias claim (15,832 R-biased vs 18 L-biased peaks at
# raw p<0.05) — likely the mirror of the HH17/16 L-bias methodology artefact.
# (2026-08-11: artefact confirmed, but the mechanism was the PRE-FIX
# tail-biased estimator — see pipeline/_lr_binom.py; it could only call the
# deeper (R) side at HH25.  Corrected re-run in
# figures/data/h4_summary_HH25.tsv (2026-08 cohort: 10L/15R): 17,677 L /
# 15,610 R at raw p<0.05, count-level emp_p L=0.001 / R=0.016 — read per
# as pseudobulk composition structure, not direction.)
#
# Pseudobulk side imbalance at HH25 INFERENCE_POOL (Ventral spl., Syringeal,
# Chondrogenic LPM): 12L / 21R, total 33 pseudobulks across 8 genotypes (all
# paired).  ("Chondrogenic LPM" is a retired HH25 label — since the
# 2026-07-02 relabel the pool is Ventral spl. + Syringeal; see
# figures/_shared/subcluster_map.py POOLS.)  Permutation-shuffle the side
# labels 1,000 times and ask: how often does the (R-biased - L-biased)
# count under the null exceed the observed?  Expected outcome: empirical
# p ~ 0.4 (artefact confirmed).
#
# Expected wall: ~30-90 min (HH25 has more pseudobulks than HH17/16 + same
# peak count ~283K).

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

# Use floor=50 (primary) counts; STAGE=HH25 triggers POOLS["HH25"] lookup
# in the script for INFERENCE_POOL membership (Ventral spl., Syringeal,
# Chondrogenic LPM — retired label; the live POOLS lookup now resolves to
# Ventral spl. + Syringeal per figures/_shared/subcluster_map.py).  Outputs
# to default OUT_DIR (perm_test) with stage-suffixed filenames.
N_PERM=${N_PERM:-1000} STAGE=${STAGE:-HH25} python3 pipeline/h4_permutation_test.py

date
echo "=== DONE ==="
