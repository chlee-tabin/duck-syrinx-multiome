#!/bin/bash
# Genome-wide syrinx-population pseudobulk RAW-COUNT export per (donor x side)
# at E09/E11, summed from the CellBender filtered h5s (post-purge route —
# replaces run_cr9_syrinx_counts.sh, whose substage RDS inputs were purged;
# that wrapper is retired and not distributed).
# One-off, small: 3 h5 loads + grouped rowsums. priority queue.
#   sbatch --parsable run_cr9_syrinx_counts_h5.sh
#SBATCH -J cr9_syrinx_counts_h5
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/syrinx_counts_h5_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/syrinx_counts_h5_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# NOTE: these MUST sit below the #SBATCH block — sbatch stops parsing directives
# at the first executable line.
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER DUCK_SCRATCH_ROOT

echo "HOST=$(hostname) START=$(date)"
module purge; module load gcc/14.2.0 conda/miniforge3
conda activate vireo || { echo "conda activate vireo FAILED"; exit 1; }
python3 ${DUCK_SCRATCH_ROOT}/pipeline/scrna_e6e9e11/cr9_syrinx_pseudobulk_from_h5.py
RC=$?
echo "RC=$RC END=$(date)"
exit $RC
