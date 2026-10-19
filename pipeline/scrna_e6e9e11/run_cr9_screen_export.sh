#!/bin/bash
# Full WNT/BMP/TGF-b screen (S8/S9/S10 gene panels) per-cell re-export from cr9
# substage objects. One-off, small: load 3 RDS + extract ~50 genes + write. priority queue.
#   sbatch --parsable run_cr9_screen_export.sh
#SBATCH -J cr9_screen_export
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/screen_export_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/screen_export_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# NOTE: these MUST sit below the #SBATCH block — sbatch stops parsing directives
# at the first executable line.
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module purge; module load gcc/14.2.0 R/4.4.2
# Seurat v5 lives in the personal lib; --vanilla skips .Renviron (which points
# R_LIBS_USER at ~/R/library_dropest, no Seurat), so export the right lib here.
export R_LIBS_USER="$HOME/R/4.4.2"
# payload repointed 2026-08-26 (review fix): execute the REPO copy synced by
# the repo-to-scratch sync, not the hand-staged old_duck_cr9/analysis tree
Rscript --vanilla ${DUCK_SCRATCH_ROOT}/pipeline/scrna_e6e9e11/cr9_screen_export.R
RC=$?
echo "RC=$RC END=$(date)"
exit $RC
