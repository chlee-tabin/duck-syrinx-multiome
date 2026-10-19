#!/bin/bash
# cr9 Fig-6a substrate: per-cell 6-ER log-norm + global cluster + demux from the master.
#   sbatch --parsable run_cr9_er_survey.sh
#SBATCH -J cr9_er
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/er_survey_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/er_survey_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# NOTE: these MUST sit below the #SBATCH block — sbatch stops parsing directives
# at the first executable line.
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module purge; module load gcc/14.2.0 R/4.4.2
# Seurat lives in the user lib; Rscript --vanilla skips .Renviron, so set it
# explicitly (the R script also prepends this via .libPaths as the primary guard).
export R_LIBS_USER="$HOME/R/4.4.2"
# payload repointed 2026-08-26 (review fix): execute the REPO copy synced by
# the repo-to-scratch sync, not the hand-staged old_duck_cr9/analysis tree
Rscript --vanilla ${DUCK_SCRATCH_ROOT}/pipeline/scrna_e6e9e11/cr9_export_er_survey.R
RC=$?
echo "RC=$RC END=$(date)"
exit $RC
