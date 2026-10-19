#!/bin/bash
# cr9 compendium annotation — U-method UEG + Wilcoxon + signed-lineage scoring on
# the cr9 substage per-stage objects.  Same script as cr6; SUBSTAGE_DIR points it
# at the cr9 substage dir.  Writes substage_compendium_clusterlabels.tsv +
# substage_E*_percell_compendium.tsv.gz under old_duck_cr9/substage/.
#   sbatch --parsable run_cr9_compendium.sh   # after cr9 substage objects exist
#SBATCH -J cr9_compendium
#SBATCH -p short
#SBATCH -t 0-02:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/compendium_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/compendium_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# NOTE: these MUST sit below the #SBATCH block — sbatch stops parsing directives
# at the first executable line.
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module purge; module load gcc/14.2.0 R/4.4.2
export SUBSTAGE_DIR=${DUCK_SCRATCH_USER}/old_duck_cr9/substage
# payload repointed 2026-08-26 (review fix): execute the REPO copy synced by
# the repo-to-scratch sync, not the hand-staged old_duck_cr9/analysis tree
Rscript --vanilla ${DUCK_SCRATCH_ROOT}/pipeline/scrna_e6e9e11/substage_compendium_annotate.R
RC=$?
echo "RC=$RC END=$(date)"
exit $RC
