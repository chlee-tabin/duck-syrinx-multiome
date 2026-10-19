#!/bin/bash
# cr9 substage subclustering — per-stage CP10K+ccreg(+sex.score)+Harmony+Louvain
# on the cr9-CB master.  Writes old_duck_cr9/substage/ (separate from cr6).
#   sbatch --parsable run_cr9_substage.sh   # after route_b master is built
#SBATCH -J cr9_substage
#SBATCH -p short
#SBATCH -t 0-08:00
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/substage_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/substage_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# NOTE: these MUST sit below the #SBATCH block — sbatch stops parsing directives
# at the first executable line.
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module purge; module load gcc/14.2.0 R/4.4.2
# payload repointed 2026-08-26 (review fix): execute the REPO copy synced by
# the repo-to-scratch sync, not the hand-staged old_duck_cr9/analysis tree
Rscript --vanilla ${DUCK_SCRATCH_ROOT}/pipeline/scrna_e6e9e11/substage_subcluster_cr9.R
RC=$?
echo "RC=$RC END=$(date)"
exit $RC
