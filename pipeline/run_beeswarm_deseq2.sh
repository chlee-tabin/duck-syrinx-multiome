#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:15
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/beeswarm_deseq2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/beeswarm_deseq2_%j.err
#SBATCH -J beeswarm_deseq2
# Beeswarm DESeq2 stats (deps-light; DESeq2 not available in the local env).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 R/4.4.2
Rscript figures/figS03/beeswarm_deseq2_stats.R || { echo "FATAL: beeswarm DESeq2 failed" >&2; exit 2; }
echo "COMPLETE_OK"
