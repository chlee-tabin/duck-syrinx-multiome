#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:15
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3_cfc1_deseq2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3_cfc1_deseq2_%j.err
#SBATCH -J g3_cfc1_deseq2

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate jupyter
# DESeq2 lives in /home/cl266/.conda/envs/jupyter/lib/R/library/DESeq2.
# The jupyter env's R is the only one with DESeq2 installed; the
# bare R/4.4.2 module + ~/R/4.4.2 user lib does NOT have DESeq2.

echo "=== g3 CFC1 DESeq2 (unpaired+sex) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
Rscript pipeline/g3_cfc1_deseq2.R "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
