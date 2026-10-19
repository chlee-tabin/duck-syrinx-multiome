#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/step07c_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/step07c_%j.err

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cd ${DUCK_SCRATCH_ROOT}
python FC_08656/20260329_step07c_pseudobulk_pca_subsets.py
