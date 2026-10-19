#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/pilot3_af_correlation/slurm_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/pilot3_af_correlation/slurm_%j.err
#SBATCH -J pilot3_af

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

mkdir -p ${DUCK_SCRATCH_ROOT}/pilot3_af_correlation

python ${DUCK_SCRATCH_ROOT}/pilot3_af_correlation.py
