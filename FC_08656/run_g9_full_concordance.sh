#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:30
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/g9_full_concordance_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/g9_full_concordance_%j.err
#SBATCH -J g9_conc

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

mkdir -p ${DUCK_SCRATCH_ROOT}/logs

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

SCRATCH=${DUCK_SCRATCH_ROOT}
python3 "${SCRATCH}/20260330_g9_full_concordance.py"
