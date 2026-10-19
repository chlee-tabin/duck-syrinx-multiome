#!/bin/bash
#SBATCH -J hh25_pcc
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/hh25_pcc_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/hh25_pcc_%j.err
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "host=$(hostname) start=$(date)"
python3 ${DUCK_SCRATCH_ROOT}/scripts/export_hh25_percell_counts.py
echo "PCC_DONE end=$(date)"
