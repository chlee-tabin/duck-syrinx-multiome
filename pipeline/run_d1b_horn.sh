#!/bin/bash
#SBATCH -J d1b_horn
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/d1b_horn_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/d1b_horn_%j.err
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


echo "=== d1b Horn's PC Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

REPO="${DUCK_SCRATCH_ROOT}"
cd "$REPO"

python pipeline/d1b_horn_pc_analysis.py

echo "--- Running annotation sweep ---"
python pipeline/d2e_annotate_resolution_sweep.py

echo "End: $(date)"
