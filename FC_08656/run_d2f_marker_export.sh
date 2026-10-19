#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:30
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/d2f_marker_export_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/d2f_marker_export_%j.err
#SBATCH -J d2f_markers

# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

REPO="${DUCK_SCRATCH_ROOT}"
cd "$REPO"

python3 pipeline/d2f_canonical_marker_export.py
