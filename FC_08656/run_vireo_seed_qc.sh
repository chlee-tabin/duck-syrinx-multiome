#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/vireo_seed_qc_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/vireo_seed_qc_%j.err
#SBATCH -J seed_qc

# Vireo seed QC: rank seeds, compare to baseline, compute GT distance.
# Run AFTER all vireo seeded rerun jobs have completed.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

mkdir -p ${DUCK_SCRATCH_ROOT}/logs

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

SCRATCH=${DUCK_SCRATCH_ROOT}
python3 "$SCRATCH/20260330_vireo_seed_qc.py"
