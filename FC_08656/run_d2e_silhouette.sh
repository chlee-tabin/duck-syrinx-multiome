#!/bin/bash
#SBATCH -J d2e_sweep
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/d2e_sweep_%j.log
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/d2e_sweep_%j.log

# Resolution sweep: Leiden at res=[0.2, 0.3, 0.4, 0.5]
# Computes silhouette scores + U-method markers at each resolution.
# Output: figures/data/d2e_resolution_sweep_summary.tsv (decision table)

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

SCRATCH="${DUCK_SCRATCH_ROOT}"
mkdir -p "${SCRATCH}/logs" "${SCRATCH}/phase2/figures"

echo "=== D2e: Resolution sweep with U-method ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
python3 "${SCRATCH}/pipeline/d2e_resolution_umethod.py"
echo "Done: $(date)"
