#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/step02e_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/step02e_%j.err
#SBATCH -J step02e_preCB

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


echo "=== Step 2e: Pre-CB Genotype Filter ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

export NUMBA_CACHE_DIR="/tmp/numba_cache_$$"
mkdir -p "$NUMBA_CACHE_DIR"

cd ${DUCK_SCRATCH_ROOT}

# Step 1: Produce updated genotype map (lightweight, ~30s)
echo ""
echo "--- Producing cb_genotype_map_v2 ---"
python3 20260326_35_produce_cb_genotype_map_v2.py

# Step 2: Apply genotype filter to pre-CB h5ad (heavy, ~30min)
echo ""
echo "--- Applying genotype filter to pre-CB h5ad ---"
python3 20260326_step02e_preCB_genotype.py

echo ""
echo "Done: $(date)"
