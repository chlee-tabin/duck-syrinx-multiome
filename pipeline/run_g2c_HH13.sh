#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-0:15
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_HH13_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_HH13_%j.err
#SBATCH -J g2c_HH13

# Run g2c_aggregated_descriptive.py at HH13 INFERENCE_POOL on ArchR peaks.
# Tests whether HH13 has the predicted L>R chromatin asymmetry
# (NODAL priming signature) — directly addresses §5.3 of hypotheses doc.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== g2c_aggregated HH13 (ArchR) ==="
date
STAGE="HH13" DA_DIR="$SCRATCH/phase2/atac/da_archr" python3 pipeline/g2c_aggregated_descriptive.py

echo ""
echo "=== DONE ==="
date
