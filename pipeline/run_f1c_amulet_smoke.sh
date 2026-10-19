#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=48G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/amulet/f1c_smoke_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/amulet/f1c_smoke_%j.err
#SBATCH -J f1c_amulet_smoke

# AMULET doublet detection — smoke test on a single ARC library.
# Defaults to ARC07 (mid-size, well-behaved by f1 QC).  Override with
# F1C_SMOKE_LIBRARY=<lib> environment variable on the sbatch command.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


LIBRARY="${F1C_SMOKE_LIBRARY:-ARC07}"
cd ${DUCK_SCRATCH_ROOT}

module purge
module load gcc/14.2.0 conda/miniforge3 java/jdk-23.0.1
conda activate snapatac2

mkdir -p ${DUCK_SCRATCH_ROOT}/phase2/atac/amulet

echo "=== f1c AMULET smoke: library=$LIBRARY ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
python3 pipeline/f1c_amulet_doublets.py --library "$LIBRARY"
echo "=== Done ==="
echo "Finished: $(date)"
