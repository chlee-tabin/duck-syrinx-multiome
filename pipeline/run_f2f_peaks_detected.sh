#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2f_peaks_detected_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2f_peaks_detected_%j.err
#SBATCH -J f2f_peaks_detected
#
# Union peaks detected per current (stage × LPM population) on the production
# ArchR peak matrix — the pool-free peak-yield QC behind figS03 B1/B2 (ms Fig S10).
# Loads one 4 GB h5ad + the f1 sidecar; writes one small TSV.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd "${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f2f peaks detected per current population ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

python3 pipeline/f2f_peaks_detected_per_population.py

echo ""
echo "=== done ==="
echo "Finished: $(date)"
