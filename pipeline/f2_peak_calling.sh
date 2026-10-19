#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 0-6:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2_peaks_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2_peaks_%j.err
#SBATCH -J f2_peaks

# Phase F Step 2: MACS3 peak calling per cell type (via SnapATAC2 built-in).
# Prerequisites: f1_coupled_atac.py complete (step_f1_combined.h5ads).
#
# Prefer using run_phase_fg.sh --from f2 --to f2 instead of this script.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
mkdir -p "$SCRATCH/phase2/atac/peaks"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f2: MACS3 Peak Calling (SnapATAC2 built-in) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

python3 pipeline/f2_peak_calling.py

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
