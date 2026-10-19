#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_lr_pileup_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_lr_pileup_%j.err
#SBATCH -J f4_lr_pileup

# f4: L/R pseudobulk fragment export + peak calling + PITX2 locus pileup.
# User reframe (2026-04-24): pool cells by (stage, snRNA-inferred L/R)
# within POOLS[stage]["inference"] to sidestep SnapATAC2's "all replicates"
# reproducibility filter.  See pipeline/f4_lr_pitx2_pileup.py docstring.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f4 L/R pileup ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
python3 pipeline/f4_lr_pitx2_pileup.py "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
