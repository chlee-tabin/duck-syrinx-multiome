#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f5_cfc1_quant_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f5_cfc1_quant_%j.err
#SBATCH -J f5_cfc1_quant

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f5 CFC1 region quantification ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
# All 4 multiome stages — HH13 is where the mRNA R>L finding
# lives (padj=0.007), and the chromatin trajectory needs the flanking
# stages (HH17/16, HH21, HH25) to distinguish "transcription tracks
# chromatin" from "chromatin-primed" or "chromatin-memory" outcomes.
python3 pipeline/f5_cfc1_region_counts.py \
    --include-stages HH13,HH21,HH25,HH17,HH16,HH17/16 \
    "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
