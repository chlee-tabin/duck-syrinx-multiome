#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f5_quant_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f5_quant_%j.err
#SBATCH -J f5_quant

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f5 PITX2 region quantification ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
# Explicit --include-stages keeps this runner self-documenting even
# though f5_pitx2_region_counts.py defaults to the same set.  HH17 +
# HH16 fold into the HH17/16 anchor via STAGE_TO_ANCHOR; "HH17/16"
# itself is also accepted in case e1b-fixed assigns the literal anchor
# string.
python3 pipeline/f5_pitx2_region_counts.py \
    --include-stages HH13,HH21,HH25,HH17,HH16,HH17/16 \
    "$@"

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
