#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_perembryo_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_perembryo_%j.err
#SBATCH -J f4_perembryo

# f4 PER-EMBRYO-HALF compartment pileup — pseudoreplication rectification for
# fig03 panel C (2026-06-17). Produces per-(stage x compartment x GENOTYPE x
# side) PITX2-locus coverage so panel_C_pileup.R can weight each embryo-HALF
# equally (mean over halves) instead of each CELL (the cell-pooled
# f4_compartment coverage let high-cell embryos dominate the TRACK SHAPE).
# The significance STARS are unaffected — they already read the rectified
# per-embryo-half g4f floor-20. --min-cells 20 matches g4f floor-20.
# Needs the f4 compose mod (--compartment-split AND --per-genotype-anchors).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
module purge; module load gcc/14.2.0 conda/miniforge3; conda activate snapatac2
cd "$SCRATCH"
echo "=== f4 per-embryo-half compartment pileup === $(date) on $(hostname)"

# (1) D/V-resolved stages: compartment x genotype x side
python3 pipeline/f4_lr_pitx2_pileup.py \
    --compartment-split \
    --per-genotype-anchors HH17/16,HH21,HH25 \
    --lr-tsv multiome \
    --firmness-also singleton_rank \
    --include-stages HH17,HH16,HH17/16,HH21,HH25 \
    --out-subdir lr_pileup_compartment_perembryo \
    --min-cells 20 \
    --skip-peaks

# (2) HH13 pool (pre-TES, no D/V split): genotype x side
python3 pipeline/f4_lr_pitx2_pileup.py \
    --per-genotype-anchors HH13 \
    --lr-tsv multiome \
    --firmness-also singleton_rank \
    --include-stages HH13 \
    --out-subdir lr_pileup_hh13_perembryo \
    --min-cells 20 \
    --skip-peaks

echo "=== Done === $(date)"
echo "--- outputs ---"
ls -la phase2/atac/lr_pileup_compartment_perembryo/*.tsv* phase2/atac/lr_pileup_hh13_perembryo/*.tsv* 2>&1 | tail
