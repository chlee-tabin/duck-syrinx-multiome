#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_hh1716_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4_hh1716_%j.err
#SBATCH -J f4_hh1716_only

# f4 HH17/16-only L/R aggregate pseudobulk pileup.
# Closes the 6→8 gap in figures/data/step_f4_pitx2_coverage.tsv.gz so the
# stage × side coverage strip in the PITX2 multi-zoom schematic includes
# the HH17/16 anchor (pre-TES, key for chromatin laterality story).
#
# Outputs to a separate subdir so the canonical lr_pileup/ stays intact.
# A downstream local Python script merges the new HH17_16_L/R columns
# into the 6-track TSV to produce the 8-track output consumed by the
# multi-zoom panels.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f4 HH17/16-only L/R pileup ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"

# --include-stages: HH17, HH16, and the bilingual "HH17/16" all fold to the
#   HH17/16 anchor via STAGE_TO_ANCHOR.
# --out-subdir lr_pileup_hh1716: dedicated dir; canonical lr_pileup/ stays
#   intact (HH13/HH21/HH25 fragments + per-base coverage from May 4 run).
# --min-cells 100: matches the canonical aggregate-mode default; passing
#   explicitly avoids a NameError in f4_lr_pitx2_pileup.py main() where the
#   default-resolution branch references per_genotype_anchors before it is
#   defined (line 356 refers to a name first bound on line 368).
python3 pipeline/f4_lr_pitx2_pileup.py \
    --include-stages HH17,HH16,HH17/16 \
    --out-subdir lr_pileup_hh1716 \
    --min-cells 100

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
