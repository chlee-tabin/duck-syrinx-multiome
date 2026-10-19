#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:30
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/donor_consolidation_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/donor_consolidation_%j.err
#SBATCH -J donor_consol

# Donor consolidation: derive unified donor mapping from evidence files.
# Reads cross-anchoring, GT distance, barcode concordance outputs.
# Produces: unified_donor_map.tsv (the ground truth for 52 embryos).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

mkdir -p ${DUCK_SCRATCH_ROOT}/logs

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

SCRATCH=${DUCK_SCRATCH_ROOT}
python3 "${SCRATCH}/pipeline/b1_donor_consolidation.py"

echo ""
echo "=== Evidence files that should be backed up to standby ==="
echo "These are the inputs to donor_consolidation.py:"
for f in \
    "$SCRATCH/libpair_g6/merged_reference/donor_mapping.tsv" \
    "$SCRATCH/unified_donor_table.tsv" \
    "$SCRATCH/concordance_matrix.tsv" \
    "$SCRATCH/gt_distance_all/within_g2_corr.tsv" \
    "$SCRATCH/genotype_distance_output/g2_pairwise_corr.tsv" \
    "$SCRATCH/arc_gex_concordance/concordance_v2.tsv" \
    "$SCRATCH/cross_anchor_output/unified_donor_map.tsv" \
    "$SCRATCH/cross_anchor_output/g5_g6_matches.tsv" \
    "$SCRATCH/universal_genotyping/identity_table/unified_donor_map.tsv" \
; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  [OK]  $SIZE  $(basename "$f")"
    else
        echo "  [--]  NOT FOUND  $(basename "$f")"
    fi
done
