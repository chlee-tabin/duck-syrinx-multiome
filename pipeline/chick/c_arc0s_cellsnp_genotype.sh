#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-3:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/c_arc0s_geno_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/c_arc0s_geno_%j.err
#SBATCH -J chk_arc0s
#
# c_arc0s — cellSNP-lite mode-1 per-cell genotyping of the 2024 ARC0S-chick
# "symmetry-breaking" lane (LEFT halves of HH15 + HH25 embryos), at the SAME
# 283,435 variant positions as the Group A (ARC01+ARC02) combined cellData.
#
# WHY the combined base (not the discovery base): the Group A donor genotypes
# (vireo_seeded/N18/seed_117/GT_donors.vireo.vcf.gz) were learned on this exact
# 283,435-position set. To project ARC0S cells onto those FIXED donors with
# `vireo --donorFile`, the ARC0S cellData rows must be position-aligned to the
# donor VCF — i.e. genotype at cellsnp/combined/cellSNP.base.vcf.gz.
#
# Output is a single-library cellData (no merge needed) -> consumed directly by
# c_arc0s_vireo_donorfile.sh.
#
# NOTE: ARC0S has significant DUCK contamination (flagged by the 2024 analyst).
# That is EXPECTED and is the point of the --donorFile projection: chick cells
# land on the 12 Group A donors; duck-contaminant cells get low-prob / unassigned.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux
REGION="$WORK/cellsnp/combined/cellSNP.base.vcf.gz"   # 283,435 variants == donor-VCF positions
BAM="$WORK/bams/ARC0S-chick/gex_possorted_bam.bam"    # staged from standby FC_08565 (41 GB)
BC="$WORK/barcodes/ARC0S-chick_raw.tsv"               # 2,956 CellRanger cells (raw CB-matching, suffix -1)
OUT="$WORK/arc0s/cellsnp"
mkdir -p "$OUT" "$WORK/logs"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== c_arc0s cellSNP mode-1 genotyping: ARC0S-chick ==="
echo "Region VCF: $REGION ($(zcat "$REGION" | grep -vc '^#') sites)"
echo "BAM:        $BAM"
echo "Barcodes:   $BC ($(wc -l < "$BC") cells)"
test -s "$BAM" || { echo "ERROR: BAM missing/empty (rsync from standby not finished?)"; exit 1; }

cellsnp-lite \
    -s "$BAM" \
    -b "$BC" \
    -R "$REGION" \
    -O "$OUT" \
    -p 8 \
    --minMAF 0 \
    --minCOUNT 0 \
    --cellTAG CB \
    --UMItag UB \
    --gzip

echo "=== c_arc0s done -> $OUT ==="
echo "Genotyped cells: $(wc -l < "$OUT/cellSNP.samples.tsv" 2>/dev/null || echo '?')"
echo "Variant rows:    $(zcat "$OUT/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#' || echo '?')"
