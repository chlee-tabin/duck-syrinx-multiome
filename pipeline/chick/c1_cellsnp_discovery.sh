#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-8:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/c1_discovery_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/c1_discovery_%j.err
#SBATCH -J chk_disc
#
# c1 — cellSNP-lite mode-2b BULK variant discovery. Mirrors the duck
# 20240729_14_cellsnp-lite-mode2b.sh, tuned for the chick pool.
#
# v2 fix (2026-06-05): the first attempt (3 BAMs, no depth cap, -p short 3h) hit the
# 3h wall at 14/34 chroms because bulk pileup over 3 deep GEX BAMs reached >1M depth
# at highly-expressed loci. Fixes:
#   - DISCOVER on ARC01+ARC02 only (both lanes carry all 12 embryos; ARC0S is an L-only
#     subset of the SAME embryos -> contributes no new donor-discriminating variants).
#     ARC0S is still GENOTYPED at the discovered sites in c2.
#   - --maxDEPTH 100000 caps the memory-hog loci (genotyping needs nowhere near 1M reads).
#   - -p short -t 0-8:00, 64G — O2 short allows <=12h; headroom to finish all chroms.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux
OUT="$WORK/cellsnp/discovery"
CHROM=$(cat "$WORK/ref/chroms.txt")
BAMS="$WORK/bams/ARC01/gex_possorted_bam.bam,$WORK/bams/ARC02/gex_possorted_bam.bam"
rm -rf "$OUT"; mkdir -p "$OUT"   # discovery output is regenerated; clear any partial run

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== c1 cellSNP mode-2b discovery (v2) ==="
echo "BAMs: $BAMS"
echo "CHROM: $CHROM"
cellsnp-lite --version 2>&1 | head -1 || true

cellsnp-lite \
    -s "$BAMS" \
    -O "$OUT" \
    -p 8 \
    --minMAF 0.05 \
    --minCOUNT 100 \
    --maxDEPTH 100000 \
    --cellTAG None \
    --UMItag UB \
    --gzip \
    --chrom "$CHROM"

echo "=== c1 done ==="
echo "Discovered variants: $(zcat "$OUT/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#')"
