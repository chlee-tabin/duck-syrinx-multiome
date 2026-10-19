#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-3:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s2_geno_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s2_geno_%j.err
#SBATCH -J g0s_geno
#
# g0s2 — GEX0S-chick cellSNP mode-1 genotyping at the discovered sites (CellRanger barcodes).
# Single library => the mode-1 output IS the vireo cellData (no cross-library merge). We copy
# it into cellsnp/combined/ (the conventional vireo input path) and re-suffix the barcodes with
# the lane tag '-0S' so per-cell joins onto the sex/QC tables are consistent with the rest of
# the chick pipeline (which suffixes cells by lane).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

WORK=${DUCK_SCRATCH_USER}/chick_demux/gex0s_demux
REGION="$WORK/cellsnp/discovery/cellSNP.base.vcf.gz"
BAM="$WORK/bams/GEX0S-chick_GEX/possorted_genome_bam.bam"
BC="$WORK/barcodes/GEX0S-chick_GEX_raw.tsv"
OUT="$WORK/cellsnp/GEX0S-chick_GEX"; mkdir -p "$OUT"
COMB="$WORK/cellsnp/combined"; mkdir -p "$COMB"

module purge; module load gcc/14.2.0 conda/miniforge3; conda activate vireo
echo "=== g0s2 genotyping GEX0S-chick ($(wc -l < "$BC") cells) at $(zcat "$REGION"|grep -vc '^#') sites ==="
cellsnp-lite -s "$BAM" -b "$BC" -R "$REGION" -O "$OUT" -p 8 \
    --minMAF 0 --minCOUNT 0 --cellTAG CB --UMItag UB --gzip

echo "=== assembling single-library cellData -> $COMB (lane suffix -0S) ==="
# Copy the sparse matrices + base VCF unchanged; only re-suffix the samples list.
for f in cellSNP.tag.AD.mtx.gz cellSNP.tag.DP.mtx.gz cellSNP.tag.OTH.mtx.gz cellSNP.base.vcf.gz; do
  [ -f "$OUT/$f" ] && cp -f "$OUT/$f" "$COMB/$f"
done
awk 'NF{print $1"-0S"}' "$OUT/cellSNP.samples.tsv" > "$COMB/cellSNP.samples.tsv"

echo "=== g0s2 done ==="
echo "cellData cells:   $(wc -l < "$COMB/cellSNP.samples.tsv")"
echo "cellData variants:$(zcat "$COMB/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#')"
