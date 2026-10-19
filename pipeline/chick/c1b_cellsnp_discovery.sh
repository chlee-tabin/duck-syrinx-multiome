#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-8:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux_B/logs/c1b_discovery_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux_B/logs/c1b_discovery_%j.err
#SBATCH -J chkB_disc
#
# c1b — Group B (GEX09/GEX10 GFP/WT syrinx) cellSNP mode-2b bulk discovery.
# Same recipe as the Group A v2 fix (depth-capped, 2 BAMs, short/8h/64G). Both lanes
# carry the same 12 embryos (L/R halves), so discovering on both captures all donors.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux_B
OUT="$WORK/cellsnp/discovery"
CHROM=$(cat "$WORK/ref/chroms.txt")
BAMS="$WORK/bams/GEX09_GEX/possorted_genome_bam.bam,$WORK/bams/GEX10_GEX/possorted_genome_bam.bam"
rm -rf "$OUT"; mkdir -p "$OUT"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== c1b cellSNP mode-2b discovery (Group B) ==="
echo "BAMs: $BAMS"; echo "CHROM: $CHROM"
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

echo "=== c1b done ==="
echo "Discovered variants: $(zcat "$OUT/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#')"
