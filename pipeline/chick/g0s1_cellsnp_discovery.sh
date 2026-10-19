#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-8:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s1_disc_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s1_disc_%j.err
#SBATCH -J g0s_disc
#
# g0s1 — GEX0S-chick cellSNP mode-2b de-novo bulk variant discovery.
# Single GEX library (the LEFT-half symmetry-breaking lane). Same recipe as the chick
# Group A/B discovery (c1/c1b): depth-capped genome-wide mode-2b, default minMAF 0.05 /
# minCOUNT 100. The discovered base is shared across all cells of this one lane, so c2
# genotyping rows align and no cross-library merge is needed (single library).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux/gex0s_demux
OUT="$WORK/cellsnp/discovery"
CHROM=$(cat "$WORK/ref/chroms.txt")
BAM="$WORK/bams/GEX0S-chick_GEX/possorted_genome_bam.bam"
rm -rf "$OUT"; mkdir -p "$OUT"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== g0s1 cellSNP mode-2b discovery (GEX0S-chick): minMAF 0.05 minCOUNT 100 ==="
echo "BAM: $BAM"; echo "CHROM: $CHROM"
cellsnp-lite --version 2>&1 | head -1 || true
test -s "$BAM" || { echo "ERROR: BAM missing/empty (g0s0 rsync not finished?)"; exit 1; }

cellsnp-lite \
    -s "$BAM" \
    -O "$OUT" \
    -p 8 \
    --minMAF 0.05 \
    --minCOUNT 100 \
    --maxDEPTH 100000 \
    --cellTAG None \
    --UMItag UB \
    --gzip \
    --chrom "$CHROM"

echo "=== g0s1 done ==="
echo "Discovered variants: $(zcat "$OUT/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#')"
