#!/bin/bash
#SBATCH -p medium
#SBATCH -t 1-12:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux_B/logs/c1c_disc_fine_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux_B/logs/c1c_disc_fine_%j.err
#SBATCH -J chkB_discF
#
# c1c — Group B FINER cellSNP mode-2b discovery for sibling / under-captured-embryo recovery.
# The default discovery (c1b: --minMAF 0.05 --minCOUNT 100 -> 274K variants) merges close
# M/F sibling pairs and loses 2 under-captured embryos into the doublet pile (resolution
# bottleneck, not cell count). Lowering minMAF captures rarer alleles private to 1 embryo
# (a het allele private to 1 of 12 embryos sits at ~4% pool MAF; under-captured embryos
# even lower), and lowering minCOUNT keeps variants covered by fewer reads. Writes to a
# SEPARATE path (discovery_fine) — the canonical discovery/ is untouched.
# Downstream after this: c2b genotype on the fine VCF -> c2c merge -> c3b vireo sweep ->
# c5g sex scan, all pointed at the fine variant set.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux_B
OUT="$WORK/cellsnp/discovery_fine"
CHROM=$(cat "$WORK/ref/chroms.txt")
BAMS="$WORK/bams/GEX09_GEX/possorted_genome_bam.bam,$WORK/bams/GEX10_GEX/possorted_genome_bam.bam"
rm -rf "$OUT"; mkdir -p "$OUT"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== c1c FINER cellSNP mode-2b discovery (Group B): minMAF 0.02 minCOUNT 50 ==="
echo "BAMs: $BAMS"
cellsnp-lite --version 2>&1 | head -1 || true

cellsnp-lite \
    -s "$BAMS" \
    -O "$OUT" \
    -p 8 \
    --minMAF 0.02 \
    --minCOUNT 50 \
    --maxDEPTH 100000 \
    --cellTAG None \
    --UMItag UB \
    --gzip \
    --chrom "$CHROM"

echo "=== c1c done ==="
echo "Fine variants: $(zcat "$OUT/cellSNP.base.vcf.gz" 2>/dev/null | grep -vc '^#')   (vs 274K at minMAF 0.05)"
