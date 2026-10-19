#!/bin/bash
# Phase A2: cellSNP variant calling (per genotype group)
#
# This step was run outside this repo on O2. Documenting for reproducibility.
# Software: cellsnp-lite 1.2.2
#
# Variant discovery strategy depends on pool composition:
#
# G6 (Groups 5/6, all Mallard):
#   Universal variant set (665K positions) — works because pool is genetically homogeneous.
#   cellsnp-lite -s <BAM> -b <BARCODES> -O <OUT> \
#     -R /path/to/universal_variants_665K.vcf.gz \
#     -p 8 --minMAF 0.01 --minCOUNT 20 --UMItag Auto
#
# G2 (KC + Mallard mix):
#   Library-specific variant discovery — required because KC/Mallard divergence
#   means universal Mallard variants have no discriminative power for KC donors.
#   cellsnp-lite -s <BAM> -b <BARCODES> -O <OUT> \
#     -R /path/to/universal_variants_665K.vcf.gz \
#     -p 8 --minMAF 0.01 --minCOUNT 20 --UMItag Auto
#
# G3 (unknown pool):
#   De novo variant discovery (mode 2b) — the 665K universal set has ZERO
#   discriminative power for this pool and caused vireo to shatter donors.
#   cellsnp-lite -s <BAM> -b <BARCODES> -O <OUT> \
#     --geneAnnoFile /path/to/genes.gtf \
#     -p 8 --minMAF 0.1 --minCOUNT 100 --UMItag Auto
#   Result: 173K pool-specific variants.
#
# G9 (all Mallard, PILOT4B):
#   Same universal variant set as G6 (dp500 filtering applied post-vireo).
#
# Outputs on scratch:
#   /n/scratch/users/c/cl266/FC_08656/<GROUP>_cellsnp/
#     cellSNP.base.vcf.gz, cellSNP.samples.tsv, cellSNP.tag.AD.mtx, etc.

echo "This is a documentation stub. cellSNP was run outside this pipeline."
echo "See comments in this script for the exact commands and variant strategies."
exit 0
