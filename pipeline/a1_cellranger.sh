#!/bin/bash
# Phase A1: cellranger count (12 libraries)
#
# This step was run outside this repo on O2. Documenting for reproducibility.
# Reference genome: Mallard duck zju1 (GCF_015476345.1)
# Software: cellranger 7.1.0 (GEX-only); the multiome counts of record are cellranger-arc 2.1.0
# (pipeline/rerun_arc_v210.sh)
#
# Each library was processed independently:
#   cellranger count --id=<LIB>_GEX \
#     --transcriptome=/path/to/zju1_cellranger_ref \
#     --fastqs=/path/to/fastqs/<LIB> \
#     --localcores=8 --localmem=64
#
# Libraries (12 total, 4 genotype groups):
#   Group 2 (G2): GEX03-1A, GEX03-1B, GEX03-2   (KC + Mallard mix)
#   Group 3 (G3): GEX04, GEX0S-duck               (unknown pool)
#   Group 5 (G5): GEX05, GEX06, GEX07, GEX08      (all Mallard, 36 donors)
#   Group 5 (G5): PILOTGEX3, PILOTGEX4A, PILOTGEX4B (PILOT libraries)
#
# Outputs on scratch:
#   /n/scratch/users/c/cl266/FC_08656/<LIB>_GEX/outs/
#     filtered_feature_bc_matrix/  (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#     raw_feature_bc_matrix/
#     possorted_genome_bam.bam
#
# CellBender (ambient RNA correction) was run on all 12 libraries:
#   cellbender remove-background --input raw_feature_bc_matrix.h5 \
#     --output <LIB>_cellbender.h5 --epochs 150 --fpr 0.01
#   Outputs: /n/scratch/users/c/cl266/FC_08656/phase2/cellbender/<LIB>_cellbender_filtered.h5

echo "This is a documentation stub. cellranger was run outside this pipeline."
echo "See comments in this script for the exact commands and paths."
exit 0
