#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-3:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux_B/logs/c2b_geno_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux_B/logs/c2b_geno_%A_%a.err
#SBATCH -J chkB_geno
#SBATCH --array=0-1
#
# c2b — Group B cellSNP mode-1 genotyping at discovered sites, one GEX library per task.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

WORK=${DUCK_SCRATCH_USER}/chick_demux_B
REGION="$WORK/cellsnp/discovery/cellSNP.base.vcf.gz"
LIBS=(GEX09_GEX GEX10_GEX)
L=${LIBS[$SLURM_ARRAY_TASK_ID]}
BAM="$WORK/bams/$L/possorted_genome_bam.bam"
BC="$WORK/barcodes/${L}_raw.tsv"
OUT="$WORK/cellsnp/$L"; mkdir -p "$OUT"

module purge; module load gcc/14.2.0 conda/miniforge3; conda activate vireo
echo "=== c2b genotyping $L ($(wc -l < "$BC") cells) at $(zcat "$REGION"|grep -vc '^#') sites ==="
cellsnp-lite -s "$BAM" -b "$BC" -R "$REGION" -O "$OUT" -p 8 \
    --minMAF 0 --minCOUNT 0 --cellTAG CB --UMItag UB --gzip
echo "=== c2b done: $L ==="
