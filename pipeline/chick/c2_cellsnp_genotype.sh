#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/c2_geno_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/c2_geno_%A_%a.err
#SBATCH -J chk_geno
#SBATCH --array=0-1
#
# c2 — cellSNP-lite mode-1 per-cell genotyping at the discovered sites, ONE library
# per array task (raw CellRanger barcodes match the raw BAM CB tags — no suffixing
# needed here; suffixing happens at the c2b merge so cell IDs are unique).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux
REGION="$WORK/cellsnp/discovery/cellSNP.base.vcf.gz"

LIBS=(ARC01 ARC02 ARC0S-chick)
L=${LIBS[$SLURM_ARRAY_TASK_ID]}
BAM="$WORK/bams/$L/gex_possorted_bam.bam"
BC="$WORK/barcodes/${L}_raw.tsv"
OUT="$WORK/cellsnp/$L"
mkdir -p "$OUT"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

echo "=== c2 cellSNP mode-1 genotyping: $L ==="
echo "Region VCF: $REGION ($(zcat "$REGION" | grep -vc '^#') sites)"
echo "Barcodes:   $BC ($(wc -l < "$BC") cells)"

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

echo "=== c2 done: $L -> $OUT ==="
