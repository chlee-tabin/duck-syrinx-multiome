#!/bin/bash
# c0 — stage chick Group-A inputs (run on a cluster transfer node; NO modules needed).
# Read-only on standby; writes only to the fresh scratch workspace.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_STANDBY_DATA_ROOT DUCK_REFERENCE_ROOT DUCK_SCRATCH_USER


STANDBY=${DUCK_STANDBY_DATA_ROOT}/FC_08565
REF=${DUCK_REFERENCE_ROOT}/crarc_galgal6a
WORK=${DUCK_SCRATCH_USER}/chick_demux

mkdir -p "$WORK"/{bams,barcodes,ref,cellsnp,vireo_seeded,logs}

# Libraries (Group A) and barcode suffixes used for the merged cellData
LIBS=(ARC01 ARC02 ARC0S-chick)
SUFS=(1 2 S)

for i in "${!LIBS[@]}"; do
  L=${LIBS[$i]}; s=${SUFS[$i]}
  src="$STANDBY/$L/outs"
  dst="$WORK/bams/$L"; mkdir -p "$dst"
  echo "=== staging $L (suffix -$s) ==="
  rsync -a "$src/gex_possorted_bam.bam" "$src/gex_possorted_bam.bam.bai" "$dst/"
  # CellRanger filtered barcodes (the cell-calling we key off, per project rule)
  zcat "$src/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$WORK/barcodes/${L}_raw.tsv"
  echo "  $L: $(wc -l < "$WORK/barcodes/${L}_raw.tsv") CellRanger cells"
done

# Reference: genome.fa (+.fai) for any bcftools fallback; gtf optional for gene-restricted modes
rsync -a "$REF/fasta/genome.fa" "$REF/fasta/genome.fa.fai" "$WORK/ref/"
rsync -a "$REF/genes/genes.gtf.gz" "$WORK/ref/" 2>/dev/null || echo "  (genes.gtf.gz not staged — only needed for gene-restricted mode)"

# Chick galGal6a CHROMS (UCSC naming): chr1-28, chr30-33, chrZ, chrW; EXCLUDE chrM (mito).
# chrZ/chrW retained — also enable ZW sex determination downstream.
echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr30,chr31,chr32,chr33,chrZ,chrW" > "$WORK/ref/chroms.txt"

echo "=== c0 done -> $WORK ==="
echo "BAMs:"; ls -la "$WORK"/bams/*/*.bam
echo "Barcodes:"; wc -l "$WORK"/barcodes/*_raw.tsv
