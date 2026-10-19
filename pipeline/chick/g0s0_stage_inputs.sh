#!/bin/bash
# g0s0 — stage GEX0S-chick (Group A symmetry-breaking LEFT lane) inputs (run on a cluster transfer node; NO modules).
#
# GEX0S-chick_GEX (FC_08565) is the GEX/RNA side of the ARC0S symmetry-breaking lane: the
# design-known LEFT halves of the 8 HH15-16/HH25 chick embryos (Group A). It is aligned to
# cr7.reference.gg6a.ncbi.EGFP (galGal6a + NCBI annot + EGFP transgene) — the SAME reference
# as Group B GEX09/GEX10, so it is batch-matched to that demux. (Group A is WT-only; the EGFP
# feature is in the reference but no GFP line was used, so no GFP+ donors are expected.)
#
# Demux off the GEX library's existing possorted_genome_bam.bam + filtered_feature_bc_matrix
# (the established duck/chick process — NO cellranger rerun needed for the GEX modality).
# Single library => one lane => the LEFT halves of <=8 embryos; expect N≈8 in vireo.
# Read-only on standby; writes only to scratch.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_STANDBY_DATA_ROOT DUCK_SCRATCH_USER


STANDBY=${DUCK_STANDBY_DATA_ROOT}/FC_08565
WORK=${DUCK_SCRATCH_USER}/chick_demux/gex0s_demux

mkdir -p "$WORK"/{bams,barcodes,matrix,ref,cellsnp,vireo_seeded,logs}

L=GEX0S-chick_GEX
src="$STANDBY/$L/outs"
dst="$WORK/bams/$L"; mkdir -p "$dst" "$WORK/matrix/$L"

echo "=== staging $L (single lane; suffix -0S) ==="
rsync -a "$src/possorted_genome_bam.bam" "$dst/"
if [ -f "$src/possorted_genome_bam.bam.bai" ]; then
  rsync -a "$src/possorted_genome_bam.bam.bai" "$dst/"
else
  echo "  WARN: no .bai on standby (will need samtools index before cellSNP)"
fi
# CellRanger filtered barcodes (the cell-calling we key off, per project rule)
zcat "$src/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$WORK/barcodes/${L}_raw.tsv"
echo "  $L: $(wc -l < "$WORK/barcodes/${L}_raw.tsv") CellRanger cells"
# filtered matrix (for the sex W/Z gate, mirrors c5c/cA_sex_scan)
rsync -a "$src/filtered_feature_bc_matrix/"{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz} "$WORK/matrix/$L/"

# Reuse the W/Z gene panels already staged for Group A (compute can't see standby).
for g in chick_W_genes.tsv chick_Z_genes.tsv; do
  if [ -f "${DUCK_SCRATCH_USER}/chick_demux/$g" ]; then
    cp -n "${DUCK_SCRATCH_USER}/chick_demux/$g" "$WORK/$g"
    echo "  reused $g ($(wc -l < "$WORK/$g") genes)"
  else
    echo "  WARN: $g not found in chick_demux/ — stage it before the sex gate (g0s4)"
  fi
done

# galGal6a UCSC CHROMS (same as Group A/B): chr1-28, chr30-33, chrZ, chrW; EXCLUDE chrM + EGFP.
# chrZ/chrW retained for ZW sex determination downstream.
echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr30,chr31,chr32,chr33,chrZ,chrW" > "$WORK/ref/chroms.txt"

echo "=== g0s0 done -> $WORK ==="
ls -la "$WORK"/bams/*/*.bam
wc -l "$WORK"/barcodes/*_raw.tsv
