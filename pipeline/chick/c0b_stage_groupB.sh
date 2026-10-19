#!/bin/bash
# c0b — stage chick Group-B (ARC09/10 GFP/WT syrinx) inputs (run on a cluster transfer node; NO modules).
# Group B uses the dedicated GEX libraries GEX09_GEX / GEX10_GEX (FC_08656), aligned to
# cr7.reference.gg6a.ncbi.EGFP (galGal6a + NCBI annot + EGFP transgene). Consistent with
# the duck process (genotype off the GEX library's possorted_genome_bam.bam).
# DESIGN: each embryo's L half -> one lane, R half -> the other (GEX09 vs GEX10 == side),
# 12 embryos (6 GFP + 6 WT, 4 per stage x HH13/15/25). EGFP is in the reference => GFP+
# embryos quantifiable directly. Read-only on standby; writes only to scratch.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_STANDBY_DATA_ROOT DUCK_REFERENCE_ROOT DUCK_SCRATCH_USER


STANDBY=${DUCK_STANDBY_DATA_ROOT}/FC_08656
REF=${DUCK_REFERENCE_ROOT}/cr7.reference.gg6a.ncbi.EGFP
WORK=${DUCK_SCRATCH_USER}/chick_demux_B

mkdir -p "$WORK"/{bams,barcodes,ref,cellsnp,vireo_seeded,logs,pipeline}

LIBS=(GEX09_GEX GEX10_GEX)
SUFS=(9 10)   # suffix encodes lane == anatomical side (cross-lane = same embryos, opposite halves)

for i in "${!LIBS[@]}"; do
  L=${LIBS[$i]}; s=${SUFS[$i]}
  src="$STANDBY/$L/outs"
  dst="$WORK/bams/$L"; mkdir -p "$dst"
  echo "=== staging $L (lane suffix -$s) ==="
  rsync -a "$src/possorted_genome_bam.bam" "$dst/"
  [ -f "$src/possorted_genome_bam.bam.bai" ] && rsync -a "$src/possorted_genome_bam.bam.bai" "$dst/" || echo "  WARN: no .bai (will need samtools index)"
  zcat "$src/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$WORK/barcodes/${L}_raw.tsv"
  echo "  $L: $(wc -l < "$WORK/barcodes/${L}_raw.tsv") CellRanger cells"
done

# Reference fasta (+.fai) for any bcftools fallback; cellSNP de-novo needs only BAMs.
rsync -a "$REF/fasta/genome.fa" "$WORK/ref/"
rsync -a "$REF/fasta/genome.fa.fai" "$WORK/ref/"

# Same galGal6a UCSC CHROMS as Group A (verified identical naming); exclude chrM, chrUn_*, EGFP.
echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr30,chr31,chr32,chr33,chrZ,chrW" > "$WORK/ref/chroms.txt"

echo "=== c0b done -> $WORK ==="
ls -la "$WORK"/bams/*/*.bam
wc -l "$WORK"/barcodes/*_raw.tsv
