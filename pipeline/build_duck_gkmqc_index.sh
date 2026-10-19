#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/references/gkmqc_idx_zju1_build_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/references/gkmqc_idx_zju1_build_%j.err
#SBATCH -J gkmqc_idx_zju1

# One-shot gkmQC null-sequence index build for duck ZJU1.0.
#
# Splits the cellranger-ARC genome FASTA into per-chromosome files,
# tars them in UCSC chromFa convention, then runs `gkmqc.py buildidx`.
# Output under $SCRATCH/references/gkmqc_idx_zju1/ is consumed by
# pipeline/f2b_gkmqc_rescue.py via the GKMQC_IDX_DIR env var.
#
# Prerequisites:
#   - conda env `gkmqc` with the upstream tool installed
#     (https://github.com/Dongwon-Lee/gkmQC — conda + make + make install).
#   - Reference genome FASTA at the crarc location below (already present
#     since the cellranger-arc v2.1.0 rerun).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
# Canonical cellranger-arc ZJU1.0 genome FASTA, built 2024-05-26 by
# /n/groups/tabin/chlee/singlecell/reference/20240314_01_build_arc_duck_pitx2.sh.
# That mkref run took `GCF_015476345.1_ZJU1.0_genomic.fna` (the NCBI
# RefSeq file) as its --fasta input with no sequence-side modification;
# the "_PITX2" suffix in the reference name refers to a GTF-only patch
# (3'UTR of PITX2 extended by 300 bp).  We confirmed byte-for-byte
# identity via md5sum on 2026-04-24:
#
#   ca4caae89d8c43ef54cddb6fa88a20cd  /n/groups/.../crarc_zju1_PITX2/fasta/genome.fa
#   ca4caae89d8c43ef54cddb6fa88a20cd  NCBI download of GCF_015476345.1
#
# Point at /n/groups/ directly — persistent storage, no duplication.
GENOME_FA="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa"
WORK_DIR="$SCRATCH/references/_gkmqc_build"
# gkmqc.py buildidx writes into $GKMQC_SRC/data/$GENOME_ID/ by
# convention — that's also where `gkmqc.py evaluate -g $GENOME_ID`
# looks at runtime.  Don't relocate the artifacts; just note the path
# so f2b_gkmqc_rescue.py consumers know where it lives.
GKMQC_SRC="$SCRATCH/references/_gkmqc_src"
GENOME_ID="zju1"
IDX_DIR="$GKMQC_SRC/data/$GENOME_ID"
N_THREADS=10

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

module purge
module load gcc/14.2.0 conda/miniforge3 samtools
conda activate gkmqc

echo "=== gkmQC index build: duck ZJU1.0 ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Genome FA: $GENOME_FA"
echo "Output dir: $OUT_DIR"

if ! command -v gkmqc.py &>/dev/null; then
    echo "ERROR: gkmqc.py not on PATH — is the conda env `gkmqc` set up?"
    exit 1
fi

if [ ! -f "$GENOME_FA" ]; then
    echo "ERROR: genome FA missing at $GENOME_FA"
    exit 1
fi
if [ ! -f "${GENOME_FA}.fai" ]; then
    echo "Indexing genome FASTA..."
    samtools faidx "$GENOME_FA"
fi

echo ""
echo "--- Splitting into per-chromosome FASTAs ---"
CHROM_DIR="$WORK_DIR/chromFa"
mkdir -p "$CHROM_DIR"
rm -f "$CHROM_DIR"/*.fa 2>/dev/null || true

while IFS=$'\t' read -r CHROM LEN _REST; do
    OUT="$CHROM_DIR/${CHROM}.fa"
    samtools faidx "$GENOME_FA" "$CHROM" > "$OUT"
    echo "  $CHROM ($LEN bp) → $(basename "$OUT")"
done < "${GENOME_FA}.fai"

echo ""
echo "--- Packaging chromFa.tar.gz ---"
TARBALL="$WORK_DIR/${GENOME_ID}.chromFa.tar.gz"
(cd "$CHROM_DIR" && tar -czf "$TARBALL" ./*.fa)
echo "  Tarball: $TARBALL ($(du -h "$TARBALL" | cut -f1))"

echo ""
echo "--- Running gkmqc.py buildidx ---"
gkmqc.py buildidx \
    -i "$TARBALL" \
    -g "$GENOME_ID" \
    -@ "$N_THREADS"

echo ""
echo "--- Index location ---"
echo "  $IDX_DIR"
du -sh "$IDX_DIR"/* 2>/dev/null | head -10

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
