#!/bin/bash
# rederive_cellranger_from_geo.sh — rebuild cellranger{,-arc} `outs/` from the
# GEO-deposited possorted BAMs (GSE336327/GSE336329). This is the public
# re-derivation path for the QC-panel inputs the record does NOT deposit
# (summary.csv / per_barcode_metrics.csv / metrics_summary.csv; ordering:
# docs/geo_tier2_build_order.md).
#
#   GEX:  rederive_cellranger_from_geo.sh gex <LIB> <GEX_BAM> <TRANSCRIPTOME_REF> <WORKDIR>
#   ARC:  rederive_cellranger_from_geo.sh arc <LIB> <GEX_BAM> <ATAC_BAM> <ARC_REF> <WORKDIR>
#
# Versions of record (Methods; pipeline/a1_cellranger.sh; pipeline/rerun_arc_v210.sh):
#   cellranger 7.1.0 (GEX; chick Group-B against the pipeline/chick/refbuild
#   reference), cellranger-arc 2.1.0 (multiome). bamtofastq is 10x's tool, bundled
#   under cellranger*/lib/bin — override discovery with BAMTOFASTQ=/path.
# 10x possorted BAMs retain ALL reads (including unmapped, e.g. EGFP), so the
# round-trip reproduces the original FASTQs; identical tool version + reference
# then reproduces `outs/` including the metrics files.
#
# Run inside an sbatch you size; measured on the original run:
#   bamtofastq leg alone: -c 4 --mem 8G;  cellranger{,-arc} count: ~ -c 8
#   --mem 64G -t 1-00:00 per library (arc libraries take ~a day).
set -euo pipefail

MODE=${1:?usage: rederive_cellranger_from_geo.sh gex|arc <LIB> <BAM(s)> <REF> <WORKDIR>}
LIB=${2:?library name (e.g. ARC05, GEX05, GEX09)}
BTF=${BAMTOFASTQ:-$(command -v bamtofastq || true)}
[ -n "$BTF" ] || { echo "ERROR: bamtofastq not on PATH; set BAMTOFASTQ=/path/to/bamtofastq" >&2; exit 1; }

b2f() {  # <bam> <dest-dir> — refuses to overwrite (partial reruns must clean up first)
    local bam=$1 dest=$2
    [ -f "$bam" ] || { echo "ERROR: BAM not found: $bam" >&2; exit 1; }
    [ -e "$dest" ] && { echo "ERROR: $dest exists; remove it to re-run" >&2; exit 1; }
    "$BTF" --nthreads="${B2F_THREADS:-4}" "$bam" "$dest"
}

fq_dirs() {  # bamtofastq nests one subdir per source flowcell/run — return all, comma-joined
    local root=$1 dirs
    dirs=$(find "$root" -mindepth 1 -maxdepth 1 -type d | sort | paste -sd, -)
    [ -n "$dirs" ] || { echo "ERROR: bamtofastq produced no FASTQ dirs under $root" >&2; exit 1; }
    echo "$dirs"
}

case "$MODE" in
  gex)
    BAM=${3:?gex possorted BAM}; REF=${4:?transcriptome reference dir}; WORK=${5:?workdir}
    mkdir -p "$WORK"; cd "$WORK"
    b2f "$BAM" "$WORK/fastq_${LIB}"
    cellranger count --id="${LIB}_GEX" --transcriptome="$REF" \
        --fastqs="$(fq_dirs "$WORK/fastq_${LIB}")" --sample=bamtofastq \
        --localcores="${CR_CORES:-8}" --localmem="${CR_MEM_GB:-64}"
    echo "DONE gex ${LIB}: $WORK/${LIB}_GEX/outs (metrics_summary.csv inside)"
    ;;
  arc)
    GEXBAM=${3:?gex possorted BAM}; ATACBAM=${4:?atac possorted BAM}
    REF=${5:?cellranger-arc reference dir}; WORK=${6:?workdir}
    mkdir -p "$WORK"; cd "$WORK"
    b2f "$GEXBAM"  "$WORK/fastq_${LIB}_gex"
    b2f "$ATACBAM" "$WORK/fastq_${LIB}_atac"
    CSV="$WORK/${LIB}_libraries.csv"
    {
        echo "fastqs,sample,library_type"
        find "$WORK/fastq_${LIB}_gex"  -mindepth 1 -maxdepth 1 -type d | sort | \
            while read -r d; do echo "$d,bamtofastq,Gene Expression"; done
        find "$WORK/fastq_${LIB}_atac" -mindepth 1 -maxdepth 1 -type d | sort | \
            while read -r d; do echo "$d,bamtofastq,Chromatin Accessibility"; done
    } > "$CSV"
    [ "$(wc -l < "$CSV")" -ge 3 ] || { echo "ERROR: libraries CSV came out empty: $CSV" >&2; exit 1; }
    cellranger-arc count --id="${LIB}" --reference="$REF" --libraries="$CSV" \
        --localcores="${CR_CORES:-8}" --localmem="${CR_MEM_GB:-64}"
    echo "DONE arc ${LIB}: $WORK/${LIB}/outs (summary.csv + per_barcode_metrics.csv inside)"
    ;;
  *) echo "ERROR: mode must be gex or arc (got '$MODE')" >&2; exit 1 ;;
esac
