#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-2:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3h_bmp7_chick_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3h_bmp7_chick_%j.err
#SBATCH -J bmp7_g3h_chick

# g3h — CHICK lift + PITX-motif scan of the duck BMP7 locus.
#
# Question: Table S3 reports three left-open distal elements at the duck BMP7
# locus bearing PITX2-family motifs (candidate PITX2 enhancers). Do their CHICK
# orthologs carry the motifs, and does the chick BMP7 promoter?
#
# Reuses pipeline/g3g_hh21_chick_lift_scan.R unchanged apart from its new
# TSS_ANCHOR env hook (dist annotation only) — same Smith-Waterman Class I/II
# lift, same fixed-background JASPAR panel (PITX2 = MA1547.1, + PITX1/PITX3 via
# the family and the SMAD/FOX/GATA controls), so duck-vs-chick presence stays
# comparable with g3c/g3e/g3g at the PITX2 locus.
#
# Interpretation guard: the duck background rate in this window is 38/119 =
# 31.9% PITX2-motif+ per 501-bp peak, so single-window presence/absence means
# little on its own. The BED carries 16 stratified background windows precisely
# so the chick claim-element rate can be read against a chick background.
#
# Writes: figures/data/g3g_bmp7_chick_motifs.tsv, g3g_bmp7_lift.tsv

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 samtools/1.21 R/4.4.2
echo "=== g3h: chick lift + PITX-motif scan of the duck BMP7 locus ==="; date

DUCK_GTF="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"
CHICK_GTF="${DUCK_REFERENCE_ROOT}/crarc_galgal6a/genes/genes.gtf.gz"
CHICK_FA="${DUCK_REFERENCE_ROOT}/crarc_galgal6a/fasta/genome.fa"
for f in "$DUCK_GTF" "$CHICK_GTF" "$CHICK_FA"; do
  [ -f "$f" ] || { echo "ERROR: missing reference file: $f"; exit 1; }
done

BED_SRC="figures/data/g3h_bmp7_motif_windows.bed"
[ -f "$BED_SRC" ] || { echo "ERROR: WIN_BED missing on scratch: $BED_SRC"; exit 1; }

# --- duck BMP7 gene record: strand + TSS (for the promoter window + dist anchor)
# The duck GTF is Gnomon (gene_id/gene attributes, NO gene_name); the chick GTF
# carries gene_name. Match the quoted symbol in any of them, semicolon-anchored
# so BMP7 does not also match a BMP7-prefixed symbol.
read -r DUCK_CHR DUCK_START DUCK_END DUCK_STRAND < <(
  zcat "$DUCK_GTF" | awk -F'\t' '$3=="gene" && $9 ~ /"BMP7";/ {
    print $1, $4, $5, $7; exit }')
[ -n "${DUCK_STRAND:-}" ] || { echo "ERROR: BMP7 not found in duck GTF"; exit 1; }
if [ "$DUCK_STRAND" = "+" ]; then DUCK_TSS=$DUCK_START; else DUCK_TSS=$DUCK_END; fi
echo "duck BMP7: ${DUCK_CHR}:${DUCK_START}-${DUCK_END} (${DUCK_STRAND}) TSS=${DUCK_TSS}"

# --- chick BMP7 gene record: defines the region Smith-Waterman searches
read -r CHICK_CHR CHICK_START CHICK_END CHICK_STRAND < <(
  zcat "$CHICK_GTF" | awk -F'\t' '$3=="gene" && $9 ~ /"BMP7";/ {
    print $1, $4, $5, $7; exit }')
[ -n "${CHICK_STRAND:-}" ] || { echo "ERROR: BMP7 not found in chick GTF"; exit 1; }
if [ "$CHICK_STRAND" = "+" ]; then CHICK_TSS=$CHICK_START; else CHICK_TSS=$CHICK_END; fi
echo "chick BMP7: ${CHICK_CHR}:${CHICK_START}-${CHICK_END} (${CHICK_STRAND}) TSS=${CHICK_TSS}"

# chick search region: gene body +- 150 kb, wide enough for the duck +79 kb
# distal element even with avian indel drift. Clamp the start at 1.
FLANK=150000
REG_START=$(( CHICK_START - FLANK )); [ "$REG_START" -lt 1 ] && REG_START=1
REG_END=$(( CHICK_END + FLANK ))
CHICK_PROX_FA="chick_bmp7_region.fa"
samtools faidx "$CHICK_FA" "${CHICK_CHR}:${REG_START}-${REG_END}" > "$CHICK_PROX_FA"
echo "chick region: ${CHICK_CHR}:${REG_START}-${REG_END} -> $CHICK_PROX_FA"

# --- append the duck BMP7 promoter window (TSS +-250 bp = the same 501 bp the
# ArchR peaks use, so promoter and distal windows are scanned at equal width)
WIN_BED_RUN="g3h_bmp7_motif_windows_run.bed"
cp "$BED_SRC" "$WIN_BED_RUN"
printf '%s\t%d\t%d\t%s\t0\t%s\n' "$DUCK_CHR" $(( DUCK_TSS - 251 )) $(( DUCK_TSS + 250 )) \
  "BMP7_promoter|prom" "$DUCK_STRAND" >> "$WIN_BED_RUN"
echo "windows: $(grep -vc '^#' "$WIN_BED_RUN")"

export WIN_BED="$WIN_BED_RUN"
export OUT_TAG=bmp7
export KEEP_CLASS="Lopen_motif,bg,prom"
export CHICK_PROX_FA
export CHICK_PROX_START="$REG_START"
export TSS_ANCHOR="$DUCK_TSS"
export MOTIF_P_CUTOFF=1e-3      # same cutoff as g3f/g3g

Rscript pipeline/g3g_hh21_chick_lift_scan.R

echo "=== DONE g3h ==="; date
