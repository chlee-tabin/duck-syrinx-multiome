#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_trans_motif_seq_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_trans_motif_seq_%j.err
#SBATCH -J hh21_trans_motif_seq
#
# TSS-anchored SEQUENCE promoters (-1000/+500, strand-aware, merged per gene) for EVERY
# expressed HH21 Ventral-spl. gene — independent of ATAC, unlike run_hh21_trans_motif.sh.
#   prep (python) cdw GTF -> transcript TSS -> clamped windows -> merge per gene -> BED
#   scan (R)      the SAME shared scanner, re-pointed via SCAN_BED / SCAN_OUT / SCAN_ID_COL
# The GLM + DA correlation run LOCALLY afterwards (hh21_trans_motif_seq.py test).

set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH" || { echo "MARKER_RC=13 (cd $SCRATCH failed)"; exit 13; }
export SCRATCH_ROOT="$SCRATCH"
export CDW="${CDW:-$DUCK_REFERENCE_ROOT}"
export GTF="${GTF:-$CDW/crarc_zju1_PITX2/genes/genes.gtf.gz}"
export FAI="${FAI:-$CDW/crarc_zju1_PITX2/fasta/genome.fa.fai}"
export DATA_DIR="$SCRATCH/figures/data"
TX="${TMPDIR:-/tmp}/hh21_seq_tx_$$.gtf"
echo "=== hh21_trans_motif_seq ==="; date

for f in "$GTF" "$FAI" "$DATA_DIR/h4_peak_empirical_p_HH21.tsv.gz"; do
  [ -f "$f" ] || { echo "MARKER_RC=10 (missing: $f)"; exit 10; }
done

echo "--- extracting transcript lines from the cdw GTF ---"
zcat "$GTF" | awk -F'\t' '$3=="transcript"' > "$TX" || { echo "MARKER_RC=12"; exit 12; }
echo "GTF=$GTF FAI=$FAI"
echo "transcript lines: $(wc -l < "$TX")"
# awk exits 0 on zero matches, so an assembly/GTF mismatch would sail past pipefail
[ -s "$TX" ] || { echo "MARKER_RC=14 (no transcript lines from $GTF)"; exit 14; }

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2 || { echo "MARKER_RC=13 (conda activate failed)"; exit 13; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "--- prep ---"
python figures/figS03/hh21_trans_motif_seq.py prep "$TX" "$FAI"
rc=$?; rm -f "$TX"
[ $rc -eq 0 ] || { echo "MARKER_RC=$rc (prep failed)"; exit $rc; }

module purge
module load gcc/14.2.0 R/4.4.2
# two promoter variants: merged-per-gene (length-offset) and single-most-5'-TSS (equal length)
for V in "hh21_seq:merged" "hh21_seq1:single"; do
  TAG="${V%%:*}"; NAME="${V##*:}"
  echo "--- scan ($NAME) ---"
  SCAN_BED="$DATA_DIR/${TAG}_promoters.bed" \
  SCAN_OUT="$DATA_DIR/${TAG}_promoter_motifs.tsv" \
  SCAN_ID_COL="gene_id" \
  Rscript figures/figS03/hh21_trans_motif_scan.R
  rc=$?
  [ $rc -eq 0 ] || { echo "MARKER_RC=$rc (scan $NAME failed)"; exit $rc; }
done

echo "MARKER_RC=0"
echo "=== DONE ==="; date
