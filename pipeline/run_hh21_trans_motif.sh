#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:40
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_trans_motif_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_trans_motif_%j.err
#SBATCH -J hh21_trans_motif
#
# Trans-motif follow-up: does PITX2-motif content distinguish the
# chromatin-FLAT promoters of HH21 L/R DEGs from bilateral promoters?
#   prep (python) cdw GTF -> transcript TSS -> containing ArchR peak -> response/bilateral
#   scan (R)      getSeq + the g3f 12-motif JASPAR panel + raw TAATCC, per promoter peak
# The three-way signed-z split + GLM run LOCALLY afterwards (hh21_trans_motif.py test).

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
export DATA_DIR="$SCRATCH/figures/data"
TX="${TMPDIR:-/tmp}/hh21_trans_tx_$$.gtf"
echo "=== hh21_trans_motif ==="; date

[ -f "$GTF" ] || { echo "MARKER_RC=10 (GTF missing: $GTF)"; exit 10; }
[ -f "$DATA_DIR/h4_peak_empirical_p_HH21.tsv.gz" ] || { echo "MARKER_RC=11 (h4 HH21 missing)"; exit 11; }

echo "--- extracting transcript lines from the cdw GTF ---"
zcat "$GTF" | awk -F'\t' '$3=="transcript"' > "$TX" || { echo "MARKER_RC=12"; exit 12; }
echo "GTF=$GTF"
echo "transcript lines: $(wc -l < "$TX")"
# awk exits 0 on zero matches, so an assembly/GTF mismatch would sail past pipefail
[ -s "$TX" ] || { echo "MARKER_RC=14 (no transcript lines from $GTF)"; exit 14; }

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2 || { echo "MARKER_RC=13 (conda activate failed)"; exit 13; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "--- prep ---"
python figures/figS03/hh21_trans_motif.py prep "$TX"
rc=$?; rm -f "$TX"
[ $rc -eq 0 ] || { echo "MARKER_RC=$rc (prep failed)"; exit $rc; }

module purge
module load gcc/14.2.0 R/4.4.2
echo "--- scan ---"
Rscript figures/figS03/hh21_trans_motif_scan.R
rc=$?
[ $rc -eq 0 ] || { echo "MARKER_RC=$rc (scan failed)"; exit $rc; }

echo "MARKER_RC=0"
echo "=== DONE ==="; date
