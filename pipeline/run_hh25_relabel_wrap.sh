#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-2:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_wrap_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_wrap_%j.err
#SBATCH -J hh25_relabel_wrap
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)
#
# HH25 relabel follow-up wrap-up: (1) h6 paired HH25 with the INFERENCE_POOL fix, and
# (2) g4f on the floor20 pseudobulk (the file fig03 panels B/C actually read).
# Both off the fresh da_archr pseudobulks g1 already wrote (SSOT relabel).
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export LOCUS_TSV="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"

# ── h6 paired DESeq2 HH25 — pool comes from POOLS[STAGE] inside h6 itself
#    (2026-08-11: h6 is SSOT-derived + set-equality-validated; a hardcoded
#    override here would rot silently the next time the pool changes) ──
echo "################ h6 paired DESeq2 HH25 ################"; date
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
STAGE="HH25" \
    OUT_TSV="$DA_DIR/h6_paired_deseq2_HH25_balanced.tsv" \
    OUT_SUMMARY="$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" \
    python3 pipeline/h6_paired_deseq2.py || { echo "FATAL: h6 failed (its outputs are synced to figures/data)" >&2; exit 4; }

# ── g4f on the FLOOR20 pseudobulk = the fig03 panel B/C input (g4f_lr_long_floor20) ──
echo "################ g4f PITX2-locus dynamics FLOOR20 (fig03 input) ################"; date
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 N_CPUS=8
export COUNTS_TSV="$DA_DIR/pseudobulk_counts_floor20.tsv.gz"
export META_TSV="$DA_DIR/pseudobulk_metadata_floor20.tsv"
export OUT_DIR="$DA_DIR/pitx2_locus_dynamics_floor20"
mkdir -p "$OUT_DIR"
for f in "$COUNTS_TSV" "$META_TSV" "$LOCUS_TSV"; do
    [ -f "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
python3 pipeline/g4f_pitx2_locus_dynamics.py || { echo "FATAL: g4f floor20 failed" >&2; exit 3; }

echo "################ WRAP DONE ################"; date
echo "=== h6 summary ==="; cat "$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" 2>/dev/null || echo "(no h6 summary)"
echo "=== g4f floor20 c-promoter HH25 ventral (sanity) ==="
grep "3971182" "$OUT_DIR/g4f_lr_long.tsv" 2>/dev/null | grep "HH25" || echo "(c-promoter row not found)"
echo "COMPLETE_OK"
