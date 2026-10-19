#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 0-6:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_chain_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_chain_%j.err
#SBATCH -J hh25_relabel_chain
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)
#
# HH25 relabel follow-up — HH25 fig03 ATAC confirmation after the d2g relabel
# (cl8 Mesothelium->Ventral spl., cl3->Syringeal, cl5->Ventral spl.,
#  cl10->Venous-pole spl.).  g1 now RE-DERIVES d2g_subcluster from the SSOT
# (patched), so the relabel flows into the pseudobulks WITHOUT an f1 rebuild.
#
# One consistent output dir (da_archr, currently empty -> no clobber of the
# Jun-17 da/).  g4f primary (floor50) = fig03 PITX2c-promoter headline;
# h6/h7 default to POOLS["HH25"]["inference"] from the synced subcluster_map.
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
# fig03 uses the ArchR 283K-peak matrix (DA_DIR=da_archr), NOT the SnapATAC2
# 25K matrix.  f2d copies d2g_cluster (leiden) into this obs, which the g1
# SSOT-relabel patch remaps from.
export F2_H5AD="$SCRATCH/phase2/atac/step_f2_peak_matrix_archr.h5ad"
export G1_LR_TSV="$SCRATCH/figures/data/lr_inference_multiome.tsv"
export LOCUS_TSV="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
mkdir -p "$DA_DIR"

for f in "$F2_H5AD" "$G1_LR_TSV" "$LOCUS_TSV"; do
    [ -f "$f" ] || { echo "FATAL: required input missing: $f" >&2; exit 1; }
done
grep -q "SSOT relabel via sidecar" pipeline/g1_pseudobulk_atac.py \
    || { echo "FATAL: g1 on scratch is NOT the SSOT-relabel patch" >&2; exit 1; }
[ -f "$SCRATCH/phase2/atac/step_f1_coupled_atac_metadata.tsv" ] \
    || { echo "FATAL: f1 metadata sidecar (leiden) missing on scratch" >&2; exit 1; }

echo "################ g1 pseudobulk (SSOT relabel) ################"; date
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
python3 pipeline/g1_pseudobulk_atac.py || { echo "FATAL: g1 failed" >&2; exit 2; }
[ -f "$DA_DIR/pseudobulk_counts.tsv.gz" ] \
    || { echo "FATAL: g1 wrote no pseudobulk_counts.tsv.gz" >&2; exit 2; }
echo "--- HH25 per-subcluster pseudobulk summary (expect cl8 cells now Ventral spl.) ---"
grep -E "HH25" "$DA_DIR/pseudobulk_summary.tsv" 2>/dev/null || true

echo "################ g4f PITX2-locus dynamics (fig03 headline, floor50) ################"; date
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 N_CPUS=8
export OUT_DIR="$DA_DIR/pitx2_locus_dynamics"
mkdir -p "$OUT_DIR"
python3 pipeline/g4f_pitx2_locus_dynamics.py || { echo "FATAL: g4f failed" >&2; exit 3; }

echo "################ h6 paired DESeq2 HH25 (INFERENCE_POOL fix) ################"; date
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
# 2026-08-11: h6 now derives the pool from POOLS[STAGE] itself (SSOT,
# set-equality validated) — no override needed, and a stale one would exit.
STAGE="HH25" \
    OUT_TSV="$DA_DIR/h6_paired_deseq2_HH25_balanced.tsv" \
    OUT_SUMMARY="$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" \
    python3 pipeline/h6_paired_deseq2.py || { echo "FATAL: h6 failed (its outputs are synced to figures/data)" >&2; exit 4; }

echo "################ h7 unpaired DESeq2 HH25 (sensitivity, default POOLS) ################"; date
STAGE="HH25" python3 pipeline/h7_unpaired_deseq2.py \
    || { if [ "${ALLOW_H7_FAILURE:-0}" = "1" ]; then echo "WARN: h7 failed (ALLOW_H7_FAILURE=1: sensitivity arm skipped)" >&2; else echo "FATAL: h7 failed (set ALLOW_H7_FAILURE=1 to downgrade this sensitivity arm)" >&2; exit 5; fi; }

echo "################ CHAIN DONE ################"; date
echo "--- key outputs ---"
ls -la "$DA_DIR/pseudobulk_counts.tsv.gz" \
       "$OUT_DIR/g4f_lr_grid.tsv" \
       "$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" 2>&1 | head
