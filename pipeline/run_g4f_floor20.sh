#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-3:00
#SBATCH -p priority
# -c 8 + N_CPUS=8 (below): parallelize pydeseq2 dispersion/LFC across peaks
# (joblib) — the genome-wide fit is the slow tail.  BLAS single-threaded to
# avoid oversubscription (joblib does the parallelism).
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/pitx2_locus_dyn_floor20_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/pitx2_locus_dyn_floor20_%j.err
#SBATCH -J pitx2_dyn_floor20
#
# g4f PITX2-locus dynamics at the FLOOR=20 pseudobulk set, to expand
# panel_D_lr_genetrack from 4 → 7 testable L/R stripes (adds HH13 ventral +
# dorsal HH17/16 + dorsal HH25).  Reads the floor-20 counts+metadata from the
# fresh da_archr_floorsens dir (g1 floor-sens regen, job 43946718); the
# regenerated floor-50 substrate was validated byte-identical to canonical
# (160/160 pseudobulks, n_cells corr 1.0), so this is a clean floor sensitivity
# variant — NOT a substrate change.  Writes to a separate OUT_DIR; the canonical
# floor-50 g4f outputs (figures/data/g4f_*) are untouched.
# mem=48G (canonical floor-50 used 32G for 160 pbs; floor-20 has 265 → headroom).
#
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH"
export SCRATCH_ROOT="$SCRATCH"
# 2026-08-11: repointed off the purged da_archr_floorsens — the unified
# relabel-campaign g1 writes all floors (+ per-floor metadata) into da_archr,
# and the campaign runner itself runs this floor20 leg; this standalone
# wrapper remains for a floor20-only rerun.
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export COUNTS_TSV="$DA_DIR/pseudobulk_counts_floor20.tsv.gz"
export META_TSV="$DA_DIR/pseudobulk_metadata_floor20.tsv"
export LOCUS_TSV="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
export OUT_DIR="$DA_DIR/pitx2_locus_dynamics_floor20"
mkdir -p "$OUT_DIR"

for f in "$COUNTS_TSV" "$META_TSV" "$LOCUS_TSV"; do
    [ -f "$f" ] || { echo "FATAL: required input missing: $f" >&2; exit 1; }
done

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export N_CPUS=8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== g4f PITX2-locus dynamics FLOOR=20 @ $(date) ==="
echo "COUNTS_TSV=$COUNTS_TSV"
echo "META_TSV=$META_TSV"
echo "OUT_DIR=$OUT_DIR"
python3 pipeline/g4f_pitx2_locus_dynamics.py
rc=$?
echo "=== g4f floor20 exit rc=$rc @ $(date) ==="
echo "OUT_DIR listing:"; ls -la "$OUT_DIR" 2>&1 | head -25
exit $rc
