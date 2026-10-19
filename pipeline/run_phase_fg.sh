#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -t 0-8:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/phase_fg_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/phase_fg_%j.err
#SBATCH -J duck_atac

# Phase F (ATAC integration) + Phase G (differential accessibility) pipeline.
# Prerequisites: Phase 1 complete (step04_preCB_clustered.h5ad exists).
#
# Usage:
#   sbatch pipeline/run_phase_fg.sh                  # full F+G
#   sbatch pipeline/run_phase_fg.sh --from f1         # f1 onward
#   sbatch pipeline/run_phase_fg.sh --from g1 --to g2l # g1+g2l only
#   sbatch pipeline/run_phase_fg.sh --from f3 --to f3  # single step
#
# Resource guide:
#   Step      Memory   Time    Notes
#   f1        64G      2:00    Fragment loading + barcode matching
#   f2        48G      6:00    MACS3 peak calling per cell type
#   f2b       32G      8:00    gkmQC rescue for sub-500 pools (needs gkmqc env)
#   f3        32G      1:00    Spectral embedding
#   g1        32G      0:30    Pseudobulk aggregation (d2g subclusters x side)
#   g2l       16G      1:00    HH17/16 LPM paired PyDESeq2 (poscounts)
#   full      64G      12:00   (defaults)
#
# Conda env setup (run ONCE before first use):
#   module load gcc/14.2.0 conda/miniforge3
#   conda create -n snapatac2 python=3.11 -y
#   conda activate snapatac2
#   pip install snapatac2 pydeseq2 scikit-learn statsmodels scanpy

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
PHASE2="$SCRATCH/phase2"
ATAC_DIR="$PHASE2/atac"
# Pin the d2n L/R sidecar (g1 honors G1_LR_TSV; SLURM --export=ALL would otherwise let
# a stale G1_LR_TSV from the submit shell silently force the legacy d2l axis).
# This is g1's own default, made explicit. For an intentional d2l build, run g1 directly
# with G1_LR_TSV=.../d2l_lr_inference.tsv (not via this canonical runner).
export G1_LR_TSV="$SCRATCH/figures/data/lr_inference_multiome.tsv"
mkdir -p "$ATAC_DIR/peaks" "$ATAC_DIR/da" "$ATAC_DIR/motifs"
cd "$SCRATCH"

# ── Phase definitions ──────────────────────────────────────────────────
VALID_PHASES="f1 f2 f2b f3 g1 g2l"

validate_phase() {
    local phase="$1" label="$2"
    local found=false
    for p in $VALID_PHASES; do [ "$p" = "$phase" ] && found=true; done
    if ! $found; then
        echo "Error: Invalid $label phase '$phase'. Valid: $VALID_PHASES" >&2
        exit 1
    fi
}

# ── Parse arguments ────────────────────────────────────────────────────
START="f1"
END="g2l"

while [ $# -gt 0 ]; do
    case "$1" in
        --from)
            [ -z "${2:-}" ] && echo "Error: --from requires a phase" >&2 && exit 1
            START="$2"; shift 2 ;;
        --from=*)
            START="${1#--from=}"; shift ;;
        --to)
            [ -z "${2:-}" ] && echo "Error: --to requires a phase" >&2 && exit 1
            END="$2"; shift 2 ;;
        --to=*)
            END="${1#--to=}"; shift ;;
        *)
            START="$1"; shift ;;
    esac
done

validate_phase "$START" "start"
validate_phase "$END" "end"

should_run() {
    local phase="$1"
    local in_range=false
    for p in $VALID_PHASES; do
        [ "$p" = "$START" ] && in_range=true
        if [ "$p" = "$phase" ] && $in_range; then return 0; fi
        [ "$p" = "$END" ] && in_range=false
    done
    return 1
}

# ── Environment setup ──────────────────────────────────────────────────
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# Verify tools
python3 -c "import snapatac2; print(f'SnapATAC2 {snapatac2.__version__}')" || \
    { echo "ERROR: snapatac2 not installed. See conda setup instructions above."; exit 1; }
python3 -c "import pydeseq2; print(f'PyDESeq2 {pydeseq2.__version__}')" || \
    { echo "ERROR: pydeseq2 not installed (needed for g2l)"; exit 1; }

echo ""
echo "=== Phase F+G Pipeline (from: $START, to: $END) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

# ── Phase F: ATAC Integration ─────────────────────────────────────────
if should_run f1; then
    echo ""
    echo "--- Phase F1: ATAC QC + Coupled Barcode Matching ---"
    # f1 reads the COMMITTED d2g sidecars from figures/data (synced from the
    # repo), NOT d2g's live phase2/figures output — so after a d2g rerun the
    # labels here are stale until the derived-data round-trip (sync_derived →
    # git → sync_to_scratch) completes.  f1 hard-errors on a missing or
    # id-mismatched sidecar, but a same-partition re-annotation is invisible;
    # warn the operator the same way run_pipeline.sh does for h1.
    echo "NOTE: f1 uses COMMITTED d2g sidecars (figures/data). If d2g was"
    echo "  re-run more recently than the last derived-data sync, STOP and"
    echo "  complete the round-trip first — f1 would bake stale labels into"
    echo "  every f4/f5/f6/g1 consumer."
    python3 pipeline/f1_coupled_atac.py
else
    echo "--- Phase F1: SKIPPED ---"
fi

if should_run f2; then
    echo ""
    echo "--- Phase F2: MACS3 Peak Calling (SnapATAC2 built-in) ---"
    python3 pipeline/f2_peak_calling.py
else
    echo "--- Phase F2: SKIPPED ---"
fi

if should_run f2b; then
    echo ""
    echo "--- Phase F2b: gkmQC rescue (sub-500 pools) ---"
    python3 pipeline/f2b_gkmqc_rescue.py
else
    echo "--- Phase F2b: SKIPPED ---"
fi

if should_run f3; then
    echo ""
    echo "--- Phase F3: ATAC Embedding (TF-IDF + LSI + UMAP) ---"
    python3 pipeline/f3_atac_embed.py
else
    echo "--- Phase F3: SKIPPED ---"
fi

# ── Phase G: Differential Accessibility ────────────────────────────────
if should_run g1; then
    echo ""
    echo "--- Phase G1: Pseudobulk ATAC Aggregation ---"
    python3 pipeline/g1_pseudobulk_atac.py
else
    echo "--- Phase G1: SKIPPED ---"
fi

if should_run g2l; then
    echo ""
    echo "--- Phase G2l: HH17/16 LPM inference-pool DESeq2 ---"
    python3 pipeline/g2_lpm_hh1716_deseq2.py
else
    echo "--- Phase G2l: SKIPPED ---"
fi

echo ""
echo "=== Pipeline Complete (from: $START, to: $END) ==="
echo "Finished: $(date)"
