#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/atac_chain_archr_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/atac_chain_archr_%j.err
#SBATCH -J atac_hh1716_chain_archr

# Re-run with ArchR's 283K peak set (instead of the SnapATAC2 25K set).
#
# Step 0: rebuild cell × peak matrix using ArchR peaks (f2d)
# Step 1: g1 pseudobulk on the ArchR peak matrix (env F2_H5AD)
# Step 2: g2l HH17/16 inference pool DESeq2
# Step 3: g4 motif enrichment
# Step 4: g5 wave-2 anticipation
# Step 5: g3 PITX2 locus deep-dive
#
# Outputs land at separate paths (suffix _archr) so SnapATAC2-track
# results stay intact for comparison.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
# Pin the d2n L/R sidecar so a stale G1_LR_TSV from the submit shell can't silently
# force g1 onto the legacy d2l axis under SLURM --export=ALL (g1's own default).
export G1_LR_TSV="$SCRATCH/figures/data/lr_inference_multiome.tsv"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== atac_hh1716_chain_archr ==="
date
echo "Host: $(hostname)"

# Pre-flight checks
for f in phase2/atac/step_f1_combined.h5ads \
         phase2/atac_archr/peaks/union_peaks.bed; do
    if [ ! -f "$f" ]; then
        echo "FATAL: required input missing: $f" >&2
        exit 1
    fi
done

# ── f2d: build the ArchR cell × peak matrix ─────────────────────────
echo ""
echo "--- f2d: cell × ArchR-peak matrix ---"
python3 pipeline/f2d_archr_peak_matrix.py

# Point downstream at the new matrix + ArchR peak BED via env vars
export F2_H5AD="$SCRATCH/phase2/atac/step_f2_peak_matrix_archr.h5ad"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export PEAKS_BED="$SCRATCH/phase2/atac_archr/peaks/union_peaks.bed"
mkdir -p "$DA_DIR"

echo "F2_H5AD=$F2_H5AD"
echo "DA_DIR=$DA_DIR"
echo "PEAKS_BED=$PEAKS_BED"

# ── g1: pseudobulk on ArchR peaks ───────────────────────────────────
echo ""
echo "--- g1: pseudobulk count matrix (ArchR peaks) ---"
python3 pipeline/g1_pseudobulk_atac.py

# ── g2l: HH17/16 LPM inference-pool DESeq2 ──────────────────────────
echo ""
echo "--- g2l: HH17/16 LPM inference-pool DESeq2 (ArchR peaks) ---"
python3 pipeline/g2_lpm_hh1716_deseq2.py

# Switch module for R-side step
module unload conda/miniforge3
module load R/4.4.2

# ── g4: motif enrichment ────────────────────────────────────────────
echo ""
echo "--- g4_hh1716: motif enrichment (ArchR peaks) ---"
Rscript pipeline/g4_hh1716_motif_enrichment.R

# Back to conda for Python steps
module unload R/4.4.2
module load conda/miniforge3
conda activate snapatac2

# ── g5: HH21 wave-2 anticipation ────────────────────────────────────
echo ""
echo "--- g5: HH21 wave-2 anticipation cross-modal ---"
python3 pipeline/g5_hh21_wave2_anticipation.py

echo ""
echo "=== atac_hh1716_chain_archr DONE ==="
date
