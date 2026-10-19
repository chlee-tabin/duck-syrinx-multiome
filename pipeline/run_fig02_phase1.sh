#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/fig02_phase1_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/fig02_phase1_%j.err
#SBATCH -J fig02_p1

# Figure 2 restructure — Phase 1: Data preparation
# All steps are independent and lightweight.
#
# Steps:
#   1A: Export per-genotype metrics (pct_ribo, cell cycle, etc.)
#   1B: Parse cellranger metrics_summary.csv
#   1D: Ribosomal gene orthology check
#
# Usage:
#   sbatch pipeline/run_fig02_phase1.sh           # all steps
#   sbatch pipeline/run_fig02_phase1.sh 1a        # single step
#   sbatch pipeline/run_fig02_phase1.sh 1b
#   sbatch pipeline/run_fig02_phase1.sh 1d

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

STEP="${1:-all}"

echo "=== Figure 2 Phase 1: Data Preparation ==="
echo "Date: $(date)"
echo "Step: $STEP"
echo ""

if [[ "$STEP" == "all" || "$STEP" == "1a" ]]; then
    echo "--- 1A: Export per-genotype metrics ---"
    python3 pipeline/d2b_export_genotype_metrics.py
    echo ""
fi

if [[ "$STEP" == "all" || "$STEP" == "1b" ]]; then
    echo "--- 1B: Parse cellranger metrics ---"
    python3 pipeline/d2c_parse_cellranger_metrics.py
    echo ""
fi

if [[ "$STEP" == "all" || "$STEP" == "1d" ]]; then
    echo "--- 1D: Ribosomal gene check ---"
    python3 FC_08656/20260401_ribo_gene_check.py
    echo ""
fi

echo "=== Phase 1 Complete ==="
echo "Finished: $(date)"
