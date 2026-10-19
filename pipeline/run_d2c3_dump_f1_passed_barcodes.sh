#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-0:30
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/d2c3_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/d2c3_%j.err
#SBATCH -J d2c3_dump_f1_barcodes

# Chained job: d2c3 reads each per-library f1 h5ad (backed-read only
# touches obs_names; fast) and writes step_f1_passed_barcodes.tsv;
# d2c2_sample then re-runs to dual-reservoir-sample per (lib × qc_set)
# using that barcode TSV.  Doing both in one job keeps the derived
# data set atomic — the cluster submission wrapper syncs both TSVs back together.
# 8G / 30min is comfortable headroom (d2c3 ~minutes, d2c2_sample ~10min
# for the 12 × 88-MB per_barcode_metrics streams).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


cd ${DUCK_SCRATCH_ROOT}

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== d2c3 dump f1-passed barcodes ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
python3 pipeline/d2c3_dump_f1_passed_barcodes.py
echo "=== d2c3 complete ==="
echo ""
echo "=== d2c2_sample (dual reservoir using new barcode TSV) ==="
python3 pipeline/d2c2_sample_arc_per_cell_metrics.py
echo "=== d2c2_sample complete ==="
echo "Finished: $(date)"
