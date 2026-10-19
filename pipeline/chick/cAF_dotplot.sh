#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:30
#SBATCH -c 4
#SBATCH --mem 48G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/cAF_dotplot_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/cAF_dotplot_%j.err
#SBATCH -J cAF_dotplot
#
# cAF — export chick atlas dotplot over the DUCK fig01 marker panel.
#   sbatch cAF_dotplot.sh <umap_tsv_gz> <out_tsv>
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "Host: $(hostname); python: $(command -v python3)"
# CellBender counts are the committed default-of-record (matches cAE/cAG); override to matrix for raw.
export CHICK_MATRIX_SUBDIR="${CHICK_MATRIX_SUBDIR:-matrix_cb}"
echo "CHICK_MATRIX_SUBDIR=$CHICK_MATRIX_SUBDIR"
UMAP="$1"; OUT="$2"
PY=${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/cAF_duck_panel_dotplot.py
python3 "$PY" --umap "$UMAP" --out "$OUT"
echo "=== done: $(date) ==="
