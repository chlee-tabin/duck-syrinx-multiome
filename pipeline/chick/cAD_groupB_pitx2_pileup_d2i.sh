#!/bin/bash
#SBATCH -p priority
#SBATCH -t 0-0:10
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/cAD_pileup_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/cAD_pileup_%j.err
#SBATCH -J cAD_pileup
#
# Chick PITX2-locus pileup — FAITHFUL duck d2i/fig02D port (stage x population).
# Needs samtools (view/depth/faidx) + the vireo conda env (scanpy/anndata/harmonypy). Mirrors cZ/cAC.
#
#   stage the repo on scratch first: this .sh + the .py in ${DUCK_SCRATCH_ROOT}/pipeline/chick/
#   and chick_stage_alignment.tsv in ${DUCK_SCRATCH_ROOT}/figures/data/chick_control/
#   sbatch ${DUCK_SCRATCH_ROOT}/pipeline/chick/cAD_groupB_pitx2_pileup_d2i.sh
#   (DUCK_SCRATCH_ROOT default: /n/scratch/users/c/cl266/FC_08656. CAD_PY overrides the .py path.)
#   (another user: export DUCK_SCRATCH_USER=... and pass --output=/--error= by hand — SLURM
#    parses the #SBATCH -o/-e literals above before any shell runs.)
#
# Outputs go to CAD_OUT_DIR (default: chick_demux_B/vireo_seeded/seed_qc_cAD_harmonyfix/), NOT the
# published seed_qc/ (the python refuses seed_qc/ itself).
#
# Sizing: -c 1 (single-threaded python; BLAS/OMP pinned; samtools unthreaded). Measured on job
# 56118016 (2026-10-13): MaxRSS 5.4 GB, Elapsed 1m56s -> --mem 8G, -t 10 min.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
: "${CAD_OUT_DIR:=${DUCK_SCRATCH_USER}/chick_demux_B/vireo_seeded/seed_qc_cAD_harmonyfix}"
export CAD_OUT_DIR
# Default = the staging SSOT in the scratch copy of the repo (figures/data/).
: "${CAD_STAGE_TSV:=${DUCK_SCRATCH_ROOT}/figures/data/chick_control/chick_stage_alignment.tsv}"
export CAD_STAGE_TSV
[ -f "$CAD_STAGE_TSV" ] || { echo "ERROR: CAD_STAGE_TSV not found: $CAD_STAGE_TSV" >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# the committed Group B split-exclusion list (the scratch pipeline/chick/ copy has no repo tree);
# pass e.g. CHICK_SPLIT_EXCLUDED_TSV=<synced repo>/figures/data/chick_control/groupB_split_excluded_cells.tsv
: "${CHICK_SPLIT_EXCLUDED_TSV:?set CHICK_SPLIT_EXCLUDED_TSV to the committed groupB_split_excluded_cells.tsv}"
[ -s "$CHICK_SPLIT_EXCLUDED_TSV" ] || { echo "missing: $CHICK_SPLIT_EXCLUDED_TSV" >&2; exit 1; }
export CHICK_SPLIT_EXCLUDED_TSV

mkdir -p "${DUCK_SCRATCH_USER}/chick_demux/logs" "$CAD_OUT_DIR"
module purge
module load gcc/14.2.0 samtools conda/miniforge3
conda activate vireo
echo "Host: $(hostname); samtools: $(command -v samtools); CAD_OUT_DIR=$CAD_OUT_DIR"
# CellBender counts are the committed default-of-record (matches cAE/cAG/cAF); override to matrix for raw.
export CHICK_MATRIX_SUBDIR="${CHICK_MATRIX_SUBDIR:-matrix_cb}"
echo "CHICK_MATRIX_SUBDIR=$CHICK_MATRIX_SUBDIR"
# Default = the scratch copy of the repo (the .py is self-contained: no
# _chick_* helper imports, all data paths absolute). Refuse a stale pre-B6 copy.
: "${CAD_PY:=${DUCK_SCRATCH_ROOT}/pipeline/chick/cAD_groupB_pitx2_pileup_d2i.py}"
[ -f "$CAD_PY" ] || { echo "ERROR: CAD_PY not found: $CAD_PY" >&2; exit 1; }
grep -q "def harmony_lane" "$CAD_PY" || { echo "ERROR: $CAD_PY is a pre-B6 copy (no harmony_lane)" >&2; exit 1; }
echo "CAD_PY=$CAD_PY"
python3 "$CAD_PY"
echo "=== done: $(date) ==="
