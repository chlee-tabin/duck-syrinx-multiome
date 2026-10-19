#!/bin/bash
#SBATCH -p priority
#SBATCH -t 0-0:10
#SBATCH -c 1
#SBATCH --mem 14G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/cb_restage_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/cb_restage_%j.err
#SBATCH -J cb_restage
# Re-stage CellBender-corrected counts into matrix_cb/ (cb_restage.py), in the vireo env
# (h5py + scanpy; the h5 read is cross-checked between the two). CHICK_RESTAGE_LIBS selects libraries (default: all four).
# Sizing measured: job 56118015 MaxRSS 9.8 GB, 4m58s for GEX09+GEX10 (2026-10-13).
# Another user: export DUCK_SCRATCH_USER and pass --output=/--error= by hand (SLURM reads #SBATCH first).
set -euo pipefail
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
# Default = the repo copy staged on scratch; refuse a stale pre-cross-check copy.
: "${CB_RESTAGE_PY:=${DUCK_SCRATCH_ROOT}/pipeline/chick/cb_restage.py}"
[ -f "$CB_RESTAGE_PY" ] || { echo "ERROR: CB_RESTAGE_PY not found: $CB_RESTAGE_PY" >&2; exit 1; }
grep -q "MIN_CB_OVERLAP" "$CB_RESTAGE_PY" || { echo "ERROR: $CB_RESTAGE_PY is a stale copy (no MIN_CB_OVERLAP guard)" >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "Host: $(hostname); CHICK_RESTAGE_LIBS=${CHICK_RESTAGE_LIBS:-all}"
echo "CB_RESTAGE_PY=$CB_RESTAGE_PY"
python3 "$CB_RESTAGE_PY"
echo "=== done: $(date) ==="
