#!/bin/bash
# cb_cellbender.sbatch — CellBender remove-background on the 4 chick atlas libraries, to close
# the methods-parity gap (the duck is CellBender-corrected; the chick was not) and test the
# suspected GEX10-library ambient PITX2 floor. Port of a CellBender 0.3.2
# recipe used in an earlier mouse retina project.
#
# SLURM ARRAY: one GPU task per library. Idempotent (skips libs whose *_filtered.h5 exists).
# Group A (ARC01/02) are multiome (GEX+Peaks) → --exclude-feature-types Peaks so CellBender models
# GEX ambient only (matches the GEX-only GEX09/10 and the snRNA atlas usage).
#
# Stage raw matrices to scratch first (compute nodes can't read /n/standby), then:
#   copy them to scratch ; sbatch --parsable cb_cellbender.sbatch
#
#SBATCH -J cb_chick
#SBATCH -p gpu_requeue
#SBATCH --gres=gpu:1
#SBATCH -t 0-8:00
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH --array=0-3
#SBATCH -o /n/scratch/users/c/cl266/chick_cellbender/logs/cb_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/chick_cellbender/logs/cb_%A_%a.err
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


ENV=${DUCK_SCRATCH_USER}/conda_envs/cellbender
ROOT=${DUCK_SCRATCH_USER}/chick_cellbender
OUTDIR="$ROOT/out"; mkdir -p "$OUTDIR" "$ROOT/logs"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate "$ENV"
export NUMBA_CACHE_DIR="/tmp/numba_cache_$$"

# Fail-fast: correct env + GPU, or abort before wasting the slot.
python3 - <<'PY'
import pyro, torch, cellbender
assert pyro.__version__ == "1.9.0", f"pyro {pyro.__version__} != 1.9.0 (1.9.1 crashes at checkpoint save)"
assert torch.cuda.is_available(), "no CUDA/GPU visible"
print(f"env OK: pyro {pyro.__version__} | torch {torch.__version__} | cellbender {cellbender.__version__} | GPU {torch.cuda.get_device_name(0)}")
PY

# name:expected_cells (cellranger 'Estimated Number of Cells'). GEX09/10 = standalone cr7 GEX;
# ARC01/02 = cellranger-arc v2.1.0 GEX cell calls (the v210 the atlas used).
LIBS=( "GEX09:25100" "GEX10:12622" "ARC01:4337" "ARC02:5732" )
IFS=':' read -r LIB EC <<< "${LIBS[$SLURM_ARRAY_TASK_ID]}"
RAW="$ROOT/raw/$LIB/raw_feature_bc_matrix.h5"
OUT="$OUTDIR/${LIB}_cellbender"
echo "[$LIB] expected_cells=$EC raw=$RAW host=$(hostname) start=$(date)"

if [ -f "${OUT}_filtered.h5" ]; then echo "[$LIB] already done; skipping"; exit 0; fi
[ -f "$RAW" ] || { echo "[$LIB] MISSING raw: $RAW"; exit 1; }

TD=$((EC * 2)); if [ "$TD" -gt 100000 ]; then TD=100000; fi
EXTRA=""; case "$LIB" in ARC*) EXTRA="--exclude-feature-types Peaks";; esac

# Per-library CWD: CellBender writes its checkpoint to a fixed name in CWD; array tasks sharing one
# CWD would clobber each other's ckpt (corrupting gpu_requeue preemption-resume). --output is absolute.
WORKDIR="$OUTDIR/$LIB"; mkdir -p "$WORKDIR"; cd "$WORKDIR"

# Run via cb_wrapper.py: in-process monkeypatch of the CellBender 0.3.2 _subset_coo zero-size bug
# (empty posterior chunks from the ARC GEX-subset). Does NOT modify the shared env.
# shellcheck disable=SC2086
python "$ROOT/cb_wrapper.py" remove-background \
    --input "$RAW" \
    --output "${OUT}.h5" \
    --cuda \
    --expected-cells "$EC" \
    --total-droplets-included "$TD" \
    --epochs 150 \
    --fpr 0.01 \
    --learning-rate 0.00005 $EXTRA

if [ -f "${OUT}_filtered.h5" ]; then echo "[$LIB] DONE $(date)"; else echo "[$LIB] FAILED — no _filtered.h5"; exit 1; fi
