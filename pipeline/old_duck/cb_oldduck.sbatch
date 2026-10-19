#!/bin/bash
# cb_oldduck.sbatch — CellBender remove-background on the 3 old-duck (FC_07344/FC_07361) E6/E9/E11
# 3' scRNA libraries (D1/D2/D3), to close the methods-parity gap with the multiome (which IS
# CellBender-corrected) and test whether ambient cleaning sharpens the E06 de-novo L/R contrast.
# Faithful port of pipeline/chick/cb_cellbender.sbatch (retina 0.3.2 recipe). These are STANDALONE
# 3' GEX libraries (NOT ARC) so NO --exclude-feature-types.
#
# SLURM ARRAY: one GPU task per library. Idempotent (skips libs whose *_filtered.h5 exists).
# Stage raw matrices to scratch first (done), then: sbatch --parsable cb_oldduck.sbatch
#
#SBATCH -J cb_oldduck
#SBATCH -p gpu_requeue
#SBATCH --gres=gpu:1
#SBATCH -t 0-8:00
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH --array=0-2
#SBATCH -o /n/scratch/users/c/cl266/old_duck_review/cellbender/logs/cb_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_review/cellbender/logs/cb_%A_%a.err
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


ENV=${DUCK_SCRATCH_USER}/conda_envs/cellbender
ROOT=${DUCK_SCRATCH_USER}/old_duck_review/cellbender
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

# name:expected_cells (cellranger 'Estimated Number of Cells'), per metrics_summary.csv.
# D1=FC_07361/D1_zju1 (E6+E9), D2=FC_07344/D2_zju1 (E6+E9), D3=FC_07344/D3_zju1 (E11).
LIBS=( "D1:28742" "D2:15005" "D3:4028" )
IFS=':' read -r LIB EC <<< "${LIBS[$SLURM_ARRAY_TASK_ID]}"
RAW="$ROOT/raw/$LIB/raw_feature_bc_matrix.h5"
OUT="$OUTDIR/${LIB}_cellbender"
echo "[$LIB] expected_cells=$EC raw=$RAW host=$(hostname) start=$(date)"

if [ -f "${OUT}_filtered.h5" ]; then echo "[$LIB] already done; skipping"; exit 0; fi
[ -f "$RAW" ] || { echo "[$LIB] MISSING raw: $RAW"; exit 1; }

TD=$((EC * 2)); if [ "$TD" -gt 100000 ]; then TD=100000; fi

# Per-library CWD: CellBender writes its checkpoint to a fixed name in CWD; array tasks sharing one
# CWD would clobber each other's ckpt (corrupting gpu_requeue preemption-resume). --output is absolute.
WORKDIR="$OUTDIR/$LIB"; mkdir -p "$WORKDIR"; cd "$WORKDIR"

# Run via cb_wrapper.py (in-process _subset_coo zero-size monkeypatch; harmless for full GEX libs).
python "$ROOT/cb_wrapper.py" remove-background \
    --input "$RAW" \
    --output "${OUT}.h5" \
    --cuda \
    --expected-cells "$EC" \
    --total-droplets-included "$TD" \
    --epochs 150 \
    --fpr 0.01 \
    --learning-rate 0.00005

if [ -f "${OUT}_filtered.h5" ]; then echo "[$LIB] DONE $(date)"; else echo "[$LIB] FAILED — no _filtered.h5"; exit 1; fi
