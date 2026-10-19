#!/bin/bash
# Adapted from the lab's CellBender environment builder. Paths are re-rooted onto
# DUCK_SCRATCH_USER / DUCK_SCRATCH_ROOT (Tier B, docs/audit_hardcoded_paths.md); logs +
# the patched checkpoint default to this repo's scratch mirror. The duck and chick
# CellBender runs depend on this env build.
#SBATCH -J cb_build
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/cb_build_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/cb_build_%j.err
#
# Step 0 — rebuild the scratch-purged CellBender conda env, REPRODUCIBLY.
# Pins are load-bearing (after 6 failed GPU runs):
#   python 3.10 | torch 2.2.2+cu121 | cellbender 0.3.2 (@v0.3.2 tag) | pyro-ppl 1.9.0
# Build runs on a CPU node (short); the GPU smoke test is run as a separate job.
#
# CRITICAL (learned 2026-06-01): pyro 1.9.0 does NOT by itself prevent the weakref
# pickle crash under torch 2.2 — torch.save(model)/torch.save(scheduler) both raise
# "cannot pickle 'weakref.ReferenceType'", aborting the checkpoint -> no ckpt.tar.gz
# -> the posterior step fails its existence assert -> NO _filtered.h5. The fix is the
# checkpoint.py patch applied below. The patched file is vendored next to this script
# (scripts/cellbender/) and is read from the repo copy on scratch;
# override with PATCHED_CK=<path>.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"

ENV="${DUCK_SCRATCH_USER}/conda_envs/cellbender"
LOGD="${DUCK_SCRATCH_ROOT}/logs"
SENTINEL="$LOGD/cb_build.done"
PATCHED_CK="${PATCHED_CK:-${DUCK_SCRATCH_ROOT}/scripts/cellbender/cellbender_checkpoint_v0.3.2_patched.py}"
mkdir -p "$LOGD"
rm -f "$SENTINEL"

echo "=== CellBender env build ==="
echo "host=$(hostname) start=$(date)"

module purge
module load gcc/14.2.0 conda/miniforge3

# Keep conda/pip caches off the 100 GB home quota (a documented past trap)
export CONDA_PKGS_DIRS="${DUCK_SCRATCH_USER}/conda_envs/.condapkgs"
export PIP_CACHE_DIR="${DUCK_SCRATCH_USER}/conda_envs/.pipcache"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"

# Clean slate: remove any scratch-purged husk so `conda create` won't refuse
rm -rf "$ENV"

set +u
conda create -y -p "$ENV" python=3.10
conda activate "$ENV"
set -u
echo "python: $(python -V 2>&1) at $(which python)"

# Install order is load-bearing:
#   1) torch 2.2.2 + torchvision 0.17.2 on cu121 (forward-compat with O2's newer driver)
#   2) cellbender pinned to the v0.3.2 git tag (PyPI only ever had 0.3.0)
#   3) pyro-ppl 1.9.0 LAST — force-downgrade off cellbender's 1.9.1 dep
echo "--- [1/3] torch 2.2.2+cu121 ---"
pip install torch==2.2.2 torchvision==0.17.2 --index-url https://download.pytorch.org/whl/cu121
echo "--- [2/3] cellbender @ v0.3.2 ---"
pip install 'cellbender @ git+https://github.com/broadinstitute/CellBender.git@v0.3.2'
echo "--- [3/3] pyro-ppl 1.9.0 (pinned LAST) ---"
pip install pyro-ppl==1.9.0

# ---- apply the weakref-pickle checkpoint patch (REQUIRED) ----
CK="$ENV/lib/python3.10/site-packages/cellbender/remove_background/checkpoint.py"
PRISTINE_MD5=d1c22d4e6ecff300826c1a004727eec2
got=$(md5sum "$CK" | awk '{print $1}')
if grep -q '__cellbender_unpicklable_placeholder__' "$CK"; then
  echo "checkpoint.py already patched"
elif [ "$got" = "$PRISTINE_MD5" ]; then
  [ -f "$CK.orig" ] || cp "$CK" "$CK.orig"
  cp "$PATCHED_CK" "$CK"
  echo "applied checkpoint.py patch (pristine backup at $CK.orig)"
else
  echo "ERROR: installed checkpoint.py md5 $got != pristine v0.3.2 $PRISTINE_MD5"
  echo "       refusing to drop in the patched file blindly (upstream may have changed)."
  exit 1
fi

echo "=== validate imports + pins + patch ==="
python - <<'PY'
import inspect
import pyro, torch, cellbender
import cellbender.remove_background.checkpoint as ck
print(f"pyro {pyro.__version__} | torch {torch.__version__} | cellbender {cellbender.__version__}")
assert pyro.__version__ == "1.9.0", f"pyro {pyro.__version__} != 1.9.0 (1.9.1 crashes at checkpoint save)"
assert torch.__version__.startswith("2.2.2"), f"torch {torch.__version__} != 2.2.2"
assert "cu121" in torch.__version__, f"torch not cu121 (a CPU wheel leaked in): {torch.__version__}"
assert "__cellbender_unpicklable_placeholder__" in inspect.getsource(ck), "checkpoint patch NOT applied"
if cellbender.__version__ != "0.3.2":
    print(f"WARN: cellbender {cellbender.__version__} != 0.3.2 (pinned @v0.3.2 — verify)")
print("VERSIONS + PATCH OK")
PY

cellbender remove-background --help >/dev/null 2>&1 && echo "cellbender CLI OK" || { echo "cellbender CLI FAILED"; exit 1; }

touch "$SENTINEL"
echo "=== BUILD DONE $(date) ==="
