#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-4:00
#SBATCH -c 4
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/cAE_atlas_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/cAE_atlas_%j.err
#SBATCH -J cAE_atlas
#
# cAE whole-embryo chick atlas v2 (duck-d1 parity + d2b annotation + multi-res + DE + h5ad).
#   sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAE_atlas.sh
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "=== ensure scikit-misc.loess (seurat_v3_paper HVG; required — no fallback) ==="
# Probe the COMPILED extension (skmisc.loess), not just the package — the numpy-ABI
# mismatch only fires when the C ext loads. vireo has numpy 1.x (dtype size 88); the
# default scikit-misc wheel is built for numpy 2.x (96) -> pin 0.1.4 (numpy-1.x ABI),
# --no-deps so vireo's numpy is never touched (it's the demux env, must not break).
python3 -c "from skmisc.loess import loess" 2>/dev/null && echo "skmisc.loess OK" || {
  echo "skmisc missing/ABI-broken; reinstalling scikit-misc==0.1.4 (numpy-1.x ABI, --no-deps)..."
  pip install --quiet --force-reinstall --no-deps "scikit-misc==0.1.4" 2>&1 | tail -3 || true
  python3 -c "from skmisc.loess import loess" 2>/dev/null && echo "skmisc.loess OK after reinstall" \
    || { echo "skmisc STILL BROKEN -> aborting (seurat_v3_paper HVG is required)"; exit 1; }
}
echo "=== cAE_full_atlas.py ==="
# CellBender-corrected counts are the committed default-of-record (cb_restage builds matrix_cb/);
# override to CHICK_MATRIX_SUBDIR=matrix for the raw-count comparison atlas.
export CHICK_MATRIX_SUBDIR="${CHICK_MATRIX_SUBDIR:-matrix_cb}"
echo "CHICK_MATRIX_SUBDIR=$CHICK_MATRIX_SUBDIR"
python3 ${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/cAE_full_atlas.py
rc=$?
echo "=== done (rc=$rc) ==="
exit $rc
