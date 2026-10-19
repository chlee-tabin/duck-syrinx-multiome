#!/bin/bash
# cb_patch_subset_coo.sh — fix the CellBender 0.3.2 `_subset_coo` zero-size bug.
#
# CellBender `remove_background/estimation.py::_subset_coo` builds `sp.coo_matrix((data, (row,col)))`
# without `shape=`. When the boolean `logic` mask selects ZERO entries (an empty posterior chunk —
# which the chick ARC GEX-subset produces after --exclude-feature-types Peaks, because some nuclei
# have ~no GEX counts), scipy raises "cannot infer dimensions from zero sized index arrays" and the
# whole run dies in compute_denoised_counts. Passing the explicit shape fixes it; for non-empty
# chunks the inferred shape == coo.shape, so this is backward-compatible (retina/duck unaffected).
#
# Idempotent + backed up. Run once on the shared scratch env.
set -e
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

ENV=${DUCK_SCRATCH_USER}/conda_envs/cellbender
F=$ENV/lib/python3.10/site-packages/cellbender/remove_background/estimation.py
[ -f "$F" ] || { echo "ABORT: $F not found"; exit 1; }
[ -e "$F.bak_pre_shapefix" ] || cp "$F" "$F.bak_pre_shapefix"   # keep the first pristine copy
"$ENV/bin/python" - "$F" <<'PY'
import sys
f = sys.argv[1]
s = open(f).read()
old = "return sp.coo_matrix((coo.data[logic], (coo.row[logic], coo.col[logic])))"
new = "return sp.coo_matrix((coo.data[logic], (coo.row[logic], coo.col[logic])), shape=coo.shape)"
if new in s and old not in s:
    print("already patched; no change")
else:
    assert s.count(old) == 1, f"expected exactly 1 occurrence of the buggy line, found {s.count(old)}"
    open(f, "w").write(s.replace(old, new))
    print("PATCHED _subset_coo: added shape=coo.shape")
PY
echo "--- verify ---"
grep -n "shape=coo.shape" "$F"
