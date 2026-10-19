#!/usr/bin/env python3
"""Run CellBender with an IN-PROCESS fix for the 0.3.2 `_subset_coo` zero-size bug.

CellBender `remove_background/estimation.py::_subset_coo` builds `sp.coo_matrix((data,(row,col)))`
without `shape=`. When the boolean mask selects ZERO entries (an empty posterior chunk — produced by
the chick ARC GEX-subset after --exclude-feature-types Peaks, since some nuclei have ~no GEX counts),
scipy raises "cannot infer dimensions from zero sized index arrays" and the run dies in
compute_denoised_counts. The fix is to pass the explicit shape; for non-empty chunks the inferred
shape == coo.shape, so it is behavior-preserving.

This wrapper monkeypatches the function **at runtime, in this process only** — it does NOT modify the
shared conda env file (which retina/duck also use). Usage is identical to the CLI:
  python cb_wrapper.py remove-background --input ... --output ... [flags]
"""
import sys
import scipy.sparse as sp
import cellbender.remove_background.estimation as _est


def _subset_coo_fixed(coo, logic):
    return sp.coo_matrix((coo.data[logic], (coo.row[logic], coo.col[logic])), shape=coo.shape)


_est._subset_coo = _subset_coo_fixed
# belt-and-suspenders: patch any module that imported the name directly
for _modname in ("cellbender.remove_background.posterior",):
    try:
        _m = __import__(_modname, fromlist=["_subset_coo"])
        if hasattr(_m, "_subset_coo"):
            _m._subset_coo = _subset_coo_fixed
    except ImportError:   # module absent in this CellBender build: nothing to patch
        pass

print("cb_wrapper: monkeypatched estimation._subset_coo (shape=coo.shape) in-process", flush=True)

from cellbender.base_cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
