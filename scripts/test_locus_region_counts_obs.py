#!/usr/bin/env python3
"""Self-check for pipeline/_locus_region_counts.py::load_obs_tsv (--from-cached-matrix).

A barcode repeated across libraries is legal (f1 stacks per-library h5ads without
prefixing obs_names); a barcode repeated WITHIN one library must raise.
snapatac2 is stubbed: the loader does not use it.
"""
import os
import sys
import tempfile
import types

sys.modules.setdefault("snapatac2", types.ModuleType("snapatac2"))
for _d in ("pipeline", "figures"):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", _d))
from _locus_region_counts import load_obs_tsv  # noqa: E402

HEAD = "barcode\tunified_genotype_id\tarc_library\td2g_stage\td2g_subcluster\n"


def write(rows):
    fh = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    fh.write(HEAD + "".join(rows))
    fh.close()
    return fh.name


ok = write(["AAAC-1\tg1\tARC05\tHH21\t3\n", "AAAC-1\tg2\tARC06\tHH21\t3\n", "TTTG-1\t\tARC06\t\t\n"])
barcodes, obs = load_obs_tsv(ok)
assert barcodes == ["AAAC-1", "AAAC-1", "TTTG-1"], barcodes
assert obs["arc_library"] == ["ARC05", "ARC06", "ARC06"] and obs["unified_genotype_id"][2] is None, obs

bad = write(["AAAC-1\tg1\tARC05\tHH21\t3\n", "AAAC-1\tg2\tARC05\tHH21\t3\n"])
try:
    load_obs_tsv(bad)
    raise AssertionError("duplicate (arc_library, barcode) did not raise")
except ValueError as e:
    assert "duplicated (arc_library, barcode)" in str(e), e
for p in (ok, bad):
    os.unlink(p)
print("load_obs_tsv self-check: PASS")
