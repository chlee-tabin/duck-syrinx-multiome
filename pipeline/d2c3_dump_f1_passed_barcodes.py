#!/usr/bin/env python3
"""Phase D2c3: Dump per-library f1-passed barcode lists.

Reads each ARC library's f1-output AnnData (``<scratch>/phase2/atac/
per_library/<arc>.h5ad``), which f1_coupled_atac.py has already
subset to cells passing TSS ≥ 2 & fragments ≥ 1000.  ``obs_names``
therefore IS the f1-passed barcode set.  We write one row per
(arc_library, barcode) to a small TSV that downstream consumers can
join against cellranger-arc per-barcode data.

This unblocks the figS03 panel A "f1-passed" overlay in A2-A6 (dodged
cellranger-arc vs f1 violins, mirroring A1's structure).
``pipeline/d2c2_sample_arc_per_cell_metrics.py`` consumes the output
to do dual reservoir sampling.

Reads (per library):
  <scratch>/phase2/atac/per_library/<arc_name>.h5ad   (snapatac2-backed)

Writes:
  <FIGDIR>/step_f1_passed_barcodes.tsv
    columns: arc_library, barcode

Run via sbatch:
  sbatch pipeline/run_d2c3_dump_f1_passed_barcodes.sh
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
import sys
import time

import snapatac2 as snap


SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
LIB_DIR = f"{ATAC_DIR}/per_library"
FIGDIR = f"{SCRATCH}/phase2/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Single source of truth for the 12-ARC-library list.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
from _shared.library_naming import ARC_LIB_ORDER  # noqa: E402


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase D2c3: dump f1-passed barcodes ===")

    out = f"{FIGDIR}/step_f1_passed_barcodes.tsv"
    missing = []
    total = 0
    with open(out, "w") as f:
        f.write("arc_library\tbarcode\n")
        for arc in ARC_LIB_ORDER:
            h5_path = os.path.join(LIB_DIR, f"{arc}.h5ad")
            if not os.path.exists(h5_path):
                ts(f"  WARNING: {h5_path} not found, skipping {arc}")
                missing.append(arc)
                continue

            # Backed-read is cheap when we only touch obs_names.
            ad = snap.read(h5_path, backed="r")
            try:
                barcodes = list(ad.obs_names)
            finally:
                ad.close()
            for bc in barcodes:
                f.write(f"{arc}\t{bc}\n")
            total += len(barcodes)
            ts(f"  {arc}: {len(barcodes):,} f1-passed barcodes")

    n_done = len(ARC_LIB_ORDER) - len(missing)
    ts(f"  Wrote {total:,} (arc_library, barcode) rows for {n_done}/"
       f"{len(ARC_LIB_ORDER)} libraries → {out}"
       + (f"  (missing: {', '.join(missing)})" if missing else ""))
    ts("=== Phase D2c3 Complete ===")


if __name__ == "__main__":
    main()
