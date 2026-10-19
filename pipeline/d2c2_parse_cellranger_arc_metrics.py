#!/usr/bin/env python3
"""Phase D2c2: Parse cellranger-arc summary.csv across all 12 ARC libraries.

ATAC-side analog of d2c_parse_cellranger_metrics.py.  Feeds figS03
panel A (per-library snATAC QC) the way d2c feeds figS01 panel B
(per-library snRNA QC).

cellranger-arc reports ATAC metrics in `outs/summary.csv` (single-row
CSV with both ATAC- and GEX-prefixed columns from the multiome run).
We extract the ATAC half; the GEX half is already covered by d2c
against the standalone GEX runs.

Note on "saturation": cellranger-arc does NOT report a single
"Sequencing Saturation" column for ATAC.  The closest proxy is
`ATAC Percent duplicates`, which rises as library complexity is
exhausted; (1 - pct_duplicates) is the unique-fragment rate.

Uses stdlib only (csv module) so it runs on a cluster transfer node
without conda / pandas.

Input:  <SCRATCH>/<LIB>/outs/summary.csv  (× 12 libraries)
Output: <FIGDIR>/step01_cellranger_arc_metrics.tsv

Run on a node that sees the cellranger-arc outs/ (file parsing only):
    python3 pipeline/d2c2_parse_cellranger_arc_metrics.py
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import csv
import os
import sys
import time

# Single source of truth for the 12-ARC-library list.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
from _shared.library_naming import ARC_LIB_ORDER  # noqa: E402


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Library name → cellranger-arc output directory on scratch (post-v2.1.0
# rerun; docs/cellranger_arc_v210_rerun.md, not distributed).  Mirrors d2c's GEX
# mapping so per-library bars in figS03 panel A align 1:1 with figS01
# panel B.  Library names happen to match their on-disk directories, so
# this is an identity map sourced from the canonical ARC_LIB_ORDER.
LIB_DIRS = {lib: lib for lib in ARC_LIB_ORDER}

# (cellranger-arc header) → (output column name).  Order is preserved.
METRICS = [
    ("Estimated number of cells",                                  "est_cells"),
    ("ATAC Sequenced read pairs",                                  "atac_sequenced_pairs"),
    ("ATAC Mean raw read pairs per cell",                          "atac_mean_reads_per_cell"),
    ("ATAC Median high-quality fragments per cell",                "atac_median_hq_fragments_per_cell"),
    ("ATAC Fraction of high-quality fragments in cells",           "atac_frac_hq_fragments_in_cells"),
    ("ATAC Fraction of high-quality fragments overlapping peaks",  "atac_frac_hq_fragments_in_peaks"),
    ("ATAC Fraction of high-quality fragments overlapping TSS",    "atac_frac_hq_fragments_in_tss"),
    ("ATAC TSS enrichment score",                                  "atac_tss_enrichment"),
    ("ATAC Percent duplicates",                                    "atac_pct_duplicates"),
    ("ATAC Number of peaks",                                       "atac_n_peaks_arc"),
    ("ATAC Fraction of genome in peaks",                           "atac_frac_genome_in_peaks"),
    ("ATAC Fraction of transposition events in peaks in cells",   "atac_frac_tn5_in_peaks_in_cells"),
]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_cellranger_value(s):
    """Parse a cellranger metric string (handles commas, percent signs)."""
    if s is None:
        return ""
    s = s.strip()
    if not s:
        return ""
    if s.endswith("%"):
        try:
            return float(s.rstrip("%").replace(",", "")) / 100.0
        except ValueError:
            return s
    s2 = s.replace(",", "")
    try:
        return float(s2)
    except ValueError:
        return s


def read_summary(csv_path):
    """Read single-row cellranger-arc summary.csv → {header: value} dict.

    Fails loud on a multi-row summary (e.g., a future dual-species or
    multi-target variant of cellranger-arc) — silently dropping extra
    rows would mask a schema change.
    """
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        row = next(reader)
        extras = list(reader)
    if extras:
        raise RuntimeError(
            f"read_summary: expected single-row summary at {csv_path}, "
            f"found {1 + len(extras)} rows.  cellranger-arc schema may have "
            "changed; the parser needs an audit before this is safe to use."
        )
    return {h: v for h, v in zip(header, row)}


def _lookup_metric(summary, header):
    """Look up a cellranger-arc metric header in a summary dict.

    Tries exact match first, then case- and whitespace-insensitive
    match (mirroring ``d2c_parse_cellranger_metrics.py:88-95``).  Raises
    ``KeyError`` listing available columns if no match — a cellranger-arc
    version bump renaming a column would otherwise produce a TSV full
    of empty cells that downstream figures render as missing bars.
    """
    if header in summary:
        return summary[header]
    target = header.lower().strip()
    for k in summary:
        if k.lower().strip() == target:
            return summary[k]
    raise KeyError(
        f"_lookup_metric: cellranger-arc summary missing metric "
        f"{header!r}; available columns = {sorted(summary.keys())}"
    )


def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return "" if v is None else str(v)


def main():
    ts("=== Phase D2c2: Parse cellranger-arc metrics ===")

    out_rows = []
    missing = []
    for lib_name, lib_dir in LIB_DIRS.items():
        csv_path = os.path.join(SCRATCH, lib_dir, "outs", "summary.csv")
        if not os.path.exists(csv_path):
            ts(f"  WARNING: {csv_path} not found, skipping {lib_name}")
            missing.append(lib_name)
            continue

        summary = read_summary(csv_path)
        ts(f"  {lib_name}: {len(summary)} metrics in summary.csv")

        row = {"library": lib_name}
        for header, colname in METRICS:
            row[colname] = parse_cellranger_value(_lookup_metric(summary, header))
        out_rows.append(row)

    ts(f"  Parsed {len(out_rows)} libraries"
       + (f" (missing: {', '.join(missing)})" if missing else ""))
    for row in out_rows:
        est = row.get("est_cells")
        dup = row.get("atac_pct_duplicates")
        tsse = row.get("atac_tss_enrichment")
        est_s = f"{est:,.0f}" if isinstance(est, float) else "?"
        dup_s = f"{dup:.1%}" if isinstance(dup, float) else "?"
        tsse_s = f"{tsse:.2f}" if isinstance(tsse, float) else "?"
        ts(f"    {row['library']}: {est_s} cells, dup={dup_s}, TSSe={tsse_s}")

    out = f"{FIGDIR}/step01_cellranger_arc_metrics.tsv"
    cols = ["library"] + [c for _, c in METRICS]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in out_rows:
            f.write("\t".join(fmt(row.get(c)) for c in cols) + "\n")
    ts(f"  Saved: {out}")
    ts("=== Phase D2c2 Complete ===")


if __name__ == "__main__":
    main()
