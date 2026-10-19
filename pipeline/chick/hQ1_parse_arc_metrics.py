#!/usr/bin/env python3
"""hQ1 — parse cellranger-arc summary.csv across the staged chick ARC libraries.

Chick analog of the duck ``d2c2_parse_cellranger_arc_metrics.py``.  Feeds the
chick snATAC QC panel (figS04_supp panel E) its per-library summary metrics
(E1 dark bars = est_cells; plus library-level TSS-enrichment / FRIP / duplicate
fractions reported in summary.csv).

Self-contained (stdlib csv only, no repo import) so it runs from
``chick_demux/pipeline`` on a cluster transfer node without conda.  Libraries are
discovered by globbing the staged tree; the figure re-orders rows by its inline
ARC_LIB_ORDER (figures/figS04_supp/generate_libqc_atac.py).

Input:  $CHICK_QC_BASE/<lib>/outs/summary.csv   (staged by hQ0)
Output: $CHICK_QC_BASE/figdata/chick_atac_metrics.tsv
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import csv
import glob
import os
import time

QC_BASE = os.environ.get("CHICK_QC_BASE",
                         f"{SCRATCH_USER}/chick_demux/atac_qc")
FIGDIR = os.path.join(QC_BASE, "figdata")
os.makedirs(FIGDIR, exist_ok=True)

# (cellranger-arc header) -> (output column name).  Same set as the duck parser
# so the chick metrics TSV is schema-compatible with the figS03 reader.
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


def discover_libs():
    """Each staged library == a <QC_BASE>/<lib>/outs dir."""
    return [os.path.basename(os.path.dirname(o))
            for o in sorted(glob.glob(os.path.join(QC_BASE, "*", "outs")))]


def parse_cellranger_value(s):
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
    try:
        return float(s.replace(",", ""))
    except ValueError:
        return s


def read_summary(csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        row = next(reader)
        extras = list(reader)
    if extras:
        raise RuntimeError(
            f"read_summary: expected single-row summary at {csv_path}, found "
            f"{1 + len(extras)} rows.  cellranger-arc schema may have changed."
        )
    return {h: v for h, v in zip(header, row)}


def _lookup_metric(summary, header):
    if header in summary:
        return summary[header]
    target = header.lower().strip()
    for k in summary:
        if k.lower().strip() == target:
            return summary[k]
    raise KeyError(
        f"_lookup_metric: cellranger-arc summary missing metric {header!r}; "
        f"available columns = {sorted(summary.keys())}"
    )


def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return "" if v is None else str(v)


def main():
    ts("=== hQ1: parse chick cellranger-arc metrics ===")
    libs = discover_libs()
    ts(f"  discovered {len(libs)} staged libraries: {libs}")
    out_rows = []
    for lib_name in libs:
        csv_path = os.path.join(QC_BASE, lib_name, "outs", "summary.csv")
        if not os.path.exists(csv_path):
            ts(f"  WARNING: {csv_path} not found, skipping {lib_name}")
            continue
        summary = read_summary(csv_path)
        row = {"library": lib_name}
        for header, colname in METRICS:
            row[colname] = parse_cellranger_value(_lookup_metric(summary, header))
        out_rows.append(row)
        est = row.get("est_cells")
        tsse = row.get("atac_tss_enrichment")
        est_s = f"{est:,.0f}" if isinstance(est, float) else "?"
        tsse_s = f"{tsse:.2f}" if isinstance(tsse, float) else "?"
        ts(f"    {lib_name}: {est_s} cells, TSSe={tsse_s}")

    out = os.path.join(FIGDIR, "chick_atac_metrics.tsv")
    cols = ["library"] + [c for _, c in METRICS]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in out_rows:
            f.write("\t".join(fmt(row.get(c)) for c in cols) + "\n")
    ts(f"  Saved {len(out_rows)} rows -> {out}")
    ts("=== hQ1 complete ===")


if __name__ == "__main__":
    main()
