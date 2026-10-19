#!/usr/bin/env python3
"""hQ2 — subsample per-cell ATAC QC metrics from cellranger-arc
per_barcode_metrics.csv for the chick snATAC QC panel (figS04_supp E2-E6).

Chick analog of the duck ``d2c2_sample_arc_per_cell_metrics.py``.  Streams each
chick ARC library's per_barcode_metrics.csv, filters ``is_cell``, and
deterministically reservoir-samples up to N cells per (library x qc_set):

  qc_set = "arc" : all cellranger-arc called cells
  qc_set = "f1"  : the subset that ALSO passes the snATAC QC gate
                   (TSS enrichment >= 2 & fragments >= 1000), read from
                   chick_atac_f1_passed_barcodes.tsv (produced by hQ3).  When
                   that file is absent, only the "arc" reservoir is emitted.

Per-cell metrics (mirroring the duck sampler exactly):
  atac_raw_reads, atac_fragments,
  atac_dup_rate  = atac_dup_reads / atac_raw_reads,
  atac_tss_frac  = atac_TSS_fragments / atac_fragments,
  atac_frip      = atac_peak_region_fragments / atac_fragments

Stdlib only (csv + random); runs on a cluster transfer node without conda.

Input:  $CHICK_QC_BASE/<lib>/outs/per_barcode_metrics.csv (staged by hQ0)
        $CHICK_QC_BASE/figdata/chick_atac_f1_passed_barcodes.tsv (hQ3, optional)
Output: $CHICK_QC_BASE/figdata/chick_atac_per_cell_sample.tsv
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import csv
import glob
import os
import random
import time

QC_BASE = os.environ.get("CHICK_QC_BASE",
                         f"{SCRATCH_USER}/chick_demux/atac_qc")
FIGDIR = os.path.join(QC_BASE, "figdata")
os.makedirs(FIGDIR, exist_ok=True)


def discover_libs():
    """Each staged library == a <QC_BASE>/<lib>/outs dir (self-contained;
    the figure re-orders by the canonical lib SoT)."""
    return [os.path.basename(os.path.dirname(o))
            for o in sorted(glob.glob(os.path.join(QC_BASE, "*", "outs")))]

N_PER_LIB = 200
SEED = 42

REQUIRED_PER_BARCODE_COLS = (
    "barcode", "is_cell",
    "atac_raw_reads", "atac_fragments", "atac_dup_reads",
    "atac_TSS_fragments", "atac_peak_region_fragments",
)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Reservoir:
    """Algorithm R reservoir sampler — single-pass, O(N) time, O(k) memory."""

    def __init__(self, k, rng):
        self.k = k
        self.rng = rng
        self.n_seen = 0
        self.items = []

    def add(self, item):
        if self.n_seen < self.k:
            self.items.append(item)
        else:
            j = self.rng.randint(0, self.n_seen)
            if j < self.k:
                self.items[j] = item
        self.n_seen += 1


def _load_f1_barcode_sets(figdir):
    """Load chick_atac_f1_passed_barcodes.tsv -> {arc_library: set(barcodes)}.
    Returns None if absent (arc-only fallback).  Raises on a header-only file
    (hQ3 ran but read no per-library QC h5ad — silently emitting arc-only there
    would mask the upstream failure)."""
    path = os.path.join(figdir, "chick_atac_f1_passed_barcodes.tsv")
    if not os.path.exists(path):
        return None
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["arc_library", "barcode"]:
            raise RuntimeError(
                f"_load_f1_barcode_sets: unexpected header {header} in {path}; "
                f"expected (arc_library, barcode)"
            )
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out.setdefault(parts[0], set()).add(parts[1])
    if not out:
        raise RuntimeError(
            f"_load_f1_barcode_sets: {path} exists but has zero data rows — "
            "hQ3 likely produced no per-library f1 barcodes.  Re-run hQ3."
        )
    return out


def _is_cell_truthy(value):
    """cellranger-arc 2.0/2.1 emit 0/1; accept True/False too."""
    return str(value or "").strip().lower() in ("1", "true")


def stream_cells(csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in REQUIRED_PER_BARCODE_COLS
                   if c not in (reader.fieldnames or [])]
        if missing:
            raise KeyError(
                f"stream_cells: per_barcode_metrics.csv at {csv_path} is "
                f"missing required columns {missing}; available "
                f"{reader.fieldnames}.  cellranger-arc schema may have changed."
            )
        for row in reader:
            if not _is_cell_truthy(row.get("is_cell")):
                continue
            try:
                yield (
                    row["barcode"],
                    int(row["atac_raw_reads"]),
                    int(row["atac_fragments"]),
                    int(row["atac_dup_reads"]),
                    int(row["atac_TSS_fragments"]),
                    int(row["atac_peak_region_fragments"]),
                )
            except ValueError:
                continue


def safe_ratio(num, den):
    return (num / den) if den > 0 else float("nan")


def _emit_row(lib_name, qc_set, item):
    _, raw, frags, dup, tss_frags, peak_frags = item
    return {
        "library":        lib_name,
        "atac_raw_reads": raw,
        "atac_fragments": frags,
        "atac_dup_rate":  round(safe_ratio(dup, raw), 4),
        "atac_tss_frac":  round(safe_ratio(tss_frags, frags), 4),
        "atac_frip":      round(safe_ratio(peak_frags, frags), 4),
        "qc_set":         qc_set,
    }


def main():
    ts(f"=== hQ2: sample chick per-cell ARC metrics "
       f"(N={N_PER_LIB}/lib/qc_set, seed={SEED}) ===")

    f1_barcodes_by_lib = _load_f1_barcode_sets(FIGDIR)
    if f1_barcodes_by_lib is None:
        ts("  NOTE: chick_atac_f1_passed_barcodes.tsv not found — emitting "
           "qc_set='arc' only (run hQ3 first to unlock the f1 overlay)")
    else:
        sizes = ", ".join(f"{k}={len(v):,}" for k, v in f1_barcodes_by_lib.items())
        ts(f"  Loaded f1-passed barcode sets: {sizes}")

    rng = random.Random(SEED)
    out_rows = []
    missing = []
    for lib_name in discover_libs():
        csv_path = os.path.join(QC_BASE, lib_name, "outs", "per_barcode_metrics.csv")
        if not os.path.exists(csv_path):
            ts(f"  WARNING: {csv_path} not found, skipping {lib_name}")
            missing.append(lib_name)
            continue

        t0 = time.time()
        arc_reservoir = Reservoir(N_PER_LIB, rng)
        f1_reservoir = Reservoir(N_PER_LIB, rng) if f1_barcodes_by_lib else None
        f1_set = (f1_barcodes_by_lib or {}).get(lib_name, set())

        for item in stream_cells(csv_path):
            barcode = item[0]
            arc_reservoir.add(item)
            if f1_reservoir is not None and barcode in f1_set:
                f1_reservoir.add(item)

        for item in arc_reservoir.items:
            out_rows.append(_emit_row(lib_name, "arc", item))
        if f1_reservoir is not None:
            for item in f1_reservoir.items:
                out_rows.append(_emit_row(lib_name, "f1", item))
            ts(f"  {lib_name}: arc={arc_reservoir.n_seen:,} "
               f"({len(arc_reservoir.items)} sampled), "
               f"f1={f1_reservoir.n_seen:,} "
               f"({len(f1_reservoir.items)} sampled) in {time.time()-t0:.1f}s")
        else:
            ts(f"  {lib_name}: arc={arc_reservoir.n_seen:,} "
               f"({len(arc_reservoir.items)} sampled) in {time.time()-t0:.1f}s")

    out = os.path.join(FIGDIR, "chick_atac_per_cell_sample.tsv")
    cols = ["library", "atac_raw_reads", "atac_fragments",
            "atac_dup_rate", "atac_tss_frac", "atac_frip", "qc_set"]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in out_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    ts(f"  Wrote {len(out_rows)} rows -> {out}"
       + (f" — missing libs: {', '.join(missing)}" if missing else ""))
    ts("=== hQ2 complete ===")


if __name__ == "__main__":
    main()
