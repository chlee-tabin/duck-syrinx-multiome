#!/usr/bin/env python3
"""Phase D2c2 (companion): subsample per-cell ATAC QC metrics from
cellranger-arc per_barcode_metrics.csv for figS03 panel A violins.

cellranger-arc reports library-level means/medians in summary.csv, but
the per-barcode CSV holds per-cell distributions that drive better QC
visualization (violin plots > summary bars).  This script streams each
library's per_barcode_metrics.csv (~88 MB / library), filters
`is_cell == 1`, deterministically subsamples up to N cells per library,
and writes a compact per-cell TSV.

Per-cell metrics emitted (6 numerics + library label + qc_set):
  atac_raw_reads             # A2 violin
  atac_fragments             # A4 violin
  atac_dup_rate              # A3 violin = atac_dup_reads / atac_raw_reads
  atac_tss_frac              # A5 violin (TSS-fragment fraction, NOT
                              # the library TSS enrichment score)
  atac_frip                  # cellranger-arc per-library-peak FRIP =
                              #   atac_peak_region_fragments / atac_fragments
                              #   (retained for provenance; was the A6
                              #   stand-in before the ArchR union landed)
  atac_frip_archr            # A6 violin = FRIP at the ArchR 283,818-peak
                              #   union, joined per-barcode from f2e
                              #   (step_f2e_archr_frip_per_cell.tsv.gz);
                              #   blank when that table is absent
  qc_set                     # "arc"  = sampled from cellranger-arc
                              #          `is_cell == 1` set
                              # "f1"   = sampled from QC-passed subset
                              #          (a.k.a. "f1-passed";
                              #          TSS ≥ 2 & fragments ≥ 1000;
                              #          requires step_f1_passed_barcodes
                              #          .tsv to exist; otherwise only
                              #          the "arc" set is emitted)

The qc_set split mirrors figS03 panel A1's cellranger-arc-vs-f1-passed
dodged bars: A2-A6 render side-by-side violins so the reader can see
how each per-cell QC distribution changes when restricted to the
cells that survive downstream snATAC QC.

Stdlib only (csv + gzip + random).  Runs on a cluster transfer node without
conda.

Inputs (both on scratch):
  - <SCRATCH>/<lib>/outs/per_barcode_metrics.csv   (per-cell metrics)
  - <FIGDIR>/step_f1_passed_barcodes.tsv           (f1 qc_set, optional)
  - <FIGDIR>/step_f2e_archr_frip_per_cell.tsv.gz   (ArchR-union FRIP,
        produced by f2e_archr_frip_per_cell.py; optional — atac_frip_archr
        is left blank when absent)

Output: <FIGDIR>/step01_cellranger_arc_per_cell_sample.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import csv
import gzip
import os
import sys
import random
import time

# Single source of truth for the 12-ARC-library list.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
from _shared.library_naming import ARC_LIB_ORDER  # noqa: E402


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Library names happen to match on-disk directories; reuse the canonical
# order so any new ARC library only needs to be added in one place
# (figures/_shared/library_naming.py).
LIB_DIRS = {lib: lib for lib in ARC_LIB_ORDER}

# Subsample size per library.  200 points is plenty for a violin KDE.
# Seed fixed for determinism.
N_PER_LIB = 200
SEED = 42

# Required columns in cellranger-arc per_barcode_metrics.csv.  Mismatch
# here means a cellranger-arc version bump renamed something; raise
# rather than silently dropping every row.
REQUIRED_PER_BARCODE_COLS = (
    "barcode", "is_cell",
    "atac_raw_reads", "atac_fragments", "atac_dup_reads",
    "atac_TSS_fragments", "atac_peak_region_fragments",
)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Reservoir:
    """Algorithm R reservoir sampler — single-pass, O(N) time, O(k) memory.

    Wrapping the state in a class lets us fan one stream into multiple
    independent reservoirs (cellranger-arc set + f1-passed subset).
    """

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
    """Load step_f1_passed_barcodes.tsv (produced by d2c3) → dict of
    {arc_library: set(barcodes)}.  Returns None if the file is absent
    (back-compat fallback: only the "arc" reservoir gets sampled).
    Raises if the file exists but contains zero data rows — a
    header-only TSV means d2c3 ran but failed to read any h5ad, and
    silently producing an arc-only sample in that case would mask
    the upstream failure."""
    path = os.path.join(figdir, "step_f1_passed_barcodes.tsv")
    if not os.path.exists(path):
        return None
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["arc_library", "barcode"]:
            raise RuntimeError(
                f"_load_f1_barcode_sets: unexpected header {header} in "
                f"{path}; expected (arc_library, barcode)"
            )
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out.setdefault(parts[0], set()).add(parts[1])
    if not out:
        raise RuntimeError(
            f"_load_f1_barcode_sets: {path} exists but contains zero data "
            "rows.  This typically means d2c3_dump_f1_passed_barcodes ran "
            "but every per-library f1 h5ad was missing.  Re-run d2c3 "
            "after confirming f1 outputs at "
            "<scratch>/phase2/atac/per_library/<arc>.h5ad."
        )
    return out


def _load_archr_frip_lookup(figdir):
    """Load f2e per-cell ArchR-union FRIP → {(library, barcode): frip_str}.

    Produced by ``pipeline/f2e_archr_frip_per_cell.py`` (the union-peak FRIP
    for every cellranger-arc-called cell).  Values are kept as the
    already-4dp-rounded strings from f2e so the sample TSV reproduces them
    verbatim.  Returns ``None`` when the file is absent — back-compat: the
    sampler still emits both reservoirs, but ``atac_frip_archr`` is left
    blank and figS03 A6 will show "n/a" markers until f2e runs.  Raises if
    the file exists but is header-only (f2e ran yet produced no rows — a
    silent empty join would mask that failure, the same discipline as
    ``_load_f1_barcode_sets``)."""
    path = os.path.join(figdir, "step_f2e_archr_frip_per_cell.tsv.gz")
    if not os.path.exists(path):
        return None
    lookup = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            lib_i = header.index("library")
            bc_i = header.index("barcode")
            frip_i = header.index("frip_archr")
        except ValueError:
            raise RuntimeError(
                f"_load_archr_frip_lookup: {path} is missing one of the "
                f"library/barcode/frip_archr columns; header={header}.  "
                "Regenerate via pipeline/f2e_archr_frip_per_cell.py."
            )
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(lib_i, bc_i, frip_i):
                continue
            lookup[(parts[lib_i], parts[bc_i])] = parts[frip_i]
    if not lookup:
        raise RuntimeError(
            f"_load_archr_frip_lookup: {path} exists but contains zero data "
            "rows.  f2e_archr_frip_per_cell.py ran but emitted nothing — "
            "re-run it before sampling."
        )
    return lookup


def _is_cell_truthy(value):
    """Normalize per-barcode ``is_cell`` to a boolean.

    cellranger-arc 2.0/2.1 emits 0/1, but related 10x tools have shipped
    ``True``/``False``; accept both rather than silently coercing every
    row to non-cell on a column-format change.

    Keep in sync with ``f1c_amulet_doublets._is_cell_truthy``.
    """
    return str(value or "").strip().lower() in ("1", "true")


def stream_cells(csv_path):
    """Yield (barcode, raw_reads, fragments, dup_reads, tss_fragments,
    peak_fragments) for every is_cell-truthy row.

    Validates the CSV schema once at the top of the file (raising on
    missing required columns).  Per-row int-cast ``ValueError`` is the
    only tolerated exception — a row where one count parses as a float
    string is rare-but-real for older cellranger output.

    The barcode is yielded so the caller can route rows into qc-set-
    specific reservoirs (cellranger-arc only vs. f1-pass intersection).
    """
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in REQUIRED_PER_BARCODE_COLS
                   if c not in (reader.fieldnames or [])]
        if missing:
            raise KeyError(
                f"stream_cells: per_barcode_metrics.csv at {csv_path} is "
                f"missing required columns {missing}; available "
                f"{reader.fieldnames}.  cellranger-arc schema may have "
                "changed."
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
            except ValueError as exc:
                # a non-integer metrics row is a malformed cellranger-arc output: fail loud
                raise SystemExit(f"ERROR: non-integer per_cell_metrics row: {row}") from exc


def safe_ratio(num, den):
    """Return num/den, or NaN if den is non-positive (avoids the
    silent ``0/0 → 0`` masking that would pull the violin median down
    if degenerate cells ever sneaked past the is_cell filter)."""
    return (num / den) if den > 0 else float("nan")


def _emit_row(lib_name, qc_set, item, archr_frip_lookup):
    """Convert a (barcode, raw, frags, dup, tss_frags, peak_frags)
    tuple into the per-cell TSV row dict.  The barcode is used internally
    for qc-set routing and to join the ArchR-union FRIP, but is not
    written to the TSV.

    ``archr_frip_lookup`` is the f2e ``{(library, barcode): frip_str}``
    map (or ``None`` when f2e hasn't run).  ``atac_frip_archr`` is the
    union-peak FRIP that figS03 A6 plots; it falls back to "" when the
    lookup is absent or the barcode is missing from it."""
    barcode, raw, frags, dup, tss_frags, peak_frags = item
    if archr_frip_lookup is None:
        archr_frip = ""
    else:
        archr_frip = archr_frip_lookup.get((lib_name, barcode), "")
    return {
        "library":             lib_name,
        "atac_raw_reads":      raw,
        "atac_fragments":      frags,
        "atac_dup_rate":       round(safe_ratio(dup, raw), 4),
        "atac_tss_frac":       round(safe_ratio(tss_frags, frags), 4),
        "atac_frip":           round(safe_ratio(peak_frags, frags), 4),
        "atac_frip_archr":     archr_frip,
        "qc_set":              qc_set,
    }


def main():
    ts(f"=== Phase D2c2 sampler: per-cell ARC metrics "
       f"(N={N_PER_LIB}/lib/qc_set, seed={SEED}) ===")

    f1_barcodes_by_lib = _load_f1_barcode_sets(FIGDIR)
    if f1_barcodes_by_lib is None:
        ts("  NOTE: step_f1_passed_barcodes.tsv not found — emitting "
           "qc_set='arc' only (run d2c3_dump_f1_passed_barcodes first to "
           "unlock the f1-passed sample for figS03 A2-A6 overlay)")
    else:
        sizes = ", ".join(f"{k}={len(v):,}" for k, v in f1_barcodes_by_lib.items())
        ts(f"  Loaded f1-passed barcode sets: {sizes}")

    archr_frip_lookup = _load_archr_frip_lookup(FIGDIR)
    if archr_frip_lookup is None:
        ts("  NOTE: step_f2e_archr_frip_per_cell.tsv.gz not found — "
           "emitting atac_frip_archr blank (run f2e_archr_frip_per_cell "
           "first to populate the figS03 A6 ArchR-union FRIP violins)")
    else:
        ts(f"  Loaded ArchR-union FRIP for {len(archr_frip_lookup):,} cells")

    rng = random.Random(SEED)
    out_rows = []
    missing = []
    for lib_name, lib_dir in LIB_DIRS.items():
        csv_path = os.path.join(SCRATCH, lib_dir, "outs", "per_barcode_metrics.csv")
        if not os.path.exists(csv_path):
            ts(f"  WARNING: {csv_path} not found, skipping {lib_name}")
            missing.append(lib_name)
            continue

        t0 = time.time()
        # Two reservoirs: cellranger-arc set (all is_cell=1 rows) and
        # f1-pass subset (intersected with the f1 barcode set).  Single
        # pass over per_barcode_metrics.csv.  Both reservoirs share the
        # same seeded rng — each reservoir's Algorithm R draws are
        # statistically coupled (the shared rng advances on every
        # `add`), but the selection rule (`j < self.k`) is independent
        # per stream, so the samples remain unbiased and the seeded run
        # is fully reproducible end-to-end.
        arc_reservoir = Reservoir(N_PER_LIB, rng)
        f1_reservoir = Reservoir(N_PER_LIB, rng) if f1_barcodes_by_lib else None
        f1_set = (f1_barcodes_by_lib or {}).get(lib_name, set())

        for item in stream_cells(csv_path):
            barcode = item[0]
            arc_reservoir.add(item)
            if f1_reservoir is not None and barcode in f1_set:
                f1_reservoir.add(item)

        n_arc = arc_reservoir.n_seen
        for item in arc_reservoir.items:
            out_rows.append(_emit_row(lib_name, "arc", item, archr_frip_lookup))
        if f1_reservoir is not None:
            n_f1 = f1_reservoir.n_seen
            for item in f1_reservoir.items:
                out_rows.append(_emit_row(lib_name, "f1", item, archr_frip_lookup))
            ts(f"  {lib_name}: arc={n_arc:,} cells "
               f"({len(arc_reservoir.items)} sampled), "
               f"f1={n_f1:,} cells "
               f"({len(f1_reservoir.items)} sampled) "
               f"in {time.time()-t0:.1f}s")
        else:
            ts(f"  {lib_name}: arc={n_arc:,} cells "
               f"({len(arc_reservoir.items)} sampled) in "
               f"{time.time()-t0:.1f}s")

    # Sanity-check the ArchR-FRIP join: every sampled cell should have a
    # union-peak FRIP when the f2e table is present.  A non-trivial miss
    # rate means a barcode-format mismatch between f2e and the sampler —
    # surface it loudly rather than letting A6 quietly render "n/a"s.
    if archr_frip_lookup is not None and out_rows:
        n_missing = sum(1 for r in out_rows if r["atac_frip_archr"] == "")
        if n_missing:
            ts(f"  WARNING: {n_missing}/{len(out_rows)} sampled cells had no "
               "ArchR-union FRIP (barcode absent from f2e table) — check for "
               "barcode-format drift between f2e and d2c2_sample")
        else:
            ts(f"  ArchR-union FRIP joined for all {len(out_rows):,} "
               "sampled cells")

    out = f"{FIGDIR}/step01_cellranger_arc_per_cell_sample.tsv"
    cols = ["library", "atac_raw_reads", "atac_fragments",
            "atac_dup_rate", "atac_tss_frac", "atac_frip", "atac_frip_archr",
            "qc_set"]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in out_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    size = os.path.getsize(out)
    ts(f"  Wrote {len(out_rows)} rows × {len(cols)} cols → {out} ({size:,} bytes)"
       + (f" — missing libs: {', '.join(missing)}" if missing else ""))
    ts("=== Sampler complete ===")


if __name__ == "__main__":
    main()
