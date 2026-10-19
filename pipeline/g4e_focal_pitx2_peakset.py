#!/usr/bin/env python3
"""g4e_focal_pitx2_peakset.py — PITX2-locus focal peaks vs matched-random peaks.

ARM 2 of the right-peak BMP-SMAD probe.  The literal claim is "BMP regulates
dorsal PITX2 (a/b)," so the most specific test asks whether a direct SMAD1/5
footprint sits in the actual PITX2-locus regulatory chromatin.

This does NOT run a model — it emits a g2c-FORMAT TSV that ENCODES "PITX2-locus
peak = positive, every other ArchR peak = background," so the already-debugged
h2 (NEG_MATCH=gc, GC-matched random PEAK negatives) + h3b (degenerate SMAD core
readout) consume it UNCHANGED.  The control is therefore matched-random PEAKS
(not the genomic null), which sidesteps the CpG-promoter GC wall that broke the
h9 promoter arm.

Foreground  = the 359 ArchR peaks in the canonical *Pitx2*-TAD window
              (pitx2_locus_tad_peaks.tsv; NC_051775.1, FAM241A/C4orf32 upstream
              anchor -> PITX2 +500 kb downstream; Welsh 2015 Hi-C TAD anchor).
              This supersedes the legacy symmetric PITX2 +/-500 kb set (329 peaks,
              h0_lvr_peaks_pitx2_locus.tsv), which excluded the FAM241A boundary;
              the 329 are a strict subset (in_old_500kb column) of the 359.
Background  = every other peak in the ArchR 283,818-peak union (deduped).

Encoding (so h2 "strict" picks the focal set as positives):
  focal  peaks:  log2FC = 1.0, binom_p = 0.001, binom_padj = 0.01
  others       : log2FC = 0.0, binom_p = 1.0,   binom_padj = 1.0
h2 strict positives = log2FC>0 & binom_p<0.05 & binom_padj<0.1  -> the focal set.
h2 strict negatives = |log2FC|<0.1 & binom_p>0.5               -> the background.

Local script (inputs are committed in the repo); deterministic.  Run:
  conda run -n duckfig python3 pipeline/g4e_focal_pitx2_peakset.py
"""
from __future__ import annotations

import gzip
import os
import sys

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNION_BED = os.environ.get("UNION_BED", f"{REPO}/figures/data/archr_union_peaks.bed.gz")
# Default = the canonical 359-peak FAM241A-anchored *Pitx2*-TAD window (supersedes the
# legacy 329-peak h0_lvr_peaks_pitx2_locus.tsv; see the docstring). NOTE: the sibling
# trough engine g4h_response_locus_trough_test.py still defaults to the 329-peak h0 set
# (it reads its locus TSV directly, never via g4e — no shared consumer), so the two
# probes intentionally use different locus windows; don't assume a single canonical file.
FOCAL_TSV = os.environ.get("FOCAL_TSV", f"{REPO}/figures/data/pitx2_locus_tad_peaks.tsv")
OUT = os.environ.get("OUT", f"{REPO}/figures/data/g4e_focal_pitx2_g2c.tsv.gz")


def main() -> None:
    for p in (UNION_BED, FOCAL_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input missing: {p}")

    # ── Universe: dedupe the ArchR union to unique 501-bp intervals ──────────
    op = gzip.open if UNION_BED.endswith(".gz") else open
    universe = set()
    with op(UNION_BED, "rt") as fi:
        for line in fi:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            universe.add(f"{chrom}:{start}-{end}")
    print(f"Union peaks (deduped): {len(universe):,}")

    # ── Focal set: PITX2-locus peak_ids ─────────────────────────────────────
    focal_df = pd.read_csv(FOCAL_TSV, sep="\t")
    focal = set(focal_df["peak_id"].astype(str))
    in_universe = focal & universe
    missing = focal - universe
    print(f"Focal peaks: {len(focal):,}; in union: {len(in_universe):,}; "
          f"NOT in union: {len(missing):,}")
    if missing:
        # Coordinate-convention mismatch would silently shrink the positive set —
        # surface it loudly rather than train on a truncated focal set.
        ex = sorted(missing)[:5]
        print(f"  WARNING: {len(missing)} focal peak_id(s) absent from union "
              f"(coord mismatch?). examples: {ex}")
    if len(in_universe) < 50:
        sys.exit(f"ERROR: only {len(in_universe)} focal peaks found in union; "
                 f"too few for a gkm-SVM positive set — check peak_id conventions")

    # ── Emit g2c-format TSV ─────────────────────────────────────────────────
    rows = []
    for pid in sorted(universe):
        is_focal = pid in in_universe
        rows.append({
            "peak_id": pid,
            "in_locus": int(is_focal),
            "log2FC": 1.0 if is_focal else 0.0,
            "binom_p": 0.001 if is_focal else 1.0,
            "binom_padj": 0.01 if is_focal else 1.0,
            "abs_log2FC": 1.0 if is_focal else 0.0,
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.6g")
    # Focal peak-id list = h9 POSITIVE_PEAKS (Arm 2); the full g2c is the
    # negative universe (h9 NEG_SOURCE=peak draws GC-matched negatives from it).
    LIST_OUT = os.environ.get(
        "LIST_OUT",
        OUT.replace("_g2c.tsv.gz", "_peaks.txt").replace("_g2c.tsv", "_peaks.txt"))
    out.loc[out["in_locus"] == 1, "peak_id"].to_csv(LIST_OUT, index=False, header=False)
    print(f"Wrote {OUT}: {len(out):,} peaks "
          f"({int(out['in_locus'].sum()):,} focal positives, "
          f"{int((1-out['in_locus']).sum()):,} background)")
    print(f"Wrote {LIST_OUT}: {int(out['in_locus'].sum()):,} focal peak-ids")
    print("h9 (Arm 2): G2C_TSV=this g2c (universe), POSITIVE_PEAKS=this list, "
          "NEG_SOURCE=peak -> focal-locus peaks vs GC-matched random peaks")


if __name__ == "__main__":
    main()
