#!/usr/bin/env python3
"""Phase F Step 4b: per-base CFC1 locus coverage from existing f4
fragment BEDs.

Lightweight sibling of f4_lr_pitx2_pileup.py.  Re-uses the (stage_anchor,
side) fragment BEDs that f4 already exported for HH13/HH21/HH25 in
``phase2/atac/lr_pileup/fragments/``, and produces the per-base
coverage track at the CFC1 locus (NC_051804.1:45,101,912-45,205,798,
±50 kb around the gene body).

Why a sibling instead of extending f4?
- f4's fragment export is the expensive step (~30 GB of BEDs); this
  script only needs to scan the BEDs, which is fast.
- Mirrors the f5_cfc1 / f6_cfc1 sibling convention.

Stages covered
--------------
HH13 / HH21 / HH25 (existing aggregate-side fragment BEDs).  HH17/16
is per-genotype-only in the chromatin-canonical run
(``lr_pileup_per_genotype/``) so its L vs R aggregate at HH17/16 is
NOT computed here — defer to a separate aggregation step if the CFC1
trajectory analysis needs it.

Output
------
  figures/data/step_f4b_cfc1_coverage.tsv.gz
    Columns: chrom, pos, <group>, <group>_percell ... for each
    (stage_anchor, side) group present in lr_pileup/fragments/.

  figures/data/step_f4b_cfc1_cell_counts.tsv
    Columns: group, n_cells, frag_bed, coverage_sum_locus,
    coverage_sum_per_cell.

Run on O2
---------
  sbatch pipeline/run_f4b_cfc1_locus.sh
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
LR_PILEUP_DIR = f"{ATAC_DIR}/lr_pileup"
FRAG_DIR = f"{LR_PILEUP_DIR}/fragments"
F4_CELL_COUNTS_TSV = f"{LR_PILEUP_DIR}/step_f4_cell_counts.tsv"

OUT_TSV = f"{SCRATCH}/figures/data/step_f4b_cfc1_coverage.tsv.gz"
OUT_PROV_TSV = f"{SCRATCH}/figures/data/step_f4b_cfc1_cell_counts.tsv"

# ── Locus ─────────────────────────────────────────────────────────────────
# CFC1 gene: NC_051804.1:45,151,913-45,155,798 (1-based GTF, - strand).
# locus_50kb|CFC1: NC_051804.1:45,101,912-45,205,798 (BED 0-based).
#
# Coordinate convention here: this is a per-base coverage scan that
# uses 0-based BED indexing for array-friendly arithmetic
# (`pos = start + i`).  This is DIFFERENT from f5_cfc1 / f6_cfc1, which
# use 1-based-inclusive coords passed to SnapATAC2's region-string
# parser (see `references/cfc1_locus_regions.bed` for the canonical
# 0-based-half-open SoT and `pipeline/f5_cfc1_region_counts.py`'s
# REGIONS docblock for the convention map).  CFC1_TSS is the GTF
# 1-based TSS coordinate (also equal to the 0-based BED `end` of the
# 1-bp `tss|XM_013100073.4|CFC1` interval), and the promoter window
# `CFC1_TSS ± 300` arithmetic is symmetric so it's 0/1-based-agnostic.
CFC1_CHROM = "NC_051804.1"
CFC1_LOCUS_START = 45_101_912
CFC1_LOCUS_END   = 45_205_798
CFC1_GENE_START  = 45_151_912   # for marking the gene body in render
CFC1_GENE_END    = 45_155_798
CFC1_TSS         = 45_155_798   # minus strand → TSS at gene end
# TODO: f5/f6 region counters load coords from references/cfc1_locus_
# regions.bed via _locus_bed.load_locus_regions; this file's 0-based
# per-base scan uses a different convention (BED-native 0-based half-
# open vs SnapATAC2's 1-based inclusive) so the helper-based
# unification was deferred.  When refactored, expose a parallel
# load_locus_regions_0based or have callers convert the 1-based output,
# then drop these six constants.


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def cfc1_coverage_from_bed(bed_gz: Path, locus_size: int) -> np.ndarray:
    """Integer per-base coverage over the CFC1 ±50 kb window."""
    cov = np.zeros(locus_size, dtype=np.int64)
    opener = gzip.open if str(bed_gz).endswith(".gz") else open
    with opener(bed_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] != CFC1_CHROM:
                continue
            frag_start = int(parts[1])
            frag_end = int(parts[2])
            if frag_end <= CFC1_LOCUS_START or frag_start >= CFC1_LOCUS_END:
                continue
            s = max(frag_start, CFC1_LOCUS_START) - CFC1_LOCUS_START
            e = min(frag_end, CFC1_LOCUS_END) - CFC1_LOCUS_START
            cov[s:e] += 1
    return cov


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frag-dir", default=FRAG_DIR,
                        help=f"Directory containing <group>.bed.gz files "
                             f"from f4.  Default {FRAG_DIR}")
    parser.add_argument("--cell-counts-tsv", default=F4_CELL_COUNTS_TSV,
                        help="f4's step_f4_cell_counts.tsv (for n_cells "
                             "lookup per group).  Default reuses the "
                             "lr_pileup/ provenance file.  Pass "
                             "/dev/null to opt out of per-cell "
                             "normalization (rare).")
    parser.add_argument("--out-tsv", default=OUT_TSV,
                        help=f"Output coverage TSV path (gzip if .gz).")
    parser.add_argument("--out-prov-tsv", default=OUT_PROV_TSV,
                        help="Output provenance TSV path.")
    parser.add_argument("--groups", default="",
                        help="Comma-separated group names to process "
                             "(default: all <group>.bed.gz in frag-dir).")
    args = parser.parse_args()

    ts("=== Phase F Step 4b: CFC1 locus per-base coverage ===")
    ts(f"  frag-dir = {args.frag_dir}")
    ts(f"  out-tsv = {args.out_tsv}")
    ts(f"  CFC1 locus = {CFC1_CHROM}:{CFC1_LOCUS_START:,}-{CFC1_LOCUS_END:,}  "
       f"({CFC1_LOCUS_END - CFC1_LOCUS_START:,} bp)")

    frag_dir = Path(args.frag_dir)
    if not frag_dir.is_dir():
        sys.exit(f"ERROR: frag-dir {frag_dir} not found")

    # Discover groups
    if args.groups:
        groups = [g.strip() for g in args.groups.split(",") if g.strip()]
    else:
        groups = sorted(p.stem.replace(".bed", "")
                        for p in frag_dir.glob("*.bed.gz"))
    ts(f"  Groups: {groups}")
    if not groups:
        sys.exit(f"ERROR: no .bed.gz files found in {frag_dir}")

    # Load f4 cell counts for n_cells lookup.  Required for per-cell
    # normalization, which downstream pile-up panels rely on.  Pass
    # --cell-counts-tsv /dev/null to explicitly opt out (then per-cell
    # tracks will be NaN by design).
    cell_counts: dict[str, int] = {}
    if args.cell_counts_tsv == "/dev/null":
        ts(f"  Per-cell normalization disabled (--cell-counts-tsv /dev/null).")
    elif os.path.isfile(args.cell_counts_tsv):
        f4_prov = pd.read_csv(args.cell_counts_tsv, sep="\t")
        for _, r in f4_prov.iterrows():
            cell_counts[str(r["group"])] = int(r["n_cells"])
        ts(f"  Loaded n_cells for {len(cell_counts)} groups from "
           f"{args.cell_counts_tsv}")
    else:
        sys.exit(f"ERROR: cell counts TSV not found at "
                 f"{args.cell_counts_tsv}.  Required for per-cell "
                 f"normalization in downstream pile-up panels.  Pass "
                 f"--cell-counts-tsv /dev/null to opt out explicitly.")

    # Compute coverage per group
    locus_size = CFC1_LOCUS_END - CFC1_LOCUS_START
    cov_by_group: dict[str, np.ndarray] = {}
    for grp in groups:
        bed = frag_dir / f"{grp}.bed.gz"
        if not bed.exists():
            ts(f"  SKIP {grp}: BED not found at {bed}")
            continue
        ts(f"  Scanning {grp} ({bed.stat().st_size / 1e6:.1f} MB)...")
        t0 = time.time()
        cov = cfc1_coverage_from_bed(bed, locus_size)
        cov_by_group[grp] = cov
        n = cell_counts.get(grp, np.nan)
        ts(f"    {grp}: sum={cov.sum():,}  max={cov.max():,}  "
           f"density={cov.sum() / locus_size:.3f}/bp  "
           f"n_cells={n if not isinstance(n, float) or not np.isnan(n) else '(?)'}  "
           f"{time.time() - t0:.1f}s")

    if not cov_by_group:
        sys.exit("ERROR: no coverage produced — check frag-dir contents.")

    # Build TSV: chrom, pos, <group>, <group>_percell columns
    positions = np.arange(CFC1_LOCUS_START, CFC1_LOCUS_END, dtype=np.int64)
    out_df = pd.DataFrame({"chrom": CFC1_CHROM, "pos": positions})
    valid_groups = list(cov_by_group.keys())
    for grp in valid_groups:
        out_df[grp] = cov_by_group[grp]
        n = cell_counts.get(grp)
        if n and n > 0:
            out_df[f"{grp}_percell"] = cov_by_group[grp] / n
        else:
            out_df[f"{grp}_percell"] = np.nan

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False,
                  compression="gzip" if str(out_path).endswith(".gz") else None)
    ts(f"  Saved coverage TSV: {out_path}  "
       f"({out_path.stat().st_size / 1e6:.1f} MB)")

    # Provenance TSV
    prov = pd.DataFrame({
        "group": valid_groups,
        "n_cells": [cell_counts.get(g, np.nan) for g in valid_groups],
        "frag_bed": [str(frag_dir / f"{g}.bed.gz") for g in valid_groups],
        "coverage_sum_locus": [int(cov_by_group[g].sum()) for g in valid_groups],
        "coverage_sum_per_cell": [
            (cov_by_group[g].sum() / cell_counts[g]) if g in cell_counts else np.nan
            for g in valid_groups
        ],
    })
    prov_path = Path(args.out_prov_tsv)
    prov_path.parent.mkdir(parents=True, exist_ok=True)
    prov.to_csv(prov_path, sep="\t", index=False)
    ts(f"  Provenance: {prov_path}")

    # Quick sanity print: sum of coverage in promoter window per group.
    promoter_start = CFC1_TSS - 300
    promoter_end   = CFC1_TSS + 300
    ts("Promoter window sanity (TSS ± 300 bp):")
    p0 = promoter_start - CFC1_LOCUS_START
    p1 = promoter_end   - CFC1_LOCUS_START
    for grp in valid_groups:
        n = cell_counts.get(grp)
        sub = cov_by_group[grp][p0:p1]
        s = int(sub.sum())
        per = (s / n) if (n and n > 0) else float("nan")
        ts(f"    {grp}: promoter sum={s:,}  per-cell={per:.3f}")

    ts("=== Done ===")


if __name__ == "__main__":
    main()
