"""Shared backbone for f5 locus region-count pipelines.

`f5_pitx2_region_counts.py` and `f5_cfc1_region_counts.py` are
near-identical sibling drivers — same SnapATAC2 peak-matrix loop,
same FIRMNESS_ACCEPT, same STAGE_TO_ANCHOR / ARC_TO_D2L_LIBRARY
maps, same per-cell aggregation, same TSV schema.  They differ only
in the locus REGIONS dict, the output paths, and a few cosmetic
log strings.  This module owns the shared logic; each sibling becomes
a thin shim that calls `run_locus_region_counts(...)` with its own
constants.

Output TSV schema (long format; one row per
(genotype × library × side × stage × region) sample):
    genotype, library_d2l, side, stage, stage_raw, region,
    n_cells, n_frags, frags_per_cell
``library_d2l`` is the d2l library_m, except for a same-side duplicate (lane-loading-record
split tube: one embryo-half in two libraries), whose two libraries are
pooled into ONE sample labelled ``GEX04+GEX0S-duck`` (figures/_shared/lr_firmness.py).

Cell selection: pool membership is `POOLS[stage_anchor]["inference"]`
at runtime — i.e. the current POOLS (Splanchnic excluded)
flows in automatically when this module's caller is rerun on the
cluster.

Input modes
-----------
default            open the f1 dataset (``step_f1_combined.h5ads``) for the
                   per-cell obs, and build (or reuse) the cell x region matrix
                   from its fragments.
--from-cached-matrix OBS_TSV
                   the f1 per-library children are gone (purged from scratch),
                   so the dataset cannot be opened.  Per-cell obs come from
                   OBS_TSV (``step_f1_coupled_atac_metadata.tsv``, one row per f1
                   cell, index = barcode) and the counts from the existing cached
                   matrix ``peak_matrix_tmp``, opened read-only.  Raises unless
                   the cache exists, F5_FORCE_REBUILD is unset, the cache's
                   barcodes equal OBS_TSV's index position by position, and its
                   regions equal the requested ones in order.  Aggregation is the
                   default mode's, line for line.  Status: the committed f5
                   PITX2 table (cache built in August 2026) was NOT produced
                   by this flag.  It came from an
                   out-of-repo wrapper running the DEFAULT path with snapatac2's
                   readers monkeypatched (job log: ``[shim] read_dataset(...)``
                   then ``Reusing cached peak matrix``).  This flag is the in-repo
                   port of that wrapper's logic; on the cluster it reproduces
                   the committed table byte for byte (acceptance numbers:
                   pipeline/run_f5_pitx2_from_cache.sh, not distributed;
                   regression check: make check-locus-obs).
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from collections import Counter

import numpy as np
import pandas as pd
import snapatac2 as snap


from _shared.lr_firmness import resolve_sides, side_map  # noqa: E402  (callers put figures/ on sys.path)

# Cell-selection constant: the inferred-side firmness classes trusted here.
# HH21 G3 single-library genotypes (g48/g50/g51/g53 =
# G3_unfiltered_donor1/4/7/10) populate via `singleton_rank`.  `weak` rows are
# excluded by rule (logged); `same_side_dup` follows the shared rule (kept once,
# both libraries pooled); chromatin_-prefixed values are normalized and any
# unknown value raises (figures/_shared/lr_firmness.py).
FIRMNESS_ACCEPT = {"direct", "placement", "singleton_rank"}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_lr_map(tsv_path: str) -> dict[tuple[str, str], tuple[str, str]]:
    """(genotype, library_m) -> (side, sample library) for the usable rows of the
    multiome merger: FIRMNESS_ACCEPT inferred calls, plus same-side duplicates on
    their design side with both libraries mapped to ONE pooled sample label.
    Unknown firmness values / accepted rows without an L/R side raise."""
    lr = pd.read_csv(tsv_path, sep="\t")
    return side_map(resolve_sides(lr, FIRMNESS_ACCEPT, source=os.path.basename(tsv_path),
                                  log=ts))


# per-cell obs columns the aggregation reads
OBS_COLS = ("unified_genotype_id", "arc_library", "d2g_stage", "d2g_subcluster")


def load_obs_tsv(path: str) -> tuple[list[str], dict[str, list]]:
    """--from-cached-matrix obs source: (barcodes in file order, {col: values}) from
    the f1 coupled-ATAC metadata TSV.  Missing values become None, exactly as the
    dataset obs yields them (the aggregation drops a cell whose genotype or library
    ``is None``).  Raises on a missing column or a duplicated (arc_library, barcode)
    cell.  A barcode alone is NOT unique: f1's AnnDataSet stacks the per-library
    h5ads without prefixing obs_names, so the same 16-bp barcode recurs across
    libraries (the 2026-08-10 metadata TSV has such repeats; the first cluster run of
    this mode, 2026-10-16, raised on them).  Cell identity against the cache is
    checked positionally by check_cached_matrix."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"--from-cached-matrix obs TSV missing: {path}")
    meta = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    missing = [c for c in OBS_COLS if c not in meta.columns]
    if missing:
        raise KeyError(f"{path} lacks obs column(s) {missing}")
    barcodes = [str(b) for b in meta.index]
    cells = list(zip(meta["arc_library"].astype(str), barcodes))
    if len(set(cells)) != len(cells):
        raise ValueError(f"{path}: duplicated (arc_library, barcode) cells in the index")
    obs = {c: [None if pd.isna(v) else v for v in meta[c]] for c in OBS_COLS}
    return barcodes, obs


def check_cached_matrix(cache_barcodes: list[str], obs_barcodes: list[str],
                        cache_regions: list[str], regions: list[str],
                        cache_path: str) -> None:
    """--from-cached-matrix: the cache must hold the obs TSV's cells in the same order
    and the requested regions in the same order; anything else raises.  Only the
    barcode sequence is compared (the cache's obs_names carry no library); the library
    of each row comes from the obs TSV, whose (arc_library, barcode) pairs
    load_obs_tsv has checked are unique."""
    if cache_barcodes != obs_barcodes:
        n_diff = sum(a != b for a, b in zip(cache_barcodes, obs_barcodes))
        raise RuntimeError(
            f"{cache_path}: barcodes disagree with the obs TSV ({len(cache_barcodes):,} "
            f"vs {len(obs_barcodes):,} rows, {n_diff:,} positional differences)")
    if cache_regions != regions:
        raise RuntimeError(f"{cache_path}: regions {cache_regions} != requested {regions}")


def load_e1b_stage_map(tsv_path: str) -> dict[str, str]:
    """genotype_id -> inferred_stage from FIXED e1b sidecar."""
    if not os.path.isfile(tsv_path):
        sys.exit(f"ERROR: e1b sidecar missing at {tsv_path}.  "
                 f"Sync figures/data/e1b_stage_assignments.tsv to scratch.")
    df = pd.read_csv(tsv_path, sep="\t")
    return dict(zip(df["genotype_id"].astype(str),
                    df["inferred_stage"].astype(str)))


def run_locus_region_counts(
    *,
    label: str,
    regions: dict[str, tuple[str, int, int]],
    f1_dataset: str,
    lr_tsv: str,
    e1b_tsv: str,
    out_tsv_default: str,
    peak_matrix_tmp: str,
    include_stages_default: str,
    include_stages_help: str,
    arc_to_d2l_library: dict[str, str],
    stage_to_anchor: dict[str, str],
    pools,
) -> None:
    """Backbone for f5 locus region-count drivers.  See module docstring."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-stages",
                        default=include_stages_default,
                        help=include_stages_help)
    parser.add_argument("--out-tsv", default=out_tsv_default,
                        help=f"Output TSV path (default {out_tsv_default})")
    parser.add_argument("--from-cached-matrix", metavar="OBS_TSV", default=None,
                        help="Do not open the f1 dataset: per-cell obs from OBS_TSV "
                             "(step_f1_coupled_atac_metadata.tsv), counts from the "
                             f"existing cache {peak_matrix_tmp} opened read-only "
                             "(barcodes and regions checked; see module docstring)")
    args = parser.parse_args()

    include_stages = [s.strip() for s in args.include_stages.split(",")
                      if s.strip()]
    out_tsv = args.out_tsv

    ts(f"=== Phase F Step 5 ({label}): "
       f"per-(genotype × library) region counts ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")
    ts(f"  include_stages = {include_stages}")
    ts(f"  out_tsv = {out_tsv}")

    cached_mode = args.from_cached_matrix is not None
    if cached_mode:
        if os.environ.get("F5_FORCE_REBUILD"):
            raise RuntimeError("--from-cached-matrix with F5_FORCE_REBUILD set: the "
                               "cache cannot be rebuilt without the f1 dataset")
        if not os.path.isfile(peak_matrix_tmp):
            raise FileNotFoundError(f"--from-cached-matrix needs the cache {peak_matrix_tmp}")
    elif not os.path.isfile(f1_dataset):
        sys.exit(f"ERROR: {f1_dataset} missing")
    if not os.path.isfile(lr_tsv):
        sys.exit(f"ERROR: {lr_tsv} missing — sync figures/data/.  "
                 f"Run pipeline/d2o_merge_multiome_inference.py first.")

    if cached_mode:
        ts(f"Loading per-cell obs from {args.from_cached_matrix} (--from-cached-matrix)...")
        obs_barcodes, obs = load_obs_tsv(args.from_cached_matrix)
        ts(f"  {len(obs_barcodes):,} cells")
    else:
        ts("Loading f1 dataset...")
        dataset = snap.read_dataset(f1_dataset)
        ts(f"  {dataset.n_obs:,} cells")

    ts("Loading L/R map (multiome merger)...")
    lr_map = load_lr_map(lr_tsv)
    ts(f"  {len(lr_map):,} (genotype, library) entries")

    ts("Loading FIXED e1b stage map...")
    stage_map = load_e1b_stage_map(e1b_tsv)
    ts(f"  {len(stage_map):,} genotype -> stage entries")

    # ── Build cell × region peak matrix ──────────────────────────────────
    peak_list = [f"{c}:{s}-{e}" for c, s, e in regions.values()]
    region_names = list(regions.keys())
    ts(f"Building cell × region matrix for {len(peak_list)} regions: "
       f"{peak_list}")

    # make_peak_matrix iterates fragments — ~20 min cold pass.  Cache the
    # output h5ad so reruns finish in seconds.  F5_FORCE_REBUILD=1 to
    # bypass.
    if cached_mode:
        ts(f"Reading cached peak matrix {peak_matrix_tmp} (read-only)")
        peak_matrix = snap.read(peak_matrix_tmp, backed="r")
        check_cached_matrix([str(b) for b in peak_matrix.obs_names], obs_barcodes,
                            [str(v) for v in peak_matrix.var_names], peak_list,
                            peak_matrix_tmp)
        ts("  cache barcodes == obs TSV rows (positional); regions == requested")
    elif (os.path.exists(peak_matrix_tmp)
            and not os.environ.get("F5_FORCE_REBUILD")):
        ts(f"Reusing cached peak matrix at {peak_matrix_tmp}")
        peak_matrix = snap.read(peak_matrix_tmp)
    else:
        if os.path.exists(peak_matrix_tmp):
            os.remove(peak_matrix_tmp)
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, use_rep=peak_list, file=peak_matrix_tmp,
        )
    ts(f"  Peak matrix: {peak_matrix.shape[0]:,} cells × "
       f"{peak_matrix.shape[1]:,} regions")

    import scipy.sparse as sp
    X_raw = peak_matrix.X[:]
    if sp.issparse(X_raw):
        X = X_raw.toarray()
    else:
        X = np.asarray(X_raw)
    if X.ndim != 2:
        sys.exit(f"ERROR: unexpected peak_matrix.X[:] shape {X.shape}")
    if X.shape != (peak_matrix.shape[0], peak_matrix.shape[1]):
        if X.shape == (peak_matrix.shape[1], peak_matrix.shape[0]):
            X = X.T
        else:
            sys.exit(f"ERROR: peak_matrix.X shape {X.shape} != "
                     f"({peak_matrix.shape[0]}, {peak_matrix.shape[1]})")
    ts(f"  Counts: total fragments per region = "
       f"{[int(X[:, k].sum()) for k in range(X.shape[1])]}")

    # ── Pull obs columns (dataset.obs, or the obs TSV in cached mode) ───
    if not cached_mode:
        obs = {c: list(dataset.obs[c]) for c in OBS_COLS}
    genotype = obs["unified_genotype_id"]
    arc_lib = obs["arc_library"]
    stage = obs["d2g_stage"]
    subcluster = obs["d2g_subcluster"]

    # ── Aggregate per (genotype, d2l_library, side, stage_anchor) ────────
    drops: Counter = Counter()
    groups: dict[tuple[str, str, str, str, str], list[int]] = {}
    for i, (g, alib, s_raw, sc) in enumerate(
            zip(genotype, arc_lib, stage, subcluster)):
        if g is None or alib is None:
            drops["uncoupled_or_unassigned"] += 1
            continue
        # Re-derive stage from FIXED e1b (NOT cell.d2g_stage).
        fixed_stage = stage_map.get(str(g))
        if fixed_stage is None:
            drops["genotype_not_in_e1b"] += 1
            continue
        if fixed_stage not in include_stages:
            drops["stage_not_included"] += 1
            continue
        d2l_lib = arc_to_d2l_library.get(str(alib))
        if d2l_lib is None:
            drops["unknown_arc_library"] += 1
            continue
        # Strict lookups: an unmapped stage or an anchor missing from POOLS would
        # otherwise empty the pool and file every cell as not_in_inference_pool,
        # indistinguishable from a real exclusion.
        if fixed_stage not in stage_to_anchor:
            raise KeyError(f"stage {fixed_stage!r} (genotype {g}) has no entry in "
                           f"STAGE_TO_ANCHOR ({sorted(stage_to_anchor)})")
        anchor = stage_to_anchor[fixed_stage]
        if anchor not in pools or "inference" not in pools[anchor]:
            raise KeyError(f"stage anchor {anchor!r} has no POOLS inference list "
                           f"(POOLS keys: {sorted(pools)})")
        pool_defs = pools[anchor]["inference"]
        if sc not in pool_defs:
            drops["not_in_inference_pool"] += 1
            continue
        lr = lr_map.get((str(g), d2l_lib))
        if lr is None:
            # no usable side: weak / pair-less rows (excluded by the logged
            # firmness rule in load_lr_map) or no sidecar row at all
            drops["no_usable_side"] += 1
            continue
        side, sample_lib = lr
        groups.setdefault(
            (str(g), sample_lib, side, anchor, fixed_stage), []
        ).append(i)
    ts(f"  Drops: {dict(drops)}")
    ts(f"  Kept groups: {len(groups)}")

    rows: list[dict] = []
    for (g, lib, side, anchor, st_raw), idxs in groups.items():
        n_cells = len(idxs)
        sums = X[idxs].sum(axis=0)
        for k, region in enumerate(region_names):
            rows.append({
                "genotype": g,
                "library_d2l": lib,
                "side": side,
                "stage": anchor,
                "stage_raw": st_raw,
                "region": region,
                "n_cells": n_cells,
                "n_frags": int(sums[k]),
                "frags_per_cell": float(sums[k]) / n_cells,
            })
    out_df = pd.DataFrame(rows)
    out_df = out_df.sort_values(
        ["region", "stage", "side", "genotype", "library_d2l"]
    ).reset_index(drop=True)

    out_df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")
    if not os.path.exists(out_tsv):
        sys.exit(f"ERROR: to_csv reported success but {out_tsv} doesn't exist")
    sz = os.path.getsize(out_tsv)
    ts(f"  Saved: {out_tsv}  ({len(out_df):,} rows, {sz:,} bytes)")

    summary = (out_df.groupby(["stage", "side", "region"])
               .agg(n_samples=("genotype", "nunique"),
                    total_cells=("n_cells", "sum"),
                    total_frags=("n_frags", "sum"),
                    median_frags_per_cell=("frags_per_cell", "median"))
               .reset_index())
    ts("=== Per-stage summary ===")
    for _, r in summary.iterrows():
        ts(f"  {r['stage']:<8} {r['side']} {r['region']:<18}: "
           f"{int(r['n_samples'])} samples, "
           f"{int(r['total_cells']):,} cells, "
           f"{int(r['total_frags']):,} frags, "
           f"median frags/cell={r['median_frags_per_cell']:.4f}")

    # NOTE: do NOT delete PEAK_MATRIX_TMP — it's the cache for the next
    # run.  F5_FORCE_REBUILD=1 to bypass on a given run.
    ts("=== Done ===")
