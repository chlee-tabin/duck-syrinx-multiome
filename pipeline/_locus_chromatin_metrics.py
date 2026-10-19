"""Shared backbone for f6 locus chromatin-metric pipelines.

`f6_chromatin_metrics.py` (PITX2) and `f6_cfc1_chromatin_metrics.py`
(CFC1) are near-identical sibling drivers — same SnapATAC2 peak-matrix
loop, same e1b stage re-derivation, same per-(genotype × library ×
stage) aggregation.  They differ only in the locus REGIONS dict, the
output paths, the column friendly-alias names (e.g. `pitx2c_open_frac`
vs `cfc1_open_frac`), and a few cosmetic log strings.  This module
owns the shared logic; each sibling becomes a thin shim that calls
`run_locus_chromatin_metrics(...)` with its own constants.

Output TSV schema (one row per (genotype × library × stage) sample):
    genotype, library, stage, n_cells,
    <friendly aliases for primary `frags_per_cell` per region>,
    <friendly aliases for `n_*_frags` per region>,
    <friendly aliases for binary `open_frac` per region>,
    <friendly aliases for `n_*_open` per region>

Each sibling supplies a `column_aliases` dict that maps the canonical
output names (`<REGION>_frags_per_cell`, `n_<REGION>_frags`,
`<REGION>_open_frac`, `n_<REGION>_open`) to the locus-specific friendly
column names (`pitx2c_frags_per_cell` / `cfc1_frags_per_cell` / ...).
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


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_e1b_stage_map(tsv_path: str) -> dict[str, str]:
    """genotype_id -> inferred_stage from FIXED e1b sidecar."""
    if not os.path.isfile(tsv_path):
        sys.exit(f"ERROR: e1b sidecar missing at {tsv_path}.  "
                 f"Sync figures/data/e1b_stage_assignments.tsv to scratch.")
    df = pd.read_csv(tsv_path, sep="\t")
    return dict(zip(df["genotype_id"].astype(str),
                    df["inferred_stage"].astype(str)))


def run_locus_chromatin_metrics(
    *,
    label: str,
    regions: dict[str, tuple[str, int, int]],
    column_aliases: dict[str, str],
    keep_cols: list[str],
    summary_columns: tuple[str, str],
    f1_dataset: str,
    e1b_tsv: str,
    out_tsv_default: str,
    peak_matrix_tmp: str,
    arc_to_d2l_library: dict[str, str],
    stage_to_anchor: dict[str, str],
    pools,
) -> None:
    """Backbone for f6 locus chromatin-metric drivers.  See module docstring.

    Parameters specific to each sibling:
      label             : log-string label ("PITX2" / "CFC1").
      regions           : {region_name: (chrom, start, end)}.
      column_aliases    : map of canonical column name → friendly alias.
                           e.g. {"PITX2c_promoter_frags_per_cell":
                                  "pitx2c_frags_per_cell", ...}.
      keep_cols         : output column ordering for the final TSV.
      summary_columns   : 2-tuple of (frags_per_cell column name,
                           open_frac column name) used in the per-stage
                           summary print.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-tsv", default=out_tsv_default,
                        help=f"Output TSV path (default {out_tsv_default})")
    args = parser.parse_args()

    ts(f"=== Phase F Step 6 ({label}): per-(genotype x library x stage) "
       f"chromatin metrics across all stages ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")
    ts(f"  out_tsv = {args.out_tsv}")

    if not os.path.isfile(f1_dataset):
        sys.exit(f"ERROR: {f1_dataset} missing")

    ts("Loading f1 dataset...")
    dataset = snap.read_dataset(f1_dataset)
    ts(f"  {dataset.n_obs:,} cells")

    ts("Loading FIXED e1b stage map...")
    stage_map = load_e1b_stage_map(e1b_tsv)
    ts(f"  {len(stage_map):,} genotype -> stage entries")

    # --- Re-derive cell stage_anchor from FIXED e1b ----------------------
    ts("Re-deriving stage_anchor per cell from FIXED e1b...")
    genotype = list(dataset.obs["unified_genotype_id"])
    arc_lib = list(dataset.obs["arc_library"])
    subcluster = list(dataset.obs["d2g_subcluster"])

    inference_pool_by_anchor = {
        anchor: set(pools[anchor]["inference"]) for anchor in pools
    }

    drops: Counter = Counter()
    cell_groups: dict[tuple[str, str, str], list[int]] = {}
    for i, (g, alib, sc) in enumerate(zip(genotype, arc_lib, subcluster)):
        if g is None or alib is None:
            drops["uncoupled_or_unassigned"] += 1
            continue
        fixed_stage = stage_map.get(str(g))
        if fixed_stage is None:
            drops["genotype_not_in_e1b"] += 1
            continue
        anchor = stage_to_anchor.get(fixed_stage)
        if anchor is None:
            drops[f"anchor_{fixed_stage}"] += 1
            continue
        inference_pool = inference_pool_by_anchor.get(anchor, set())
        if sc not in inference_pool:
            drops[f"pool_{sc}"] += 1
            continue
        d2l_lib = arc_to_d2l_library.get(str(alib))
        if d2l_lib is None:
            drops[f"unknown_arc_library_{alib}"] += 1
            continue
        cell_groups.setdefault((str(g), d2l_lib, anchor), []).append(i)

    ts(f"  Drops: {dict(drops)}")
    ts(f"  Kept (genotype, library, stage_anchor) groups: {len(cell_groups)}")
    if not cell_groups:
        sys.exit("ERROR: no cells passed the anchor + inference "
                 "pool filter.  Check fixed e1b output and POOLS definition.")
    n_cells_total = sum(len(v) for v in cell_groups.values())
    ts(f"  Total cells in pool: {n_cells_total:,}")

    # --- Build cell x region peak matrix ---------------------------------
    peak_list = [f"{c}:{s}-{e}" for c, s, e in regions.values()]
    region_names = list(regions.keys())
    ts(f"Building cell x region matrix for {len(peak_list)} regions: "
       f"{peak_list}")

    if (os.path.exists(peak_matrix_tmp)
            and not os.environ.get("F6_FORCE_REBUILD")):
        ts(f"Reusing cached peak matrix at {peak_matrix_tmp}")
        peak_matrix = snap.read(peak_matrix_tmp)
    else:
        if os.path.exists(peak_matrix_tmp):
            os.remove(peak_matrix_tmp)
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, use_rep=peak_list, file=peak_matrix_tmp,
        )
    ts(f"  Peak matrix: {peak_matrix.shape[0]:,} cells x "
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

    region_col_idx = {name: i for i, name in enumerate(region_names)}
    ts(f"  Region -> column index (use_rep order): {region_col_idx}")

    # --- Per-(genotype, library, stage) pseudobulk aggregation ------------
    rows: list[dict] = []
    for (g, lib, stage), idxs in cell_groups.items():
        n_cells = len(idxs)
        record: dict = {
            "genotype":  g,
            "library":   lib,
            "stage":     stage,
            "n_cells":   n_cells,
        }
        for region_name in region_names:
            col = region_col_idx[region_name]
            slice_col = X[idxs, col]
            n_total_frags = int(slice_col.sum())
            n_open = int((slice_col > 0).sum())
            record[f"n_{region_name}_frags"] = n_total_frags
            record[f"{region_name}_frags_per_cell"] = n_total_frags / n_cells
            record[f"n_{region_name}_open"] = n_open
            record[f"{region_name}_open_frac"] = n_open / n_cells
        rows.append(record)

    out_df = pd.DataFrame(rows)
    # Apply friendly-alias renames (sibling-specific).  Each sibling
    # provides a column_aliases dict mapping canonical column names to
    # locus-specific friendly aliases (e.g. PITX2c_promoter_frags_per_cell
    # → pitx2c_frags_per_cell, CFC1_promoter_frags_per_cell →
    # cfc1_frags_per_cell).  CFC1 columns include "cfc1_genebody_*"
    # namespacing to avoid shadowing the PITX2 "n_genebody_*" columns
    # if anyone naive-concats both files.
    for canonical, friendly in column_aliases.items():
        out_df[friendly] = out_df[canonical]

    out_df = (out_df[keep_cols]
              .sort_values(["stage", "genotype", "library"])
              .reset_index(drop=True))

    out_df.to_csv(args.out_tsv, sep="\t", index=False, float_format="%.6f")
    if not os.path.exists(args.out_tsv):
        sys.exit(f"ERROR: to_csv reported success but {args.out_tsv} "
                 f"doesn't exist")
    sz = os.path.getsize(args.out_tsv)
    ts(f"  Saved: {args.out_tsv}  ({len(out_df):,} rows, {sz:,} bytes)")

    summary_fpc, summary_open = summary_columns
    ts("Summary by (stage, library):")
    for (stage, lib), grp in out_df.groupby(["stage", "library"]):
        ts(f"  {stage:8s} {lib:14s}  n_genotypes={len(grp):2d}  "
           f"n_cells_total={int(grp['n_cells'].sum()):,}  "
           f"median_{summary_fpc}={grp[summary_fpc].median():.3f}  "
           f"median_{summary_open}={grp[summary_open].median():.3f}")

    ts("Done.")
