#!/usr/bin/env python3
"""Phase H Step 1: per-(genotype × library_m × side × subcluster × stage)
gene-level UMI pseudobulks.

Counterpart to ``g1_pseudobulk_atac.py`` for the snRNA-seq DEG track.
Aggregates the master clustered h5ad into one count vector per pseudobulk;
those pseudobulks feed h2 DESeq2 (R) for L vs R contrasts per
(stage × subcluster) scope.

Pseudobulk key
--------------
Tuple ``(unified_genotype_id, library_m, inferred_side, d2g_subcluster,
stage_anchor)`` with library-collapsed technical replicates:

  * GEX03-1A + GEX03-1B → GEX03-1
  * PILOTGEX4A + PILOTGEX4B → PILOTGEX4

Size floor MIN_CELLS = 30 (project standard; matches f4/d2l/g1 discipline).
The h2 R script applies the floor — h1 emits all (genotype × library_m ×
side × subcluster × stage) buckets including under-floor ones; metadata
carries n_cells per bucket so h2 can filter at multiple thresholds for
sensitivity analysis.

Stage anchor (FIXED e1b)
------------------------
Stage anchor is re-derived per cell from the FIXED e1b sidecar
(``figures/data/e1b_stage_assignments.tsv``), not from
``cell.d2g_stage`` which carries the OLD pre-2026-04-26 e1b values.
Cells where OLD ≠ FIXED stage are kept — they carry the OLD-stage's
d2g_subcluster annotation (since d2g hasn't been re-run with FIXED
e1b yet) but are tagged in metadata as ``stage_mismatch=True`` so h2
can optionally filter them.

Side inference
--------------
Per pseudobulk, side is taken from ``lr_inference_multiome.tsv``
(chromatin-canonical at HH17/16, snRNA at other stages) joined on
(unified_genotype_id, library_m).  Firmness tiers accepted:
``direct``, ``placement``, ``singleton_rank`` (and the chromatin_*
prefixed equivalents — normalized via prefix strip).  ``weak``
firmness rows have no clear L/R call and are excluded.

Inference-pool restriction
--------------------------
Cells must have d2g_subcluster ∈ POOLS[stage_anchor]["inference"]
(per ``figures/_shared/subcluster_map.py``).  Excludes Skeletogenic LPM /
Sclerotome / Dermomyotome / paraxial / ambiguous / Contamination lineages.

Input
-----
  phase2/step04_preCB_clustered.h5ad           (master snRNA h5ad)
  figures/data/e1b_stage_assignments.tsv        (FIXED e1b)
  figures/data/lr_inference_multiome.tsv        (post pass-2 d2o output)
  references/d2g_subcluster_annotations.tsv     (population label dictionary)
  figures/_shared/subcluster_map.py             (POOLS)

Output
------
  figures/data/h1_pseudobulk_snrna_counts.tsv.gz   (genes × samples; gzipped TSV)
  figures/data/h1_pseudobulk_snrna_metadata.tsv    (one row per sample)

Run on O2
---------
  sbatch pipeline/run_pipeline.sh --from h1 --to h1
(runs all three h1 pseudobulk variants in one job: LPM + paraxial + excluded.)
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
from _same_side_duplicates import SAME_SIDE_DUP, pool_duplicate_libraries  # noqa: E402


import gzip
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# Resolve the figures/ helper dir relative to this file (sync script
# rsyncs figures/_shared/ to scratch).
sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"),
)
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.lr_firmness import normalize_firmness as check_firmness_alphabet  # noqa: E402

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
# step02_preCB_clean.h5ad = CellBender-corrected RAW integer counts,
# all ~25K genes (NOT HVG-subsetted, NOT log-normalized).  This is what
# DESeq2 expects.  step04_preCB_clustered.h5ad has only 3000 HVGs and
# log-normalized .X — wrong substrate for pseudobulk DEG.
H5AD = f"{PHASE2}/step02_preCB_clean.h5ad"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
LR_TSV = f"{SCRATCH}/figures/data/lr_inference_multiome.tsv"
ANNOT_TSV = f"{SCRATCH}/references/d2g_subcluster_annotations.tsv"

# Per-stage ccreg d2g LPM subclustering TSVs, matched to the annotation TSV.
# Do not join default-variant Leiden IDs to ccreg annotations.
# Master h5ad doesn't carry
# d2g_stage / d2g_subcluster columns — those live in these per-stage
# sidecars keyed by (genotype, library, barcode).  d2g uses the LPM
# subset only, so non-LPM cells are absent from these TSVs and get
# excluded from h1 by design (consistent with g1 ATAC discipline).
D2G_TSVS = {
    "HH13":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}

# Canonical leiden resolution per stage (per references/d2g_subcluster_annotations.tsv)
LEIDEN_COL_PER_STAGE = {
    "HH13":    "leiden_0.2",
    "HH17/16": "leiden_0.3",
    "HH21":    "leiden_0.3",
    "HH25":    "leiden_0.3",
}

OUT_COUNTS = f"{SCRATCH}/figures/data/h1_pseudobulk_snrna_counts.tsv.gz"
OUT_META = f"{SCRATCH}/figures/data/h1_pseudobulk_snrna_metadata.tsv"

# Library tech-rep collapse — matches g1/f5 conventions
LIBRARY_COLLAPSE = {
    "GEX03-1A":   "GEX03-1",
    "GEX03-1B":   "GEX03-1",
    "PILOTGEX4A": "PILOTGEX4",
    "PILOTGEX4B": "PILOTGEX4",
}

# Stage canonicalization — handle HH16/HH17 → HH17/16 anchor.  All
# other stages are passthrough.
STAGE_TO_ANCHOR = {
    "HH13":    "HH13",
    "HH16":    "HH17/16",
    "HH17":    "HH17/16",
    "HH17/16": "HH17/16",
    "HH21":    "HH21",
    "HH25":    "HH25",
}

# Firmness tiers accepted (post chromatin_ prefix strip).
# same_side_dup = lane-loading-record split tube (one half in two libraries):
# a usable design-fixed side; its two libraries are pooled into ONE
# pseudobulk below (pool_duplicate_libraries), like the tech-rep collapse.
ACCEPTED_FIRMNESS = {"direct", "placement", "singleton_rank", SAME_SIDE_DUP}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def collapse_library(raw: str) -> str:
    if pd.isna(raw):
        return raw
    return LIBRARY_COLLAPSE.get(raw, raw)


def normalize_firmness(firm: str) -> str:
    if pd.isna(firm):
        return firm
    return firm.removeprefix("chromatin_")


def load_e1b_stage_map(tsv_path: str) -> dict[str, str]:
    """genotype_id -> inferred_stage (FIXED e1b)."""
    df = pd.read_csv(tsv_path, sep="\t")
    return dict(zip(df["genotype_id"].astype(str),
                    df["inferred_stage"].astype(str)))


def load_lr_map(tsv_path: str) -> dict[tuple[str, str], tuple[str, str]]:
    """(genotype, library_m) -> (side, firmness_normalized)."""
    df = pd.read_csv(tsv_path, sep="\t")
    # unknown firmness values raise, naming them (shared alphabet incl. same_side_dup)
    check_firmness_alphabet(df["firmness"], os.path.basename(tsv_path))
    lr = {}
    for r in df.itertuples():
        side = getattr(r, "inferred_side", None)
        firm = getattr(r, "firmness", None)
        if pd.isna(side) or pd.isna(firm):
            continue
        firm_norm = normalize_firmness(str(firm))
        lr[(str(r.genotype), str(r.library_m))] = (str(side), firm_norm)
    return lr


def load_d2g_subcluster_map() -> dict:
    """Build per-cell (genotype, library, barcode) -> (anchor, population)
    map by joining the per-stage d2g LPM TSVs with the population
    annotation in references/d2g_subcluster_annotations.tsv.

    Returns: dict[(str, str, str), (str, str)] keyed on
    (genotype, library, barcode), value (stage_anchor, population).
    """
    annot = pd.read_csv(ANNOT_TSV, sep="\t")
    out: dict[tuple[str, str, str], tuple[str, str]] = {}
    for anchor, tsv_path in D2G_TSVS.items():
        if not os.path.isfile(tsv_path):
            ts(f"  ⚠ skipping anchor {anchor}: missing {tsv_path}")
            continue
        leiden_col = LEIDEN_COL_PER_STAGE[anchor]
        # Pull just the columns we need
        df = pd.read_csv(tsv_path, sep="\t", compression="gzip",
                         usecols=["genotype", "library", "barcode",
                                  leiden_col])
        # Build leiden -> population for THIS anchor + leiden_col
        annot_sub = annot[(annot["stage"] == anchor) &
                          (annot["leiden_col"] == leiden_col)]
        leiden_to_pop = dict(zip(annot_sub["cluster"].astype(int),
                                 annot_sub["population"]))
        df["population"] = df[leiden_col].astype(int).map(leiden_to_pop)
        # The d2g TSVs sometimes carry rows whose leiden cluster is
        # absent from the annotation table (rare ambiguous clusters).
        # Such rows are dropped here — population NaN means "no
        # canonical annotation" and they'd fail the inference-pool
        # filter downstream anyway.
        df = df.dropna(subset=["population"])
        ts(f"  d2g[{anchor}]: {len(df):,} cells with annotated population "
           f"(leiden_col={leiden_col})")
        for r in df.itertuples():
            out[(str(r.genotype), str(r.library), str(r.barcode))] = (
                anchor, str(r.population))
    return out


def main() -> None:
    ts("=== Phase H Step 1: snRNA-seq gene-level pseudobulk ===")
    ts(f"  H5AD: {H5AD}")
    ts(f"  E1B:  {E1B_TSV}")
    ts(f"  LR:   {LR_TSV}")

    # ── 1. Load h5ad ────────────────────────────────────────────────────
    if not os.path.isfile(H5AD):
        sys.exit(f"ERROR: {H5AD} missing")
    ts("Loading master h5ad ...")
    adata = sc.read_h5ad(H5AD)
    ts(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Required obs columns.  The master snRNA h5ad uses "library".
    # Stage and subcluster come from the per-stage d2g LPM sidecars
    # (NOT obs.d2g_stage / obs.d2g_subcluster — those don't exist in
    # the master h5ad), joined by (genotype, library, barcode).
    for col in ("unified_genotype_id", "library", "scDblFinder_class"):
        if col not in adata.obs.columns:
            sys.exit(f"ERROR: obs missing column '{col}'.  "
                     f"Available: {sorted(adata.obs.columns.tolist())}")

    # ── 2. Drop doublets ────────────────────────────────────────────────
    n_pre = adata.n_obs
    keep = adata.obs["scDblFinder_class"].astype(str) != "doublet"
    adata = adata[keep].copy()
    ts(f"  Dropped {n_pre - adata.n_obs:,} doublets; "
       f"{adata.n_obs:,} singlets remaining")

    # ── 3. Load FIXED e1b stage map + LR map ────────────────────────────
    ts("Loading FIXED e1b stage map ...")
    stage_map = load_e1b_stage_map(E1B_TSV)
    ts(f"  {len(stage_map):,} genotype → inferred_stage entries")

    ts("Loading multiome LR map ...")
    lr_map = load_lr_map(LR_TSV)
    ts(f"  {len(lr_map):,} (genotype, library_m) → (side, firmness) entries")

    ts("Loading d2g per-stage subcluster map (LPM-only) ...")
    d2g_map = load_d2g_subcluster_map()
    ts(f"  {len(d2g_map):,} (genotype, library, barcode) -> "
       f"(d2g_anchor, population) entries")

    # ── 4. Build per-cell metadata: stage_anchor (FIXED), library_m, side,
    #       firmness, subcluster, in_pool flag ──────────────────────────
    ts("Building per-cell metadata ...")
    obs = adata.obs.copy()

    # The h5ad obs index is the cell barcode (with potential library
    # suffix).  Verify and extract the bare barcode for d2g lookup.
    # Standard scanpy convention: index is the obs name; need to also
    # check whether barcodes are stored in a column.
    if "barcode" in obs.columns:
        obs["_bc"] = obs["barcode"].astype(str)
    else:
        obs["_bc"] = obs.index.astype(str)
    obs["library_m"] = obs["library"].map(collapse_library)

    # FIXED stage anchor (per genotype, not per cell — all cells of one
    # genotype share the same e1b inferred_stage).
    obs["fixed_stage"] = obs["unified_genotype_id"].astype(str).map(stage_map)
    obs["stage_anchor"] = obs["fixed_stage"].map(
        lambda s: STAGE_TO_ANCHOR[s] if pd.notna(s) else s)
    n_no_e1b = obs["stage_anchor"].isna().sum()
    if n_no_e1b:
        ts(f"  ⚠ {n_no_e1b:,} cells with no FIXED e1b stage (genotype "
           f"missing from sidecar).  Tagging stage_anchor=NaN — these "
           f"cells will be excluded.")

    # Side + firmness from multiome LR map
    obs["_lr"] = obs.apply(
        lambda r: lr_map.get((str(r["unified_genotype_id"]),
                              str(r["library_m"]))), axis=1)
    obs["inferred_side"] = obs["_lr"].map(
        lambda x: x[0] if x is not None else np.nan)
    obs["firmness_norm"] = obs["_lr"].map(
        lambda x: x[1] if x is not None else np.nan)
    obs.drop(columns=["_lr"], inplace=True)
    n_no_lr = obs["inferred_side"].isna().sum()
    ts(f"  {n_no_lr:,} cells with no multiome LR entry")

    # d2g subcluster lookup: returns (d2g_anchor, population) per cell.
    # The d2g_anchor is the OLD pre-2026-04-26-e1b-fix stage; the
    # population is the canonical d2g annotation we want.  We use the
    # FIXED stage (above) but the OLD subcluster identity (since d2g
    # hasn't been re-run with the FIXED e1b yet — same convention as
    # f6_chromatin_metrics_hh1716.py).
    obs["_d2g"] = obs.apply(
        lambda r: d2g_map.get((str(r["unified_genotype_id"]),
                               str(r["library"]),
                               str(r["_bc"]))), axis=1)
    obs["d2g_anchor"] = obs["_d2g"].map(
        lambda x: x[0] if x is not None else np.nan)
    obs["d2g_subcluster"] = obs["_d2g"].map(
        lambda x: x[1] if x is not None else np.nan)
    obs.drop(columns=["_d2g"], inplace=True)
    n_in_d2g = obs["d2g_subcluster"].notna().sum()
    ts(f"  Cells with d2g LPM annotation: {n_in_d2g:,} "
       f"(LPM-only — non-LPM cells absent by design)")

    # OLD vs FIXED stage mismatch flag (informational)
    obs["stage_mismatch"] = (obs["stage_anchor"].astype(str) !=
                             obs["d2g_anchor"].astype(str)) & \
                            obs["d2g_anchor"].notna()
    n_mismatch = obs["stage_mismatch"].sum()
    ts(f"  Stage mismatch (OLD d2g_anchor ≠ FIXED e1b): {n_mismatch:,} cells")
    if n_mismatch > 0:
        breakdown = (obs[obs["stage_mismatch"]]
                     .groupby(["unified_genotype_id"])
                     .size().sort_values(ascending=False).head(10))
        ts(f"    top mismatched genotypes (cell counts):\n      "
           f"{breakdown.to_dict()}")

    # ── 5. Inference-pool filter ────────────────────────────────────────
    ts("Applying inference-pool filter (POOLS) ...")

    def in_pool(row) -> bool:
        anchor = row["stage_anchor"]
        sc_name = row["d2g_subcluster"]
        if pd.isna(anchor) or pd.isna(sc_name):
            return False
        pool = POOLS[anchor]["inference"]
        return sc_name in pool

    obs["in_pool"] = obs.apply(in_pool, axis=1)
    ts(f"  In-pool cells: {int(obs['in_pool'].sum()):,} "
       f"(of {len(obs):,} total = {100*obs['in_pool'].mean():.1f}%)")

    # ── 6. Apply firmness gate ──────────────────────────────────────────
    obs["pass_firm"] = obs["firmness_norm"].isin(ACCEPTED_FIRMNESS)
    ts(f"  Firmness pass ({ACCEPTED_FIRMNESS}): "
       f"{int(obs['pass_firm'].sum()):,} cells")

    # ── 7. Final cell-pass mask ─────────────────────────────────────────
    keep_cell = (obs["in_pool"] &
                 obs["pass_firm"] &
                 obs["inferred_side"].isin(["L", "R"]))
    n_keep = int(keep_cell.sum())
    ts(f"  Cells entering pseudobulks: {n_keep:,}")
    if n_keep == 0:
        sys.exit("ERROR: no cells passed all gates.  Check inputs.")

    obs_keep = obs[keep_cell].copy()
    obs_keep = pool_duplicate_libraries(obs_keep, LR_TSV, log=ts)

    # ── 8. Aggregate pseudobulks ────────────────────────────────────────
    ts("Aggregating pseudobulk counts per "
       "(genotype × library_m × side × subcluster × stage) ...")

    bucket_keys = ["unified_genotype_id", "library_m", "inferred_side",
                   "d2g_subcluster", "stage_anchor"]
    obs_keep["_bucket"] = obs_keep[bucket_keys].astype(str).agg("|".join, axis=1)

    # Index alignment between adata.X (sparse) and obs_keep
    X = adata[obs_keep.index].X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    # Sum per bucket via group-wise sparse multiplication
    buckets = obs_keep["_bucket"].astype("category")
    bucket_codes = buckets.cat.codes.values
    n_buckets = len(buckets.cat.categories)
    n_genes = adata.n_vars

    ts(f"  {n_buckets} buckets × {n_genes:,} genes")
    incidence = sp.csr_matrix(
        (np.ones(len(bucket_codes), dtype=np.float32),
         (bucket_codes, np.arange(len(bucket_codes)))),
        shape=(n_buckets, len(bucket_codes)))
    bulk = (incidence @ X).toarray().astype(np.int64)  # [buckets × genes]

    # ── 9. Write counts (genes × buckets) gzipped TSV ───────────────────
    ts(f"Writing counts to {OUT_COUNTS} ...")
    counts_df = pd.DataFrame(bulk.T,
                             index=adata.var_names,
                             columns=list(buckets.cat.categories))
    counts_df.index.name = "gene_id"
    counts_df.to_csv(OUT_COUNTS, sep="\t", compression="gzip",
                     float_format="%d")
    ts(f"  Wrote {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples "
       f"({os.path.getsize(OUT_COUNTS):,} bytes)")

    # ── 10. Write metadata (one row per bucket) ─────────────────────────
    ts(f"Writing metadata to {OUT_META} ...")
    meta_rows = []
    for bucket_id in buckets.cat.categories:
        sub = obs_keep[obs_keep["_bucket"] == bucket_id]
        first = sub.iloc[0]
        n_cells = len(sub)
        n_mismatch_b = int(sub["stage_mismatch"].sum())
        meta_rows.append({
            "sample_id": bucket_id,
            "genotype": first["unified_genotype_id"],
            "library_m": first["library_m"],
            "side": first["inferred_side"],
            "firmness_norm": first["firmness_norm"],
            "subcluster": first["d2g_subcluster"],
            "stage_anchor": first["stage_anchor"],
            "n_cells": n_cells,
            "n_stage_mismatch_cells": n_mismatch_b,
        })
    meta_df = pd.DataFrame(meta_rows).sort_values(
        ["stage_anchor", "subcluster", "side", "genotype", "library_m"]
    ).reset_index(drop=True)
    meta_df.to_csv(OUT_META, sep="\t", index=False)
    ts(f"  Wrote {len(meta_df)} sample rows")

    # Per-(stage, subcluster) summary at MIN=30
    print("\nPer-(stage, subcluster) summary at MIN_CELLS=30:")
    summary = (meta_df[meta_df["n_cells"] >= 30]
               .groupby(["stage_anchor", "subcluster", "side"])
               .agg(n_buckets=("sample_id", "count"),
                    total_cells=("n_cells", "sum"))
               .reset_index())
    print(summary.to_string(index=False))

    # Pool design (collapse subcluster) summary at MIN=30
    print("\nPool design (genotype × library_m × side × stage):")
    pool_grid = (obs_keep.groupby(["unified_genotype_id", "library_m",
                                    "inferred_side", "stage_anchor"])
                          .size().reset_index(name="n_cells"))
    pool_summary = (pool_grid[pool_grid["n_cells"] >= 30]
                    .groupby(["stage_anchor", "inferred_side"])
                    .agg(n_buckets=("n_cells", "count"),
                         total_cells=("n_cells", "sum"))
                    .reset_index())
    print(pool_summary.to_string(index=False))

    ts("=== Done ===")


if __name__ == "__main__":
    main()
