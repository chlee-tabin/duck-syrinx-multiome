#!/usr/bin/env python3
"""Phase H Step 1 (excluded-pool variant): per-(genotype × library_m × side ×
d2g_subcluster × stage) gene-level UMI pseudobulks restricted to the
**excluded** pool subclusters (per ``figures/_shared/subcluster_map.py``
POOLS[stage]["excluded"]).

Sister to:
  * ``h1_pseudobulk_snrna.py``  (LPM inference pool — the production substrate)

This variant runs L/R DEG on the **excluded** subclusters carried inside
the d2g LPM-restricted subset — i.e. ``Venous-pole spl.``, ``Mesothelium``,
``Skeletogenic LPM`` — to investigate whether L/R DEG signal originally seen in
the OLD "Splanchnic" catchall (lumping pre-2026-05-04 ``Ventral spl.`` +
``Venous-pole spl.``) is in fact specific to the now-excluded
``Venous-pole spl.`` (off-axis cardiac-pole introgression) cells.

Motivation
----------
Pre-2026-05-04 the "Splanchnic" label conflated:
  - the laterality-relevant ventral splanchnic LPM (now ``Ventral spl.``,
    kept in the inference pool) and
  - off-axis venous-pole splanchnic introgression (now ``Venous-pole spl.``,
    excluded as off-axis).

LHX9 came up Yan-direction-OPPOSITE in pre-merge "Splanchnic HH21"
(lfc=−5.27, padj=9e-8).  After the rename + reclassification, h1
restricts to the inference pool only and LHX9 effectively disappears at
HH21 in inference-pool scopes.  Hypothesis: post-split, LHX9 R>L is
cleanly localized to ``Venous-pole spl.`` and not present in
``Mesothelium`` or ``Skeletogenic LPM``.

Pseudobulk key
--------------
Same as h1_pseudobulk_snrna.py: ``(unified_genotype_id, library_m,
inferred_side, d2g_subcluster, stage_anchor)`` with library-collapsed
technical replicates (GEX03-1A/B → GEX03-1; PILOTGEX4A/B → PILOTGEX4).

Cluster substrate
-----------------
Same as h1_pseudobulk_snrna.py — uses the d2g per-stage LPM subclustering
TSVs (``figures/data/d2g_lpm_subclustering_HH*_ccreg.tsv.gz``) joined with
``references/d2g_subcluster_annotations.tsv`` (the d2g population label
SoT).  Cells in the master h5ad whose (genotype, library, barcode) is
absent from these d2g sidecars are non-LPM (filtered out by d2g's LPM
selection) and excluded by design.

Side inheritance
----------------
Identical to h1_pseudobulk_snrna.py: side is
inherited from the per-(genotype × library_m) sidecar
``figures/data/lr_inference_multiome.tsv``.  No d2l rerun on excluded
cells — each library half = one embryo side, so excluded-pool cells
share the side assignment of the ``library_m`` they originate from.

Stage anchor (FIXED e1b)
------------------------
Identical to h1_pseudobulk_snrna.py — re-derived per cell from the FIXED
e1b sidecar (``figures/data/e1b_stage_assignments.tsv``); cells whose
OLD d2g_anchor differs from the FIXED stage are kept and tagged
``stage_mismatch=True``.

Excluded-pool filter
--------------------
Cells must have ``d2g_subcluster ∈ POOLS[stage_anchor]["excluded"]``
(excludes inference + bilateral pools, plus `ambiguous` / paraxial
labels which never appear in the d2g LPM subset anyway).  At each
stage the concrete excluded labels are typically a subset of:
``Skeletogenic LPM``, ``Mesothelium``, ``Venous-pole spl.``, ``Contamination``
(Sclerotome / Dermomyotome / ambiguous are biologically absent from the d2g
LPM subset for the most part — they're flagged optional in
``EXCLUDED_OPTIONAL_PER_STAGE`` upstream).

Inputs
------
  phase2/step02_preCB_clean.h5ad                (master snRNA raw counts)
  figures/data/e1b_stage_assignments.tsv         (FIXED e1b)
  figures/data/lr_inference_multiome.tsv         (post pass-2 d2o output)
  references/d2g_subcluster_annotations.tsv      (d2g population SoT)
  figures/data/d2g_lpm_subclustering_HH*_ccreg.tsv.gz  (per-stage d2g LPM cells)
  figures/_shared/subcluster_map.py              (POOLS)

Outputs
-------
  figures/data/h1_excluded_pseudobulk_counts.tsv.gz   (genes × samples)
  figures/data/h1_excluded_pseudobulk_metadata.tsv    (one row per sample)

Run on O2
---------
  sbatch pipeline/run_pipeline.sh --from h1 --to h1
(runs all three h1 pseudobulk variants in one job; this is the excluded one.)
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
from _same_side_duplicates import SAME_SIDE_DUP, pool_duplicate_libraries  # noqa: E402


import os
import sys
import time
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
H5AD = f"{PHASE2}/step02_preCB_clean.h5ad"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
LR_TSV = f"{SCRATCH}/figures/data/lr_inference_multiome.tsv"
ANNOT_TSV = f"{SCRATCH}/references/d2g_subcluster_annotations.tsv"

# ccreg sidecars match the ccreg-derived annotation TSV.
D2G_TSVS = {
    "HH13":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}

LEIDEN_COL_PER_STAGE = {
    "HH13":    "leiden_0.2",
    "HH17/16": "leiden_0.3",
    "HH21":    "leiden_0.3",
    "HH25":    "leiden_0.3",
}

OUT_COUNTS = f"{SCRATCH}/figures/data/h1_excluded_pseudobulk_counts.tsv.gz"
OUT_META = f"{SCRATCH}/figures/data/h1_excluded_pseudobulk_metadata.tsv"

LIBRARY_COLLAPSE = {
    "GEX03-1A":   "GEX03-1",
    "GEX03-1B":   "GEX03-1",
    "PILOTGEX4A": "PILOTGEX4",
    "PILOTGEX4B": "PILOTGEX4",
}

STAGE_TO_ANCHOR = {
    "HH13":    "HH13",
    "HH16":    "HH17/16",
    "HH17":    "HH17/16",
    "HH17/16": "HH17/16",
    "HH21":    "HH21",
    "HH25":    "HH25",
}

# same_side_dup = lane-loading-record split tube (one half in two libraries):
# a usable design-fixed side; its two libraries are pooled into ONE
# pseudobulk below (pool_duplicate_libraries), like the tech-rep collapse.
ACCEPTED_FIRMNESS = {"direct", "placement", "singleton_rank", SAME_SIDE_DUP}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def collapse_library(raw):
    if pd.isna(raw):
        return raw
    return LIBRARY_COLLAPSE.get(raw, raw)


def normalize_firmness(firm):
    if pd.isna(firm):
        return firm
    return str(firm).removeprefix("chromatin_")


def load_e1b_stage_map(tsv_path: str) -> dict[str, str]:
    df = pd.read_csv(tsv_path, sep="\t")
    return dict(zip(df["genotype_id"].astype(str),
                    df["inferred_stage"].astype(str)))


def load_lr_map(tsv_path: str) -> dict[tuple[str, str], tuple[str, str]]:
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
    """(genotype, library, barcode) -> (anchor, population), restricted to
    the d2g LPM subset.  Identical to h1_pseudobulk_snrna's loader."""
    annot = pd.read_csv(ANNOT_TSV, sep="\t")
    out: dict[tuple[str, str, str], tuple[str, str]] = {}
    for anchor, tsv_path in D2G_TSVS.items():
        if not os.path.isfile(tsv_path):
            ts(f"  ⚠ skipping anchor {anchor}: missing {tsv_path}")
            continue
        leiden_col = LEIDEN_COL_PER_STAGE[anchor]
        df = pd.read_csv(tsv_path, sep="\t", compression="gzip",
                         usecols=["genotype", "library", "barcode",
                                  leiden_col])
        annot_sub = annot[(annot["stage"] == anchor) &
                          (annot["leiden_col"] == leiden_col)]
        leiden_to_pop = dict(zip(annot_sub["cluster"].astype(int),
                                 annot_sub["population"]))
        df["population"] = df[leiden_col].astype(int).map(leiden_to_pop)
        df = df.dropna(subset=["population"])
        ts(f"  d2g[{anchor}]: {len(df):,} cells with annotated population "
           f"(leiden_col={leiden_col})")
        for r in df.itertuples():
            out[(str(r.genotype), str(r.library), str(r.barcode))] = (
                anchor, str(r.population))
    return out


def main() -> None:
    ts("=== Phase H Step 1 (excluded): snRNA excluded-pool pseudobulk ===")
    ts(f"  H5AD: {H5AD}")
    ts(f"  E1B:  {E1B_TSV}")
    ts(f"  LR:   {LR_TSV}")

    # ── 1. Load h5ad ────────────────────────────────────────────────────
    if not os.path.isfile(H5AD):
        sys.exit(f"ERROR: {H5AD} missing")
    ts("Loading master h5ad ...")
    adata = sc.read_h5ad(H5AD)
    ts(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

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

    # ── 3. Load FIXED e1b stage map + LR map + d2g subcluster map ───────
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

    # ── 4. Build per-cell metadata ──────────────────────────────────────
    ts("Building per-cell metadata ...")
    obs = adata.obs.copy()
    if "barcode" in obs.columns:
        obs["_bc"] = obs["barcode"].astype(str)
    else:
        obs["_bc"] = obs.index.astype(str)
    obs["library_m"] = obs["library"].map(collapse_library)

    # FIXED stage anchor
    obs["fixed_stage"] = obs["unified_genotype_id"].astype(str).map(stage_map)
    obs["stage_anchor"] = obs["fixed_stage"].map(
        lambda s: STAGE_TO_ANCHOR[s] if pd.notna(s) else s)

    # Side + firmness
    obs["_lr"] = obs.apply(
        lambda r: lr_map.get((str(r["unified_genotype_id"]),
                              str(r["library_m"]))), axis=1)
    obs["inferred_side"] = obs["_lr"].map(
        lambda x: x[0] if x is not None else np.nan)
    obs["firmness_norm"] = obs["_lr"].map(
        lambda x: x[1] if x is not None else np.nan)
    obs.drop(columns=["_lr"], inplace=True)

    # d2g subcluster
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
    ts(f"  Cells with d2g LPM annotation: {n_in_d2g:,}")

    obs["stage_mismatch"] = (obs["stage_anchor"].astype(str) !=
                             obs["d2g_anchor"].astype(str)) & \
                            obs["d2g_anchor"].notna()
    n_mismatch = obs["stage_mismatch"].sum()
    ts(f"  Stage mismatch (OLD d2g_anchor ≠ FIXED e1b): {n_mismatch:,} cells")

    # ── 5. Excluded-pool filter ─────────────────────────────────────────
    ts("Applying EXCLUDED-pool filter (POOLS[stage]['excluded']) ...")

    def in_excluded(row) -> bool:
        anchor = row["stage_anchor"]
        sc_name = row["d2g_subcluster"]
        if pd.isna(anchor) or pd.isna(sc_name):
            return False
        excl = POOLS[anchor]["excluded"]
        return sc_name in excl

    obs["in_excluded"] = obs.apply(in_excluded, axis=1)
    ts(f"  In-excluded-pool cells: {int(obs['in_excluded'].sum()):,} "
       f"(of {len(obs):,} total = {100*obs['in_excluded'].mean():.1f}%)")

    # Pre-filter breakdown by (stage, subcluster)
    excl_breakdown = (obs[obs["in_excluded"]]
                      .groupby(["stage_anchor", "d2g_subcluster"])
                      .size().reset_index(name="n_cells")
                      .sort_values(["stage_anchor", "n_cells"],
                                   ascending=[True, False]))
    ts("  Excluded-pool breakdown (pre side/firmness gate):")
    for _, r in excl_breakdown.iterrows():
        ts(f"    {r['stage_anchor']:8s} {r['d2g_subcluster']:24s} "
           f"{r['n_cells']:6d}")

    # ── 6. Apply firmness + side gates ──────────────────────────────────
    obs["pass_firm"] = obs["firmness_norm"].isin(ACCEPTED_FIRMNESS)
    keep_cell = (obs["in_excluded"] &
                 obs["pass_firm"] &
                 obs["inferred_side"].isin(["L", "R"]))
    n_keep = int(keep_cell.sum())
    ts(f"  Cells entering pseudobulks (excl + firm + side): {n_keep:,}")
    if n_keep == 0:
        sys.exit("ERROR: no cells passed all gates.  Check inputs.")

    obs_keep = obs[keep_cell].copy()
    obs_keep = pool_duplicate_libraries(obs_keep, LR_TSV, log=ts)

    # ── 7. Aggregate pseudobulks ────────────────────────────────────────
    ts("Aggregating pseudobulk counts per "
       "(genotype × library_m × side × subcluster × stage) ...")
    bucket_keys = ["unified_genotype_id", "library_m", "inferred_side",
                   "d2g_subcluster", "stage_anchor"]
    obs_keep["_bucket"] = obs_keep[bucket_keys].astype(str).agg("|".join, axis=1)

    X = adata[obs_keep.index].X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    buckets = obs_keep["_bucket"].astype("category")
    bucket_codes = buckets.cat.codes.values
    n_buckets = len(buckets.cat.categories)
    n_genes = adata.n_vars

    ts(f"  {n_buckets} buckets × {n_genes:,} genes")
    incidence = sp.csr_matrix(
        (np.ones(len(bucket_codes), dtype=np.float32),
         (bucket_codes, np.arange(len(bucket_codes)))),
        shape=(n_buckets, len(bucket_codes)))
    bulk = (incidence @ X).toarray().astype(np.int64)

    # ── 8. Write counts ─────────────────────────────────────────────────
    ts(f"Writing counts to {OUT_COUNTS} ...")
    counts_df = pd.DataFrame(bulk.T,
                             index=adata.var_names,
                             columns=list(buckets.cat.categories))
    counts_df.index.name = "gene_id"
    counts_df.to_csv(OUT_COUNTS, sep="\t", compression="gzip",
                     float_format="%d")
    ts(f"  Wrote {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples "
       f"({os.path.getsize(OUT_COUNTS):,} bytes)")

    # ── 9. Write metadata ──────────────────────────────────────────────
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

    print("\nPer-(stage, subcluster) summary at MIN_CELLS=30:")
    summary = (meta_df[meta_df["n_cells"] >= 30]
               .groupby(["stage_anchor", "subcluster", "side"])
               .agg(n_buckets=("sample_id", "count"),
                    total_cells=("n_cells", "sum"))
               .reset_index())
    print(summary.to_string(index=False))

    ts("=== Done ===")


if __name__ == "__main__":
    main()
