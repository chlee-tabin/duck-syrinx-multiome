#!/usr/bin/env python3
"""Phase G Step 1: Pseudobulk peak counts per (genotype × d2g_subcluster × side × stage).

Plan §3.5 granularity.  Aggregates the f2 peak matrix into one count
vector per pseudobulk; those pseudobulks feed g2 DESeq2.

Pseudobulk key
--------------
Tuple ``(unified_genotype_id, d2g_subcluster, inferred_side, d2g_stage)``
with library-collapsed technical replicates:

  * GEX03-1A + GEX03-1B → GEX03-1
  * PILOTGEX4A + PILOTGEX4B → PILOTGEX4
  * ARC03-1A + ARC03-1B → ARC03-1 (ATAC-side)
  * PILOTARC4A + PILOTARC4B → PILOTARC4 (ATAC-side)

Size floor = 50 cells (plan §3.5 primary).  Sensitivity variants at
30 and 100 persisted alongside for methods-robustness supplement.

Side inference (plan §3.1)
--------------------------
Per pseudobulk, side is taken from the d2n chromatin-canonical L/R merge
``figures/data/lr_inference_multiome.tsv`` (LR_TSV): snRNA sides at
HH13/HH21/HH25 + chromatin-canonical d2n sides at HH17/16 — the canonical L/R
axis used by fig03 panel A/E, figS03, and the f4 coverage pileup.  Joined on
``(unified_genotype_id, library_m)`` through the shared firmness rule
(``figures/_shared/lr_firmness.resolve_sides``: an unknown firmness value raises;
excluded rows are logged by name).  Accepted inferred tiers: ``direct``,
``placement`` (``singleton_rank`` and ``weak`` are excluded by rule) with
``pair_constraint_ok``; ``same_side_dup`` rows are kept on their design side, and
both aliquots of a split-tube half sum into one pseudobulk (see ``load_side_lookup``).
``design/embryo_lr_design.tsv`` is already folded into the P1 placement logic,
so we do not re-query design here (d2m is SSOT for genotype → inferred_side).
Override to the legacy snRNA-only ``d2l_lr_inference.tsv`` via G1_LR_TSV.

LPM-only
--------
Pseudobulks are restricted to LPM-lineage cells (those with a
non-NA, non-``ambiguous`` d2g_subcluster).  The paraxial / epithelial
negative-control pools participate in peak calling (f2) but not in the
laterality contrast this script prepares — they are handled in a
separate negative-control comparison script (fig03 panel E).

Input
-----
  phase2/atac/step_f2_peak_matrix.h5ad
  figures/data/lr_inference_multiome.tsv  (d2n-canonical; G1_LR_TSV → d2l)

Output
------
  phase2/atac/da/pseudobulk_counts.tsv.gz        primary (floor=50)
  phase2/atac/da/pseudobulk_counts_floor20.tsv.gz   sensitivity (HH13 depth)
  phase2/atac/da/pseudobulk_counts_floor30.tsv.gz
  phase2/atac/da/pseudobulk_counts_floor100.tsv.gz
  phase2/atac/da/pseudobulk_metadata.tsv         primary metadata
  phase2/atac/da/pseudobulk_summary.tsv          per-subcluster counts

Run on O2
---------
  sbatch --mem=32G -t 0-1:00 -c 4 -p short pipeline/run_phase_fg.sh --from g1 --to g1
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import snapatac2 as snap

# Resolve the figures/ helper dir relative to this file so `GENOTYPE_META`
# (used for embryo_group metadata on pseudobulks) imports cleanly both on
# scratch and in local dev.  The sync script rsyncs figures/_shared/ to
# scratch, so the helper is reachable in either context.
sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"),
)
from _shared.style import GENOTYPE_META  # noqa: E402
from _shared.lr_firmness import resolve_sides  # noqa: E402

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
ATAC_DIR = f"{PHASE2}/atac"
# DA_DIR overridable via env var so alternative-peak-set re-runs (e.g.
# the ArchR chain at `pipeline/run_atac_hh1716_chain_archr.sh`) emit
# outputs at separate paths and don't clobber SnapATAC2-track results.
DA_DIR = os.environ.get("DA_DIR", f"{ATAC_DIR}/da")
os.makedirs(DA_DIR, exist_ok=True)

# F2_H5AD also env-overridable.
F2_H5AD = os.environ.get("F2_H5AD", f"{ATAC_DIR}/step_f2_peak_matrix.h5ad")
# Canonical L/R axis = the d2n chromatin-canonical merge (lr_inference_multiome):
# snRNA sides at HH13/HH21/HH25 + chromatin-canonical (d2n) sides at HH17/16 — the
# SAME L/R axis as fig03 panel A/E + figS03 (all d2n-canonical) and the f4 coverage
# pileup.  Identical (genotype, library_m) firm coverage to the legacy d2l-only
# sidecar; differs ONLY at HH17/16 (6 library-side flips → chromatin-canonical),
# where snRNA inference is weak.  Override to the snRNA-only d2l sidecar via G1_LR_TSV.
LR_TSV = os.environ.get("G1_LR_TSV", f"{SCRATCH}/figures/data/lr_inference_multiome.tsv")

# Plan §3.5
PRIMARY_FLOOR = 50
# floor=20 added 2026-05-12 to boost HH13 INFERENCE_POOL depth for the
# cross-stage UpSet (HH13 had only 1L+2R pseudobulks at floor=50).
SENSITIVITY_FLOORS = (20, 30, 100)

# Plan §3.1: INFERRED firmness tiers kept for primary laterality inference
# (passed to the shared rule, figures/_shared/lr_firmness.resolve_sides, which
# raises on any value outside the alphabet and logs every exclusion by name).
# same_side_dup (lane-loading-record split tube: one half in two libraries) is not
# configurable there: it is always kept on its design side.  Pseudobulks below
# group by genotype x subcluster x side x stage (NOT library), so a duplicate's
# libraries pool into one pseudobulk and, lacking an opposite half, it can only
# enter unpaired arms.
ACCEPTED_FIRMNESS = {"direct", "placement"}

# Plan §3.1: tech-rep collapse (applies to ATAC's gex_library column
# so we can join the L/R sidecar on library_m).
LIBRARY_COLLAPSE = {
    "GEX03-1A":   "GEX03-1",
    "GEX03-1B":   "GEX03-1",
    "PILOTGEX4A": "PILOTGEX4",
    "PILOTGEX4B": "PILOTGEX4",
}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def collapse_library(raw: str) -> str:
    """Fold paired tech-rep libraries; passthrough for everything else."""
    if pd.isna(raw):
        return raw
    return LIBRARY_COLLAPSE.get(raw, raw)


def load_side_lookup() -> pd.DataFrame:
    """Load the L/R sidecar (LR_TSV = d2n-canonical multiome by default) through the
    shared firmness rule: ACCEPTED_FIRMNESS inferred rows with pair_constraint_ok, plus
    same-side duplicates on their design side.  An unknown firmness value raises.

    A same-side duplicate (one half split into two tubes sequenced in two libraries:
    G3_unfiltered_donor0 / donor3 in GEX04 and GEX0S-duck) keeps BOTH aliquots, and
    because pseudobulks group by genotype x subcluster x side x stage (not library)
    the two aliquots sum into one half, as in f4, f5 and h1.

    Returns a DataFrame with columns [genotype, library_m, inferred_side,
    firmness] ready to merge on (unified_genotype_id, library_m).  The merge key
    (genotype, library_m) is unique.
    """
    if not os.path.isfile(LR_TSV):
        sys.exit(f"ERROR: L/R sidecar missing at {LR_TSV}; rerun d2l/d2o first")
    df = pd.read_csv(LR_TSV, sep="\t")
    source = os.path.basename(LR_TSV)
    resolved = resolve_sides(df, ACCEPTED_FIRMNESS, source=source,
                             require_pair_constraint=True, log=ts)
    out = resolved[["genotype", "library_m", "inferred_side", "firmness"]]
    ts(f"  L/R sidecar ({source}): {len(out):,} rows kept out of {len(df):,}")
    return out.reset_index(drop=True)


def pseudobulk_sum(X, row_idx):
    """Column-wise sum over the rows indexed by row_idx.

    Works on sparse and dense backing.  Returns a 1-D numpy array of
    length n_peaks.
    """
    sub = X[row_idx]
    if sp.issparse(sub):
        return np.asarray(sub.sum(axis=0)).ravel()
    return np.asarray(sub.sum(axis=0)).ravel()


def main() -> None:
    ts("=== Phase G Step 1: Pseudobulk "
       "(genotype × d2g_subcluster × side × stage) ===")

    ts(f"Loading L/R side lookup ({os.path.basename(LR_TSV)})...")
    side_df = load_side_lookup()

    if not os.path.isfile(F2_H5AD):
        sys.exit(f"ERROR: f2 peak matrix missing at {F2_H5AD}")

    ts("Loading f2 peak matrix...")
    adata = snap.read(F2_H5AD)
    ts(f"  {adata.n_obs:,} cells × {adata.n_vars:,} peaks")

    # SnapATAC2 2.9 quirk: adata.obs is a PyDataFrameElem (Polars-backed),
    # NOT a pandas DataFrame.  `pd.DataFrame(adata.obs)` raises
    # "DataFrame constructor not properly called!".  Pull columns
    # individually via the iterable protocol (mirrors f1.snap_obs_to_pandas).
    def _snap_obs_to_pandas(obs, columns):
        data = {}
        for col in columns:
            try:
                data[col] = list(obs[col])
            except Exception:  # noqa: BLE001
                continue
        return pd.DataFrame(data)

    obs_df = _snap_obs_to_pandas(adata.obs, [
        "arc_library", "gex_library", "unified_genotype_id",
        "cell_type", "sex", "stage",
        "d2g_stage", "d2g_cluster", "d2g_subcluster", "peak_pool",
    ])
    obs_df["cell_idx"] = np.arange(len(obs_df))

    # If gex_library missing but arc_library present, derive from the
    # canonical ARC → GEX mapping.  f1 propagates arc_library (the
    # AnnDataSet add_key) but not gex_library; f2 inherits the same
    # propagation gap.  Patching in-memory here avoids re-running f1.
    if "gex_library" not in obs_df.columns and "arc_library" in obs_df.columns:
        ARC_TO_GEX = {
            "PILOTARC3":  "PILOTGEX3",
            "PILOTARC4A": "PILOTGEX4A",
            "PILOTARC4B": "PILOTGEX4B",
            "ARC05": "GEX05", "ARC06": "GEX06", "ARC07": "GEX07",
            "ARC08": "GEX08", "ARC03-1A": "GEX03-1A", "ARC03-1B": "GEX03-1B",
            "ARC03-2": "GEX03-2", "ARC04": "GEX04",
            "ARC0S-duck": "GEX0S-duck",
        }
        obs_df["gex_library"] = obs_df["arc_library"].map(ARC_TO_GEX)
        ts(f"  Derived gex_library from arc_library ("
           f"{obs_df['gex_library'].notna().sum():,} mapped)")

    # ── Re-derive d2g_subcluster from the CANONICAL SSOT via the f1 sidecar ──
    # The d2g_subcluster label baked into the f2 peak-matrix obs is written at
    # f1-time from references/d2g_subcluster_annotations.tsv and goes STALE when
    # that TSV is updated without re-running f1 (the 2026-07 Skeletogenic/cl8
    # revamp).  The leiden id needed to re-map is NOT persisted in the peak
    # matrix, but IS in f1's per-cell metadata sidecar
    # (step_f1_coupled_atac_metadata.tsv: arc_library, cb16, d2g_stage,
    # d2g_cluster).  Join that leiden onto our cells by (arc_library, cb16) and
    # re-map leiden -> fresh label via subcluster_map.get_domain_map (f1's own
    # rule), so an annotation change flows into the pseudobulks WITHOUT a
    # fragment-level rebuild.  Fail loud on a coverage gap — never silently drop.
    from _shared.subcluster_map import get_domain_map  # noqa: E402
    META_SIDECAR = os.environ.get(
        "F1_META_TSV", f"{ATAC_DIR}/step_f1_coupled_atac_metadata.tsv")
    if not os.path.isfile(META_SIDECAR):
        sys.exit(f"ERROR: f1 metadata sidecar missing at {META_SIDECAR}; "
                 f"restore from standby (atac_backup/atac/) to re-derive labels.")
    if "arc_library" not in obs_df.columns:
        sys.exit("ERROR: obs missing 'arc_library'; cannot join the leiden "
                 "sidecar by (arc_library, cb16).")

    def _leiden_str(v):
        if pd.isna(v):
            return None
        try:
            return str(int(v))       # int / numpy-int / float-integer / "8"
        except (ValueError, TypeError):
            return str(v)

    # cb16 = the 16 bp barcode (obs_name minus the "-1" GEM suffix).
    obs_df["cb16"] = [str(b).split("-")[0] for b in list(adata.obs_names)]

    _meta = pd.read_csv(
        META_SIDECAR, sep="\t",
        usecols=lambda c: c in {"arc_library", "cb16", "d2g_stage",
                                "d2g_cluster"},
        dtype={"arc_library": str, "cb16": str, "d2g_stage": str},
    )
    _meta = _meta[_meta["d2g_cluster"].notna()
                  & _meta["d2g_stage"].notna()].copy()
    _fresh = []
    for _st, _g in _meta.groupby("d2g_stage", dropna=True):
        _dm = get_domain_map(_st)              # {str(cluster_id): population}
        _fresh.append(_g["d2g_cluster"].map(_leiden_str).map(_dm))
    _meta["fresh_label"] = pd.concat(_fresh).reindex(_meta.index)
    _lut = {
        (lib, cb): lab
        for lib, cb, lab in zip(_meta["arc_library"], _meta["cb16"],
                                _meta["fresh_label"])
        if lab is not None and not pd.isna(lab)
    }
    ts(f"  f1 sidecar: {len(_lut):,} LPM cells carry a fresh SSOT label")

    _old_label = obs_df["d2g_subcluster"].astype("object")
    obs_df["d2g_subcluster"] = [
        _lut.get((lib, cb))
        for lib, cb in zip(obs_df["arc_library"], obs_df["cb16"])
    ]
    _n_old = int(_old_label.notna().sum())
    _n_new = int(obs_df["d2g_subcluster"].notna().sum())
    _n_changed = int(
        (_old_label.fillna("∅").astype(str)
         != obs_df["d2g_subcluster"].fillna("∅").astype(str)).sum()
    )
    # CAMPAIGN-GUARD: run_atac_relabel_campaign.sh + run_hh25_relabel_chain.sh grep
    # this literal ("SSOT relabel via sidecar") as their version pre-flight —
    # do NOT reword it without updating both greps.
    ts(f"  SSOT relabel via sidecar: baked LPM={_n_old:,} fresh LPM={_n_new:,} "
       f"changed={_n_changed:,}")
    # Coverage guard: the cell set is Apr-19-ccreg-defined + stable, so the
    # fresh LPM set must cover essentially all previously-labelled LPM cells.
    # A large shortfall means the (restored) sidecar predates the current peak
    # matrix -> require a proper f1 re-run rather than a silent undercount.
    if _n_old > 0 and _n_new < 0.98 * _n_old:
        sys.exit(
            f"ERROR: SSOT relabel coverage gap — fresh LPM ({_n_new:,}) < 98% "
            f"of baked LPM ({_n_old:,}); the f1 metadata sidecar cell set does "
            f"not cover the current peak matrix. Re-run f1 to regenerate "
            f"step_f1_coupled_atac_metadata.tsv for these cells."
        )

    needed = ("unified_genotype_id", "gex_library", "d2g_subcluster",
              "d2g_stage", "stage")
    for col in needed:
        if col not in obs_df.columns:
            sys.exit(
                f"ERROR: obs missing {col!r}; "
                f"pipeline must be re-run starting from f1."
            )

    # ── Filter to LPM-lineage cells ──────────────────────────────────────
    before = len(obs_df)
    lpm_mask = (
        obs_df["d2g_subcluster"].notna()
        & (obs_df["d2g_subcluster"] != "ambiguous")
    )
    obs_df = obs_df.loc[lpm_mask].copy()
    ts(f"  LPM-lineage cells: {len(obs_df):,} / {before:,}")

    # ── Tech-rep collapse for library join ──────────────────────────────
    obs_df["library_m"] = obs_df["gex_library"].map(collapse_library)

    # ── Attach inferred_side via the L/R sidecar join ───────────────────
    # The sidecar `genotype` column aligns with unified_genotype_id (both "gNN").
    merged = obs_df.merge(
        side_df.rename(columns={"genotype": "unified_genotype_id"}),
        how="left",
        on=["unified_genotype_id", "library_m"],
    )
    n_sided = int(merged["inferred_side"].notna().sum())
    ts(f"  Cells with accepted L/R side: {n_sided:,} / {len(merged):,}")

    sided = merged[merged["inferred_side"].notna()].copy()

    # ── Prepare peak matrix ──────────────────────────────────────────────
    # SnapATAC2 AnnData may return X as a backed sparse matrix; pull a
    # view we can index.  If the file is too large to load, this will
    # need a streaming rewrite — for 30-50K ATAC cells × ~200K peaks it
    # fits in memory at pseudobulk precision.
    X = adata.X[:]
    ts(f"  Peak matrix loaded: dtype={X.dtype} "
       f"sparse={sp.issparse(X)} nnz_or_size={getattr(X, 'nnz', X.size)}")

    peak_names = list(adata.var_names)

    # ── Group + aggregate ────────────────────────────────────────────────
    group_cols = ["unified_genotype_id", "d2g_subcluster",
                  "inferred_side", "d2g_stage"]
    ts(f"  Grouping by {group_cols}...")

    all_meta = []
    count_vectors = []

    for keys, grp in sided.groupby(group_cols, dropna=False):
        (gid, subclu, side, d2g_stage) = keys
        n_cells = len(grp)
        if n_cells < min(PRIMARY_FLOOR, *SENSITIVITY_FLOORS):
            continue

        row_idx = grp["cell_idx"].to_numpy()
        counts = pseudobulk_sum(X, row_idx)

        # Stage = d2g_stage (the bucket; "HH17/16" merged). Every cell left after
        # the LPM filter above carries one; a missing stage is an error.
        if pd.isna(d2g_stage):
            raise SystemExit(f"ERROR: pseudobulk {gid}/{subclu}/{side} has no d2g_stage")
        stage_label = d2g_stage
        # Library list collapsed for this pseudobulk
        libs = sorted(grp["library_m"].dropna().unique().tolist())
        firmness_modes = sorted(grp["firmness"].dropna().unique().tolist())

        sample_id = (
            f"{gid}__{subclu}__{side}__{stage_label}".replace(" ", "_")
        )
        geno_meta = GENOTYPE_META[gid]   # KeyError on a genotype missing from the SSOT

        all_meta.append({
            "sample_id": sample_id,
            "unified_genotype_id": gid,
            "d2g_subcluster": subclu,
            "inferred_side": side,
            "d2g_stage": stage_label,
            "raw_stage": ",".join(sorted(grp["stage"].dropna().unique().tolist())),
            "libraries_collapsed": ",".join(libs),
            "firmness_modes": ",".join(firmness_modes),
            "embryo_group": geno_meta["group"],
            "n_cells": n_cells,
            "total_fragments": int(counts.sum()),
        })
        count_vectors.append(counts)

    if not all_meta:
        sys.exit("ERROR: no pseudobulks produced (all below sensitivity floors)")

    meta_df = pd.DataFrame(all_meta)
    count_mat = np.vstack(count_vectors)
    count_df = pd.DataFrame(
        count_mat, index=meta_df["sample_id"], columns=peak_names,
    )

    ts(f"  Built {len(meta_df):,} pre-filter pseudobulks "
       f"(min below floors: {int(meta_df['n_cells'].min())})")

    # ── Emit at all floors ──────────────────────────────────────────────
    # Explicit deduplicated sorted list so the iteration order is obvious
    # and every floor in the spec gets its own emission.  Two regression
    # guards added after a floor=20 silent miss in an earlier
    # run (only floor=30 + floor=100 landed despite SENSITIVITY_FLOORS
    # = (20, 30, 100)):
    #   1. ALL_FLOORS is built from a set → if PRIMARY_FLOOR happens to
    #      collide with an entry in SENSITIVITY_FLOORS the loop still
    #      visits every distinct floor exactly once (vs the old tuple-
    #      unpack, which would have produced two iterations writing the
    #      same suffix).
    #   2. Every iteration asserts the file landed on disk with non-zero
    #      bytes.  No skip-if-exists path now or ever.
    ALL_FLOORS = sorted({PRIMARY_FLOOR, *SENSITIVITY_FLOORS})
    ts(f"  Emitting at floors: {ALL_FLOORS}")
    for floor in ALL_FLOORS:
        keep = meta_df["n_cells"] >= floor
        sub_counts = count_df.loc[meta_df.loc[keep, "sample_id"]]
        suffix = "" if floor == PRIMARY_FLOOR else f"_floor{floor}"
        counts_path = f"{DA_DIR}/pseudobulk_counts{suffix}.tsv.gz"
        sub_counts.to_csv(counts_path, sep="\t", compression="gzip")
        if not os.path.isfile(counts_path) or os.path.getsize(counts_path) == 0:
            sys.exit(
                f"ERROR: floor={floor} write to {counts_path} produced "
                f"no/empty file; refusing to continue silently."
            )
        ts(f"  floor={floor:3d}: {int(keep.sum()):,} pseudobulks → "
           f"{os.path.basename(counts_path)} "
           f"({os.path.getsize(counts_path):,} bytes)")
        # Write a matched metadata file PER FLOOR (suffix mirrors the counts
        # file) so sensitivity-floor consumers (e.g. g4f at floor=20) get a
        # counts+metadata pair that agree.  Pre-2026-06-15 only the primary
        # floor's metadata was written, which silently blocked floor-20 g4f
        # runs: pseudobulk_counts_floor20.tsv.gz existed but had no matching
        # pseudobulk_metadata_floor20.tsv.  Primary path (suffix="") unchanged
        # → backward compatible.
        meta_path = f"{DA_DIR}/pseudobulk_metadata{suffix}.tsv"
        meta_df.loc[keep].to_csv(meta_path, sep="\t", index=False)
        ts(f"    metadata: {os.path.basename(meta_path)}")

    # ── Summary (primary floor) ──────────────────────────────────────────
    primary = meta_df[meta_df["n_cells"] >= PRIMARY_FLOOR].copy()
    summary = (
        primary
        .groupby(["d2g_stage", "d2g_subcluster", "inferred_side"])
        .agg(
            n_pseudobulks=("sample_id", "count"),
            n_genotypes=("unified_genotype_id", "nunique"),
            total_cells=("n_cells", "sum"),
            median_cells=("n_cells", "median"),
        )
        .reset_index()
        .sort_values(["d2g_stage", "d2g_subcluster", "inferred_side"])
    )
    summary_path = f"{DA_DIR}/pseudobulk_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    ts(f"  summary: {os.path.basename(summary_path)}")

    ts("")
    ts("=== Pseudobulk Summary (primary floor=50) ===")
    for _, row in summary.iterrows():
        ts(f"  {row['d2g_stage']:<10} {row['d2g_subcluster']:<22} "
           f"{row['inferred_side']:<2} "
           f"n_pb={int(row['n_pseudobulks']):>3} "
           f"n_geno={int(row['n_genotypes']):>3} "
           f"cells={int(row['total_cells']):>6,}")

    ts("")
    ts("=== Done ===")


if __name__ == "__main__":
    main()
