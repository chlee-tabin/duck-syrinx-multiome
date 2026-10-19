#!/usr/bin/env python3
"""Phase F Step 1: ATAC QC + coupled barcode matching.

For each ARC library:
  1. Load ATAC fragments via SnapATAC2 (pp.import_fragments, min frags=100)
  2. Compute QC: TSS enrichment, fragment count
  3. Filter: TSS >= 2, fragments >= 1000
  4. Match ATAC barcodes to Phase 1 h5ad via 16bp barcode core + library
  5. Transfer d2b metadata: unified_genotype_id, cell_type, sex, stage
  6. Transfer d2g stage-restricted subcluster labels for LPM-lineage cells
     (plan §3.3) via figures/data/d2g_lpm_subclustering_*_ccreg.tsv.gz
  7. Compute composite peak_pool = '{d2g_bucket}|{subcluster_or_negctrl}'
     used as the groupby for f2 peak calling
  8. Flag uncoupled cells for Route B

L/R side is NOT assigned here — it is deferred to g1_pseudobulk_atac.py
where per-genotype-per-library aggregation allows correct side assignment
from the design table and d2l firm tiers (matching the approach in
e2_laterality_pitx2.py).

Input:
  - cellranger-arc atac_fragments.tsv.gz per ARC library
  - Phase 1 step04_preCB_clustered.h5ad (~165K cells, 52 embryos)
  - d2g stage-restricted subclustering TSVs (ccreg variant, 4 stage files)
  - references/d2g_subcluster_annotations.tsv (via _shared.subcluster_map)

Output:
  - phase2/atac/per_library/*.h5ad  (backed h5ad per ARC library)
  - phase2/atac/step_f1_combined.h5ads  (AnnDataSet combining all libraries)
  - phase2/atac/step_f1_uncoupled_barcodes.tsv.gz  (for Route B rescue)
  - phase2/atac/step_f1_qc_summary.tsv
  - phase2/atac/step_f1_coupled_atac_metadata.tsv (cell-level metadata with
    d2g_subcluster and peak_pool columns)

Run on O2:
    sbatch --mem=64G -t 0-2:00 -c 4 -p short pipeline/run_phase_fg.sh --from f1 --to f1
    # or for a single-library smoke test:
    sbatch pipeline/run_f1_smoke.sh
"""
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import argparse
import snapatac2 as snap
import scanpy as sc
import numpy as np
import pandas as pd
import polars as pl
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "figures"))
from _shared.style import GENOTYPE_META, CLUSTER_ANNOTATIONS
from _shared.subcluster_map import get_leiden_col, get_domain_map

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
ATAC_DIR = f"{PHASE2}/atac"
LIB_DIR = f"{ATAC_DIR}/per_library"
os.makedirs(LIB_DIR, exist_ok=True)

RNA_H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"
# Canonical /n/groups reference (cdw) — never purged, reachable from
# compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror.
DUCK_GTF = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"
GENOME_FAI = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa.fai"

# d2g stage-restricted subclustering outputs (ccreg variant, matching
# scripts/annotate_d2g_via_compendium.py VARIANT).  Per-cell leiden
# assignments for LPM-lineage cells; non-LPM cells are absent from
# these TSVs and get a NULL d2g_subcluster.
# 2026-08-10: read the git-synced figures/data copy, NOT the purge-prone
# phase2/figures working dir (same migration h1 made; the 2026-06 purge
# left phase2/figures empty and f1 would have silently NULL-labeled
# every cell — see the hard error below).
D2G_TSV_DIR = f"{SCRATCH}/figures/data"
D2G_STAGE_FILES = {
    "HH13":    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}

# Plan §3.3 negative-control cell types for peak calling (paraxial at
# early stages; epithelial at later stages when paraxial has largely
# differentiated).  These are d2b-level labels, NOT d2g subclusters,
# because d2g only ran on LPM-lineage cells.  Keys are RAW stage
# labels as they appear in GENOTYPE_META — that table carries both
# the split labels ('HH17' / 'HH16') and the pre-combined form
# ('HH17/16') for different genotypes, so all three map here.
NEG_CONTROL_TYPES = {
    "HH13":    {"Sclerotome", "Dermomyotome"},
    "HH17":    {"Sclerotome", "Dermomyotome"},
    "HH16":    {"Sclerotome", "Dermomyotome"},
    "HH17/16": {"Sclerotome", "Dermomyotome"},
    "HH21":    {"Epithelial"},
    "HH25":    {"Epithelial"},
}

# ── ARC library configuration ────────────────────────────────────────────
# All fragment files rsynced from standby to scratch (compute nodes can't see standby)
ARC_LIBRARIES = {
    "PILOTARC3":  {"fragments": f"{SCRATCH}/PILOTARC3/outs/atac_fragments.tsv.gz",
                   "gex_library": "PILOTGEX3"},
    "PILOTARC4A": {"fragments": f"{SCRATCH}/PILOTARC4A/outs/atac_fragments.tsv.gz",
                   "gex_library": "PILOTGEX4A"},
    "PILOTARC4B": {"fragments": f"{SCRATCH}/PILOTARC4B/outs/atac_fragments.tsv.gz",
                   "gex_library": "PILOTGEX4B"},
    "ARC05":      {"fragments": f"{SCRATCH}/ARC05/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX05"},
    "ARC06":      {"fragments": f"{SCRATCH}/ARC06/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX06"},
    "ARC07":      {"fragments": f"{SCRATCH}/ARC07/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX07"},
    "ARC08":      {"fragments": f"{SCRATCH}/ARC08/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX08"},
    "ARC03-1A":   {"fragments": f"{SCRATCH}/ARC03-1A/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX03-1A"},
    "ARC03-1B":   {"fragments": f"{SCRATCH}/ARC03-1B/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX03-1B"},
    "ARC03-2":    {"fragments": f"{SCRATCH}/ARC03-2/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX03-2"},
    "ARC04":      {"fragments": f"{SCRATCH}/ARC04/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX04"},
    "ARC0S-duck": {"fragments": f"{SCRATCH}/ARC0S-duck/outs/atac_fragments.tsv.gz",
                   "gex_library": "GEX0S-duck"},
}

# QC thresholds (calibrated from 2024 R notebook analysis)
MIN_TSS = 2.0
MIN_FRAGS = 1000
RES = "leiden_res0.2"


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# Shared with downstream phase-F/G
# scripts that may consume AnnData/AnnDataSet.obs from SnapATAC2 2.9.
def snap_obs_to_pandas(obs, obs_names=None, columns=None):
    """Convert a SnapATAC2 PyDataFrameElem to a pandas DataFrame.

    SnapATAC2 2.9's ``AnnData.obs`` / ``AnnDataSet.obs`` is a Rust-side
    ``PyDataFrameElem``.  It rejects ``pd.DataFrame(obj)`` ("DataFrame
    constructor not properly called") and has no ``.to_pandas()``
    method.  Individual columns are iterable, so iterate per name.

    ``columns`` is the list of column names to extract.  If None, we
    try the default set used everywhere in f1; callers can pass an
    explicit list to skip unavailable columns.  ``obs_names`` is the
    iterable used as the DataFrame index; pass ``adata.obs_names`` or
    None to leave the index default (0..N-1).

    Missing columns are silently skipped (no KeyError) because
    AnnDataSet merges of heterogeneous per-library h5ads may drop
    library-specific columns.
    """
    if columns is None:
        columns = [
            "tsse", "n_fragment",
            "cb16", "arc_library", "gex_library",
            "unified_genotype_id", "cell_type", "sex", "demux_group",
            "stage", "d2g_stage", "d2g_cluster", "d2g_subcluster",
            "peak_pool",
        ]
    data = {}
    for col in columns:
        try:
            data[col] = list(obs[col])
        except Exception:  # noqa: BLE001
            continue
    idx = list(obs_names) if obs_names is not None else None
    return pd.DataFrame(data, index=idx)


def parse_chrom_sizes(fai_path):
    """Parse genome.fa.fai into {chrom: size} dict for SnapATAC2."""
    chrom_sizes = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom_sizes[parts[0]] = int(parts[1])
    ts(f"  Parsed {len(chrom_sizes)} chromosomes from {fai_path}")
    return chrom_sizes


def stage_to_d2g_bucket(raw_stage):
    """Map a raw stage label to the d2g annotations bucket.

    Accepts both the split labels (``HH17`` / ``HH16``) and the
    pre-combined label (``HH17/16``) because GENOTYPE_META already
    carries the combined form for some genotypes — the smoke test on
    ARC07 revealed that all 9,149 HH17/16-stage cells were getting
    peak_pool = None because the combined label fell through the
    membership test.  Returns None for unrecognised stages such as
    ``?`` from unstaged genotypes.
    """
    if raw_stage in ("HH17", "HH16", "HH17/16"):
        return "HH17/16"
    if raw_stage in ("HH13", "HH21", "HH25"):
        return raw_stage
    return None


def load_d2g_subcluster_lookup():
    """Build a (library, cb16) → d2g_subcluster lookup from d2g TSVs.

    Each d2g stage-restricted subclustering TSV contains one row per
    LPM-lineage cell at that stage, with `library`, `barcode`, and the
    full leiden resolution sweep.  Canonical resolution per stage comes
    from references/d2g_subcluster_annotations.tsv via
    _shared.subcluster_map.get_leiden_col.  Population names come from
    get_domain_map (same TSV, joined on cluster id).

    Returns a pandas DataFrame indexed by (library, cb16) with columns:
      d2g_stage      — bucket name ('HH13', 'HH17/16', 'HH21', 'HH25')
      d2g_cluster    — int cluster ID at the canonical resolution
      d2g_subcluster — population label (or 'ambiguous')

    A MISSING stage TSV is a HARD ERROR (FileNotFoundError): a silently
    NULL-labeled stage would corrupt every f4/f5/f6/g1 consumer — run
    the repo-to-scratch sync (figures/data) before f1.  Each
    loaded TSV's leiden ids must also EXACTLY match the annotation
    domain map for its bucket (RuntimeError otherwise): a present-but-
    different-run sidecar, or a reshuffled regen reusing the same id
    range, would otherwise map silently to wrong/NULL populations.
    """
    ts("Loading d2g stage-restricted subcluster assignments...")
    all_parts = []
    for bucket, fname in D2G_STAGE_FILES.items():
        path = os.path.join(D2G_TSV_DIR, fname)
        if not os.path.isfile(path):
            # Hard error (was a warning): a missing sidecar silently NULL-labels
            # a whole stage and every f4/f5/f6 consumer downstream — the exact
            # silent-fallback class this repo bans. Sync figures/data first.
            raise FileNotFoundError(
                f"{path} missing — {bucket} would get NULL d2g_subcluster. "
                "Copy the repository figures/data/ into the scratch working tree first.")
        leiden_col = get_leiden_col(bucket)
        domain = get_domain_map(bucket)  # {cluster_id_str: population}

        df = pd.read_csv(
            path, sep="\t", compression="gzip",
            usecols=["library", "barcode", leiden_col],
        )
        df = df.rename(columns={leiden_col: "d2g_cluster"})
        df["d2g_stage"] = bucket
        # Same-run guard: the sidecar and the annotation
        # TSV must come from the SAME d2g run.  Set equality — not mere
        # overlap — catches a reshuffled regen that reuses the id range 0..N
        # (the exact incident class behind).  Sibling convention:
        # g1_pseudobulk_atac.py's >=98%-coverage hard exit.
        sidecar_ids = set(df["d2g_cluster"].astype(str))
        domain_ids = set(domain)
        if sidecar_ids != domain_ids:
            raise RuntimeError(
                f"{bucket}: sidecar leiden ids {sorted(sidecar_ids)} != "
                f"annotation ids {sorted(domain_ids)} (only-in-sidecar="
                f"{sorted(sidecar_ids - domain_ids)}, only-in-annotation="
                f"{sorted(domain_ids - sidecar_ids)}) — the d2g sidecar and "
                f"references/d2g_subcluster_annotations.tsv are from "
                f"DIFFERENT d2g runs; re-sync or rerun d2g_validate_subclusters.")
        df["d2g_subcluster"] = df["d2g_cluster"].astype(str).map(domain)
        # cb16 = barcode core; d2g TSV barcodes carry the '-1' gem suffix
        df["cb16"] = [
            b.rsplit("-", 1)[0] if "-" in b else b
            for b in df["barcode"]
        ]
        n_total = len(df)
        n_annot = int(df["d2g_subcluster"].notna().sum())
        ts(f"  {bucket}: {n_total:,} cells ({n_annot:,} with annotation)")
        all_parts.append(df[[
            "library", "cb16", "d2g_stage", "d2g_cluster", "d2g_subcluster",
        ]])

    combined = pd.concat(all_parts, ignore_index=True)
    before = len(combined)
    combined = combined.drop_duplicates(subset=["library", "cb16"], keep="first")
    if len(combined) < before:
        ts(f"  Dropped {before - len(combined):,} duplicate (library, cb16) rows "
           f"(cross-stage overlap)")
    ts(f"  Combined: {len(combined):,} cells with d2g labels")
    return combined.set_index(["library", "cb16"])


def compute_peak_pool(stage, d2g_subcluster, cell_type, d2g_stage=None):
    """Per-cell composite label for peak calling (plan §3.3).

    LPM-lineage cells: f'{d2g_bucket}|{d2g_subcluster}' — stage-specific
    subcluster granularity.  Non-LPM cells: f'{d2g_bucket}|{cell_type}'
    when cell_type is a stage-appropriate negative control.  Otherwise
    None (excluded from peak calling).

    'ambiguous' d2g labels are excluded (cells whose compendium scores
    did not resolve to a population — these should not contribute to
    peak calls for any subcluster).

    For the LPM branch, prefer ``d2g_stage`` (set by the compendium
    scorer when a subcluster label is assigned) over the raw ``stage``
    (from RNA metadata).  ARC04 and ARC0S-duck have stage='?' for every
    coupled cell — the embryos weren't staged pre-sequencing — but the
    compendium successfully assigned 3,249 + 278 of their cells to
    HH17/16 subclusters based on marker profile.  Without this fallback
    those cells were silently dropped from peak calling (n_peak_pool=0
    in the 37810403 QC summary) despite carrying valid d2g labels.
    """
    if pd.notna(d2g_subcluster) and d2g_subcluster != "ambiguous":
        bucket = d2g_stage if pd.notna(d2g_stage) else stage_to_d2g_bucket(stage)
        if bucket is None:
            return None
        return f"{bucket}|{d2g_subcluster}"
    bucket = stage_to_d2g_bucket(stage)
    if bucket is None:
        return None
    if cell_type in NEG_CONTROL_TYPES.get(stage, set()):
        return f"{bucket}|{cell_type}"
    return None


def load_rna_metadata():
    """Load Phase 1 RNA h5ad and extract per-cell metadata for barcode matching."""
    ts("Loading Phase 1 RNA h5ad...")
    adata = sc.read_h5ad(RNA_H5AD)
    ts(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    ann_map = {str(k): v for k, v in CLUSTER_ANNOTATIONS.items()}
    adata.obs["cell_type"] = adata.obs[RES].astype(str).map(ann_map)

    cols = ["library", "unified_genotype_id", "cell_type", "sex", "demux_group"]
    available = [c for c in cols if c in adata.obs.columns]
    meta = adata.obs[available].copy()

    # Extract 16bp barcode core — strip -N or -S suffix
    meta["cb16"] = meta.index.str.replace(r"-(\d+|S)$", "", regex=True)

    ts(f"  Extracted metadata for {len(meta):,} cells")
    return meta


D2G_COLS = ["d2g_stage", "d2g_cluster", "d2g_subcluster"]


def _attach_d2g(obs_pd, d2g_lookup):
    """(Re)attach d2g subcluster labels to ``obs_pd`` in place.

    Join on (gex_library, cb16); cells absent from d2g (non-LPM lineage or
    dropped during subclustering) end up with NA.  Used by both the fresh
    import and the cache-reuse path, so a d2g rerun reaches cached
    libraries too.
    """
    if len(d2g_lookup) > 0:
        obs_keys = pd.MultiIndex.from_arrays(
            [obs_pd["gex_library"].values, obs_pd["cb16"].values],
            names=["library", "cb16"],
        )
        joined = d2g_lookup.reindex(obs_keys)
        for col in D2G_COLS:
            obs_pd[col] = joined[col].to_numpy()
    else:
        for col in D2G_COLS:
            obs_pd[col] = pd.Series([None] * len(obs_pd), dtype=object)


def _try_reuse_cached(h5ad_path, arc_name, gex_lib, d2g_lookup):
    """Probe a cached per-library h5ad for completeness.

    Returns (True, stats, uncoupled_bc, obs_pd) if peak_pool is present,
    or (False, None, None, None) otherwise.  On any read error, also
    returns the False tuple so the caller re-imports.  obs_pd is
    returned so main() can aggregate the per-library obs without
    re-opening the child h5ads — the AnnDataSet already holds an
    exclusive HDF5 handle on each child, so snap.read of the same file
    later would hit "H5Fopen: file is already open for read-only".
    """
    # BaseException, not Exception: SnapATAC2's Rust layer raises
    # pyo3_runtime.PanicException (a BaseException) on a truncated HDF5 —
    # e.g. the 2 KB ARC04.h5ad stub left by a TIMEOUT mid-write (job
    # 50248206 died on it). Treat any read failure as a cache miss.
    try:
        cached = snap.read(h5ad_path)
    except (KeyboardInterrupt, SystemExit, GeneratorExit):
        raise
    except BaseException as e:  # noqa: BLE001
        ts(f"  WARN {arc_name}: could not open cached h5ad ({e}) — re-importing")
        return False, None, None, None

    try:
        list(cached.obs["peak_pool"])
    except (KeyboardInterrupt, SystemExit, GeneratorExit):
        raise
    except BaseException as e:  # noqa: BLE001
        ts(f"  WARN {arc_name}: cached h5ad unreadable at obs['peak_pool'] "
           f"({e}) — re-importing")
        cached.close()
        return False, None, None, None

    obs_pd = snap_obs_to_pandas(cached.obs, obs_names=cached.obs_names)
    n_pass = cached.n_obs
    cached.close()

    # Re-join d2g labels from the CURRENT d2g sidecars: the cached obs
    # carries whatever d2g labels existed when the h5ad was written, so a
    # d2g rerun would otherwise never reach cached libraries.  The refreshed
    # labels flow into the in-memory aggregation (the f1 metadata sidecar);
    # like peak_pool below, the on-disk per-library obs is NOT rewritten —
    # delete the per-library caches after a d2g rerun if f4/f5 must see the
    # new labels in step_f1_combined.h5ads.
    old_cluster = (obs_pd["d2g_cluster"].copy()
                   if "d2g_cluster" in obs_pd.columns else None)
    _attach_d2g(obs_pd, d2g_lookup)
    if old_cluster is not None:
        a = pd.Series(old_cluster.to_numpy(), dtype=object)
        b = pd.Series(obs_pd["d2g_cluster"].to_numpy(), dtype=object)
        same = (a.isna() & b.isna()) | (
            a.notna() & b.notna()
            & (a.astype(str).str.replace(r"\.0$", "", regex=True)
               == b.astype(str).str.replace(r"\.0$", "", regex=True)))
        n_changed = int((~same).sum())
        if n_changed:
            ts(f"  WARN {arc_name}: {n_changed:,} cached d2g_cluster labels "
               f"differ from current d2g sidecars — using current")

    # Always recompute peak_pool from cached obs.  The column may be
    # stale if compute_peak_pool() has changed since the h5ad was
    # written — e.g. ARC04 / ARC0S-duck in job 37810403 had
    # peak_pool=all-None on disk because the old compute_peak_pool()
    # used the raw RNA ``stage`` (='?' for those embryos) instead of
    # ``d2g_stage``.  Re-derive in memory each time so downstream
    # aggregation sees the current function's output; the on-disk
    # value will drift but is not read again in this process.
    if all(k in obs_pd.columns for k in
           ("stage", "d2g_subcluster", "cell_type", "d2g_stage")):
        obs_pd["peak_pool"] = [
            compute_peak_pool(s, d, c, ds)
            for s, d, c, ds in zip(
                obs_pd["stage"].to_numpy(),
                obs_pd["d2g_subcluster"].to_numpy(),
                obs_pd["cell_type"].to_numpy(),
                obs_pd["d2g_stage"].to_numpy(),
            )
        ]

    n_coupled = int(obs_pd["unified_genotype_id"].notna().sum()) \
        if "unified_genotype_id" in obs_pd else 0
    n_uncoupled = n_pass - n_coupled
    n_pool = int(obs_pd["peak_pool"].notna().sum())
    n_d2g = int(obs_pd["d2g_subcluster"].notna().sum()) \
        if "d2g_subcluster" in obs_pd else 0

    stats = {
        "arc_library": arc_name,
        "gex_library": gex_lib,
        "n_raw_barcodes": -1,
        "n_qc_pass": n_pass,
        "n_coupled": n_coupled,
        "n_uncoupled": n_uncoupled,
        "coupled_frac": n_coupled / n_pass if n_pass > 0 else 0,
        "n_d2g_annotated": n_d2g,
        "n_peak_pool": n_pool,
        "median_tsse": float(obs_pd["tsse"].median())
            if "tsse" in obs_pd else 0.0,
        "median_fragments": float(obs_pd["n_fragment"].median())
            if "n_fragment" in obs_pd else 0.0,
    }

    if "unified_genotype_id" in obs_pd and "cb16" in obs_pd:
        uncoupled_mask = obs_pd["unified_genotype_id"].isna()
        uncoupled_bc = obs_pd.loc[
            uncoupled_mask,
            ["cb16", "arc_library", "gex_library", "tsse", "n_fragment"],
        ].copy()
    else:
        uncoupled_bc = pd.DataFrame()

    ts(f"  SKIP {arc_name}: cached h5ad complete "
       f"({n_pass:,} cells, {n_pool:,} peak_pool)")
    return True, stats, uncoupled_bc, obs_pd


def process_arc_library(arc_name, config, chrom_sizes, rna_meta, d2g_lookup):
    """Process a single ARC library: import, QC, match to RNA, attach d2g.

    Returns (stats_dict, h5ad_path, uncoupled_df) or None if missing.
    """
    frag_path = config["fragments"]
    gex_lib = config["gex_library"]

    if not os.path.isfile(frag_path):
        ts(f"  SKIP {arc_name}: fragment file not found")
        return None

    h5ad_path = os.path.join(LIB_DIR, f"{arc_name}.h5ad")

    # Idempotence: reuse a cached per-library h5ad if the writeback
    # completed (peak_pool column present).  A crash mid-writeback
    # (e.g. ARC04 in job 37796513 hit the all-None Polars quirk) leaves
    # a zombie with fragments but no obs metadata — delete it and
    # re-import so the fragments get re-imported cleanly.  The ~15 min
    # import cost dominates per-library wall time, so skipping clean
    # libraries turns a full 4-5h rerun into a ~30 min salvage.
    if os.path.isfile(h5ad_path):
        cached_reused, stats, uncoupled_bc, obs_pd = _try_reuse_cached(
            h5ad_path, arc_name, gex_lib, d2g_lookup,
        )
        if cached_reused:
            return stats, h5ad_path, uncoupled_bc, obs_pd
        ts(f"  ZOMBIE {arc_name}: cached h5ad missing peak_pool — re-importing")
        os.remove(h5ad_path)

    ts(f"  Loading {arc_name} fragments...")

    # Import fragments → backed h5ad (preserves fragment data for downstream).
    # A pre-filter of min_num_fragments=100 excludes ~95% of background
    # barcodes (typical ARC fragment file has ~700K barcodes, of which
    # only ~30K are cellranger-called cells).  The stricter QC filter
    # (MIN_TSS=2, MIN_FRAGS=1000) runs below on the pre-filtered set;
    # cells with 100-999 fragments are retained briefly so tsse() can
    # compute on a reasonable denominator before we drop them.
    adata = snap.pp.import_fragments(
        frag_path,
        chrom_sizes=chrom_sizes,
        file=h5ad_path,
        min_num_fragments=100,
        sorted_by_barcode=False,
    )
    n_raw = adata.n_obs
    ts(f"    {n_raw:,} raw barcodes (pre-filtered at min_num_fragments=100)")

    # ── QC: TSS enrichment ──
    # Parameter name varies by version: gene_anno (v2.6+) or gtf_file (older)
    try:
        snap.metrics.tsse(adata, gene_anno=DUCK_GTF)
    except TypeError:
        snap.metrics.tsse(adata, DUCK_GTF)
    ts(f"    TSS enrichment: median={adata.obs['tsse'].median():.2f}")
    ts(f"    Fragments: median={adata.obs['n_fragment'].median():,.0f}")

    # ── Apply QC filters (TSS + fragment count only) ──
    # SnapATAC2 does not compute nucleosome_signal natively; filter on TSS + frags
    mask = (adata.obs["tsse"] >= MIN_TSS) & (adata.obs["n_fragment"] >= MIN_FRAGS)
    # Use np.asarray for Polars/Pandas compatibility
    adata.subset(np.asarray(mask))
    n_pass = adata.n_obs
    pct = 100 * n_pass / n_raw if n_raw > 0 else 0
    ts(f"    {n_pass:,}/{n_raw:,} pass QC ({pct:.1f}%)")

    # ── Build a pandas working copy of obs ──
    # See snap_obs_to_pandas above for why we can't use pd.DataFrame
    # on PyDataFrameElem directly.  At this point only tsse and
    # n_fragment exist on adata.obs (added by import_fragments + tsse);
    # the rest of the columns we populate below.
    obs_pd = snap_obs_to_pandas(
        adata.obs,
        obs_names=adata.obs_names,
        columns=["tsse", "n_fragment"],
    )

    # cellranger-arc ATAC barcodes end in "-1"; strip the gem suffix.
    obs_pd["cb16"] = [
        b.rsplit("-", 1)[0] if "-" in b else b
        for b in obs_pd.index
    ]
    obs_pd["arc_library"] = arc_name
    obs_pd["gex_library"] = gex_lib

    # RNA lookup indexed by cb16 within this GEX library
    rna_lib = rna_meta[rna_meta["library"] == gex_lib].copy()
    rna_lookup = rna_lib.drop_duplicates(subset="cb16").set_index("cb16")
    rna_lookup = rna_lookup[["unified_genotype_id", "cell_type", "sex", "demux_group"]]
    for col in rna_lookup.columns:
        obs_pd[col] = obs_pd["cb16"].map(rna_lookup[col])

    n_coupled = int(obs_pd["unified_genotype_id"].notna().sum())
    n_uncoupled = n_pass - n_coupled
    if n_pass > 0:
        ts(f"    Coupled: {n_coupled:,} ({100 * n_coupled / n_pass:.1f}%)")
    else:
        ts("    Coupled: 0")
    ts(f"    Uncoupled: {n_uncoupled:,}")

    # Stage from GENOTYPE_META (raw stage; d2g bucket is derived below)
    obs_pd["stage"] = obs_pd["unified_genotype_id"].map(
        lambda x: GENOTYPE_META.get(x, {}).get("stage", None) if pd.notna(x) else None
    )

    # ── Attach d2g stage-restricted subcluster labels (plan §3.3) ──
    _attach_d2g(obs_pd, d2g_lookup)

    n_d2g = int(obs_pd["d2g_subcluster"].notna().sum())
    ts(f"    d2g-annotated: {n_d2g:,}")

    # ── Composite peak_pool label (plan §3.3) ──
    obs_pd["peak_pool"] = [
        compute_peak_pool(s, d, c, ds)
        for s, d, c, ds in zip(
            obs_pd["stage"].to_numpy(),
            obs_pd["d2g_subcluster"].to_numpy(),
            obs_pd["cell_type"].to_numpy(),
            obs_pd["d2g_stage"].to_numpy(),
        )
    ]
    n_pool = int(obs_pd["peak_pool"].notna().sum())
    ts(f"    peak_pool assigned: {n_pool:,}")

    # ── Write the new columns back onto adata.obs ──
    # SnapATAC2 2.9 builds each obs column as a strictly-typed Polars
    # Series from whatever the first non-null value is.  If the first
    # value is a string but a later value is numpy.nan, Polars raises
    # "unexpected value while building Series of type String; found
    # value of type Float64: NaN".
    #
    # pandas Series.where(notna(), None) does NOT reliably produce
    # Python None in the resulting list: for float dtypes (including
    # int columns that NaN-promoted to float on a reindex miss) the
    # None gets coerced back into numpy.nan at .tolist() time.  Do
    # the NaN -> None substitution explicitly in a list comprehension
    # so Polars sees a null instead of a stray NaN.
    #
    # d2g_cluster is semantically an integer (leiden cluster id) that
    # was promoted to float64 by the reindex; cast the valid values
    # back to int so the column ends up Int64-with-nulls in the h5ad
    # (otherwise Polars gets a Float64 column and `cluster == 0`
    # comparisons in f2/g1 would surprise).
    new_cols = [
        "cb16", "arc_library", "gex_library",
        "unified_genotype_id", "cell_type", "sex", "demux_group",
        "stage", "d2g_stage", "d2g_cluster", "d2g_subcluster",
        "peak_pool",
    ]
    # (6th Polars quirk, §8) A plain list of all-None values cannot be
    # written through ``adata.obs[col] = vals`` — the Polars backend
    # raises "Unsupported series data type: Null" when it can't infer a
    # dtype from the values.  Hand it a ``pl.Series`` with an explicit
    # dtype in that case.  Libraries with zero peak-pool assignments
    # (e.g. ARC04 in job 37796513) hit this for d2g_* + peak_pool.
    for col in new_cols:
        if col not in obs_pd.columns:
            continue
        raw = obs_pd[col].tolist()
        if col == "d2g_cluster":
            vals = [None if pd.isna(v) else int(v) for v in raw]
            dtype = pl.Int64
        else:
            vals = [None if (v is None or pd.isna(v)) else v for v in raw]
            dtype = pl.Utf8
        if all(v is None for v in vals):
            adata.obs[col] = pl.Series(name=col, values=vals, dtype=dtype)
        else:
            adata.obs[col] = vals

    # ── Collect stats (use the pandas obs copy for .median / .loc) ──
    stats = {
        "arc_library": arc_name,
        "gex_library": gex_lib,
        "n_raw_barcodes": n_raw,
        "n_qc_pass": n_pass,
        "n_coupled": n_coupled,
        "n_uncoupled": n_uncoupled,
        "coupled_frac": n_coupled / n_pass if n_pass > 0 else 0,
        "n_d2g_annotated": n_d2g,
        "n_peak_pool": n_pool,
        "median_tsse": float(obs_pd["tsse"].median()),
        "median_fragments": float(obs_pd["n_fragment"].median()),
    }

    # Uncoupled barcodes for Route B (genotype-less ATAC nuclei that
    # need reference-VCF rescue in Phase 2 strategy step 3).
    uncoupled_mask = obs_pd["unified_genotype_id"].isna()
    uncoupled_bc = obs_pd.loc[
        uncoupled_mask,
        ["cb16", "arc_library", "gex_library", "tsse", "n_fragment"],
    ].copy()

    return stats, h5ad_path, uncoupled_bc, obs_pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--library",
        default=None,
        help=(
            "Restrict to a single ARC library (smoke-test mode).  "
            "Per-library h5ad still writes to the canonical path; "
            "aggregate outputs gain a _smoke_<LIB> suffix to avoid "
            "clobbering a later full-cohort run."
        ),
    )
    args = parser.parse_args()

    ts("=== Phase F Step 1: ATAC QC + Coupled Barcode Matching ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")

    libs_to_run = ARC_LIBRARIES
    suffix = ""
    if args.library is not None:
        if args.library not in ARC_LIBRARIES:
            sys.exit(
                f"ERROR: unknown library {args.library!r}.  "
                f"Available: {sorted(ARC_LIBRARIES.keys())}"
            )
        libs_to_run = {args.library: ARC_LIBRARIES[args.library]}
        suffix = f"_smoke_{args.library}"
        ts(f"  SMOKE MODE: processing only {args.library}")

    # Parse genome chrom sizes from .fai
    chrom_sizes = parse_chrom_sizes(GENOME_FAI)

    # Load RNA metadata
    rna_meta = load_rna_metadata()

    # Load d2g stage-restricted subcluster lookup (plan §3.3)
    d2g_lookup = load_d2g_subcluster_lookup()

    # Process each ARC library → backed h5ad
    all_stats = []
    h5ad_paths = []  # (name, path) for AnnDataSet
    all_uncoupled = []
    per_lib_obs = []  # pandas obs per library — avoids re-opening h5ads

    for arc_name, config in libs_to_run.items():
        ts(f"\n--- Processing {arc_name} ---")
        result = process_arc_library(arc_name, config, chrom_sizes, rna_meta, d2g_lookup)
        if result is None:
            continue
        stats, h5ad_path, uncoupled, obs_pd = result
        all_stats.append(stats)
        h5ad_paths.append((arc_name, h5ad_path))
        all_uncoupled.append(uncoupled)
        obs_pd = obs_pd.copy()
        obs_pd["arc_library"] = arc_name
        per_lib_obs.append(obs_pd)

    if not h5ad_paths:
        ts("ERROR: No ARC libraries processed successfully")
        return

    # ── Create AnnDataSet (preserves fragment backing for downstream) ──
    ts(f"\n--- Creating AnnDataSet from {len(h5ad_paths)} libraries ---")
    dataset_path = f"{ATAC_DIR}/step_f1_combined{suffix}.h5ads"
    dataset = snap.AnnDataSet(
        adatas=h5ad_paths,
        filename=dataset_path,
        add_key="arc_library",
    )
    ts(f"  Total ATAC cells: {dataset.n_obs:,}")

    # SnapATAC2 2.9's AnnDataSet.obs only surfaces the add_key column
    # ("arc_library"); the per-library obs columns we populated above
    # are NOT automatically propagated to the merged view.  Earlier
    # attempts re-opened each child h5ad with snap.read() after the
    # AnnDataSet was constructed, but the AnnDataSet holds exclusive
    # HDF5 handles on its children and that hit "H5Fopen: file is
    # already open for read-only" (7th Polars-backend issue, §8.4).
    # Instead, each process_arc_library / _try_reuse_cached call now
    # hands back its pandas obs; concat in the order they were added
    # so the merged frame's cell order matches AnnDataSet's.
    obs_df = pd.concat(per_lib_obs, axis=0)
    ts(f"  Aggregated obs: {obs_df.shape[0]:,} rows × {obs_df.shape[1]:,} cols")

    # Propagate the columns f2 / g1 need onto AnnDataSet.obs.  The
    # default AnnDataSet only carries the add_key column ("arc_library")
    # per §8.2; snap.tl.macs3(groupby=..., replicate=...) needs these
    # names resolvable on the merged view.  The concat order of
    # per_lib_obs matches AnnDataSet's cell stacking order, so we can
    # assign directly by position.  Null columns go through pl.Series
    # with an explicit dtype to dodge §8.3.
    ts("  Propagating obs columns to AnnDataSet.obs for f2/g1...")
    for col in ("peak_pool", "unified_genotype_id", "d2g_stage",
                "d2g_subcluster", "cell_type", "stage"):
        if col not in obs_df.columns:
            continue
        vals = [
            None if (v is None or (isinstance(v, float) and pd.isna(v)))
            else v
            for v in obs_df[col].tolist()
        ]
        if all(v is None for v in vals):
            dataset.obs[col] = pl.Series(name=col, values=vals, dtype=pl.Utf8)
        else:
            dataset.obs[col] = vals
        ts(f"    {col}: {sum(1 for v in vals if v is not None):,}/{len(vals):,} non-null")

    coupled_mask = obs_df["unified_genotype_id"].notna()
    n_coupled = int(coupled_mask.sum())
    n_uncoupled = obs_df.shape[0] - n_coupled
    ts(f"  Coupled: {n_coupled:,}")
    ts(f"  Uncoupled: {n_uncoupled:,}")

    # ── Save outputs ──
    # Uncoupled barcodes for Route B
    if all_uncoupled:
        uncoupled_df = pd.concat(all_uncoupled, ignore_index=True)
        uncoupled_path = f"{ATAC_DIR}/step_f1_uncoupled_barcodes{suffix}.tsv.gz"
        uncoupled_df.to_csv(uncoupled_path, sep="\t", index=False,
                            compression="gzip")
        ts(f"  Uncoupled barcodes: {len(uncoupled_df):,} → {uncoupled_path}")

    # QC summary
    stats_df = pd.DataFrame(all_stats)
    stats_path = f"{ATAC_DIR}/step_f1_qc_summary{suffix}.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False, float_format="%.4f")
    ts(f"  QC summary: {stats_path}")

    # Export cell-level metadata as TSV for the R coverage script
    meta_cols = ["arc_library", "gex_library", "cb16",
                 "unified_genotype_id", "cell_type", "sex", "stage",
                 "d2g_stage", "d2g_cluster", "d2g_subcluster", "peak_pool"]
    available_cols = [c for c in meta_cols if c in obs_df.columns]
    meta_path = f"{ATAC_DIR}/step_f1_coupled_atac_metadata{suffix}.tsv"
    obs_df[available_cols].to_csv(meta_path, sep="\t", index=True)
    ts(f"  Cell metadata: {meta_path}")

    # Print summary table
    ts("\n=== QC Summary ===")
    for _, row in stats_df.iterrows():
        ts(f"  {row['arc_library']:>15} → {row['gex_library']:<12} "
           f"raw={row['n_raw_barcodes']:>7,} "
           f"pass={row['n_qc_pass']:>6,} "
           f"coupled={row['n_coupled']:>6,} ({row['coupled_frac']:.1%})")

    totals = stats_df[["n_raw_barcodes", "n_qc_pass", "n_coupled",
                        "n_uncoupled"]].sum()
    ts(f"\n  TOTAL: {totals['n_qc_pass']:,.0f} QC-pass, "
       f"{totals['n_coupled']:,.0f} coupled "
       f"({totals['n_coupled']/totals['n_qc_pass']:.1%}), "
       f"{totals['n_uncoupled']:,.0f} uncoupled")

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
