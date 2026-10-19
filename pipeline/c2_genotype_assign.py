#!/usr/bin/env python3
"""Phase C2: Assign unified genotype IDs and filter to confirmed embryos.

Reads the evidence-derived donor mapping (produced by b1_donor_consolidation.py),
applies it to assign unified_genotype_id, and removes noise/artifact/duplicate
genotypes. After this step, only confirmed embryos remain.

Two-priority lookup (matching step02f logic):
  1. Exact (demux_group, raw_donor_id) in unified_donor_map.tsv
  2. donor_id is already a valid unified_id (handles demux groups like
     G6_combined_hh13 that pre-resolve donor names to genotype IDs)

Idempotent: reads from c1 output, writes to a separate clean path.
Re-running always produces the same result from the same inputs.

Input:  phase2/step02_preCB_qc.h5ad  (c1 output: all QC-passed singlet cells)
        universal_genotyping/identity_table/unified_donor_map.tsv  (scratch — the
          chain SSOT written by b1/b1b; fallback to the committed KEEP-only
          references/unified_donor_map.tsv ONLY with DUCK_ALLOW_REPO_DONOR_MAP=1)
        universal_genotyping/identity_table/overmerge_barcode_overrides.tsv  (b1b;
          warns loudly when absent — barcode-level splits are then not applied)
Output: phase2/step02_preCB_clean.h5ad  (only confirmed embryos)

Run via sbatch:
    sbatch pipeline/run_pipeline.sh  (or standalone with 64G)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import pandas as pd
import numpy as np
import os
import time

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
# Donor-map resolution (review fix 2026-08-26): the demux CHAIN's own product —
# scratch identity_table, written by b1 and rewritten by b1b (EXCLUDE rows) —
# MUST win inside a pipeline run, or c2 silently ignores the map b1/b1b just
# produced (the committed 71-row references/ copy is a KEEP-only reconstruction
# with no EXCLUDE rows). The committed copy is the offline/public fallback only,
# behind an explicit opt-in (DUCK_ALLOW_REPO_DONOR_MAP=1) — a missing scratch map
# in a chain run raises instead of silently substituting the reconstruction.
_REFS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "references")
_REPO_DONOR_MAP = os.path.join(_REFS, "unified_donor_map.tsv")
_SCRATCH_DONOR_MAP = f"{SCRATCH}/universal_genotyping/identity_table/unified_donor_map.tsv"
USING_REPO_FALLBACK = False
if os.path.isfile(_SCRATCH_DONOR_MAP):
    DONOR_MAP_PATH = _SCRATCH_DONOR_MAP
elif os.environ.get("DUCK_ALLOW_REPO_DONOR_MAP") == "1":
    DONOR_MAP_PATH = _REPO_DONOR_MAP
    USING_REPO_FALLBACK = True
    print(
        "[c2] WARNING: scratch identity_table missing — using the committed\n"
        "[c2] references/unified_donor_map.tsv (DUCK_ALLOW_REPO_DONOR_MAP=1).\n"
        "[c2] It is a KEEP-only reconstruction with NO EXCLUDE rows, and the\n"
        "[c2] barcode-level overmerge overrides are scratch-only — b1b's\n"
        "[c2] exclusions/splits will NOT be applied. Offline use only.",
        flush=True,
    )
else:
    raise FileNotFoundError(
        f"scratch donor map not found: {_SCRATCH_DONOR_MAP}\n"
        "Run pipeline/b1_donor_consolidation.py (+ b1b) first — the chain's own map\n"
        "must win inside a pipeline run. For an offline run without scratch, set\n"
        "DUCK_ALLOW_REPO_DONOR_MAP=1 to use the committed references/ copy\n"
        "(KEEP-only; no EXCLUDE rows, no overmerge barcode overrides)."
    )
print(f"[c2] donor map: {DONOR_MAP_PATH}")
OVERMERGE_OVERRIDES_PATH = f"{SCRATCH}/universal_genotyping/identity_table/overmerge_barcode_overrides.tsv"

# Separate input/output paths — never overwrite the input
INPUT_PATH = f"{PHASE2}/step02_preCB_qc.h5ad"       # c1 output
OUTPUT_PATH = f"{PHASE2}/step02_preCB_clean.h5ad"    # c2 output


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_donor_map(tsv_path):
    """Load unified_donor_map.tsv and return (mapping_dict, valid_embryos_set).

    mapping_dict: (demux_group, raw_donor_id) -> unified_id  (non-EXCLUDE only)
    valid_embryos: set of unified_id values (the confirmed embryos)
    """
    if not os.path.isfile(tsv_path):
        raise FileNotFoundError(
            f"unified_donor_map.tsv not found: {tsv_path}\n"
            "Run pipeline/b1_donor_consolidation.py first to generate it."
        )
    df = pd.read_csv(tsv_path, sep="\t")
    required_cols = {"demux_group", "raw_donor_id", "unified_id"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"unified_donor_map.tsv missing columns: {missing}")

    ts(f"  Loaded {len(df)} rows from {os.path.basename(tsv_path)}")
    ts(f"  Actions: {df['action'].value_counts().to_dict()}")

    # Build mapping from non-EXCLUDE entries
    keep = df[df["unified_id"] != "EXCLUDE"]
    mapping = {
        (row["demux_group"], row["raw_donor_id"]): row["unified_id"]
        for _, row in keep.iterrows()
    }
    valid_embryos = set(keep["unified_id"].unique())

    n_excl = (df["unified_id"] == "EXCLUDE").sum()
    ts(f"  {len(valid_embryos)} confirmed embryos, {n_excl} excluded entries")
    return mapping, valid_embryos


def main():
    ts("=== Phase C2: Assign genotypes to confirmed embryos ===")

    # ---- 0. Load the single source of truth ----
    ts(f"Loading donor map: {DONOR_MAP_PATH}")
    donor_map, valid_embryos = load_donor_map(DONOR_MAP_PATH)
    ts(f"  Valid embryos: {sorted(valid_embryos)[:5]}... ({len(valid_embryos)} total)")

    # ---- 1. Load c1 output (never modified by this script) ----
    if not os.path.isfile(INPUT_PATH):
        raise FileNotFoundError(
            f"C1 output not found: {INPUT_PATH}\n"
            "Run pipeline/c1_qc_filter.py first."
        )
    ts(f"\nLoading {INPUT_PATH}...")
    adata = sc.read_h5ad(INPUT_PATH)
    ts(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 2. Map to unified genotype IDs ----
    ts("Computing unified_genotype_id from donor map...")
    # Two-priority lookup matching step02f logic:
    #   Priority 1: exact (demux_group, donor_id) in donor_map
    #   Priority 2: donor_id is already a valid unified_id (handles
    #     G6_combined_hh13 and similar demux groups that pre-resolve
    #     donor names to genotype IDs in the genotype map)
    demux_groups = adata.obs["demux_group"].astype(str).fillna("__NA__").values
    donor_ids = adata.obs["donor_id"].astype(str).fillna("__NA__").values
    unified = np.empty(len(demux_groups), dtype=object)
    n_priority1 = 0
    n_priority2 = 0
    for i in range(len(demux_groups)):
        result = donor_map.get((demux_groups[i], donor_ids[i]))
        if result is not None:
            unified[i] = result
            n_priority1 += 1
        elif donor_ids[i] in valid_embryos:
            unified[i] = donor_ids[i]
            n_priority2 += 1
        else:
            unified[i] = np.nan
    adata.obs["unified_genotype_id"] = unified

    all_ids = set(x for x in unified if isinstance(x, str))
    unmapped = sum(1 for x in unified if not isinstance(x, str))
    ts(f"  Mapped to {len(all_ids)} unique IDs, {unmapped:,} unmapped cells")
    ts(f"  Priority 1 (donor_map): {n_priority1:,}, Priority 2 (pre-resolved): {n_priority2:,}")

    # ---- 2b. Apply barcode-level overrides for over-merged donors ----
    # In repo-fallback mode the overrides are NEVER applied, even if a stale
    # scratch copy survives: mixing the committed map with scratch overrides of
    # unknown vintage is a provenance mismatch.
    if USING_REPO_FALLBACK and os.path.isfile(OVERMERGE_OVERRIDES_PATH):
        ts("\n  WARNING: stale scratch overmerge_barcode_overrides.tsv present but "
           "SKIPPED (repo-fallback mode — its vintage cannot be tied to the committed map)")
    elif os.path.isfile(OVERMERGE_OVERRIDES_PATH):
        ts(f"\nApplying over-merge barcode overrides...")
        overrides = pd.read_csv(OVERMERGE_OVERRIDES_PATH, sep="\t")
        override_map = dict(zip(overrides["barcode"], overrides["new_unified_id"]))
        ts(f"  {len(override_map)} barcode overrides loaded")

        # Apply: for each barcode in overrides, set its unified_genotype_id
        bc_index = adata.obs.index  # barcodes are the index
        override_mask = bc_index.isin(override_map)
        n_hit = override_mask.sum()
        if n_hit > 0:
            new_ids = bc_index[override_mask].map(override_map)
            # Replace EXCLUDE with NaN
            new_ids = new_ids.where(new_ids != "EXCLUDE", np.nan)
            adata.obs.loc[override_mask, "unified_genotype_id"] = new_ids.values
            n_reclassified = (new_ids.notna()).sum()
            n_excluded = (new_ids.isna()).sum()
            ts(f"  Applied: {n_reclassified:,} reclassified, {n_excluded:,} excluded")

            # Report per split ID
            recl_counts = new_ids.dropna().value_counts()
            for sid, n in recl_counts.items():
                ts(f"    {sid}: {n:,} cells")

        # Update valid_embryos to include split IDs
        split_ids = set(overrides[overrides["new_unified_id"] != "EXCLUDE"]["new_unified_id"])
        valid_embryos = valid_embryos | split_ids
        ts(f"  Valid embryos updated: {len(valid_embryos)} (added {len(split_ids)} splits)")
    else:
        ts(f"\n  WARNING: no overmerge_barcode_overrides.tsv at {OVERMERGE_OVERRIDES_PATH}")
        ts("  Over-merged donors (e.g. the g26 splits) keep donor-level assignment only;")
        ts("  b1b writes this file — rerun b1b if the chain should carry overrides.")

    noise_ids = all_ids - valid_embryos
    if noise_ids:
        ts(f"  Removing {len(noise_ids)} IDs not in valid embryos:")
        for nid in sorted(noise_ids):
            n = (adata.obs["unified_genotype_id"] == nid).sum()
            ts(f"    {nid}: {n:,} cells")

    # ---- 3. Filter to valid embryos ----
    mask = adata.obs["unified_genotype_id"].isin(valid_embryos)
    adata_clean = adata[mask].copy()
    ts(f"  After: {adata_clean.n_obs:,} cells ({(~mask).sum():,} removed)")
    n_unique = adata_clean.obs["unified_genotype_id"].nunique()
    ts(f"  Unique IDs: {n_unique}")

    if n_unique < len(valid_embryos):
        missing = valid_embryos - set(adata_clean.obs["unified_genotype_id"].unique())
        ts(f"  WARNING: {len(missing)} embryos have 0 cells: {sorted(missing)[:10]}")

    assert adata_clean.n_obs > 0, "No cells remaining after genotype filtering!"
    n_libs = adata_clean.obs["library"].nunique()
    assert n_libs == 12, (
        f"Expected 12 libraries after genotype filtering, got {n_libs}: "
        f"{sorted(adata_clean.obs['library'].unique())}")

    # ---- 4. Save to separate output path ----
    ts(f"Saving {OUTPUT_PATH}...")
    adata_clean.write_h5ad(OUTPUT_PATH)
    ts(f"  Saved: {os.path.getsize(OUTPUT_PATH)/1e9:.2f} GB")

    ts("\n=== Phase C2 Complete ===")


if __name__ == "__main__":
    main()
