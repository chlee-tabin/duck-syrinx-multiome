#!/usr/bin/env python3
"""archr_metadata_sidecar.py — barcode -> ArchR peak_pool join table (chick galGal6a).

The chick analog of the duck `step_f1_coupled_atac_metadata.tsv` consumed by
pipeline/archr_build_peaks.R.  Produced LOCALLY from already-committed inputs (cAM
per-cell populations + design + demux summary) — deterministic, seconds — then synced to
scratch for the ArchR build.  Replaces the stale scratch h0/h1 sidecar (coarse early/late
stage + the erroneous donor8a/8b split flagged in _groupB_resolved12).

Per coupled-ATAC barcode (= every cAM LPM-lineage cell; ArchR's createArrowFiles applies
its own MIN_TSS/MIN_FRAGS QC, retaining ~95% — h1 coupling check), emit:
  archr_cb            "<arc_library>#<cb_fragments>"  (ArchR cell name; build joins on this)
  cb_fragments        "<raw>-1"  (ATAC fragment barcode = cAM barcode minus the lane suffix)
  arc_library         ARC01/ARC02/ARC09/ARC10   (GEX09->ARC09, GEX10->ARC10)
  peak_pool           "<stage>|<population>"     (groupBy for addReproduciblePeakSet)
  unified_genotype_id donor (resolved-12 label)  (reproducibility="2" replicate key)
  population stage side group gfp firmness        (DA-step metadata)

cAM IS the chick d2g-compartment analog (LPM-lineage compendium-annotated cells spanning
inference + bilateral-control + excluded populations), so peak calling mirrors the duck
(all populations x stage; the 500-cell floor + the DA's inference/control selection do the
narrowing downstream — NOT this sidecar).

Group B `side` is design-alone NON-CIRCULAR (_chick_lr_orientation, the h0 rework);
Group A `side` = 'whole' (ARC01/02 genetically identical halves -> per-embryo L/R
impossible; magnitude/peak-calling only).

  conda run -n duckfig python3 pipeline/chick/archr_metadata_sidecar.py [--out PATH]
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from _groupB_split_excluded import drop_split_excluded  # noqa: E402
from _chick_lr_orientation import (  # noqa: E402
    build_class_orientation, donor_gfp_map, cell_side, FIRM_FIRMNESS)

CHICK = os.path.join(REPO, "figures", "data", "chick_control")
CAM = os.path.join(CHICK, "cAM_percell_populations.tsv.gz")
# cAL = the chick d2g subcluster SSOT (per-cluster (stage, cluster) -> population); the
# rename-proof source for the per-cell population label (see main()).
CAL = os.path.join(CHICK, "cAL_subcluster_annotations.tsv")
DESIGN = os.path.join(REPO, "design", "chick_groupB_lr_design.tsv")
SUMMARY = os.path.join(CHICK, "chick_demux_summary_n24.tsv")
DEFAULT_OUT = os.path.join(REPO, "figures", "data", "chick_archr_coupled_metadata.tsv")

# cAM library -> ARC fragment library (multiome: GEX09 lane == ARC09 lane).
LIB_TO_ARC = {"ARC01": "ARC01", "ARC02": "ARC02", "GEX09": "ARC09", "GEX10": "ARC10"}
MIN_POOL_CELLS = 500  # mirrors duck archr_build_peaks.R floor (for the sanity report)


def log(m=""):
    print(m, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    df = pd.read_csv(CAM, sep="\t")
    # Group B cS-split doublet/unassigned cells still in the pre-fix cAM -> never an ArchR replicate
    df = drop_split_excluded(df, "cAM_percell_populations", log=log)

    # ── Re-map population through the CHICK SSOT (cAL) by leiden — never trust the
    #    baked cAM name ─────────────────────────────────────────────────────────────
    # cAM carries a per-cell `population` AND its source `leiden` id.  The baked
    # `population` string goes stale if cAL is re-labelled without regenerating cAM, so
    # the rename-proof source is cAL_subcluster_annotations.tsv (the committed chick d2g
    # SSOT), keyed by (stage, cluster).  This is the chick analog of the duck g1 /
    # duck-locus leiden→get_domain_map re-map.  cAM.population == cAL
    # for every (stage, leiden) today (verified), so this is label-identical now and only
    # future-proofs against a rename.  Fail loud on an unannotated (stage, leiden).
    cal = pd.read_csv(CAL, sep="\t")
    cal_map = {(str(r.stage), int(r.cluster)): r.population for r in cal.itertuples()}

    def _cal_pop(stage, leiden):
        key = (str(stage), int(leiden))
        if key not in cal_map:
            sys.exit(f"ERROR: cAM cell (stage={stage!r}, leiden={leiden}) has no cAL "
                     f"annotation — regenerate cAL/cAM or fix the chick SSOT.")
        return cal_map[key]

    _baked = df["population"].astype("object")
    df["population"] = [_cal_pop(s, l) for s, l in zip(df["stage"], df["leiden"])]
    _n_changed = int((_baked.astype(str) != df["population"].astype(str)).sum())
    log(f"cAL re-map: {len(df):,} cells | population labels changed vs baked cAM: "
        f"{_n_changed:,}")

    design = pd.read_csv(DESIGN, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")

    cls_orient = build_class_orientation(design)
    gfp_of = donor_gfp_map(summary)

    df["arc_library"] = df["library"].map(LIB_TO_ARC)
    if df["arc_library"].isna().any():
        bad = sorted(df.loc[df["arc_library"].isna(), "library"].unique())
        sys.exit(f"ERROR: cAM libraries not in LIB_TO_ARC: {bad}")
    # cAM barcode '<raw>-1-<laneSuffix>' -> ATAC fragment barcode '<raw>-1'.
    df["cb_fragments"] = df["barcode"].astype(str).str.rsplit("-", n=1).str[0]
    df["archr_cb"] = df["arc_library"] + "#" + df["cb_fragments"]
    df["peak_pool"] = df["stage"].astype(str) + "|" + df["population"].astype(str)
    # chick donor IDs are REUSED across groups (Group-A donor10 != Group-B donor10; colliding set
    # donor4/6/8/9/10), so the ArchR replicate / pseudobulk key is GROUP-QUALIFIED. Belt-and-
    # suspenders. The prefix also makes the peak-discovery n_replicates count correct (donor10 is
    # HH15-16 in BOTH groups).
    # LIVE: the committed figures/data/chick_locus/* have been
    # regenerated and now carry "A:"/"B:" ids. Scope of what actually changed:
    #   - locus_pseudobulk_*: groupings byte-identical with or without the prefix — that path is
    #     Group-B-only (archr_pitx2_locus_pseudobulk.R keeps side in {L,R}, where the bare donor
    #     is already unique), so only the genotype LABEL changed ("donor4" -> "B:donor4").
    #   - locus_percell_counts.tsv.gz: A+B MIXED, so its `genotype` GROUPING KEY genuinely
    #     changed — 19 bare ids had been conflating 24 embryos. Consumers joining a bare donor
    #     must normalise with .split(":")[-1] (see build_isoform_data.py and the exploratory
    #     isoform_probe/{chick_atac_isoform_promoters,r4_hh25_bootstrap}.py).
    df["unified_genotype_id"] = df["group"].astype(str) + ":" + df["donor"].astype(str)

    # side / firmness
    df["gfp"] = df["donor"].astype(str).map(gfp_of)
    is_b = df["group"] == "B"
    sides = df[is_b].apply(
        lambda r: cell_side(r["stage"], r["gfp"], r["lane"], r["donor"], cls_orient),
        axis=1, result_type="expand")
    df["side"] = "whole"
    df["firmness"] = "groupA_whole_embryo"
    df.loc[is_b, "side"] = sides[0].values
    df.loc[is_b, "firmness"] = sides[1].values
    df.loc[df["group"] == "A", "gfp"] = "WT"
    # Group A is whole-embryo by design (genetically identical halves, no per-embryo L/R), so it
    # must never acquire an L/R side. Downstream relies on this: the oriented locus pseudobulk
    # filters side in {L,R} and is therefore Group-B-only, and cAN's directional control assumes
    # the same. A Group-A cell with a side would silently enter an oriented pseudobulk.
    assert (df.loc[df["group"] == "A", "side"] == "whole").all(), \
        "Group A cell assigned an L/R side -> it would leak into the Group-B oriented pseudobulk"

    out_cols = ["archr_cb", "cb_fragments", "arc_library", "peak_pool",
                "unified_genotype_id", "population", "stage", "side", "group",
                "gfp", "firmness"]
    out = df[out_cols].copy()
    # ArchR cell barcodes must be unique; a cAM cell appears once, so archr_cb is unique.
    dup = out["archr_cb"].duplicated().sum()
    if dup:
        sys.exit(f"ERROR: {dup} duplicate archr_cb — fragment-barcode collision across libs?")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    log(f"wrote {len(out):,} barcodes -> {os.path.relpath(args.out, REPO)}")

    # ── sanity report: peak_pool counts vs the 500 floor (A+B pooled) ──────────
    pc = (out.groupby("peak_pool")
          .agg(n_cells=("archr_cb", "size"),
               n_genotypes=("unified_genotype_id", "nunique")).reset_index())
    pc["clears_floor"] = (pc.n_cells >= MIN_POOL_CELLS) & (pc.n_genotypes >= 2)
    pc = pc.sort_values("n_cells", ascending=False)
    log(f"\npeak_pool counts (A+B pooled; floor {MIN_POOL_CELLS} cells & >=2 genotypes):")
    log(pc.to_string(index=False))
    log(f"\n  pools clearing floor (will call peaks): {int(pc.clears_floor.sum())}/{len(pc)}")
    log(f"  total barcodes: {len(out):,}  | groups: {out.group.value_counts().to_dict()}")
    log(f"  Group B firm-oriented (L/R) cells: "
        f"{int((is_b & df.firmness.isin(FIRM_FIRMNESS)).sum()):,}")
    log("\n(sync this TSV to scratch + run pipeline/chick/archr_build_peaks.R)")


if __name__ == "__main__":
    main()
