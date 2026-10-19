#!/usr/bin/env python3
"""h1b — coupled-cell POWER CENSUS for the chick snATAC ArchR build (STEP 0).

Decides the entire downstream ArchR design (pool granularity, whether to pool stages,
which stages are testable for L/R) BEFORE any multi-hour ArchR build — pure counting
from metadata already on disk (no ArchR, no SnapATAC2, runs in seconds in `duckfig`).

The duck production peak-calling (`pipeline/archr_build_peaks.R`) used
`groupBy = peak_pool = d2g_subcluster x d2g_stage`, `MIN_POOL_CELLS = 500`,
`reproducibility = "2"` (peak in >=2 library/Sample pseudobulks — ArchR-default sampleLabels="Sample";
genotype-level tested as a sensitivity check, focal result invariant, see
docs/archr_peak_reproducibility_decision.md, not distributed). The chick coupled-ATAC is
~10x smaller, so we must check which (population x stage) pools clear a reproducible-peak
floor and how many genotype replicates each side has before committing to a grouping.

Two power layers (per the snATAC STEP-2/3 ArchR build notes, internal; STEP 0):
  1. PEAK CALLING — enough cells per (population x stage) pool (pooled across A+B
     genotypes; peaks are genomic regions, group/side is irrelevant to discovery) to call
     reproducible peaks. Floor 500 (duck); a relaxed floor is the documented fallback.
  2. DIFFERENTIAL ACCESSIBILITY (L/R) — enough Group B genotype pseudobulk replicates per
     side (>=2, with enough cells each) to fit DESeq2 and detect an L-R difference.

L/R ORIENTATION IS DESIGN-ALONE AND NON-CIRCULAR (the chick's edge, and the rework of the
ported h0 scaffold whose old "by-elimination" branch used the PITX2-high lane — circular).
Per (stage x gfp) donor class we read ONLY `design/chick_groupB_lr_design.tsv` cardinality:
  - class-uniform (every design embryo in the class agrees on which lane is gex09=L/R, which
    includes singleton classes) -> FIRM orientation for every donor in the class.
  - non-uniform                 -> AMBIGUOUS (cannot orient a donor without knowing WHICH
                                   embryo it is; that needs an independent donor->embryo map,
                                   NOT the PITX2 signal we are trying to test).
No PITX2, no "elimination". donor9's prior "= HH15-WT-8" pin is treated as NOT firm (its
basis, "only GEX09-left early WT", is itself laterality-derived = potentially circular); it
is reported as a flagged secondary extension, never folded into the firm set.

Group A (ARC01/ARC02) is magnitude-only (genetically identical L/R halves -> per-embryo L/R
structurally impossible); it contributes to PEAK-CALLING power (more cells, same assembly)
but never to oriented L/R. Reported separately.

Inputs (all local, in repo):
  figures/data/chick_control/cAM_percell_populations.tsv.gz  per-cell {group,library,lane,donor,
                                                      stage,sex,leiden,PITX2_cp10k,population}
                                                      = the chick d2g-compartment analog
                                                      (LPM-lineage compendium-annotated cells)
  design/chick_groupB_lr_design.tsv                   per-embryo gex09/gex10 = L/R
  figures/data/chick_control/chick_demux_summary_n24.tsv     per-donor stage + gfp (Group B gfp join)
  pipeline/chick/_chick_subcluster_pools.py           CHICK_POOLS (inference/bilateral/excluded)

Optional refinement (if an ATAC per-cell QC sidecar is staged locally): pass --atac-qc TSV
with columns {barcode (or atac_barcode), library, n_fragment, tsse} to intersect with the
ArchR QC cut (MIN_FRAGS, MIN_TSS) and report POST-COUPLING counts. Without it, the census is
the GEX-side UPPER BOUND (ArchR's createArrowFiles QC only removes cells).

Outputs:
  (internal sidecars under docs/chick_control/, not distributed; no shipped figure reads them)
  h1b_orientation_map.tsv     per Group B donor -> design-alone side + firmness
  h1b_power_census.tsv        per (group x stage x population x side) census
  console: peak-calling-power view, DA-power view, and the recommended ArchR grouping.

  conda run -n duckfig python3 pipeline/chick/h1b_power_census.py [--atac-qc TSV] \
      [--min-frags 1000] [--min-tss 2.0] [--floor 500]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)  # for _chick_subcluster_pools / _chick_lr_orientation
from _chick_subcluster_pools import CHICK_POOLS, STAGES  # noqa: E402
from _chick_lr_orientation import (  # noqa: E402  — single source for non-circular L/R
    LR_EXCLUDED_DONORS, FIRM_FIRMNESS, build_class_orientation, donor_gfp_map, cell_side)

# the h1b_* OUTPUTS are internal chick sidecars under docs/chick_control/ (not distributed).
CHICK_FIG = os.path.join(REPO, "figures", "data", "chick_control")
CAM = os.path.join(CHICK_FIG, "cAM_percell_populations.tsv.gz")
DESIGN = os.path.join(REPO, "design", "chick_groupB_lr_design.tsv")
SUMMARY = os.path.join(CHICK_FIG, "chick_demux_summary_n24.tsv")
OUT_ORIENT = os.path.join(REPO, "docs", "chick_control", "h1b_orientation_map.tsv")
OUT_CENSUS = os.path.join(REPO, "docs", "chick_control", "h1b_power_census.tsv")

def log(m=""):
    print(m, flush=True)


def pool_role(stage, population):
    """Classify a (stage, population) into inference / bilateral / excluded / unclassified."""
    pools = CHICK_POOLS.get(stage)
    if pools is None:
        return "unclassified"
    if population in pools["inference"]:
        return "inference"
    if population in pools["bilateral"]:
        return "bilateral"
    if population in pools["excluded"]:
        return "excluded"
    return "unclassified"


def cells_per_genotype_stats(sub, donor_col="donor"):
    n = sub.groupby(donor_col).size()
    if len(n) == 0:
        return 0, 0, 0, 0
    return int(len(n)), int(n.min()), int(n.median()), int(n.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atac-qc", default=None,
                    help="optional ATAC per-cell QC TSV {barcode|atac_barcode, library, "
                         "n_fragment, tsse} to compute POST-COUPLING counts")
    ap.add_argument("--min-frags", type=int, default=1000, help="ArchR MIN_FRAGS (duck=1000)")
    ap.add_argument("--min-tss", type=float, default=2.0, help="ArchR MIN_TSS (duck=2.0)")
    ap.add_argument("--floor", type=int, default=500, help="MIN_POOL_CELLS peak floor (duck=500)")
    args = ap.parse_args()

    df = pd.read_csv(CAM, sep="\t")
    design = pd.read_csv(DESIGN, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    log(f"cAM cells (LPM-lineage compendium = chick d2g compartment): {len(df):,}")
    log(f"  group A {int((df.group=='A').sum()):,} | group B {int((df.group=='B').sum()):,}")

    # ── 1. design-alone orientation map (the h0 rework) ────────────────────────
    cls_orient = build_class_orientation(design)
    gfp_of = donor_gfp_map(summary)
    log("\n=== design-alone class orientation (NON-CIRCULAR; gex09 lane side) ===")
    for (stage, gfp), (side, firm) in sorted(cls_orient.items()):
        log(f"  {stage:8s} {gfp:4s}: gex09={side}  [{firm}]")

    b = df[df.group == "B"].copy()
    b["gfp_n"] = b["donor"].astype(str).map(gfp_of)
    # per-cell side from design cardinality (firm donors only; else '?') — shared SSOT
    sides = b.apply(lambda r: cell_side(r["stage"], r["gfp_n"], r["lane"], r["donor"], cls_orient),
                    axis=1, result_type="expand")
    b["side"] = sides[0]
    b["firmness"] = sides[1]

    orient_rows = (b.drop_duplicates("donor")[["donor", "stage", "gfp_n", "firmness"]]
                   .sort_values(["stage", "gfp_n", "donor"]))
    orient_rows = orient_rows.rename(columns={"gfp_n": "gfp"})
    orient_rows["gex09_side"] = orient_rows.apply(
        lambda r: cls_orient.get((r["stage"], r["gfp"]), ("?", ""))[0], axis=1)
    # per-donor lane cell counts for context
    lane_ct = b.groupby(["donor", "lane"]).size().unstack(fill_value=0)
    for lane in ("GEX09", "GEX10"):
        if lane not in lane_ct.columns:
            lane_ct[lane] = 0
    orient_rows["n_gex09"] = orient_rows["donor"].map(lane_ct["GEX09"])
    orient_rows["n_gex10"] = orient_rows["donor"].map(lane_ct["GEX10"])
    orient_rows.to_csv(OUT_ORIENT, sep="\t", index=False)
    log(f"\n=== Group B per-donor orientation (-> {os.path.relpath(OUT_ORIENT, REPO)}) ===")
    log(orient_rows.to_string(index=False))
    firm_donors = sorted(set(orient_rows.loc[orient_rows.firmness.isin(
        list(FIRM_FIRMNESS)), "donor"]))
    log(f"\nFIRM (non-circular) oriented donors [{len(firm_donors)}]: {firm_donors}")
    amb = sorted(set(orient_rows.loc[orient_rows.firmness == "ambiguous", "donor"]))
    log(f"AMBIGUOUS (need independent donor->embryo map) [{len(amb)}]: {amb}")

    # Group A: no per-embryo side
    a = df[df.group == "A"].copy()
    a["side"] = "whole"
    a["firmness"] = "groupA_whole_embryo"
    a["gfp_n"] = "WT"

    # ── 2. optional ATAC coupling (post-QC) ────────────────────────────────────
    atac_note = "GEX-side UPPER BOUND (no ATAC QC sidecar; ArchR QC will reduce counts)"
    if args.atac_qc and os.path.isfile(args.atac_qc):
        qc = pd.read_csv(args.atac_qc, sep="\t")
        bc_col = "atac_barcode" if "atac_barcode" in qc.columns else "barcode"
        qc = qc.rename(columns={bc_col: "atac_barcode"})
        qc["pass_atac"] = (qc["n_fragment"] >= args.min_frags) & (qc.get("tsse", 1e9) >= args.min_tss)
        # map cAM gex barcode '<raw>-1-<suffix>' -> atac '<raw>-1' + library
        def to_atac_key(row):
            bc = str(row["barcode"])
            # strip trailing lane suffix: '-1-9'/'-1-10'/'-1-1'/'-1-2' -> '-1'
            parts = bc.rsplit("-", 1)
            atac_bc = parts[0] if len(parts) == 2 else bc
            return atac_bc
        pass_set = set(zip(qc.loc[qc.pass_atac, "library"].astype(str),
                           qc.loc[qc.pass_atac, "atac_barcode"].astype(str)))
        # cAM stores Group B cells under their GEX lane (GEX09/GEX10) but the ATAC QC sidecar is
        # keyed by the ARC fragment library (ARC09/ARC10) — the same physical multiome lane. Without
        # this GEX->ARC normalisation every Group B cell silently fails the membership test and is
        # dropped (no oriented L/R power reported). Group A is already ARC01/ARC02 (identity).
        LIB_TO_ARC = {"ARC01": "ARC01", "ARC02": "ARC02", "GEX09": "ARC09", "GEX10": "ARC10"}
        for sub in (a, b):
            sub["atac_bc"] = sub["barcode"].astype(str).map(lambda s: s.rsplit("-", 1)[0])
            sub["atac_lib"] = sub["library"].astype(str).map(LIB_TO_ARC)
            if sub["atac_lib"].isna().any():
                bad = sorted(sub.loc[sub["atac_lib"].isna(), "library"].astype(str).unique())
                sys.exit(f"ERROR: cAM library not in LIB_TO_ARC (GEX->ARC map): {bad}")
            sub["coupled"] = [(lib, bc) in pass_set
                              for lib, bc in zip(sub["atac_lib"], sub["atac_bc"])]
        atac_note = (f"POST-COUPLING (ATAC QC: n_fragment>={args.min_frags}, tsse>={args.min_tss} "
                     f"from {os.path.basename(args.atac_qc)})")
        a = a[a["coupled"]]
        b = b[b["coupled"]]
        log(f"\n[coupling] retained A={len(a):,} B={len(b):,} after ATAC QC")
    log(f"\n[census basis] {atac_note}")

    both = pd.concat([a, b], ignore_index=True)

    # ── 3. PEAK-CALLING power: (stage x population), A+B pooled ─────────────────
    log(f"\n=== PEAK-CALLING power: (stage x population) pooled A+B, floor={args.floor}, "
        f"reproducibility needs >=2 genotypes ===")
    pk_rows = []
    for (stage, pop), sub in both.groupby(["stage", "population"]):
        ng, nmin, nmed, nmax = cells_per_genotype_stats(sub)
        pk_rows.append(dict(stage=stage, population=pop, role=pool_role(stage, pop),
                            n_cells=len(sub), n_genotypes=ng,
                            clears_floor=(len(sub) >= args.floor and ng >= 2)))
    pk = pd.DataFrame(pk_rows).sort_values(
        ["stage", "n_cells"], ascending=[True, False])
    log(pk.to_string(index=False))
    n_clear = int(pk.clears_floor.sum())
    n_clear_inf = int(((pk.role == "inference") & pk.clears_floor).sum())
    log(f"\n  pools clearing floor: {n_clear}/{len(pk)} "
        f"(inference pools clearing: {n_clear_inf})")

    # ── 4. DA (L/R) power: Group B firm-oriented, (stage x population x side) ───
    log(f"\n=== DA (L/R) power: Group B FIRM-oriented donors, per (stage x population x side) ===")
    bf = b[b.firmness.isin(list(FIRM_FIRMNESS))].copy()
    da_rows = []
    for (stage, pop, side), sub in bf.groupby(["stage", "population", "side"]):
        ng, nmin, nmed, nmax = cells_per_genotype_stats(sub)
        da_rows.append(dict(stage=stage, population=pop, role=pool_role(stage, pop),
                            side=side, n_cells=len(sub), n_genotypes=ng,
                            cells_per_geno_min=nmin, cells_per_geno_med=nmed,
                            cells_per_geno_max=nmax))
    da = pd.DataFrame(da_rows)
    if len(da):
        da = da.sort_values(["stage", "role", "population", "side"])
        # focus print: inference + bilateral roles only (the DA-relevant pools)
        da_focus = da[da.role.isin(["inference", "bilateral"])]
        log("  [inference + bilateral pools only]")
        log(da_focus.to_string(index=False) if len(da_focus) else "  (none)")
    else:
        log("  (no firm-oriented Group B cells)")

    # per (stage x side) genotype replicate counts for paired DESeq2 feasibility
    log("\n  paired-DESeq2 replicate check — firm-oriented genotypes per (stage x side):")
    rep = (bf.groupby(["stage", "side"])["donor"].nunique().rename("n_genotypes")
           .reset_index().sort_values(["stage", "side"]))
    log(rep.to_string(index=False))

    # ── 5. write combined census TSV ───────────────────────────────────────────
    census = both.groupby(["group", "stage", "population", "side"]).agg(
        n_cells=("barcode", "size"),
        n_genotypes=("donor", "nunique")).reset_index()
    cpg = (both.groupby(["group", "stage", "population", "side", "donor"]).size()
           .groupby(level=[0, 1, 2, 3]).agg(["min", "median", "max"]))
    cpg.columns = ["cells_per_geno_min", "cells_per_geno_med", "cells_per_geno_max"]
    census = census.merge(cpg.reset_index(), on=["group", "stage", "population", "side"])
    census["role"] = census.apply(lambda r: pool_role(r["stage"], r["population"]), axis=1)
    census = census.sort_values(["group", "stage", "role", "population", "side"])
    census.to_csv(OUT_CENSUS, sep="\t", index=False)
    log(f"\n=== wrote census -> {os.path.relpath(OUT_CENSUS, REPO)} ({len(census)} rows) ===")
    log("\n(STEP-0 census done; read the peak-calling + DA tables above to pick the ArchR grouping.)")


if __name__ == "__main__":
    main()
