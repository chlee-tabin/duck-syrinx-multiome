#!/usr/bin/env python3
"""export_laterality_trajectory.py — per-(stage, population) embryo-paired PITX2 L/R
across the d2g stage-restricted subclusters, for the "laterality trajectory" panel.

The biological message this exports the data for: the lateralized LPM-lineage / inference
populations (Ventral spl., Syringeal, Chondrogenic LPM, BARX1+, Early/Cardiopharyngeal LPM)
stay LEFT-biased across stages, the dorsal control stays bilateral, and **Skeletogenic LPM
(the FOXF1-low pharyngeal-cranial skeletal / pre-chondrogenic progenitor) is LEFT-biased
early (HH17/16) then resolves to bilateral by the fate stages (HH21/HH25)** — i.e. it
shares the early left-LPM program then "fizzles", while the lateralized syrinx skeletal
program continues in Chondrogenic LPM.

Self-contained (reads only committed figures/data + the SSOT); no dependence on the
exploratory old-duck vet scripts.

SIDE SOURCE = the canonical merged sidecar `figures/data/lr_inference_multiome.tsv`
(d2o: CHROMATIN-based d2n at HH17/16, snRNA d2l at HH13/HH21/HH25) — NOT raw d2l.  This
matters: at HH17/16 the snRNA d2l side is weak and fully circular (it is defined from PITX2
RNA).  The chromatin d2n side is a different modality but NOT independent of the RNA PITX2
being tested: it labels the higher-PITX2c-accessibility half "L", and PITX2c accessibility
and PITX2 RNA covary in the same halves, so the HH17/16 rows are partly circular (corrected
2026-10-14; this said "non-circular" before; see figS03/PROVENANCE.md, Side basis, and the
legend (C)).  The flat dorsal control, computed identically, checks that selection; it does
not prove it absent.

Pool roles (inference / bilateral / excluded) come from the canonical
`figures/_shared/subcluster_map.py::POOLS` — never hardcoded.

Cohort rules (2026-10-14):
  * Stage bins follow the sidecar's `stage_anchor` (HH16 / HH17 genotypes fold into
    HH17/16), not its literal `stage` string — g01, g08, g20 (stage "HH16"/"HH17")
    previously fell out of the HH17/16 bin.
  * Sides go through the shared firmness rule (`figures/_shared/lr_firmness.py`,
    paired=True): an unknown firmness value raises; a same-side duplicate (lane-loading-record
    split tube: one half in two libraries) has no opposite half and is excluded
    by name.
  * Every embryo that does not reach a row is logged BY NAME with its reason (no side
    row at this stage, one side only, < MIN_CELLS cells on a side, no PITX2 reads for
    the c-usage metric); a leiden cluster missing from the annotation SSOT raises.
  * Each (stage, population, metric) draws from its OWN seeded generator (seed =
    [0, crc32("stage|population|metric")]), so the permutation p and bootstrap CI of
    one cell no longer depend on how many draws an earlier cell consumed (the old
    single stream let an HH17/16 cohort change shift every HH21/HH25 CI).

Output: figures/data/laterality_trajectory.tsv  (per stage x population x metric).
"""
from __future__ import annotations
import sys
import zlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                      # figures/figS03 -> repo root
DATA = REPO / "figures" / "data"
ANN = REPO / "references" / "d2g_subcluster_annotations.tsv"
OUT = DATA

sys.path.insert(0, str(REPO / "figures"))
from _shared.subcluster_map import POOLS    # noqa: E402  (canonical pool membership)
from _shared.library_naming import collapse_lib  # noqa: E402
from _shared.lr_firmness import resolve_sides  # noqa: E402

MIN_CELLS = 10
N_PERM = 2000
N_BOOT = 2000   # bootstrap resamples for the median-log2FC 95% CI (figS03 forest whiskers)
SEED = 0
# side classes trusted here: every inferred call; a genotype with one library
# (singleton_rank) has no opposite half and is excluded by name below
ACCEPT = {"direct", "placement", "singleton_rank"}

# (display stage, d2g file label, leiden col) — the display stage IS the sidecar's
# stage_anchor
STAGES = [
    ("HH13",    "HH13_ccreg",        "leiden_0.2"),
    ("HH17/16", "HH1716_HH16_ccreg", "leiden_0.3"),
    ("HH21",    "HH21_ccreg",        "leiden_0.3"),
    ("HH25",    "HH25_ccreg",        "leiden_0.3"),
]


def cell_rng(stage: str, pop: str, metric: str) -> np.random.Generator:
    """One generator per (stage, population, metric): reproducible and independent of
    every other cell's draws (crc32 is stable across processes, unlike hash())."""
    return np.random.default_rng([SEED, zlib.crc32(f"{stage}|{pop}|{metric}".encode())])


def load_side_maps() -> tuple[dict, dict]:
    """(stage_anchor -> genotype -> {library_m: side}, genotype -> exclusion reason),
    through the shared firmness rule (paired: same-side duplicates excluded, logged)."""
    lr = pd.read_csv(DATA / "lr_inference_multiome.tsv", sep="\t")
    res = resolve_sides(lr, ACCEPT, source="lr_inference_multiome.tsv", paired=True)
    kept = set(zip(res["genotype"], res["library_m"]))
    why: dict = {}
    for r in lr.itertuples():
        if (r.genotype, r.library_m) not in kept:
            firm = str(r.firmness).removeprefix("chromatin_")
            why[r.genotype] = ("same-side duplicate (one half in two libraries, no opposite half)"
                               if firm == "same_side_dup" else f"firmness {firm} not accepted")
    anchor = lr.set_index(["genotype", "library_m"])["stage_anchor"]
    res["stage_anchor"] = [anchor[(g, l)] for g, l in zip(res["genotype"], res["library_m"])]
    unknown = sorted(set(res["stage_anchor"]) - {s for s, _, _ in STAGES})
    if unknown:
        raise ValueError(f"lr_inference_multiome.tsv stage_anchor outside {[s for s, _, _ in STAGES]}: {unknown}")
    maps: dict = {}
    for r in res.itertuples():
        maps.setdefault(r.stage_anchor, {}).setdefault(r.genotype, {})[r.library_m] = r.inferred_side
    for st, _, _ in STAGES:
        if st not in maps:
            raise ValueError(f"no usable lr_inference_multiome rows for stage_anchor {st}")
    return maps, why


def role_of(stage: str, pop: str) -> str:
    """Pool role of `pop` at `stage` from the POOLS SSOT.  A population in none of
    the three lists (a new or renamed compendium label) raises instead of being
    filed as 'excluded'."""
    if stage not in POOLS:
        raise KeyError(f"role_of: stage {stage!r} not in POOLS ({sorted(POOLS)})")
    p = POOLS[stage]
    roles = [r for r in ("inference", "bilateral", "excluded") if pop in p[r]]
    if len(roles) != 1:
        raise KeyError(f"role_of: {pop!r} at {stage} is in POOLS roles {roles} "
                       "(expected exactly one) — stale POOLS or a renamed population")
    return roles[0]


def per_embryo_LR(cells: pd.DataFrame, side_map: dict, why: dict, excl: list) -> pd.DataFrame:
    """Per genotype with >= MIN_CELLS cells on BOTH sides: mean PITX2 and PITX2c usage
    per half.  Every genotype that does not qualify is appended to `excl` with its reason."""
    rows = []
    for gt, g in cells.groupby("genotype"):
        if gt not in side_map:
            excl.append((gt, why.get(gt, "no sidecar row at this stage_anchor"), len(g)))
            continue
        sm = side_map[gt]
        g = g.assign(side=g["library"].map(collapse_lib).map(sm))
        unsided = g[g["side"].isna()]
        if len(unsided):
            libs = ", ".join(f"{l} ({n})" for l, n in unsided["library"].value_counts().sort_index().items())
            excl.append((gt, f"cells in a library without a side row dropped: {libs}", len(unsided)))
        L = g[g["side"] == "L"]; R = g[g["side"] == "R"]
        if len(L) == 0 or len(R) == 0:
            excl.append((gt, f"one side only (L {len(L)}, R {len(R)} cells)", len(g)))
            continue
        if len(L) < MIN_CELLS or len(R) < MIN_CELLS:
            excl.append((gt, f"< {MIN_CELLS} cells on a side (L {len(L)}, R {len(R)})", len(g)))
            continue

        def c_usage(x):
            # d2j_pitx2_isoform_percell is sparse (cells with >= 1 PITX2 read); a half
            # with no isoform-informative read has no c-usage (NaN, logged at the metric)
            tot = (x["pitx2c_umi"] + x["pitx2ab_umi"] + x["shared_umi"]).sum()
            return float(x["pitx2c_umi"].sum() / tot) if tot > 0 else np.nan

        rows.append({
            "genotype": gt,
            "pitx2_L": float(L["PITX2_cp10k"].mean()), "pitx2_R": float(R["PITX2_cp10k"].mean()),
            "cuse_L": c_usage(L), "cuse_R": c_usage(R),
        })
    return pd.DataFrame(rows, columns=["genotype", "pitx2_L", "pitx2_R", "cuse_L", "cuse_R"])


def lr_stats(df: pd.DataFrame, lcol: str, rcol: str, rng: np.random.Generator,
             excl: list) -> dict:
    nan = df[df[lcol].isna() | df[rcol].isna()]
    for gt in nan["genotype"]:
        excl.append((gt, f"{lcol[:-2]}: no PITX2 isoform-informative read in a half", 0))
    d = df.dropna(subset=[lcol, rcol])
    n = len(d)
    if n == 0:
        return {"n_emb": 0}
    L = d[lcol].values; R = d[rcol].values
    nL = int((L > R).sum()); nR = int((R > L).sum())
    diff = L - R
    nt = nL + nR
    sign_p = stats.binomtest(nL, nt, 0.5).pvalue if nt > 0 else np.nan
    # label-shuffle null: swap L/R WITHIN each embryo (one flip per embryo, applied to both
    # sides), count how often #(L>R) >= observed.  flip must be drawn ONCE and reused for Ls/Rs.
    obs = nL
    ge = 0
    for _ in range(N_PERM):
        flip = rng.integers(0, 2, n).astype(bool)
        Ls = np.where(flip, R, L); Rs = np.where(flip, L, R)
        if int((Ls > Rs).sum()) >= obs:
            ge += 1
    perm_p = (ge + 1) / (N_PERM + 1)
    # eps pseudocount avoids div-by-zero when a half has ~0 PITX2.  NOTE: for R≈0
    # populations (early stages, e.g. Cardiopharyngeal HH13) the log2FC MAGNITUDE is
    # pseudocount-dominated, not biological — the figS03 panel D caps its axis at +5
    # and marks those off-scale, and the scale-free sign_p / perm_p (not the magnitude)
    # is the robust significance basis.  The DIRECTION (sign of L−R) is eps-robust;
    # only the magnitude is eps-sensitive.
    eps = 1e-3
    lfc = np.log2((L + eps) / (R + eps))             # per-embryo log2(L/R)
    med_lfc = float(np.median(lfc))
    # bootstrap 95% CI of the median log2FC (resample embryos w/ replacement) — the
    # figS03 forest whisker.  Degenerate (lo=hi=median) when n<2 or all embryos agree.
    if n >= 2:
        bt = np.array([np.median(rng.choice(lfc, n, replace=True)) for _ in range(N_BOOT)])
        lfc_lo, lfc_hi = float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))
    else:
        lfc_lo = lfc_hi = med_lfc
    return {
        "n_emb": n, "n_L_gt_R": nL, "n_R_gt_L": nR,
        "frac_L_gt_R": round(nL / n, 4),
        "median_L_minus_R": round(float(np.median(diff)), 4),
        "median_log2FC": round(med_lfc, 3),
        "log2fc_ci_lo": round(lfc_lo, 3), "log2fc_ci_hi": round(lfc_hi, 3),
        "sign_p": round(float(sign_p), 4) if sign_p == sign_p else None,
        "perm_p": round(float(perm_p), 4),
    }


def main():
    ann = pd.read_csv(ANN, sep="\t")
    iso = pd.read_csv(DATA / "d2j_pitx2_isoform_percell.tsv.gz", sep="\t")
    if iso.duplicated(["barcode", "library"]).any():
        raise ValueError("d2j_pitx2_isoform_percell.tsv.gz has duplicate (barcode, library) rows")
    side_maps, why = load_side_maps()
    out_rows = []
    excl_rows = []
    for stage, label, leiden_col in STAGES:
        df = pd.read_csv(DATA / f"d2g_lpm_subclustering_{label}.tsv.gz", sep="\t")
        a = ann[ann["stage"] == stage]
        cl2pop = {int(c): p for c, p in zip(a["cluster"], a["population"])}
        unmapped = sorted(set(df[leiden_col].astype(int)) - set(cl2pop))
        if unmapped:
            raise KeyError(f"{stage}: {leiden_col} cluster(s) {unmapped} missing from "
                           f"{ANN.name} — stale annotation or wrong resolution")
        df["population"] = df[leiden_col].astype(int).map(cl2pop)
        # sparse isoform table: a cell absent from it has zero PITX2 isoform UMIs
        df = df.merge(iso, on=["barcode", "library"], how="left", validate="one_to_one")
        for c in ["pitx2c_umi", "pitx2ab_umi", "shared_umi"]:
            df[c] = df[c].fillna(0)
        side_map = side_maps[stage]
        for pop in sorted(df["population"].unique()):
            cells = df[df["population"] == pop]
            excl: list = []
            pe = per_embryo_LR(cells, side_map, why, excl)
            for metric, (lc, rc) in [("PITX2_mean", ("pitx2_L", "pitx2_R")),
                                     ("PITX2c_usage", ("cuse_L", "cuse_R"))]:
                mexcl: list = []
                st = lr_stats(pe, lc, rc, cell_rng(stage, pop, metric), mexcl)
                out_rows.append({"stage": stage, "population": pop,
                                 "role": role_of(stage, pop), "metric": metric,
                                 "n_cells": int(len(cells)), **st})
                excl_rows += [(stage, pop, metric, gt, why) for gt, why, _ in mexcl]
            excl_rows += [(stage, pop, "both", gt, why) for gt, why, _ in excl]
    out = pd.DataFrame(out_rows)
    dest = OUT / "laterality_trajectory.tsv"
    out.to_csv(dest, sep="\t", index=False)
    print(f"wrote {dest}  ({len(out)} rows)")
    ex = pd.DataFrame(excl_rows, columns=["stage", "population", "metric", "genotype", "reason"])
    print(f"\nnamed exclusions ({len(ex)} genotype x population x metric records; "
          f"embryos that clear every rule are the n_emb of each row):")
    for (st, pop), grp in ex.groupby(["stage", "population"], sort=False):
        print(f"  {st:8s} {pop}:")
        for r in grp.itertuples():
            print(f"      {r.genotype:22s} [{r.metric}] {r.reason}")
    # quick console view of the PITX2_mean trajectory (the panel's primary metric)
    pm = out[out.metric == "PITX2_mean"]
    print(pm[["stage", "population", "role", "n_emb", "n_L_gt_R", "n_R_gt_L",
              "frac_L_gt_R", "median_L_minus_R", "sign_p", "perm_p"]].to_string(index=False))
    print("\nLATERALITY_TRAJECTORY_EXPORT_OK")


if __name__ == "__main__":
    main()
