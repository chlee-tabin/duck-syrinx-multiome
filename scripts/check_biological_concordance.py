#!/usr/bin/env python3
"""Per-cluster biological L/R sign concordance (replaces alphabetical concordance).

The compendium scorer's `lat_n_concordant_sign` counts embryos whose
within-pair PITX2c UMI ordering matches an alphabetical library reference
(libA > libB). Under a balanced design (half genotypes have L=libA, half
have L=libB), alphabetical concordance expects 50% by chance — so high
values indicate library-ordering artefacts, not biology.

This script joins per-embryo PITX2c UMIs with the 4-metric `inferred_side`
from `figures/data/d2l_lr_inference.tsv` and reports **biological**
concordance:

    sign = +1 if L-side UMIs > R-side UMIs
           -1 if R > L
            0 if tie / unresolved

Expected behaviour:
    inference-pool clusters → near-unanimous L > R (all PITX2c strong
        embryos on the left)
    bilateral controls       → ~50/50 random split (true bilateral)
    excluded Skeletogenic LPM → mixed (no meaningful signal)

Output: `figures/data/d2l_cluster_biological_concordance.tsv` — one row
per (stage, leiden_col, cluster).

Unit and join (2026-10):
    * the per-embryo unit is the GENOTYPE (genotype = distinct embryo; the
      genotype -> embryo SSOT is d2m).  The legacy embryo key from
      b_genotype_embryo_mapping.tsv used to drop every genotype it does not
      list (all G3_unfiltered_donorN, G2_genotype_01, g26_split0/1) without a
      message.
    * sides are joined on d2l's ``stage_anchor`` (the d2g stage table the row
      came from), not on its per-genotype e1b ``stage`` -- inside the merged
      HH17/16 bin d2l labels g01/g08 "HH16" and g20 "HH17", which a literal
      ``stage == "HH17/16"`` join silently left without a side.
    * d2l writes one row per (genotype, library_m) with >= 1 inference-pool
      cell, so a (genotype, library_m) without a d2l row is allowed ONLY when it
      has no inference-pool cell at that stage (logged by name; it is an
      unresolved library, like '?').  A missing row for a library WITH
      inference-pool cells, or a d2l row with no d2g cells, means the d2l
      sidecar and the d2g tables are out of sync and raises.

Schema:
    stage, leiden_col, cluster, population,
    n_embryos_paired        embryos with 2 distinct libraries present
                            (a same-side duplicate -- firmness same_side_dup,
                            lane-loading-record group 17-B(2) = embryos 17-5/17-6: ONE left half in GEX04
                            and GEX0S-duck -- is not a paired embryo and is
                            excluded up front, logged)
    n_sides_resolved        paired embryos where both libraries have
                            non-`?` inferred_side (one L + one R)
    n_both_nonzero          n_sides_resolved AND both libraries have
                            >0 PITX2c UMIs in this cluster
    n_L_dom, n_R_dom        L>R and R>L embryo counts (of n_both_nonzero)
    n_strong_L, n_strong_R  L/R with |log2 fold| >= 1
    n_bio_concordant        max(n_L_dom, n_R_dom)
    bio_direction           'L' | 'R' (whichever side dominates; '-' if tied)
    median_log2_fold        median log2((L+0.5)/(R+0.5)) across paired embryos
                            (signed — negative = R dominant)
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "figures" / "data"
REFS = REPO / "references"
OUT_PATH = DATA / "d2l_cluster_biological_concordance.tsv"

sys.path.insert(0, str(REPO / "figures"))
from _shared.lr_firmness import normalize_firmness, same_side_duplicates  # noqa: E402
from _shared.subcluster_map import POOLS, get_pool_clusters  # noqa: E402

# Match the compendium scorer's tech-rep merge (otherwise d2j library strings
# like "GEX03-1A" won't join to inference-TSV library_m like "GEX03-1")
TECH_REP_MERGE = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1",
                  "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}

# Canonical (stage, resolution) — must match pipeline/d2g_validate_subclusters.py
CANONICAL = [
    ("HH13",    0.2, "HH13_ccreg"),
    ("HH17/16", 0.3, "HH1716_HH16_ccreg"),
    ("HH21",    0.3, "HH21_ccreg"),
    ("HH25",    0.3, "HH25_ccreg"),
]


def load_inference_sides() -> pd.DataFrame:
    inf = pd.read_csv(DATA / "d2l_lr_inference.tsv", sep="\t")
    keep = ["genotype", "library_m", "stage_anchor", "inferred_side", "firmness"]
    inf = inf[keep].copy()
    inf["firmness"] = normalize_firmness(inf["firmness"], "d2l_lr_inference.tsv")
    return inf


def same_side_duplicate_genotypes(inf: pd.DataFrame) -> set:
    """Genotypes whose two library rows hold the SAME half (same_side_dup).
    This script contrasts an embryo's L half against its R half, so they are
    excluded, by name (figures/_shared/lr_firmness.py rule: never one half of a
    paired L/R contrast)."""
    dups = same_side_duplicates(inf, "d2l_lr_inference.tsv")
    for g, (libs, side) in sorted(dups.items()):
        print(f"  excluded (same-side duplicate, one {side} half in both "
              f"{'+'.join(sorted(libs))}): {g}")
    return set(dups)


def check_side_coverage(d2g: pd.DataFrame, inf_stage: pd.DataFrame,
                        leiden_col: str, stage: str) -> None:
    """d2l rows vs d2g cells for one stage bin (rule in the module docstring).
    Raises on a sync failure; logs each library that has cells but no side."""
    cells = d2g[["genotype", "library", leiden_col]].copy()
    cells["library_m"] = cells["library"].map(lambda l: TECH_REP_MERGE.get(l, l))
    inf_clusters = set(get_pool_clusters(stage, "inference"))
    cells["inference"] = cells[leiden_col].astype(int).isin(inf_clusters)
    per_lib = (cells.groupby(["genotype", "library_m"])["inference"]
               .agg(n_cells="size", n_inference="sum").reset_index())
    d2l_keys = set(zip(inf_stage["genotype"], inf_stage["library_m"]))
    d2g_keys = set(zip(per_lib["genotype"], per_lib["library_m"]))
    no_cells = sorted(d2l_keys - d2g_keys)
    if no_cells:
        raise ValueError(f"{stage}: d2l_lr_inference.tsv rows with no cells in the d2g "
                         f"table (stale sidecar?): {no_cells}")
    no_row = per_lib[[k not in d2l_keys for k in zip(per_lib["genotype"],
                                                     per_lib["library_m"])]]
    stale = no_row[no_row["n_inference"] > 0]
    if len(stale):
        raise ValueError(f"{stage}: (genotype, library_m) with inference-pool cells but "
                         f"no d2l_lr_inference.tsv row (stale sidecar?): "
                         f"{stale.to_dict('records')}")
    for r in no_row.itertuples():
        print(f"  {stage}: no d2l side for {r.genotype} x {r.library_m} "
              f"({r.n_cells} cells, none in the inference pool) -> unresolved library")


def load_annotations() -> pd.DataFrame:
    return pd.read_csv(REFS / "d2g_subcluster_annotations.tsv", sep="\t")


def per_cluster_concordance(
    d2g_sub: pd.DataFrame,
    iso: pd.DataFrame,
    inf: pd.DataFrame,
    leiden_col: str,
    stage: str,
    cluster: int,
    dup_genotypes: set,
) -> dict:
    """Compute biological L/R concordance for one (stage, resolution, cluster)."""
    sub = d2g_sub[d2g_sub[leiden_col].astype(int) == cluster].copy()
    sub = sub[~sub["genotype"].isin(dup_genotypes)]
    if sub.empty:
        return _empty_row()

    # Attach PITX2c UMIs (barcode × raw library)
    sub = sub.merge(
        iso[["barcode", "library", "pitx2c_umi"]],
        on=["barcode", "library"], how="left",
    )
    sub["pitx2c_umi"] = sub["pitx2c_umi"].fillna(0)

    # Tech-rep merge to align with library_m used in inference TSV
    sub["library_m"] = sub["library"].map(lambda l: TECH_REP_MERGE.get(l, l))

    # Join inferred_side per (genotype, library_m) within this stage bin
    # (d2l stage_anchor; coverage checked by check_side_coverage).  A library
    # left without a side is an unresolved library (no inference-pool cells).
    inf_stage = inf[inf["stage_anchor"] == stage]
    if inf_stage.empty:
        raise ValueError(f"d2l_lr_inference.tsv has no rows with stage_anchor == {stage!r}")
    side_map = inf_stage[["genotype", "library_m", "inferred_side"]]
    if side_map.duplicated(["genotype", "library_m"]).any():
        raise ValueError(f"{stage}: duplicate (genotype, library_m) rows in "
                         f"d2l_lr_inference.tsv")
    sub = sub.merge(side_map, on=["genotype", "library_m"], how="left")

    # Per-embryo aggregation (one genotype = one embryo)
    n_paired = 0
    n_resolved = 0
    n_both_nz = 0
    l_wins = r_wins = 0
    strong_L = strong_R = 0
    log2_folds = []
    for _, esub in sub.groupby("genotype"):
        libs = esub["library_m"].dropna().unique()
        if len(libs) != 2:
            continue
        n_paired += 1
        sides = esub.groupby("library_m")["inferred_side"].first().to_dict()
        side_vals = set(sides.values())
        if side_vals != {"L", "R"}:
            # Either one lib is '?' or both are the same label (shouldn't happen
            # with pair-constraint but tolerate)
            continue
        n_resolved += 1
        l_lib = next(k for k, v in sides.items() if v == "L")
        r_lib = next(k for k, v in sides.items() if v == "R")
        l_umi = esub.loc[esub["library_m"] == l_lib, "pitx2c_umi"].sum()
        r_umi = esub.loc[esub["library_m"] == r_lib, "pitx2c_umi"].sum()
        if l_umi == 0 and r_umi == 0:
            continue
        n_both_nz += 1
        log2f = float(np.log2((l_umi + 0.5) / (r_umi + 0.5)))
        log2_folds.append(log2f)
        if log2f > 0:
            l_wins += 1
            if log2f >= 1:
                strong_L += 1
        elif log2f < 0:
            r_wins += 1
            if log2f <= -1:
                strong_R += 1
        # exact tie log2f == 0 → neither side wins

    if n_both_nz == 0:
        return _empty_row(n_paired=n_paired, n_resolved=n_resolved)

    concordant = max(l_wins, r_wins)
    if l_wins > r_wins:
        direction = "L"
    elif r_wins > l_wins:
        direction = "R"
    else:
        direction = "-"
    return {
        "n_embryos_paired": n_paired,
        "n_sides_resolved": n_resolved,
        "n_both_nonzero":   n_both_nz,
        "n_L_dom":          l_wins,
        "n_R_dom":          r_wins,
        "n_strong_L":       strong_L,
        "n_strong_R":       strong_R,
        "n_bio_concordant": concordant,
        "bio_direction":    direction,
        "median_log2_fold": round(float(np.median(log2_folds)), 3),
    }


def _empty_row(n_paired: int = 0, n_resolved: int = 0) -> dict:
    return {
        "n_embryos_paired": n_paired,
        "n_sides_resolved": n_resolved,
        "n_both_nonzero":   0,
        "n_L_dom":          0,
        "n_R_dom":          0,
        "n_strong_L":       0,
        "n_strong_R":       0,
        "n_bio_concordant": 0,
        "bio_direction":    "-",
        "median_log2_fold": float("nan"),
    }


def main() -> None:
    ann = load_annotations()
    iso = pd.read_csv(DATA / "d2j_pitx2_isoform_percell.tsv.gz", sep="\t")
    inf = load_inference_sides()
    dup_genotypes = same_side_duplicate_genotypes(inf)

    # Preflight: how many embryos have both sides resolved per stage?
    resolved = inf[inf["inferred_side"].isin(["L", "R"])]
    print("inference coverage:")
    for stage, _, _ in CANONICAL:
        s = resolved[resolved["stage_anchor"] == stage]
        embs = s.groupby("genotype")["inferred_side"].nunique()
        n_both = (embs == 2).sum()
        print(f"  {stage:8s}  {n_both:3d} genotypes with L+R resolved "
              f"(of {s['genotype'].nunique():3d} with ≥1 resolved side)")

    rows = []
    for stage, res, label in CANONICAL:
        leiden_col = f"leiden_{res}"
        d2g_path = DATA / f"d2g_lpm_subclustering_{label}.tsv.gz"
        if not d2g_path.exists():   # never drop a stage from the table
            raise FileNotFoundError(f"{stage}: missing {d2g_path}")
        d2g = pd.read_csv(d2g_path, sep="\t")
        check_side_coverage(d2g, inf[inf["stage_anchor"] == stage], leiden_col, stage)

        ann_stage = ann[ann["stage"] == stage].sort_values("cluster")
        for _, ar in ann_stage.iterrows():
            cl = int(ar["cluster"])
            row = per_cluster_concordance(
                d2g, iso, inf, leiden_col, stage, cl, dup_genotypes,
            )
            row.update({
                "stage":       stage,
                "leiden_col":  leiden_col,
                "cluster":     cl,
                "population":  ar["population"],
                "n_cells":     int(ar["n_cells"]),
            })
            rows.append(row)

    out = pd.DataFrame(rows)
    cols = ["stage", "leiden_col", "cluster", "population", "n_cells",
            "n_embryos_paired", "n_sides_resolved", "n_both_nonzero",
            "n_L_dom", "n_R_dom", "n_strong_L", "n_strong_R",
            "n_bio_concordant", "bio_direction", "median_log2_fold"]
    out = out[cols]
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nWrote {OUT_PATH} ({len(out)} rows)")

    # Summarise by pool class for a quick sanity check at the console
    print("\nBiological direction by pool (inference = laterality-showing "
          "LPM-lineage):")
    def _class(stage, pop):
        p = POOLS[stage]
        if pop in p["inference"]: return "inference"
        if pop in p["bilateral"]: return "bilateral"
        return "excluded"
    out["klass"] = [_class(s, p) for s, p in zip(out["stage"], out["population"])]
    summary = (out[out["n_both_nonzero"] >= 3]
               .groupby("klass")
               .agg(n_clusters=("cluster", "size"),
                    n_L_direction=("bio_direction", lambda s: (s == "L").sum()),
                    n_R_direction=("bio_direction", lambda s: (s == "R").sum()),
                    mean_L_frac=("n_L_dom",
                        lambda s: np.mean(s / (s + out.loc[s.index, "n_R_dom"])
                                          .replace(0, np.nan)).round(3)),
                    median_log2=("median_log2_fold", "median")))
    print(summary.to_string())


if __name__ == "__main__":
    main()
