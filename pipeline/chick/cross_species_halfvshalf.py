#!/usr/bin/env python3
"""(i) CROSS-SPECIES orientation-free half-vs-half at the PITX2c promoter — duck vs chick,
SAME metric, SAME null (is the chick "balanced" a real negative where the duck is
off-diagonal, or is the metric just underpowered for both?).

Both species' embryos are bisected into two physical-half libraries/lanes (duck: g1b
top/bot library = L/R; chick: ARC01/02, ARC09/10). The orientation-free half-vs-half asks,
per embryo, whether its two halves DIVERGE (off the y=x diagonal) or are BALANCED (on it) in
PITX2c-promoter accessibility — discarding the L/R label, so the divergence MAGNITUDE is
non-circular for BOTH species (the physical split is real; only the label is inferred).

Inputs (chick schema, same columns): cell(=lane#barcode), stage, population, side, group,
genotype, nFrags, <locus peak cols>.
  chick: figures/data/chick_locus/{locus_percell_counts.tsv.gz, locus_peak_coords.tsv}
  duck : figures/data/duck_locus/{...}        (pipeline/chick/duck_pitx2_locus_percell_extract.py)

Metric = depth-robust FRiP-at-promoter = sum(insertions in the POOLED in_promoter_c peaks) /
sum(nFrags). Pooling all in_promoter_c peaks = the (ii)-validated promoter-region readout
(chick: both shoulders; duck: the TSS-covering peak). Within-embryo cell-shuffle null
calibrates the divergence expected from sampling alone.

  conda run -n duckfig python3 pipeline/chick/cross_species_halfvshalf.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))            # pipeline/chick
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)                                     # _chick_subcluster_pools is a sibling
sys.path.insert(0, os.path.join(REPO, "figures"))
from _chick_subcluster_pools import CHICK_POOLS  # noqa: E402
from _groupB_split_excluded import drop_split_excluded, excluded_archr_cbs  # noqa: E402
from _shared.subcluster_map import POOLS as DUCK_POOLS  # noqa: E402

EPS = 1e-6
FACET_VENTRAL = "Ventral splanchnic pool"
FACET_DORSAL = "Dorsal sub-epi."
DATA = os.path.join(REPO, "figures", "data")
# diagnostics home: docs/chick_control/ (not distributed; no shipped figure reads these
# diagnostics). The directory is created on demand.
XS = os.path.join(REPO, "docs", "chick_control")
os.makedirs(XS, exist_ok=True)


def lane_of(cell):
    return str(cell).split("#", 1)[0]


def facet_of(stage, population, pools):
    """Fold (stage, population) → Ventral/Dorsal facet via a POOLS dict (duck or chick).

    inference → Ventral splanchnic pool (HH13 inference also folds here for the
    cross-species ventral comparison); bilateral → Dorsal sub-epi.; else dropped.
    Robust to `population` being either an annotated NAME or (fallback) anything not in POOLS
    → None (fail-closed). Asserts at call sites if nothing folds.
    """
    p = pools.get(stage)
    if p is None:
        return None
    if population in p.get("inference", []):
        return FACET_VENTRAL
    if population in p.get("bilateral", []):
        return FACET_DORSAL
    return None


def pooled_promoter_col(pc, coords):
    """Add a 'PITX2c_prom_pooled' column = sum of all in_promoter_c peak columns."""
    prom = coords.loc[coords.in_promoter_c, "peak_id"].tolist()
    prom = [p for p in prom if p in pc.columns]
    if not prom:
        raise SystemExit("no in_promoter_c peaks present in the per-cell matrix")
    pc = pc.copy()
    pc["PITX2c_prom_pooled"] = pc[prom].sum(axis=1).to_numpy(float)
    return pc, prom


def halfvshalf(pc, peak, facet, pools, min_cells=5, k_null=2000, seed=0, metric="ratio"):
    """Per-embryo half-vs-half FRiP-at-peak within one facet — orientation-free.

    Identical logic to the chick generate.py orientation_free_halfvshalf, parametrised by the
    species POOLS so the duck and chick run through the SAME code path.

    metric: "ratio" → |log2(FRiP_a/FRiP_b)| (scale-free, but EXPLODES when one half ≈ 0 →
    zero-inflation under small counts); "diff" → |FRiP_a − FRiP_b| (depth-normalised
    DIFFERENCE — bounded, self-gates small counts: few insertions → both FRiP tiny → ΔFRiP≈0).
    The `abs_log2_div`/`null_median_div` columns carry whichever metric was chosen (name kept
    for backward-compat with the plot code).
    """
    rng = np.random.default_rng(seed)
    keep = [facet_of(s, p, pools) == facet for s, p in zip(pc.stage, pc.population)]
    sub = pc[keep].copy()
    sub["lane"] = sub["cell"].map(lane_of)

    def _stat(a, b):
        return abs(a - b) if metric == "diff" else abs(np.log2((a + EPS) / (b + EPS)))

    rows = []
    for (grp, geno, stage), d in sub.groupby(["group", "genotype", "stage"], observed=True):
        lanes = sorted(d.lane.unique())
        if len(lanes) != 2:
            continue
        la, lb = lanes
        da, db = d[d.lane == la], d[d.lane == lb]
        if len(da) < min_cells or len(db) < min_cells:
            continue
        cnt = d[peak].to_numpy(float)
        nf = d["nFrags"].to_numpy(float)
        na = len(da)
        acc_a = da[peak].sum() / max(da["nFrags"].sum(), EPS)
        acc_b = db[peak].sum() / max(db["nFrags"].sum(), EPS)
        obs = _stat(acc_a, acc_b)
        idx = np.arange(len(d))
        nulls = np.empty(k_null)
        for j in range(k_null):
            pp = rng.permutation(idx)
            ia, ib = pp[:na], pp[na:]
            aa = cnt[ia].sum() / max(nf[ia].sum(), EPS)
            bb = cnt[ib].sum() / max(nf[ib].sum(), EPS)
            nulls[j] = _stat(aa, bb)
        emp_p = (np.sum(nulls >= obs) + 1) / (k_null + 1)
        rows.append({"facet": facet, "group": grp, "genotype": geno, "stage": stage,
                     "n_a": na, "n_b": len(db), "acc_a": acc_a, "acc_b": acc_b,
                     "abs_log2_div": obs, "null_median_div": float(np.median(nulls)),
                     "emp_p": emp_p})
    return pd.DataFrame(rows)


def load_species(name, pools):
    d = os.path.join(DATA, f"{name}_locus")
    coords = pd.read_csv(os.path.join(d, "locus_peak_coords.tsv"), sep="\t")
    pc = pd.read_csv(os.path.join(d, "locus_percell_counts.tsv.gz"), sep="\t")
    if name == "chick":
        # drop the Group B cS-split doublet/unassigned nuclei still in the committed pre-fix
        # locus table (ArchR ids ARC09#/ARC10#…; 0 hits once the ArchR chain is rerun)
        pc = drop_split_excluded(pc, "chick_locus_percell_counts", barcode_col="cell",
                                 group_col="group", ids=excluded_archr_cbs())
    pc, prom = pooled_promoter_col(pc, coords)
    # diagnostic: what populations are present, and do they fold?
    present = sorted(set(zip(pc.stage.astype(str), pc.population.astype(str))))
    folded = {(s, p): facet_of(s, p, pools) for s, p in present}
    n_fold = sum(v is not None for v in folded.values())
    print(f"[{name}] {len(pc):,} cells | in_promoter_c peaks pooled: {prom}")
    print(f"[{name}] (stage,population)→facet ({n_fold}/{len(present)} fold):")
    for (s, p), f in folded.items():
        print(f"    {s:10s} {p:28s} -> {f}")
    if n_fold == 0:
        raise SystemExit(f"[{name}] NOTHING folded — population values are likely cluster IDs, "
                         f"not names. Inspect and add an ID→name map.")
    return pc, coords


def main():
    smark = {"HH13": "o", "HH15-16": "s", "HH16-17": "s", "HH17/16": "s",
             "HH21": "D", "HH25": "^"}
    species = [("duck", DUCK_POOLS), ("chick", CHICK_POOLS)]
    facets = [FACET_VENTRAL, FACET_DORSAL]

    loaded = {}
    for name, pools in species:
        loaded[name] = (load_species(name, pools), pools)   # a missing extract fails the run

    all_rows = []
    fig, axes = plt.subplots(len(species), len(facets),
                             figsize=(5.0 * len(facets), 4.6 * len(species)), squeeze=False)
    for ri, (name, pools) in enumerate(species):
        (pc, _coords), _ = loaded[name]
        for ci, fac in enumerate(facets):
            ax = axes[ri][ci]
            s = halfvshalf(pc, "PITX2c_prom_pooled", fac, pools, min_cells=5)
            if len(s):
                s = s.assign(species=name)
                all_rows.append(s)
            if not len(s):
                ax.text(0.5, 0.5, f"{name}\n{fac}\nno embryos (≥5/half)", ha="center",
                        va="center", transform=ax.transAxes, fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
                continue
            mx = max(s[["acc_a", "acc_b"]].to_numpy().max(), EPS) * 1.15
            ax.plot([0, mx], [0, mx], color="0.6", lw=0.9, zorder=1)
            for _, r in s.iterrows():
                ax.scatter(r.acc_a, r.acc_b, marker=smark.get(r.stage, "o"), s=46,
                           facecolor="k", edgecolor="k", linewidth=0.8, zorder=3)
            ax.set_xlim(0, mx); ax.set_ylim(0, mx); ax.set_aspect("equal")
            obs, null = s.abs_log2_div.median(), s.null_median_div.median()
            ndiv = int((s.emp_p < 0.05).sum())
            ax.set_title(f"{name.upper()} | {fac}\nn={len(s)} | obs {obs:.2f} vs null {null:.2f} "
                         f"| div {ndiv}/{len(s)}", fontsize=7.5)
            if ci == 0:
                ax.set_ylabel(f"{name}\nhalf B FRiP")
            if ri == len(species) - 1:
                ax.set_xlabel("half A FRiP")
            ax.spines[["top", "right"]].set_visible(False)
    hh = [plt.Line2D([], [], marker=m, color="k", ls="", mfc="k", label=st)
          for st, m in [("HH13", "o"), ("HH15-16/17", "s"), ("HH21", "D"), ("HH25", "^")]]
    axes[0][-1].legend(handles=hh, fontsize=6, frameon=False, loc="lower right")
    fig.suptitle("Cross-species orientation-free half-vs-half at the PITX2c promoter "
                 "(y=x = balanced) — duck (top) diverges, chick (bottom) balanced?", fontsize=10)
    fig.tight_layout()
    out = os.path.join(XS, "xspecies_halfvshalf.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print(f"wrote {os.path.relpath(out, REPO)}")

    if all_rows:
        df = pd.concat(all_rows, ignore_index=True)
        tsv = os.path.join(XS, "cross_species_halfvshalf.tsv")
        df.to_csv(tsv, sep="\t", index=False)
        print(f"wrote {os.path.relpath(tsv, REPO)}")
        print("\n=== summary (median obs vs null | divergent) ===")
        g = (df.groupby(["species", "facet"])
             .apply(lambda d: pd.Series({
                 "n": len(d), "obs": d.abs_log2_div.median(),
                 "null": d.null_median_div.median(),
                 "div": int((d.emp_p < 0.05).sum())}))
             .reset_index())
        print(g.to_string(index=False))


if __name__ == "__main__":
    main()
