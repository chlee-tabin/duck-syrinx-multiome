#!/usr/bin/env python3
"""E6->E9 responding population (label-free) + HH25 multiome back-check.

Question (2026-08-22): what distinguishes E6 (no proliferative
asymmetry, P-S1 primary) from E9 (left cycling excess)? Approach:
skip E06 L/R inference entirely (the phase-3 de-novo attempt found a latent-axis
signature but no confirmable identity; analysis not distributed); instead identify the responding population in
the old scRNA E06->E09 series WITHOUT side labels, then back-check in the HH25
multiome (where sides are known) for a left/right compositional bias of that
population's precursor state.

All inputs are committed tables; no cluster access, no new clustering.
  - figures/data/scrna_e6e9e11/cr9_er_survey_percell.tsv.gz
      (per-cell substage_population x stage x donor; singlets)
  - figures/data/scrna_e6e9e11/esr1_cellcycle_lr_tidy.tsv
      (figS07 panel-F export: per donor x side x population cycling fraction,
       REAL dissected sides at E09/E11)
  - figures/data/m5_cycling_lr_hh25_stats.tsv (the P-S1 record: HH25 back-check
       rows — Syringeal share L-R and within-Syringeal cycling)

Writes figures/data/m5_responding_population_e6e9.tsv.
Record + caveats: figures/figS16/STATS.md section 6c (not distributed).
"""
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "figures" / "data"
OLD = DATA / "scrna_e6e9e11"
MIN_CELLS_LANE = 10          # per donor x side x population, matches figS07 panel F
MIN_CELLS_DONOR = 5          # a donor "carries" a population

rows = []

# ---- 1. composition: which populations exist at each old-scRNA stage ----
s = pd.read_csv(OLD / "cr9_er_survey_percell.tsv.gz", sep="\t")
s = s[(s.scDblFinder_class == "singlet") & s.substage_population.notna()
      & (s.substage_population != "NA")]
comp = pd.crosstab(s.substage_population, s.stage)
donors = (s.groupby(["stage", "substage_population", "donor_id"]).size()
            .reset_index(name="n"))
donors = donors[donors.n >= MIN_CELLS_DONOR]
dcnt = donors.groupby(["stage", "substage_population"]).donor_id.nunique()
for (stage, pop), nd in dcnt.items():
    rows.append({"block": "composition", "stage": stage, "population": pop,
                 "estimate": int(comp.loc[pop, stage]), "unit": "cells",
                 "n": int(nd), "n_unit": f"donors with >={MIN_CELLS_DONOR} cells",
                 "p": None, "direction": None, "notes": ""})

for pop in ["Syringeal chondrocyte", "Syringeal fibroblast"]:
    if comp.loc[pop, "E06"] != 0:
        raise AssertionError(f"{pop} now annotated at E06 ({comp.loc[pop,'E06']} cells) "
                             "— the emergence claim in STATS.md §6c must be re-verified")
    if dcnt.get(("E09", pop), 0) != 8:
        raise AssertionError(f"{pop} not in 8/8 donors at E09 — re-verify §6c")
if comp.loc["Ventral spl.", "E06"] == 0:
    raise AssertionError("Ventral spl. absent at E06 — annotation vintage drift")

# ---- 2. cycling L-R per population at E09/E11 (real dissected sides) ----
t = pd.read_csv(OLD / "esr1_cellcycle_lr_tidy.tsv", sep="\t")
t = t[t.stage != "stage"]                      # embedded duplicate header row
t["cycling"] = t.cycling.astype(float)
t["n_cells"] = t.n_cells.astype(int)
for stage in ["E09", "E11"]:
    for lane in ["Syringeal chondrocyte", "Syringeal fibroblast", "Dorsal sub-epi."]:
        d = t[(t.stage == stage) & (t.lane == lane) & (t.n_cells >= MIN_CELLS_LANE)]
        w = d.pivot_table(index=["unit", "sex"], columns="side", values="cycling").dropna()
        if len(w) < 4:
            rows.append({"block": "cycling_LR", "stage": stage, "population": lane,
                         "estimate": None, "unit": "cycling fraction (L - R)",
                         "n": int(len(w)), "n_unit": "donor pairs", "p": None,
                         "direction": None, "notes": "fewer than 4 pairs — not tested"})
            continue
        diff = w["L"] - w["R"]
        p = wilcoxon(diff, alternative="two-sided", mode="exact").pvalue
        sexes = {sx: diff[[i for i in diff.index if i[1] == sx]] for sx in ("M", "F")}
        rows.append({"block": "cycling_LR", "stage": stage, "population": lane,
                     "estimate": round(float(diff.mean()), 4),
                     "unit": "cycling fraction (L - R, mean of donor diffs)",
                     "n": int(len(w)), "n_unit": "donor pairs (dissected sides)",
                     "p": round(float(p), 4),
                     "direction": f"{int((diff > 0).sum())}/{len(w)} L>R",
                     "notes": "; ".join(
                         f"{sx}: mean {d2.mean():+.3f}, {int((d2 > 0).sum())}/{len(d2)} L>R"
                         for sx, d2 in sexes.items() if len(d2))})

cyc = {(r["stage"], r["population"]): r for r in rows if r["block"] == "cycling_LR"}
e09c = cyc[("E09", "Syringeal chondrocyte")]
if not (e09c["p"] <= 0.05 and e09c["direction"].startswith(f"{e09c['n']}/")):
    raise AssertionError(f"E09 Chondrogenic cycling L-bias no longer unanimous+significant "
                         f"({e09c}) — re-verify §6c")
for lane in ["Syringeal fibroblast", "Dorsal sub-epi."]:
    if cyc[("E09", lane)]["p"] <= 0.05:
        raise AssertionError(f"E09 {lane} cycling now significant — §6c specificity "
                             "claim must be re-verified")

# ---- 3. HH25 multiome back-check (rows from the committed P-S1 record) ----
m = pd.read_csv(DATA / "m5_cycling_lr_hh25_stats.tsv", sep="\t")
for arm, gloss in [("falsification_check",
                    "HH25 chondro-precursor (Syringeal) SHARE of pool, L-R — no "
                    "compositional bias of the precursor state at E6"),
                   ("S2_Syringeal",
                    "HH25 cycling WITHIN the chondro-precursor (Syringeal), L-R")]:
    r = m[m.arm == arm]
    if len(r) != 1:
        raise KeyError(f"{arm} row missing/duplicated in m5_cycling_lr_hh25_stats.tsv")
    r = r.iloc[0]
    rows.append({"block": "hh25_backcheck", "stage": "HH25", "population": "Syringeal",
                 "estimate": round(float(r.estimate), 4), "unit": r.unit,
                 "n": int(r.n_embryos), "n_unit": "embryos",
                 "p": round(float(r.p), 4), "direction": r.direction_count,
                 "notes": gloss})

# ---- 4. direct HH25 x E09 alignment ----
# The subway map (manuscript S14C) draws only ADJACENT stage columns, so multiome
# HH25 only ever gets edges to E06 there; the committed all-vs-all MetaNeighbor
# AUROC matrix answers the direct HH25 -> E09 question.
aur = pd.read_csv(OLD / "xstage_auroc_matrix.tsv", sep="\t", index_col=0)
hh25 = [l for l in aur.index if l.startswith("HH25|")]
e09l = [l for l in aur.index if l.startswith("E09|")]
blk = aur.loc[hh25, e09l]
for r in blk.index:
    best = blk.loc[r].idxmax()
    recip = blk.loc[:, best].idxmax() == r
    rows.append({"block": "xstage_direct_HH25xE09", "stage": "HH25xE09",
                 "population": r.split("|", 1)[1],
                 "estimate": round(float(blk.loc[r, best]), 3),
                 "unit": f"MetaNeighbor AUROC vs best E09 match ({best.split('|',1)[1]})",
                 "n": None, "n_unit": None, "p": None,
                 "direction": "reciprocal-best" if recip else "NOT reciprocal",
                 "notes": "; ".join(f"{c.split('|',1)[1]} {blk.loc[r,c]:.3f}"
                                    for c in blk.columns)})
if not (blk.loc["HH25|Ventral spl."].idxmax() == "E09|Syringeal chondrocyte"
        and blk.loc[:, "E09|Syringeal chondrocyte"].idxmax() == "HH25|Ventral spl."):
    raise AssertionError("HH25 Ventral spl. <-> E09 Syringeal chondrocyte no longer "
                         "reciprocal-best — §6c antecedent claim must be re-verified")
syr_best = blk.loc["HH25|Syringeal"].idxmax()
if blk.loc[:, syr_best].idxmax() == "HH25|Syringeal":
    raise AssertionError("HH25 Syringeal now has a reciprocal-best E09 partner "
                         f"({syr_best}) — §6c orphan statement must be re-verified")

out = pd.DataFrame(rows)
dest = DATA / "m5_responding_population_e6e9.tsv"
out.to_csv(dest, sep="\t", index=False)
print(out.to_string(index=False))
print(f"\nwrote {dest}; all §6c assertions hold")
