#!/usr/bin/env python3
"""Local analysis of the g4f PITX2-locus dynamics outputs (depth-controlled).

Answers the central Part-B/C questions from the ferried g4f tables + the peak
annotation:
  (1) Is the HH17/16 ventral L-lean FOCAL to PITX2 or a locus-wide/regional offset?
      -> median log2FC + frac_L per (compartment, stage, region).
  (2) Per-peak L/R trajectory of the focal elements (PITX2 a/b/c promoters, gene
      body, Playrr, FAM241A-desert) across stages.
  (3) Stage CPM trajectory (depth-normalized) of those elements per side.

Reads figures/data/* + figures/data/pitx2_locus_tad_peaks.tsv. Local.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = os.path.join(REPO, "figures/data")
PEAKS = os.path.join(REPO, "figures/data/pitx2_locus_tad_peaks.tsv")
PADJ = 0.1

pk = pd.read_csv(PEAKS, sep="\t")
lr = pd.read_csv(os.path.join(F, "g4f_lr_long.tsv"), sep="\t")
# g4f_lr_long now carries a `design` column (unpaired/paired over EMBRYO-HALVES);
# this analysis uses the primary UNPAIRED laterality (all halves), else rows double.
if "design" in lr.columns:
    lr = lr[lr["design"] == "unpaired"].copy()
lr = lr.merge(pk[["peak_id", "region", "dist_genebody", "dist_PITX2c_TSS", "in_old_500kb"]],
              on="peak_id", how="left")

# Coarser zone for the focal-vs-regional decomposition
def zone(r):
    reg, d = r["region"], r["dist_genebody"]
    if reg in ("PITX2_core", "PITX2_NAT_3UTR"):
        return "1_PITX2_core"
    if abs(d) <= 50_000:
        return "2_PITX2_prox(<=50kb)"
    if reg == "PLAYRR_LOC106016802":
        return "3_PLAYRR"
    if reg in ("FAM241A_C4orf32",) or d <= -500_000:
        return "4_FAM241A/farUp"
    if reg == "PANCR_LOC106017192":
        return "5_PANCR_3prime"
    return "6_desert_up" if d < 0 else "7_desert_down"

lr["zone"] = lr.apply(zone, axis=1)

VENT = ["HH17/16", "HH21", "HH25"]
print("=" * 78)
print("(1) FOCAL-vs-REGIONAL: median L/R log2FC by ZONE × stage  (+ = Left)")
print("    VENTRAL compartment (the powered one across 3 stages)")
print("=" * 78)
v = lr[lr["compartment"] == "ventral"]
for metric, fmt in (("median_log2FC", "%+.2f"), ("frac_L", "%.2f"), ("n", "%d")):
    print(f"\n  -- {metric} --")
    rows = []
    for z, g in v.groupby("zone"):
        rec = {"zone": z}
        for st in VENT:
            gg = g[g["stage"] == st]
            if metric == "median_log2FC":
                rec[st] = (fmt % gg["log2FoldChange"].median()) if len(gg) else "."
            elif metric == "frac_L":
                rec[st] = (fmt % (gg["log2FoldChange"] > 0).mean()) if len(gg) else "."
            else:
                rec[st] = fmt % len(gg)
        rows.append(rec)
    print(pd.DataFrame(rows).to_string(index=False))

print("\n" + "=" * 78)
print("(2) FOCAL PEAK L/R log2FC trajectory (ventral), key annotated elements")
print("=" * 78)
TSS = {"PITX2a": 3_961_387, "PITX2b": 3_964_253, "PITX2c": 3_971_385}
# label peaks of interest
def label_peak(pid, start, end, region, dgb):
    for nm, t in TSS.items():
        if start <= t <= end:
            return f"{nm}_promoter"
    if region == "PITX2_core":
        return "PITX2_genebody"
    if region == "PLAYRR_LOC106016802":
        return "Playrr"
    if region == "FAM241A_C4orf32":
        return "FAM241A"
    if region == "PANCR_LOC106017192":
        return "PANCR(LOC106017192)"
    return None

pk2 = pk.copy()
pk2["label"] = pk2.apply(lambda r: label_peak(r["peak_id"], r["start"], r["end"],
                                              r["region"], r["dist_genebody"]), axis=1)
focal = pk2[pk2["label"].notna()][["peak_id", "label"]]
lrf = lr[lr["compartment"] == "ventral"].merge(focal, on="peak_id")  # lr already has dist_PITX2c_TSS
vt = lrf.pivot_table(index=["label", "dist_PITX2c_TSS", "peak_id"], columns="stage",
                     values="log2FoldChange", aggfunc="first")
sig = lrf.pivot_table(index=["label", "dist_PITX2c_TSS", "peak_id"], columns="stage",
                      values="padj_locus", aggfunc="first")
out = vt[VENT].copy()
for st in VENT:
    out[st] = out[st].map(lambda x: f"{x:+.2f}" if pd.notna(x) else ".")
    out[st] = out[st] + np.where(sig[st].reindex(out.index) < PADJ, "*", "")
print(out.sort_values("dist_PITX2c_TSS").to_string())

print("\n  (Playrr / FAM241A / PITX2c in DORSAL HH21 too:)")
dt = lr[lr["compartment"] == "dorsal"].merge(focal, on="peak_id")
if len(dt):
    dd = dt[dt["label"].isin(["Playrr", "FAM241A", "PITX2c_promoter", "PITX2a_promoter",
                              "PITX2b_promoter", "PITX2_genebody"])][
        ["label", "stage", "log2FoldChange", "padj_locus"]]
    print(dd.sort_values(["label", "stage"]).to_string(index=False))

# how many focal-core peaks stay L across all stages vs the desert
print("\n" + "=" * 78)
print("(3) PERSISTENCE: median |log2FC| & frac_L for PITX2-core vs far-desert, ventral")
print("=" * 78)
for z in ("1_PITX2_core", "2_PITX2_prox(<=50kb)", "4_FAM241A/farUp", "6_desert_up"):
    g = v[v["zone"] == z]
    line = [f"{z:22s}"]
    for st in VENT:
        gg = g[g["stage"] == st]
        line.append(f"{st}: med={gg['log2FoldChange'].median():+.2f} fracL={ (gg['log2FoldChange']>0).mean():.2f} n={len(gg)}")
    print("  " + " | ".join(line))
