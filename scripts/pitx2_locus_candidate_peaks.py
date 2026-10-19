#!/usr/bin/env python3
"""Part-C candidate-peak shortlist for the PITX2 Pitx2-TAD locus.

Ranks the 359 peaks by depth-robust biological interest from the g4f outputs:
  - persistence: # ventral stages with locus-FDR-sig L (padj_locus<0.1 & log2FC>0)
  - effect: mean ventral L/R log2FC
  - stage-dynamic: late/early log2FC (ventral HH17/16->HH21 opening) from g4f_stage_long
  - dorsal: HH21 log2FC (bilateral check)
Writes figures/data/pitx2_locus_candidate_peaks.tsv (sorted), classifies each peak,
and prints the shortlist.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, "figures/data")
PADJ = 0.1
VENT = ["HH17/16", "HH21", "HH25"]
TSS = {"PITX2a": 3_961_387, "PITX2b": 3_964_253, "PITX2c": 3_971_385}

pk = pd.read_csv(os.path.join(D, "pitx2_locus_tad_peaks.tsv"), sep="\t")
lr = pd.read_csv(os.path.join(D, "g4f_lr_long.tsv"), sep="\t")
st = pd.read_csv(os.path.join(D, "g4f_stage_long.tsv"), sep="\t")

# g4f_lr_long now carries a `design` column (unpaired/paired over EMBRYO-HALVES);
# the candidate shortlist uses the primary UNPAIRED laterality (all halves).
if "design" in lr.columns:
    lr = lr[lr["design"] == "unpaired"].copy()

vlr = lr[lr["compartment"] == "ventral"]
sigL = vlr[(vlr["log2FoldChange"] > 0) & (vlr["padj_locus"] < PADJ)]
nsig = sigL.groupby("peak_id")["stage"].nunique().rename("n_vent_stages_sigL")
meanL = vlr.groupby("peak_id")["log2FoldChange"].mean().rename("mean_vent_log2FC")
dor = (lr[(lr["compartment"] == "dorsal") & (lr["stage"] == "HH21")]
       .set_index("peak_id")["log2FoldChange"].rename("dorsal_HH21_log2FC"))
# ventral opening HH17/16->HH21
opn = (st[(st["compartment"] == "ventral") & (st["stage_pair"] == "HH17/16->HH21")]
       .set_index("peak_id")["log2FoldChange"].rename("vent_open_HH1716_HH21_log2FC"))

def nearest_iso(r):
    best, bd = None, 1e9
    for nm, t in TSS.items():
        d = (r["start"] + r["end"]) // 2 - t
        if abs(d) < bd:
            best, bd = nm, abs(d)
    return best

df = pk.copy()
df["nearest_isoform"] = df.apply(nearest_iso, axis=1)
df = (df.merge(nsig, on="peak_id", how="left").merge(meanL, on="peak_id", how="left")
        .merge(dor, on="peak_id", how="left").merge(opn, on="peak_id", how="left"))
df["n_vent_stages_sigL"] = df["n_vent_stages_sigL"].fillna(0).astype(int)

def classify(r):
    d = r["dist_PITX2c_TSS"]
    if r["n_vent_stages_sigL"] >= 2 and -600 <= d <= 13000:
        return "1_PERSISTENT_focal_L"          # the core laterality signal
    if r["n_vent_stages_sigL"] >= 1 and -600 <= d <= 13000:
        return "2_focal_L_partial"
    if r["region"] == "PLAYRR_LOC106016802":
        return "5_Playrr(no clean L/R)"
    if r["region"] == "FAM241A_C4orf32":
        return "6_FAM241A(RNA-R is artefact)"
    if r["dist_genebody"] <= -100000 and r["mean_vent_log2FC"] < -0.15:
        return "4_desert_R_drift"
    if abs(d) <= 50000:
        return "3_PITX2_prox"
    return "7_other"

df["class"] = df.apply(classify, axis=1)
cols = ["peak_id", "region", "nearest_isoform", "dist_PITX2c_TSS", "dist_genebody",
        "class", "n_vent_stages_sigL", "mean_vent_log2FC", "dorsal_HH21_log2FC",
        "vent_open_HH1716_HH21_log2FC"]
df = df[cols].sort_values(["class", "n_vent_stages_sigL", "mean_vent_log2FC"],
                          ascending=[True, False, False])
out = os.path.join(D, "pitx2_locus_candidate_peaks.tsv")
df.round(3).to_csv(out, sep="\t", index=False)
print(f"-> wrote {out}")
print("\n=== class tally ===")
print(df["class"].value_counts().sort_index().to_string())
print("\n=== TOP candidates (class 1-2, persistent/partial focal L) ===")
top = df[df["class"].str.startswith(("1_", "2_"))].copy()
print(top.round(2).to_string(index=False))
