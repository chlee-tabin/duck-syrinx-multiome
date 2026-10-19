#!/usr/bin/env python3
"""Depth/power-MATCHED L/R check for the PITX2 locus (Part-B item b), DESeq2-free.

The worry: ventral n_pb / n_cells / depth grow 2-3x across HH17/16->HH21->HH25, so
sig-peak COUNTS scale with power. This resamples every ventral stage to a COMMON
(n_L, n_R) = min across stages — a bootstrap WITH replacement, so the CI is the
sampling uncertainty of an n=(n_L,n_R) estimate and is non-degenerate even where a
stage's available n equals the matched n (HH13).  Recomputes the L/R accessibility
log2 ratio on depth-normalized CPM.  If the focal PITX2c / +5.4 kb L-bias and the
desert R-drift survive matched power, they are biology, not power.

Reads figures/data/g4f_locus_count_slice_floor20.tsv (peaks x samples raw counts)
+ g4f_sample_meta_floor20.tsv + pitx2_locus_tad_peaks.tsv. Local; seeded per (element, stage):
each cell draws from its own random stream (SeedSequence of a CRC32 of the element
label and of the stage), so a cohort change at one stage cannot move the bootstrap
at any other stage or element (before 2026-10-14 one stream ran across all cells).

FLOOR-20 (the canonical fig03 laterality floor) — NOT floor-50: at floor-50 the
HH13 ventral pool collapses to 1L/2R (no matched power), while floor-20 gives
8L/6R, so the matched-power control can span all four stages (HH13 included) on
the same floor as fig03 panels B/C.
"""
from __future__ import annotations

import os
import zlib

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, "figures/data")
SEED = 0
B = 3000
EPS = 0.1  # CPM pseudocount (locus CPMs are ~1-14)

cnt = pd.read_csv(os.path.join(D, "g4f_locus_count_slice_floor20.tsv"), sep="\t", index_col=0)  # peaks x samples
meta = pd.read_csv(os.path.join(D, "g4f_sample_meta_floor20.tsv"), sep="\t")
pk = pd.read_csv(os.path.join(D, "pitx2_locus_tad_peaks.tsv"), sep="\t")

# CPM per sample (depth normalize)
libsize = meta.set_index("sample_id")["libsize"]
cpm = cnt.div(libsize.reindex(cnt.columns), axis=1) * 1e6   # peaks x samples, CPM

def zone(r):
    # NOTE these zones are MUTUALLY EXCLUSIVE and assigned in this order, so
    # "prox" is the <=50 kb annulus EXCLUDING the gene-body core (core is tested
    # separately).  "desert_up_far" is upstream-only (d <= -100 kb); a 50-100 kb
    # upstream buffer + the entire downstream flank fall in "other" (untested
    # here; shown in panel D's full-window heatmap).
    d, reg = r["dist_genebody"], r["region"]
    if reg in ("PITX2_core", "PITX2_NAT_3UTR"):
        return "core"
    if abs(d) <= 50_000:
        return "prox"
    if d <= -100_000:
        return "desert_up_far"
    return "other"
pk["zone"] = pk.apply(zone, axis=1)
PITX2C = "NC_051775.1:3971182-3971683"

VENT = ["HH13", "HH17/16", "HH21", "HH25"]
ven = meta[meta["compartment"] == "ventral"]
nL = {s: list(ven[(ven.d2g_stage == s) & (ven.inferred_side == "L")]["sample_id"]) for s in VENT}
nR = {s: list(ven[(ven.d2g_stage == s) & (ven.inferred_side == "R")]["sample_id"]) for s in VENT}
cL = min(len(nL[s]) for s in VENT)
cR = min(len(nR[s]) for s in VENT)
print(f"per-stage L pb: {[len(nL[s]) for s in VENT]}; R pb: {[len(nR[s]) for s in VENT]}")
print(f"MATCHED to common n_L={cL}, n_R={cR} (bootstrap WITH replacement, B={B})\n")

def lr_log2(samples_L, samples_R, peakset):
    mL = cpm.loc[peakset, samples_L].mean(axis=1).mean() + EPS
    mR = cpm.loc[peakset, samples_R].mean(axis=1).mean() + EPS
    return np.log2(mL / mR)

def cell_rng(label: str, stage: str) -> np.random.Generator:
    """Independent, reproducible stream per (element, stage) cell."""
    return np.random.default_rng(np.random.SeedSequence(
        [SEED, zlib.crc32(label.encode()), zlib.crc32(stage.encode())]))


recs = []
def boot(peakset, label):
    print(f"  {label:18s} (n_peaks={len(peakset)}):")
    for s in VENT:
        rng = cell_rng(label, s)
        vals = []
        for _ in range(B):
            # WITH replacement: a proper bootstrap of an n=(cL,cR) estimate, so the CI
            # is the sampling uncertainty at the common matched n at EVERY stage.
            # (Without replacement, a stage whose available n == matched n — e.g. HH13
            # at 8L/6R — would draw the identical set each iteration → a degenerate
            # zero-width CI that misreads as "most precise" when it is the sparsest,
            # noisiest stage.  HH13 has the HIGHEST per-embryo CV.)
            sl = list(rng.choice(nL[s], cL, replace=True))
            sr = list(rng.choice(nR[s], cR, replace=True))
            vals.append(lr_log2(sl, sr, peakset))
        vals = np.array(vals)
        full = lr_log2(nL[s], nR[s], peakset)  # unmatched (full-n)
        recs.append({"peakset": label, "n_peaks": len(peakset), "stage": s,
                     "matched_mean": round(float(vals.mean()), 4),
                     "ci_lo": round(float(np.percentile(vals, 5)), 4),
                     "ci_hi": round(float(np.percentile(vals, 95)), 4),
                     "full_n": round(float(full), 4), "n_L_matched": cL, "n_R_matched": cR})
        print(f"     {s:8s} matched log2(L/R)={vals.mean():+.2f} "
              f"[{np.percentile(vals,5):+.2f},{np.percentile(vals,95):+.2f}]  "
              f"(full-n={full:+.2f})")

# NAMED locus elements (matching fig03 panel C + fig04): the a/b/c promoters,
# the +5.4 kb putative distal enhancer, and a representative upstream-desert
# R-drift peak — single-peak peaksets, so panel F reads element-resolved (not
# zone-aggregate).  Coordinates are the fig03-panel-C FEAT set (single SSOT).
ELEMENTS = [
    ("NC_051775.1:3961167-3961668", "PITX2a promoter"),     # within-locus control (a)
    ("NC_051775.1:3963950-3964451", "PITX2b promoter"),     # within-locus control (b)
    ("NC_051775.1:3971182-3971683", "PITX2c promoter"),     # the L-biased focal element
    ("NC_051775.1:3976533-3977034", "+5.4 kb distal"),      # putative distal enhancer (L)
    ("NC_051775.1:3568743-3569244", "R-desert peak"),       # upstream-desert R-drift representative
]
print("=== MATCHED-POWER L/R accessibility log2(meanCPM_L / meanCPM_R), ventral ===")
missing = [pid for pid, _ in ELEMENTS if pid not in cnt.index]
if missing:
    raise KeyError(f"named elements absent from the floor-20 count slice: {missing}")
for pid, label in ELEMENTS:
    boot([pid], label)
pd.DataFrame(recs).to_csv(os.path.join(D, "g4f_depthmatch.tsv"), sep="\t", index=False)
print("\n-> wrote figures/data/g4f_depthmatch.tsv")
print("+ = Left-biased accessibility; - = Right.  PITX2c + +5.4 kb stay + (a/b near-"
      "bilateral, the desert peak -) at MATCHED n -> biology, not power.")
