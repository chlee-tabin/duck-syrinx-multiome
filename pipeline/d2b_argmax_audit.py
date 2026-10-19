#!/usr/bin/env python3
"""AUDIT (one-off): does the DUCK d2b GREEDY assignment differ from per-cluster ARGMAX?

Prompted by the chick figS04 finding that d2b's greedy one-cluster-per-celltype (specific
types first) MIS-labels the chick — it assigned the main LPM blob to "Splanchnic" because the
chick doesn't separate LPM from Splanchnic, forcing "LPM" onto a leftover. The duck has
distinct LPM (GPC6/TSHZ2) and Splanchnic (NFIA/ADAMTS19) clusters, so it SHOULD be safe — but
the markers + greedy order are identical, so verify rather than assume.

Replicates d2b's EXACT scoring (frac-expressing per cluster x celltype on .raw) + the greedy
assignment, and also computes plain argmax, then prints both side by side and flags mismatches.

Reads: $DUCK_SCRATCH_ROOT/phase2/step04_preCB_clustered.h5ad (leiden_res0.2)
Prints: per-cluster argmax(score) vs greedy + MISMATCH flags (stdout only; no writes)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np, scanpy as sc

H5AD = f"{SCRATCH_ROOT}/phase2/step04_preCB_clustered.h5ad"
RES = "leiden_res0.2"
CELLTYPE_MARKERS = {
    "RBC": ["HBZ", "HBA1", "LOC101798111", "FTH1"],
    "Notochord": ["COL8A1", "SHH", "KIAA1217"],
    "Cardiomyocyte": ["RYR2", "RBM20", "GATA5"],
    "Otic placode": ["AGR3", "EYA1", "PAX2"],
    "Endothelial": ["FLI1", "EGFL7", "PTPRM", "PECAM1"],
    "Pharyngeal endoderm": ["ADAMTS17", "TFAP2E", "FOXI2"],
    "Ventral neural prog.": ["SPON1", "RFX4", "ARX"],
    "Dorsal neural prog.": ["LMX1B", "SLIT2", "DMD"],
    "Diff. neuron": ["BRSK2", "MYT1L", "KALRN", "ELAVL4"],
    "Dermomyotome": ["CDH4", "TMTC2", "PAX3", "SPOCK3"],
    "Epithelial": ["EPCAM", "KIAA1217", "PTPRF"],
    "Neural crest": ["COL25A1", "TFAP2B"],
    "Neural progenitor": ["PRTG", "NPAS3", "PTPRN2", "SOX2"],
    "Sclerotome": ["PAX1", "PAX9", "FOXC2", "MEOX1"],
    "Splanchnic mesenchyme": ["FOXF1", "GATA4", "ZFPM2", "ISL1"],
    "LPM mesenchyme": ["KIF26B", "CACNA1C", "TBX5", "TBX4", "NFIA"],
}


def log(m): print(m, flush=True)


adata = sc.read_h5ad(H5AD)
log(f"[load] {adata.n_obs} cells; raw={adata.raw is not None}")
raw = adata.raw if adata.raw is not None else adata
rv = raw.var_names
clusters = sorted(adata.obs[RES].unique(), key=lambda x: int(x))

# scores: mean fraction-expressing across a type's present markers (d2b exact)
import pandas as pd
score = pd.DataFrame(index=[str(c) for c in clusters], columns=list(CELLTYPE_MARKERS), dtype=float)
for ct, mk in CELLTYPE_MARKERS.items():
    pres = [m for m in mk if m in rv]
    if not pres:
        log(f"  WARNING no markers for {ct}"); continue
    mi = [np.where(rv == m)[0][0] for m in pres]
    for cl in clusters:
        idx = np.where((adata.obs[RES] == cl).values)[0]
        x = raw.X[idx][:, mi]
        if hasattr(x, "toarray"): x = x.toarray()
        score.loc[str(cl), ct] = float((x > 0).mean(axis=0).mean())

# greedy (d2b): specific types first, one cluster per celltype, <0.01 floor, leftover Unassigned
greedy, remaining = {}, set(str(c) for c in clusters)
for ct in CELLTYPE_MARKERS:
    cand = score.loc[list(remaining), ct].sort_values(ascending=False) if remaining else pd.Series(dtype=float)
    if len(cand) == 0 or cand.iloc[0] < 0.01: continue
    greedy[cand.index[0]] = ct; remaining.discard(cand.index[0])
for cl in remaining: greedy[cl] = f"Unassigned_{cl}"

# argmax (per-cluster best type) + compare
log(f"\n{'cl':>3} {'n':>6}  {'ARGMAX(score)':<28} {'GREEDY':<24} flag")
nmis = 0
for cl in clusters:
    c = str(cl); n = int((adata.obs[RES] == cl).sum())
    am = score.loc[c].astype(float).idxmax(); ams = score.loc[c, am]
    am2 = score.loc[c].astype(float).sort_values(ascending=False).index[1]
    gd = greedy[c]
    mism = "" if am == gd else "  <-- MISMATCH"
    if am != gd: nmis += 1
    log(f"{c:>3} {n:>6}  {am+'('+format(ams,'.2f')+')':<28} {gd:<24} {mism}  [argmax2={am2}({score.loc[c,am2]:.2f})]")
log(f"\n=== {nmis} / {len(clusters)} clusters where argmax != greedy ===")
log("If 0: duck greedy == argmax (safe). If >0: inspect those clusters (greedy artifact like chick cl0).")
