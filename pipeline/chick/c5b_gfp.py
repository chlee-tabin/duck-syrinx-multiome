#!/usr/bin/env python3
"""c5b — Group B GFP validation + per-donor GFP/WT call (orthogonal disambiguator).

Per-cell EGFP UMI (from the CellRanger matrix; EGFP is in the gg6a.EGFP reference) ->
per-donor GFP+ fraction. Every vireo donor must be uniformly GFP+ or WT (expect 6/6);
a donor with an intermediate GFP+ fraction is a GFP/WT merger. GFP-status then anchors
donor->embryo assignment (c6b), together with lane(=side) and stage.

Prereq: stage the GEX09/GEX10 filtered_feature_bc_matrix to WORK/matrix/<LIB>/
  (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz) — small (~MB).
  conda run -n vireo python3 pipeline/chick/c5b_gfp.py [WORK_DIR]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
LIBS = [("GEX09_GEX", "9"), ("GEX10_GEX", "10")]
OUT = f"{WORK}/vireo_seeded/seed_qc"; os.makedirs(OUT, exist_ok=True)
GFP_POS_UMI = 1   # a cell is EGFP+ if EGFP UMI >= this

# ── per-cell EGFP UMI from each library's CellRanger matrix ──────────────────
egfp = {}
for L, suf in LIBS:
    mdir = f"{WORK}/matrix/{L}"
    if not os.path.isdir(mdir):
        sys.exit(f"ERROR: {mdir} not found — stage filtered_feature_bc_matrix for {L} first")
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{mdir}/features.tsv.gz", "rt")]
    bcs = [l.strip() for l in gzip.open(f"{mdir}/barcodes.tsv.gz", "rt")]
    egi = [i for i, f in enumerate(feats)
           if (len(f) > 1 and f[1].strip().upper() == "EGFP") or f[0].strip().upper() == "EGFP"]
    M = sio.mmread(f"{mdir}/matrix.mtx.gz").tocsr()       # features x cells
    row = np.asarray(M[egi[0]].todense()).ravel() if egi else np.zeros(M.shape[1])
    if not egi:
        print(f"WARN: no EGFP feature in {L}")
    for bc, v in zip(bcs, row):
        egfp[f"{bc}-{suf}"] = float(v)

# ── best seed (expected 12) ──────────────────────────────────────────────────
rep = pd.read_csv(f"{OUT}/seed_qc_report.tsv", sep="\t")
pick = rep[(rep.N == 12) & (rep.n_donors == 12)]
if pick.empty: pick = rep[rep.n_donors == 12]
if pick.empty: pick = rep
bdir = pick.sort_values("logLik", ascending=False).iloc[0]["dir"]
did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t")
did["egfp"] = did["cell"].map(egfp).fillna(0.0)

# ── per-donor GFP call ───────────────────────────────────────────────────────
rows = []
for d in sorted(set(did.donor_id) - {"doublet", "unassigned"}):
    g = did.loc[did.donor_id == d, "egfp"].values
    rows.append({"donor": d, "n_cells": len(g), "mean_EGFP": round(float(g.mean()), 3),
                 "frac_EGFP_pos": round(float((g >= GFP_POS_UMI).mean()), 3)})
sdf = pd.DataFrame(rows).sort_values("frac_EGFP_pos").reset_index(drop=True)
# EGFP transgene is MOSAIC/lowly expressed -> GFP+ donors sit at frac~0.2-0.4 over a
# ~0.07 ambient WT baseline. Split at the largest natural gap among well-covered
# donors (>=100 cells), not a fixed cutoff.
real = sdf[sdf.n_cells >= 100].sort_values("frac_EGFP_pos")
if len(real) >= 3:
    fr = real["frac_EGFP_pos"].values
    gi = int(np.diff(fr).argmax())
    thresh = float((fr[gi] + fr[gi + 1]) / 2)
else:
    thresh = 0.15
sdf["call"] = [("GFP+" if f > thresh else "WT") + ("" if n >= 100 else " (low-cell)")
               for f, n in zip(sdf["frac_EGFP_pos"], sdf["n_cells"])]
sdf.to_csv(f"{OUT}/per_donor_gfp.tsv", sep="\t", index=False)
ngfp = sum(c == "GFP+" for c in sdf["call"]); nwt = sum(c == "WT" for c in sdf["call"])
print(f"[Group B GFP] data-driven split frac_EGFP_pos>{thresh:.3f} (EGFP mosaic/low-expressed); expect ~6/6 of 12:")
print(sdf.sort_values("frac_EGFP_pos", ascending=False).to_string(index=False))
print(f"\nGFP+={ngfp}  WT={nwt}  (well-covered; '(low-cell)' = under-captured/phantom)")
print(f"Wrote {OUT}/per_donor_gfp.tsv")
