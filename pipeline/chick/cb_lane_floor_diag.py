#!/usr/bin/env python3
"""Diagnose the GEX09-vs-GEX10 PITX2 lane floor: depth/detection-sensitivity vs empty-droplet ambient.

Numbers only (no big outputs) -- read from the job log.
  PART 1: per-lane cell depth (median UMI, genes/cell) + PITX2 (+frac, cp10k) on the pre-CB filtered cells.
  PART 2: PITX2+ fraction by shared UMI-depth bin -- is the GEX10>GEX09 gap depth-explained?
  PART 3: ambient PITX2 in empty/background droplets (raw_feature_bc_matrix.h5) -- is it truly ambient?

Reads pre-CB filtered matrices (chick_demux_B/matrix/{GEX09_GEX,GEX10_GEX}) and the staged raw
all-droplet h5 (chick_cellbender/raw/{GEX09,GEX10}/raw_feature_bc_matrix.h5).
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys
import numpy as np
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chick_matrix_io import load_filtered as _load_filtered   # shared cellranger-matrix loader

WORKB = f"{SCRATCH_USER}/chick_demux_B"
RAWROOT = f"{SCRATCH_USER}/chick_cellbender/raw"
FILT = [("GEX09_GEX", "GEX09"), ("GEX10_GEX", "GEX10")]


def load_filtered(lib):
    """(M genes x cells, sym) from WORKB/matrix/{lib} — pinned pre-CB (matches cb_libqc_export)."""
    return _load_filtered(WORKB, lib, "matrix")


print("=== PART 1: pre-CB filtered cells -- depth + PITX2 by lane ===", flush=True)
store = {}
for lib, name in FILT:
    M, sym = load_filtered(lib)
    tot = np.asarray(M.sum(0)).ravel()
    ng = np.asarray((M > 0).sum(0)).ravel()
    pj = np.where(sym == "PITX2")[0]
    pit = np.asarray(M[pj[0]].todense()).ravel() if len(pj) else np.zeros(M.shape[1])
    cp = np.where(tot > 0, pit / tot * 1e4, 0.0)
    med_among = float(np.median(cp[pit > 0])) if (pit > 0).any() else 0.0
    print(f"[{name}] n={M.shape[1]}  medUMI={np.median(tot):.0f}  meanUMI={tot.mean():.0f}  "
          f"medGenes={np.median(ng):.0f}  PITX2+frac={(pit > 0).mean():.4f}  "
          f"PITX2_cp10k_mean={cp.mean():.3f}  PITX2_cp10k_median(among+)={med_among:.3f}", flush=True)
    store[name] = (tot, pit)

print("\n=== PART 2: PITX2+ fraction by shared UMI-depth bin (depth-controlled) ===", flush=True)
allt = np.concatenate([store[n][0] for _, n in FILT])
edges = np.quantile(allt, [0, .2, .4, .6, .8, 1.0])
print("shared UMI quintile edges:", [int(e) for e in edges], flush=True)
for lib, name in FILT:
    tot, pit = store[name]
    out = []
    for i in range(len(edges) - 1):
        hi = edges[i + 1]
        m = (tot >= edges[i]) & ((tot < hi) if i < len(edges) - 2 else (tot <= hi))
        out.append(f"{(pit[m] > 0).mean():.3f}(n{int(m.sum())})" if m.sum() else "na")
    print(f"[{name}] PITX2+ by UMI-bin: {out}", flush=True)

print("\n=== PART 3: empty/background-droplet ambient PITX2 (raw_feature_bc_matrix.h5) ===", flush=True)
for name in ["GEX09", "GEX10"]:
    p = f"{RAWROOT}/{name}/raw_feature_bc_matrix.h5"
    if not os.path.exists(p):
        print(f"[{name}] raw h5 MISSING: {p}", flush=True)
        continue
    a = sc.read_10x_h5(p)
    a.var_names_make_unique()
    X = a.X.tocsr()
    tot = np.asarray(X.sum(1)).ravel()
    has = "PITX2" in a.var_names
    pcol = np.asarray(X[:, list(a.var_names).index("PITX2")].todense()).ravel() if has else np.zeros(X.shape[0])
    for lo, hi in [(5, 100), (10, 50), (2, 20)]:
        em = (tot >= lo) & (tot <= hi)
        at = float(X[em].sum())
        pt = float(pcol[em].sum())
        ppm = pt / at * 1e6 if at else 0.0
        posf = float((pcol[em] > 0).mean()) if em.any() else 0.0
        print(f"[{name}] empty[{lo:>3}-{hi:<3}] n_drops={int(em.sum()):>7}  ambient_PITX2_ppm={ppm:7.2f}  "
              f"PITX2+drops={posf:.5f}", flush=True)

print("\n=== interpretation guide ===", flush=True)
print("DEPTH effect      if: GEX10 medUMI >> GEX09  AND  PART2 gap closes within matched bins.", flush=True)
print("TRUE AMBIENT      if: GEX10 ambient_PITX2_ppm >> GEX09 in PART3 (then CB underremoved it).", flush=True)
print("PER-CELL/CAPTURE  if: PART2 gap persists within bins AND PART3 ambient ~equal.", flush=True)
