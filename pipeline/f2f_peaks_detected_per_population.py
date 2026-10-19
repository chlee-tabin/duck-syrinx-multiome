#!/usr/bin/env python3
"""f2f: union peaks detected per CURRENT (stage × LPM population) — figS03 panel B1/B2.

The production ArchR union (283,818 peaks) was CALLED per April-2026 peak-calling
pool, so per-pool yields carry retired group labels and cannot be re-attributed
to today's clusters without re-calling peaks (which would change the union every
downstream result is defined on).  This step recounts a pool-free yield on the
SAME union: for every current d2g population × stage (the Figure 4B vocabulary),
the number of union peaks detected in that group's nuclei.  Every threshold is a
column, so the panel's choice never needs a rerun:

  n_peaks_any      >= 1 fragment in the group's pseudobulk
  n_peaks_frac01   >= 1 fragment in >= 1% of the group's nuclei (min 3)   <- B1/B2
  n_peaks_frac05   >= 1 fragment in >= 5% of the group's nuclei (min 3)
  n_peaks_cpm1     pseudobulk CPM >= 1

Population per nucleus is re-derived from the f1 sidecar's (d2g_stage, d2g_cluster)
via figures/_shared/subcluster_map.get_domain_map — the same recipe as g1 — never
from the name baked into the matrix.  Fails loud on an unmapped (stage, leiden)
pair, a matrix that is not the 283,818-peak union, or a join that matches nothing.

Input : phase2/atac/step_f2_peak_matrix_archr.h5ad      (F2_H5AD)
        phase2/atac/step_f1_coupled_atac_metadata.tsv   (F1_META_TSV)
Output: phase2/atac/step_f2f_peaks_detected_per_population.tsv
Env   : snapatac2 conda env (snapatac2, scipy, pandas); ~64 GB, one core.
"""
import math
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import snapatac2 as snap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import SCRATCH_ROOT, REPO_ROOT  # noqa: E402
sys.path.insert(0, str(REPO_ROOT / "figures"))
from _shared.subcluster_map import get_domain_map  # noqa: E402

ATAC_DIR = f"{SCRATCH_ROOT}/phase2/atac"
H5AD = os.environ.get("F2_H5AD", f"{ATAC_DIR}/step_f2_peak_matrix_archr.h5ad")
META = os.environ.get("F1_META_TSV", f"{ATAC_DIR}/step_f1_coupled_atac_metadata.tsv")
OUT = f"{ATAC_DIR}/step_f2f_peaks_detected_per_population.tsv"
N_UNION = 283_818


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def leiden_str(v):
    if pd.isna(v):
        return None
    try:
        return str(int(v))          # int / numpy-int / float-integer / "8"
    except (ValueError, TypeError):
        return str(v)


def main() -> None:
    for p in (H5AD, META):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing input {p}")
    ts(f"reading {H5AD}")
    # SnapATAC2 wrote the h5ad (its own HDF5 codec — plain anndata cannot decode X here);
    # read it back with snapatac2 exactly as g1 does.  obs is Polars-backed: pull columns
    # one by one.  adata.X[:] materialises the backed sparse matrix.
    adata = snap.read(H5AD)
    if adata.n_vars != N_UNION:
        sys.exit(f"ERROR: matrix has {adata.n_vars:,} peaks — not the {N_UNION:,}-peak production union")
    try:
        lib = [str(v) for v in adata.obs["arc_library"]]
    except Exception as e:  # noqa: BLE001
        sys.exit(f"ERROR: obs lacks 'arc_library' ({e}); cannot join the sidecar by (arc_library, cb16)")
    cb16 = [str(b).split("-")[0] for b in list(adata.obs_names)]
    X = adata.X[:]
    if not sp.isspmatrix_csr(X):
        X = sp.csr_matrix(X)
    n_obs, n_vars = X.shape
    if len(lib) != n_obs:
        sys.exit(f"ERROR: obs rows {len(lib):,} != matrix rows {n_obs:,}")
    ts(f"  {n_obs:,} nuclei × {n_vars:,} peaks; nnz {X.nnz:,}")

    ts(f"reading sidecar {META}")
    meta = pd.read_csv(META, sep="\t",
                       usecols=lambda c: c in {"arc_library", "cb16", "d2g_stage", "d2g_cluster"},
                       dtype={"arc_library": str, "cb16": str, "d2g_stage": str})
    need = {"arc_library", "cb16", "d2g_stage", "d2g_cluster"}
    if not need <= set(meta.columns):
        sys.exit(f"ERROR: sidecar lacks column(s) {sorted(need - set(meta.columns))}")
    meta = meta[meta["d2g_cluster"].notna() & meta["d2g_stage"].notna()].copy()
    if meta.duplicated(["arc_library", "cb16"]).any():
        sys.exit("ERROR: duplicate (arc_library, cb16) rows in the sidecar — join key not unique")
    parts = []
    for st, g in meta.groupby("d2g_stage"):
        dm = get_domain_map(st)                       # {str(cluster_id): population}
        parts.append(g["d2g_cluster"].map(leiden_str).map(dm))
    meta["population"] = pd.concat(parts).reindex(meta.index)
    un = meta.loc[meta["population"].isna(), ["d2g_stage", "d2g_cluster"]]
    if len(un):
        pairs = sorted({(s, leiden_str(c)) for s, c in zip(un["d2g_stage"], un["d2g_cluster"])})
        sys.exit(f"ERROR: sidecar (stage, leiden) pairs with no SSOT entry: {pairs}")
    lut = {(l, c): (s, p) for l, c, s, p in
           zip(meta["arc_library"], meta["cb16"], meta["d2g_stage"], meta["population"])}
    keys = [lut.get((l, c)) for l, c in zip(lib, cb16)]
    n_lpm = sum(k is not None for k in keys)
    ts(f"  {n_lpm:,} matrix nuclei carry a current LPM population; "
       f"{n_obs - n_lpm:,} lie outside the d2g LPM selection")
    if n_lpm == 0:
        sys.exit("ERROR: no matrix nucleus joins the sidecar — check obs_names / arc_library")

    key_arr = np.array([f"{k[0]}|{k[1]}" if k else "" for k in keys])
    rows = []
    for st, pop in sorted({k for k in keys if k is not None}):
        idx = np.flatnonzero(key_arr == f"{st}|{pop}")
        n = int(len(idx))
        Xg = X[idx]
        pb = np.asarray(Xg.sum(axis=0)).ravel()
        det = np.asarray((Xg > 0).sum(axis=0)).ravel()
        tot = float(pb.sum())
        cpm = pb / tot * 1e6 if tot > 0 else pb
        r = dict(stage=st, population=pop, n_nuclei=n, n_fragments_in_peaks=int(tot),
                 n_peaks_any=int((pb > 0).sum()),
                 n_peaks_frac01=int((det >= max(3, math.ceil(0.01 * n))).sum()),
                 n_peaks_frac05=int((det >= max(3, math.ceil(0.05 * n))).sum()),
                 n_peaks_cpm1=int((cpm >= 1).sum()),
                 n_union_peaks=N_UNION)
        rows.append(r)
        ts(f"  {st:8s} {pop:24s} n={n:6,d}  any={r['n_peaks_any']:7,d}  "
           f"1%={r['n_peaks_frac01']:7,d}  5%={r['n_peaks_frac05']:7,d}  cpm1={r['n_peaks_cpm1']:7,d}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    ts(f"wrote {OUT} ({len(out)} stage × population groups)")


if __name__ == "__main__":
    main()
