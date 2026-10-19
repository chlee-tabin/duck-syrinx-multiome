#!/usr/bin/env python3
"""cb_restage — re-stage CellBender-corrected counts into matrix_cb/ dirs the atlas (cAE) reads.

Mirrors the duck `c1_qc_filter.py` pre-CB pattern: keep the SAME cell set the atlas used (the
cellranger filtered barcodes), swap in CellBender-corrected counts where CellBender kept the cell,
and RESCUE (pre-CB counts) any atlas cell CellBender dropped — so cell-calling is unchanged and the
barcode-keyed demux/stage/sex/population labels stay aligned. CellBender output is non-integer →
clip(>=0) + round to int32 before writing (downstream tools expect ints).

Writes `{work}/matrix_cb/{library}/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}` (GEX-only;
ARC Peaks were excluded from CellBender's model but are passed through in its output h5, so
load_cb_h5 reads only the Gene Expression features). Originals in `{work}/matrix/{library}/` are NOT
touched. Then run cAE with CHICK_MATRIX_SUBDIR=matrix_cb to re-cluster on corrected counts.

Hard self-checks (this is the silent-failure-risk step — fail loud, never silently drop cells):
  - output barcode set == the atlas (pre-CB) barcode set, exactly
  - every output count is a non-negative integer
  - CellBender-cell overlap with the atlas set is reported; rescued (CB-dropped) cells counted
  - per-library global removal fraction reported (sanity; >>~10% or <0 = inspect)

Run (needs h5py + scanpy, e.g. the vireo env):
  conda run -n duckfig python3 pipeline/chick/cb_restage.py     # local, if h5s pulled
  # or on cluster: CHICK_CB_OUT=<scratch>/chick_cellbender/out python3 ...
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os
import gzip
import shutil
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

CB_OUT = os.environ.get("CHICK_CB_OUT", f"{SCRATCH_USER}/chick_cellbender/out")
WORK_A = os.environ.get("CHICK_WORK_A", f"{SCRATCH_USER}/chick_demux")
WORK_B = os.environ.get("CHICK_WORK_B", f"{SCRATCH_USER}/chick_demux_B")
# Hard floor on atlas-cell overlap with the CellBender output (measured 2026-10-13: GEX09 99.3%,
# GEX10 99.7%); below it the pre-CB "rescue" would silently replace corrected counts wholesale.
MIN_CB_OVERLAP = 0.95
SUBDIR_SRC = "matrix"
SUBDIR_DST = "matrix_cb"

# (work, atlas-library-name, cellbender-output-prefix)
LIBS = [
    (WORK_A, "ARC01", "ARC01"),
    (WORK_A, "ARC02", "ARC02"),
    (WORK_B, "GEX09_GEX", "GEX09"),
    (WORK_B, "GEX10_GEX", "GEX10"),
]


def log(m):
    print(m, flush=True)


def load_cb_h5(path):
    """Load a CellBender _filtered.h5 -> (cells x genes CSR, gene_ids array, barcodes array).

    Gene Expression features only. Reads the 10x-format h5 with h5py, cross-checked against
    scanpy. Returns RAW (non-integer) counts.
    """
    # CellBender's *_filtered.h5 is a 10x-format count matrix (group "matrix": data / indices /
    # indptr / shape / barcodes / features). For the ARC (multiome) libraries CellBender RETAINS
    # the excluded Peaks features in its output, while scanpy's read_10x_h5 (gex_only=True, the
    # default) returns Gene Expression only — so the direct h5py read is restricted to the
    # "Gene Expression" feature_type before the cross-check (like with like); any remaining
    # disagreement raises (no fallback). A file without feature_type, or with no GEX features,
    # raises too.
    import h5py
    import scanpy as sc
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        if "features" not in g or "feature_type" not in g["features"]:
            raise KeyError(f"{path}: matrix/features/feature_type missing — cannot select the "
                           f"Gene Expression features of the CellBender matrix")
        n_feat, n_cells = (int(x) for x in g["shape"][()])
        X_all = sp.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]),
                              shape=(n_feat, n_cells)).T.tocsr()
        bcs_d = [b.decode() for b in g["barcodes"][()]]
        fid = np.array([b.decode() for b in g["features"]["id"][()]])
        ftype = np.array([b.decode() for b in g["features"]["feature_type"][()]])
    is_gex = ftype == "Gene Expression"
    if not is_gex.any():
        raise ValueError(f"{path}: no 'Gene Expression' features (feature types {sorted(set(ftype))})")
    Xd = X_all[:, np.where(is_gex)[0]]
    gid_d = list(fid[is_gex])
    other = pd.Series(ftype[~is_gex]).value_counts().to_dict()
    a = sc.read_10x_h5(path)                                 # gex_only=True: GEX features only
    Xs = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    if (Xs.shape != Xd.shape or (Xs != Xd).nnz != 0
            or list(a.obs_names) != bcs_d or list(a.var["gene_ids"]) != gid_d):
        raise RuntimeError(f"{path}: h5py and scanpy reads of the CellBender GEX matrix disagree")
    log(f"    read {path}: {n_cells} cells x {int(is_gex.sum())} GEX features (h5py == scanpy); "
        f"non-GEX features ignored: {other or 'none'}")
    a.var_names_make_unique()
    gid = (a.var["gene_ids"].to_numpy() if "gene_ids" in a.var.columns
           else a.var["gene_id"].to_numpy() if "gene_id" in a.var.columns
           else a.var_names.to_numpy())
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    return X, np.asarray(gid).astype(str), np.asarray(a.obs_names).astype(str)


def read_features_barcodes(srcdir):
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{srcdir}/features.tsv.gz", "rt") if l.strip()]
    bcs = [l.strip() for l in gzip.open(f"{srcdir}/barcodes.tsv.gz", "rt") if l.strip()]
    feats = np.array(feats, dtype=object)
    gid = np.array([f[0] for f in feats]).astype(str)
    gsym = np.array([f[1] if len(f) > 1 else f[0] for f in feats]).astype(str)
    gtype = np.array([f[2] if len(f) > 2 else "Gene Expression" for f in feats]).astype(str)
    return gid, gsym, gtype, np.array(bcs).astype(str)


def write_mtx_dir(dstdir, mat_genes_x_cells, gid, gsym, barcodes):
    os.makedirs(dstdir, exist_ok=True)
    # features.tsv.gz: id<TAB>symbol<TAB>Gene Expression
    with gzip.open(f"{dstdir}/features.tsv.gz", "wt") as fh:
        for i, s in zip(gid, gsym):
            fh.write(f"{i}\t{s}\tGene Expression\n")
    with gzip.open(f"{dstdir}/barcodes.tsv.gz", "wt") as fh:
        for b in barcodes:
            fh.write(f"{b}\n")
    tmp = f"{dstdir}/matrix.mtx"
    sio.mmwrite(tmp, sp.csr_matrix(mat_genes_x_cells).astype(np.int32), field="integer")
    with open(tmp, "rb") as fi, gzip.open(f"{dstdir}/matrix.mtx.gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(tmp)


def restage_one(work, lib, cb):
    src = f"{work}/{SUBDIR_SRC}/{lib}"
    dst = f"{work}/{SUBDIR_DST}/{lib}"
    log(f"\n=== {lib} (CB={cb}) ===")
    gid_o, gsym_o, gtype_o, bcs_o = read_features_barcodes(src)
    gex = gtype_o == "Gene Expression"
    gid_g, gsym_g = gid_o[gex], gsym_o[gex]                  # GEX gene order the atlas uses
    Morig = sio.mmread(f"{src}/matrix.mtx.gz").tocsr()[gex]   # GEX genes x cells (pre-CB)
    log(f"  atlas (pre-CB): {len(bcs_o)} cells, {gex.sum()} GEX genes")

    Xcb, gid_cb, bcs_cb = load_cb_h5(f"{CB_OUT}/{cb}_cellbender_filtered.h5")  # cells x genes
    log(f"  cellbender:     {Xcb.shape[0]} cells, {Xcb.shape[1]} genes")

    # gene alignment by gene_id (robust to symbol-dup renaming); every atlas GEX gene must be in CB
    cb_gpos = {g: j for j, g in enumerate(gid_cb)}
    miss = [g for g in gid_g if g not in cb_gpos]
    assert not miss, f"{lib}: {len(miss)} atlas GEX genes absent from CB output (e.g. {miss[:3]})"
    gcol = np.array([cb_gpos[g] for g in gid_g])             # CB col index per atlas GEX gene
    Xcb_g = Xcb[:, gcol]                                     # cells x atlas-GEX-genes, CB order

    cb_bpos = {b: i for i, b in enumerate(bcs_cb)}
    in_cb = np.array([b in cb_bpos for b in bcs_o])
    if in_cb.mean() < MIN_CB_OVERLAP:
        raise ValueError(
            f"{lib}: only {int(in_cb.sum())}/{len(bcs_o)} ({100*in_cb.mean():.1f}%) atlas cells found in "
            f"{CB_OUT}/{cb}_cellbender_filtered.h5 (< {MIN_CB_OVERLAP:.0%}; barcode-format mismatch?); "
            f"e.g. atlas {list(bcs_o[:3])} vs CB {list(bcs_cb[:3])}")
    log(f"  atlas cells in CB output: {int(in_cb.sum())}/{len(bcs_o)} "
        f"({100*in_cb.mean():.1f}%); rescuing {int((~in_cb).sum())} CB-dropped with pre-CB counts")

    # build atlas-cells x GEX-genes: CB-corrected where present, else pre-CB (rescue)
    n_cells, n_genes = len(bcs_o), len(gid_g)
    out = sp.lil_matrix((n_cells, n_genes), dtype=np.float64)
    cb_rows = np.array([cb_bpos[b] for b in bcs_o[in_cb]])
    out[np.where(in_cb)[0]] = Xcb_g[cb_rows]
    if (~in_cb).any():                                       # rescue: pre-CB counts (genes x cells -> cells x genes)
        out[np.where(~in_cb)[0]] = Morig[:, np.where(~in_cb)[0]].T
    out = out.tocsr()
    out.data = np.rint(np.clip(out.data, 0, None))           # CB non-integer -> clip>=0 + round

    # removal fraction (CB cells only): how much total count CellBender removed vs pre-CB
    tot_pre = float(Morig[:, np.where(in_cb)[0]].sum())
    tot_cb = float(out[np.where(in_cb)[0]].sum())
    frac = 1 - tot_cb / tot_pre if tot_pre else float("nan")
    log(f"  global removal fraction (CB cells): {frac:.3%}")

    # hard checks
    assert out.shape == (n_cells, n_genes), f"{lib}: shape {out.shape} != ({n_cells},{n_genes})"
    assert np.all(out.data >= 0) and np.allclose(out.data, np.rint(out.data)), f"{lib}: non-int/neg counts"
    write_mtx_dir(dst, out.T.tocsr(), gid_g, gsym_g, bcs_o)   # write genes x cells
    # verify round-trip: re-read dims + barcodes match exactly
    _, _, _, bcs_chk = read_features_barcodes(dst)
    assert list(bcs_chk) == list(bcs_o), f"{lib}: written barcodes != atlas barcodes"
    log(f"  wrote {dst}  ({n_genes} genes x {n_cells} cells, int)")
    return dict(lib=lib, cells=n_cells, genes=n_genes, in_cb=int(in_cb.sum()),
                rescued=int((~in_cb).sum()), removal_frac=round(frac, 4))


def main():
    # Optional subset, e.g. CHICK_RESTAGE_LIBS=GEX09_GEX,GEX10_GEX (unknown names raise).
    want = [x for x in (os.environ.get("CHICK_RESTAGE_LIBS") or "").split(",") if x]
    known = {l for _, l, _ in LIBS}
    bad = sorted(set(want) - known)
    if bad:
        raise ValueError(f"CHICK_RESTAGE_LIBS has unknown libraries {bad}; known {sorted(known)}")
    libs = [t for t in LIBS if not want or t[1] in want]
    rows = [restage_one(w, l, c) for (w, l, c) in libs]
    log("\n=== re-stage summary ===")
    log(pd.DataFrame(rows).to_string(index=False))
    log("\nNext: re-run cAE with CHICK_MATRIX_SUBDIR=matrix_cb, then cAL -> cAM -> cAN.")


if __name__ == "__main__":
    main()
