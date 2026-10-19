"""Shared cellranger-matrix loader for the chick pipeline (ARC v2.1.0 blank-line tolerant).

The same filtered-matrix read was copy-pasted across cAE/cAG/cAF (atlas) and cb_libqc_export /
cb_lane_floor_diag (diagnostics). This is the single implementation. The ARC v2.1.0
`features.tsv.gz` can carry a trailing blank line that `sc.read_10x_mtx` chokes on, so the reader
filters blank lines explicitly.

One low-level reader, two thin shapes on top:
    read_matrix()   -> (M genes x cells CSR, sym, typ, barcodes)      [raw components]
    load_anndata()  -> GEX-only AnnData (cells x genes), barcodes suffixed   [atlas scripts]
    load_filtered() -> (M_gex genes x cells, sym_gex)                 [cb_* diagnostics]

`subdir` is an EXPLICIT argument (not a buried env read) so the matrix/ vs matrix_cb choice is
visible at each call site: 'matrix' = raw cellranger filtered; 'matrix_cb' = CellBender-corrected
(re-staged by cb_restage.py). The cb_* diagnostics deliberately pin 'matrix' so their pre-CB
numbers match each other.
"""
from __future__ import annotations

import gzip

import numpy as np
import scipy.io as sio
import scipy.sparse as sp


def read_matrix(work: str, library: str, subdir: str = "matrix"):
    """Load `{work}/{subdir}/{library}` -> (M, sym, typ, barcodes).

    M = genes x cells CSR (all feature types); sym/typ = per-feature symbol/type arrays;
    barcodes = list. Blank-line tolerant (ARC v2.1.0 quirk).
    """
    md = f"{work}/{subdir}/{library}"
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
    sym = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    typ = np.array([f[2] if len(f) > 2 else "Gene Expression" for f in feats])
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    return M, sym, typ, bcs


def load_anndata(work: str, library: str, suffix: str, subdir: str = "matrix"):
    """GEX-only AnnData (cells x genes) with barcodes `f'{bc}-{suffix}'`; var_names made unique.

    Caller adds any obs columns (e.g. `library`) and logging. Matches the old cAE/cAG/cAF
    `load_library` body exactly (those set `obs['library']` afterward; cAF does not).
    """
    # lazy on purpose: keeps read_matrix()/load_filtered() importable in an env without anndata
    # (only this AnnData projection needs it); the cluster atlas env always has anndata.
    import anndata as ad
    import pandas as pd
    M, sym, typ, bcs = read_matrix(work, library, subdir)
    gex = typ == "Gene Expression"
    a = ad.AnnData(sp.csr_matrix(M[gex].T),
                   obs=pd.DataFrame(index=[f"{bc}-{suffix}" for bc in bcs]),
                   var=pd.DataFrame(index=sym[gex]))
    a.var_names_make_unique()
    return a


def load_filtered(work: str, library: str, subdir: str = "matrix"):
    """(M_gex genes x cells CSR, sym_gex) — the cb_* diagnostic shape."""
    M, sym, typ, _ = read_matrix(work, library, subdir)
    gex = typ == "Gene Expression"
    return M[gex], sym[gex]
