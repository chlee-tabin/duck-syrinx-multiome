#!/usr/bin/env python3
"""Regenerate the Fig S8 (figures/fig01_supp) panel A/D library-QC tables for all 12
duck snRNA libraries.

Panel A reads   step01_library_summary.tsv      (uncorrected cellranger nuclei)
                step01_cb_library_summary.tsv   (CellBender nuclei)
Panel D reads   step02_mito_by_library.tsv.gz   (per-nucleus pct_mito, uncorrected)
                step02_cb_mito_by_library.tsv.gz (per-nucleus pct_mito, CellBender)

The committed tables (2026-03-22) were written by the archived scripts
FC_08656/archive_superseded/20260320_step01_load_libraries.py (+ 20260320_step02d_filter_merge.py; not distributed)
and 20260322_step01_cb_load.py (+ 20260322_step02d_cb_filter.py) for 10 libraries, before
the G3 pool (GEX04, GEX0S-duck) had CellBender. This script reproduces that recipe:

  * per-nucleus metrics: sc.pp.calculate_qc_metrics(qc_vars=[mito, ribo, rbc],
    percent_top=None, log1p=True), renamed n_genes / n_counts / pct_mito / ...
  * uncorrected: cellranger outs/filtered_feature_bc_matrix (read_10x_mtx, gene_symbols)
    + metrics_summary.csv columns merged into the library summary
  * CellBender: <LIB>_cellbender_filtered.h5 (read_10x_h5 + var_names_make_unique);
    nuclei with 0 detected genes after ambient removal are dropped (as in step01_cb)
  * library summary: groupby(library) n_cells / median_genes / median_counts /
    median_pct_mito / mean_pct_mito / median_pct_rbc, round(2)
  * mito table: library, pct_mito, n_genes, n_counts, filter_stage, index {cb}_{library}

Deliberate differences from the 2026-03 recipe (both documented in
figures/fig01_supp/PROVENANCE.md):
  1. Gene lists. When the committed tables were exported (2026-03-22) BOTH archived
     scripts used chicken (galGal6) names: mito "MT-ND2", "MT-CO1", "COII", ... (matches
     only 10 of the 13 duck ZJU1.0 mitochondrial genes) and RBC "HBZ", "HBBR", ... (matches
     none). The uncorrected script was fixed on 2026-03-27, after the export.
     Both arms here use the 13-gene mito list of pipeline/c1_qc_filter.py and the ZJU1.0
     hemoglobin LOC list. Both arms ALSO compute the legacy values (pct_mito_legacy10,
     pct_rbc_legacy) into validation sidecars, so the 10 committed libraries can be
     checked value-for-value.
  2. Only pre_filter rows are written to the mito tables. The 2026-03 tables also held
     post_filter rows (after scDblFinder + the 2026-03 genotype maps); those maps no longer
     exist, and panel D reads pre_filter rows only.

CellBender input source: the 10 original-library CellBender h5 files were purged with
scratch; the merged object they were loaded into, phase2/step01_cb_raw_merged.h5ad
(standby phase2_backup_20260512, 10 libraries, CellBender counts unmodified in X), is the
source for those 10. GEX04 / GEX0S-duck come from their CellBender h5 files
(phase2_backup_20260512/cellbender/). Do NOT use phase2/cellbender_as_mtx: that copy was
clipped and rounded to int32.

Usage (SLURM, see pipeline/s8ad_regen_libqc_tables.sh):
    python3 s8ad_regen_libqc_tables.py --part unc [--libs GEX0S-duck]
    python3 s8ad_regen_libqc_tables.py --part cb  [--libs GEX0S-duck]
"""
import argparse
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md

LIBRARIES = [
    "GEX03-1A", "GEX03-1B", "GEX03-2", "GEX04", "GEX0S-duck",
    "GEX05", "GEX06", "GEX07", "GEX08",
    "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B",
]
CB_FROM_H5 = {"GEX04", "GEX0S-duck"}  # the rest come from step01_cb_raw_merged.h5ad

# 13 duck (ZJU1.0) mitochondrial genes -- same as pipeline/c1_qc_filter.py MITO_GENES
MITO_GENES = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
              "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"]
# duck (ZJU1.0) hemoglobin genes, as in the archived uncorrected script after its 2026-03-27 fix
# (median_pct_rbc summary column only; not plotted)
RBC_GENES = ["LOC101798290", "LOC113839608", "LOC113844783", "LOC101798111",
             "LOC101792412", "LOC101800520", "LOC101800713"]
# the 2026-03-22 lists (galGal6 names); used ONLY for validation
MITO_GENES_LEGACY = ["ND1", "MT-ND2", "MT-CO1", "COII", "ATP8", "ATP6",
                     "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"]
RBC_GENES_LEGACY = ["HBZ", "HBBR", "HBA1", "HBE", "HBBA"]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def barcodes_to_index(adata, lib):
    cb = [bc.rsplit("-", 1)[0] for bc in adata.obs_names]
    if {len(c) for c in cb} != {16}:
        raise ValueError(f"{lib}: unexpected barcode lengths")
    if len(set(cb)) != len(cb):
        raise ValueError(f"{lib}: duplicate barcodes")
    adata.obs_names = [f"{c}_{lib}" for c in cb]
    adata.obs["library"] = lib
    adata.obs["cb"] = cb


def qc_obs(adata):
    """calculate_qc_metrics as in the archived step01 scripts; returns obs.

    pct_mito / pct_rbc use the duck lists; pct_mito_legacy10 / pct_rbc_legacy use the
    2026-03-22 lists (validation only)."""
    genes = set(adata.var_names)
    missing = [g for g in MITO_GENES if g not in genes]
    if missing:
        raise ValueError(f"mito genes missing from var_names: {missing}")
    rbc_missing = [g for g in RBC_GENES if g not in genes]
    if rbc_missing:
        raise ValueError(f"RBC genes missing from var_names: {rbc_missing}")
    legacy_m = [g for g in MITO_GENES_LEGACY if g in genes]
    legacy_r = [g for g in RBC_GENES_LEGACY if g in genes]
    ts(f"  legacy lists match mito {len(legacy_m)}/13 {legacy_m}, rbc {len(legacy_r)}/5")
    adata.var["mito"] = adata.var_names.isin(MITO_GENES)
    adata.var["ribo"] = adata.var_names.str.startswith("RPS") | adata.var_names.str.startswith("RPL")
    adata.var["rbc"] = adata.var_names.isin(RBC_GENES)
    adata.var["mito_legacy10"] = adata.var_names.isin(legacy_m)
    adata.var["rbc_legacy"] = adata.var_names.isin(legacy_r)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "ribo", "rbc", "mito_legacy10", "rbc_legacy"],
                               percent_top=None, log1p=True, inplace=True)
    adata.obs.rename(columns={
        "n_genes_by_counts": "n_genes", "total_counts": "n_counts",
        "pct_counts_mito": "pct_mito", "pct_counts_ribo": "pct_ribo",
        "pct_counts_rbc": "pct_rbc", "log1p_total_counts": "log1p_n_counts",
        "pct_counts_mito_legacy10": "pct_mito_legacy10",
        "pct_counts_rbc_legacy": "pct_rbc_legacy",
    }, inplace=True)
    return adata.obs.copy()


def per_lib_summary(obs, mito_col="pct_mito", rbc_col="pct_rbc"):
    return obs.groupby("library").agg(
        n_cells=("n_genes", "size"),
        median_genes=("n_genes", "median"),
        median_counts=("n_counts", "median"),
        median_pct_mito=(mito_col, "median"),
        mean_pct_mito=(mito_col, "mean"),
        median_pct_rbc=(rbc_col, "median"),
    ).round(2)


def parse_metrics_summary(csv_path):
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(csv_path)
    df = pd.read_csv(csv_path)
    out = {}
    for col in df.columns:
        val = df[col].iloc[0]
        if isinstance(val, str):
            try:
                val = float(val.replace(",", "").replace("%", ""))
            except ValueError:
                pass
        out[col] = val
    return out


def mito_table(obs):
    t = obs[["library", "pct_mito", "n_genes", "n_counts"]].copy()
    t["filter_stage"] = "pre_filter"
    return t


def write_validation(obs, outdir, arm, tag):
    """Sidecars (not committed to figures/data): summary + per-nucleus legacy values."""
    per_lib_summary(obs, "pct_mito_legacy10", "pct_rbc_legacy").reset_index().to_csv(
        f"{outdir}/validation_{arm}_library_summary_legacy{tag}.tsv", sep="\t", index=False)
    obs[["library", "pct_mito", "pct_mito_legacy10", "n_genes", "n_counts"]].to_csv(
        f"{outdir}/validation_{arm}_percell{tag}.tsv.gz", sep="\t", index=True, compression="gzip")


def run_unc(libs, indir, outdir, tag):
    obs_all, metrics = [], []
    ref_genes = None
    for lib in libs:
        mtx = f"{indir}/gex/{lib}_GEX/outs/filtered_feature_bc_matrix"
        ts(f"{lib}: {mtx}")
        a = sc.read_10x_mtx(mtx, var_names="gene_symbols", cache=False)
        barcodes_to_index(a, lib)
        if ref_genes is None:
            ref_genes = list(a.var_names)
        elif list(a.var_names) != ref_genes:
            raise ValueError(f"{lib}: gene list differs from {libs[0]}")
        if "feature_types" in a.var.columns and set(a.var["feature_types"]) != {"Gene Expression"}:
            raise ValueError(f"{lib}: non-GEX features")
        o = qc_obs(a)
        if not ((o["n_genes"] > 0).all() and (o["n_counts"] > 0).all()):
            raise ValueError(f"{lib}: zero-count barcodes in cellranger filtered matrix")
        ts(f"  {lib}: {len(o):,} nuclei, pct_mito dtype {o['pct_mito'].dtype}")
        obs_all.append(o)
        m = parse_metrics_summary(f"{indir}/gex/{lib}_GEX/outs/metrics_summary.csv")
        m["library"] = lib
        metrics.append(m)
        del a
    obs = pd.concat(obs_all)
    summary = per_lib_summary(obs).reset_index().merge(pd.DataFrame(metrics), on="library", how="left")
    summary.to_csv(f"{outdir}/step01_library_summary{tag}.tsv", sep="\t", index=False)
    mito_table(obs).to_csv(f"{outdir}/step02_mito_by_library{tag}.tsv.gz", sep="\t",
                           index=True, compression="gzip")
    write_validation(obs, outdir, "unc", tag)
    print(summary.iloc[:, :7].to_string(index=False))


def run_cb(libs, indir, outdir, tag):
    parts = []
    from_h5ad = [l for l in libs if l not in CB_FROM_H5]
    if from_h5ad:
        p = f"{indir}/step01_cb_raw_merged.h5ad"
        ts(f"loading {p}")
        a = sc.read_h5ad(p)
        ts(f"  {a.n_obs:,} x {a.n_vars:,}; X dtype {a.X.dtype}")
        print("  libraries in h5ad:", a.obs["library"].value_counts().sort_index().to_dict())
        stored = a.obs[["n_genes", "n_counts", "pct_mito"]].copy()
        if set(a.obs["library"]) != set(from_h5ad):
            a = a[a.obs["library"].isin(from_h5ad)].copy()
        for lib in from_h5ad:
            if (a.obs["library"] == lib).sum() == 0:
                raise ValueError(f"{lib}: absent from {p}")
        keep = ["library", "cb"]
        a.obs = a.obs[keep].copy()
        a.var = a.var[[]].copy()
        o = qc_obs(a)
        # the h5ad is already zero-gene filtered by step01_cb; recomputed metrics must
        # equal the ones step01_cb stored (legacy list)
        missing = o.index.difference(stored.index)
        if len(missing):
            raise ValueError(f"{len(missing):,} recomputed barcodes absent from stored step01_cb obs")
        st = stored.loc[o.index]
        d_g = int((st["n_genes"].values != o["n_genes"].values).sum())
        d_c = float(np.abs(st["n_counts"].values - o["n_counts"].values).max())
        d_m = float(np.abs(st["pct_mito"].values - o["pct_mito_legacy10"].values).max())
        ts(f"  vs stored step01_cb obs: n_genes mismatches {d_g}, max|dn_counts| {d_c:.3g}, "
           f"max|dpct_mito_legacy10| {d_m:.3g}")
        # fail loud: a wrong/modified h5ad must not silently write regenerated QC tables
        # (tolerances absorb float32 rounding; the committed run logged 0 mismatches)
        if d_g != 0 or d_c > 1e-3 or d_m > 1e-4:
            raise ValueError(
                f"recomputed QC != stored step01_cb obs (n_genes mismatches {d_g}, "
                f"max|dn_counts| {d_c:.3g}, max|dpct_mito| {d_m:.3g}) — wrong/modified h5ad?")
        if (o["n_genes"] == 0).any():
            raise ValueError("zero-gene nuclei in step01_cb_raw_merged.h5ad (expected pre-filtered)")
        parts.append(o)
        del a
    for lib in [l for l in libs if l in CB_FROM_H5]:
        p = f"{indir}/cellbender/{lib}_cellbender_filtered.h5"
        ts(f"{lib}: {p}")
        a = sc.read_10x_h5(p)
        a.var_names_make_unique()
        barcodes_to_index(a, lib)
        o = qc_obs(a)
        n0 = int((o["n_genes"] == 0).sum())
        ts(f"  {lib}: {len(o):,} barcodes, dropping {n0:,} with 0 genes")
        parts.append(o[o["n_genes"] > 0])
        del a
    obs = pd.concat(parts)
    obs = obs.loc[obs["library"].map({l: i for i, l in enumerate(LIBRARIES)}).sort_values(kind="stable").index]
    summary = per_lib_summary(obs).reset_index()
    summary["cellbender_h5"] = summary["library"] + "_cellbender_filtered.h5"
    summary.to_csv(f"{outdir}/step01_cb_library_summary{tag}.tsv", sep="\t", index=False)
    mito_table(obs).to_csv(f"{outdir}/step02_cb_mito_by_library{tag}.tsv.gz", sep="\t",
                           index=True, compression="gzip")
    write_validation(obs, outdir, "cb", tag)
    print(summary.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["unc", "cb"], required=True)
    ap.add_argument("--libs", default=",".join(LIBRARIES))
    ap.add_argument("--workdir", default=f"{SCRATCH_ROOT}/s8ad_regen_20261013")
    ap.add_argument("--tag", default="", help="suffix for output names (e.g. _pilot)")
    args = ap.parse_args()
    libs = args.libs.split(",")
    bad = [l for l in libs if l not in LIBRARIES]
    if bad:
        raise ValueError(f"unknown libraries {bad}")
    libs = [l for l in LIBRARIES if l in libs]
    indir, outdir = f"{args.workdir}/inputs", f"{args.workdir}/out"
    os.makedirs(outdir, exist_ok=True)
    ts(f"part={args.part} libs={libs} scanpy {sc.__version__} anndata {ad.__version__}")
    (run_unc if args.part == "unc" else run_cb)(libs, indir, outdir, args.tag)
    ts("done")


if __name__ == "__main__":
    sys.exit(main())
