#!/usr/bin/env python3
"""Stage-restricted mesenchyme subclustering.

Hypothesis: multi-stage data has stage as the dominant PC, obscuring
spatial axes (D/V, radial) that Yan et al. (2025) resolved in their
single-timepoint chick foregut data. By restricting to one stage and
re-computing PCA/clustering, we may recover D/V-aligned subclusters.

Approach:
  1. Load h5ad, filter to LPM mesenchyme (res=0.2 clusters 0,3) at target stage
  2. Re-compute HVGs on subset
  3. PCA → Harmony (batch = library) → UMAP → Leiden at multiple resolutions
  4. Export: UMAP coords, cluster assignments, key marker expression

Usage:
  python d2g_stage_restricted_subclustering.py [--stage HH21] [--stage HH25]
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os

# ── Pin BLAS/OpenMP threading BEFORE the numpy-pulling imports (matches d1) ──
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import scanpy as sc
import pandas as pd
import numpy as np
import random
import sys, argparse

# ── Global random seeds + single-threaded scanpy (matches d1) ──
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
sc.settings.n_jobs = 1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _d2g_lpm_selection import (  # noqa: E402
    LPM_ANNOTATIONS, load_cluster_annotations, fold_hh17_hh16,
)

H5AD = f"{SCRATCH_ROOT}/phase2/step04_preCB_clustered.h5ad"
STAGES_FILE = f"{SCRATCH_ROOT}/phase2/figures/e1b_stage_assignments.tsv"
OUTDIR = f"{SCRATCH_ROOT}/phase2/figures"

# Key markers to export
MARKERS = [
    # D/V axis (Yan et al.)
    "FOXF1", "PITX2", "FOXF2", "LOC113843492",  # NKX6-1
    "ISL1", "TBX5", "TBX4", "ZFPM2", "BMP4",
    # Somatic/splanchnic
    "PRRX1", "TWIST1", "NR2F2", "HAND1", "HAND2", "IRX3",
    # Radial
    "TBX1", "FOXC2", "PTCH1",
    # Other
    "BARX1", "GATA4", "COL1A1", "ACTA2",
    "WT1", "TCF21", "ALDH1A2",
    # Paraxial controls
    "PAX1", "PAX9", "PAX3", "PAX7",
]

# LPM mesenchyme selection lives in _d2g_lpm_selection.LPM_ANNOTATIONS
# (imported above).  Names come from d2b auto-annotator — the single
# source of truth.  Select by annotation name rather than cluster ID so
# downstream code is robust to seed-driven cluster renumbering.

# Column names in step04_preCB_clustered.h5ad
GENO_COL = "unified_genotype_id"
RES02 = "leiden_res0.2"
RES05 = "leiden_res0.5"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", nargs="+", default=["HH21"],
                        help="Stage(s) to subset. Default: HH21")
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--resolutions", nargs="+", type=float,
                        default=[0.1, 0.2, 0.3, 0.5, 0.8, 1.0])
    parser.add_argument("--regress-cc", action="store_true",
                        help="Experimental: also regress out S_score + G2M_score "
                             "(d1 does this globally; d2g default does NOT). "
                             "Output filename suffixed _ccreg so default + ccreg "
                             "variants can coexist side-by-side.")
    parser.add_argument("--outdir", default=None,
                        help="Override output directory (default: production "
                             f"{OUTDIR}). Use for experimental reruns so "
                             "production sidecars are never overwritten.")
    args = parser.parse_args()
    outdir = args.outdir or OUTDIR
    os.makedirs(outdir, exist_ok=True)

    stage_label = "_".join(args.stage).replace("/", "")
    if args.regress_cc:
        stage_label += "_ccreg"
    print(f"=== Stage-restricted subclustering: {args.stage} (regress_cc={args.regress_cc}) ===")

    # Load
    print(f"Loading: {H5AD}")
    adata = sc.read_h5ad(H5AD)
    print(f"  Total: {adata.n_obs:,} cells, {adata.n_vars:,} genes (.X)")
    if adata.raw is not None:
        print(f"  raw: {adata.raw.n_vars:,} genes")

    # Stage assignments (C1 input assertion: schema + row count)
    stages_df = pd.read_csv(STAGES_FILE, sep="\t")
    _E1B_COLS = ["genotype_id", "design_stage", "sex",
                 "inferred_stage", "PC1", "PC1_adj"]
    _missing = set(_E1B_COLS) - set(stages_df.columns)
    if _missing:
        raise RuntimeError(
            f"{STAGES_FILE}: missing columns {sorted(_missing)}.  "
            f"Re-run pipeline/e1b_stage_inference.py — the scratch copy "
            f"is likely stale (predates the sex-adjusted staging fix)."
        )
    if len(stages_df) != 52:
        raise RuntimeError(
            f"{STAGES_FILE}: expected 52 rows, got {len(stages_df)}.  "
            f"Re-run pipeline/e1b_stage_inference.py."
        )
    stage_map = dict(zip(stages_df["genotype_id"], stages_df["inferred_stage"]))
    adata.obs["stage"] = adata.obs[GENO_COL].map(stage_map)
    # HH17 genotypes (17-Md-13, 17-KC-1, 17-A..D) are from the same
    # developmental window as HH17/16 and should be subclustered together.
    fold_hh17_hh16(adata)

    # Filter to LPM mesenchyme at target stage(s) by annotation name.
    cl_to_name = load_cluster_annotations()
    celltype = adata.obs[RES02].astype(str).map(cl_to_name)
    lpm_mask = celltype.isin(LPM_ANNOTATIONS)
    stage_mask = adata.obs["stage"].isin(args.stage)
    mask = lpm_mask & stage_mask

    # Use .raw (all genes, CP10K log-normalized) for subclustering
    # step04's .X has only ~3K HVGs (scaled) — not suitable for de novo HVG
    if adata.raw is not None:
        print("  Using adata.raw for subclustering (all genes, CP10K)")
        obs_idx = np.where(mask.values)[0]
        raw_X = adata.raw.X[obs_idx]
        sub = sc.AnnData(
            X=raw_X.copy() if hasattr(raw_X, "copy") else raw_X,
            obs=adata.obs.iloc[obs_idx].copy(),
            var=adata.raw.var.copy(),
        )
        sub.obs_names = adata.obs_names[obs_idx]
    else:
        sub = adata[mask].copy()
    print(f"  LPM at {args.stage}: {sub.n_obs:,} cells, {sub.n_vars:,} genes from "
          f"{sub.obs[GENO_COL].nunique()} genotypes, "
          f"{sub.obs['library'].nunique()} libraries")

    # Drop libraries with < 2 cells (causes ZeroDivisionError in HVG)
    lib_counts = sub.obs["library"].value_counts()
    small_libs = lib_counts[lib_counts < 2].index.tolist()
    if small_libs:
        print(f"  Dropping {len(small_libs)} libraries with <2 cells: {small_libs}")
        sub = sub[~sub.obs["library"].isin(small_libs)].copy()

    # Re-compute HVGs on subset (using log-normalized data from .raw)
    print(f"  Computing {args.n_hvg} HVGs on subset...")
    sc.pp.highly_variable_genes(sub, n_top_genes=args.n_hvg, batch_key="library")
    print(f"  HVGs: {sub.var['highly_variable'].sum()}")

    # Save CP10K expression before scaling (for marker export)
    sub.raw = sub.copy()

    # Subset to HVGs before scaling (matching d1 pipeline)
    sub = sub[:, sub.var["highly_variable"]].copy()

    # Scale before PCA (matching d1 pipeline: regress_out + scale).
    # Default d2g regresses only pct_mito; cell cycle is retained so
    # proliferating cells can cluster apart from resting cells — this is
    # the default 2026-04-xx behavior.  With --regress-cc we additionally
    # regress S_score + G2M_score (matching d1's global behavior), letting
    # us test whether cell-cycle-based clusters are a collider artifact.
    regress_keys = []
    if "pct_mito" in sub.obs.columns:
        regress_keys.append("pct_mito")
    if args.regress_cc:
        for k in ("S_score", "G2M_score"):
            if k in sub.obs.columns and sub.obs[k].abs().sum() > 0:
                regress_keys.append(k)
            else:
                print(f"  WARNING: --regress-cc requested but {k} missing/zero; skipping")
    if regress_keys:
        print(f"  Regressing out: {regress_keys}")
        sc.pp.regress_out(sub, regress_keys)
    sc.pp.scale(sub, max_value=10)

    # PCA
    print(f"  PCA ({args.n_pcs} PCs)...")
    sc.tl.pca(sub, n_comps=args.n_pcs)

    # Harmony batch correction (use harmonypy directly, matching d1 pipeline).
    # harmonypy already defaults to random_state=0; the pre-fix nondeterminism
    # came from multithreaded BLAS/OpenMP jitter amplified by the n_init=1
    # kmeans, which made the embedding — and every leiden label downstream of
    # it — non-deterministic.  The fix pins threads (top of file) AND makes the
    # seed explicit (SEED=42, i.e. a different init than the implicit 0).
    print("  Harmony correction (batch=library)...")
    import harmonypy
    ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, "library",
                               max_iter_harmony=30, random_state=SEED)
    Z = ho.Z_corr
    if Z.shape[0] == sub.n_obs:
        sub.obsm["X_pca_harmony"] = Z
    else:
        sub.obsm["X_pca_harmony"] = Z.T
    use_rep = "X_pca_harmony"

    # Neighbors + UMAP
    print("  Computing neighbors and UMAP...")
    sc.pp.neighbors(sub, use_rep=use_rep, n_pcs=args.n_pcs)
    sc.tl.umap(sub, random_state=42)

    # Leiden at multiple resolutions
    for res in args.resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(sub, resolution=res, key_added=key, random_state=42)
        n_cl = sub.obs[key].nunique()
        print(f"  Leiden res={res}: {n_cl} clusters")

    # Extract marker expression from .raw (CP10K log-normalized, all genes)
    expr_source = sub.raw if sub.raw is not None else sub
    available = [m for m in MARKERS if m in expr_source.var_names]
    missing = [m for m in MARKERS if m not in expr_source.var_names]
    print(f"  Markers: {len(available)}/{len(MARKERS)} available")
    if missing:
        print(f"  Missing: {missing}")

    # Build output DataFrame
    out = sub.obs[[GENO_COL, "library", RES02, RES05, "stage"]].copy()
    out = out.rename(columns={GENO_COL: "genotype", RES02: "cluster_02", RES05: "cluster_05"})
    out["barcode"] = sub.obs_names

    # QC + cell cycle columns (investigative: test whether unregressed
    # cell cycle at d2g creates proliferating-bias clusters)
    for col in ("scDblFinder_class", "scDblFinder_score",
                "S_score", "G2M_score", "phase", "pct_mito", "pct_ribo"):
        if col in sub.obs.columns:
            out[col] = sub.obs[col].values

    # UMAP
    out["UMAP1"] = sub.obsm["X_umap"][:, 0]
    out["UMAP2"] = sub.obsm["X_umap"][:, 1]

    # Leiden columns
    for res in args.resolutions:
        key = f"leiden_{res}"
        out[key] = sub.obs[key].values

    # Marker expression (from CP10K log-normalized .raw)
    for marker in available:
        vals = expr_source[:, marker].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().ravel()
        else:
            vals = np.asarray(vals).ravel()
        name = marker.replace("LOC113843492", "NKX6_1")
        out[f"{name}_cp10k"] = vals

    # Save
    outfile = os.path.join(outdir, f"d2g_lpm_subclustering_{stage_label}.tsv.gz")
    out.to_csv(outfile, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {outfile}  shape={out.shape}")

    # Quick marker audit per cluster at res=0.3
    res_key = "leiden_0.3"
    print(f"\n=== Marker audit at res=0.3 ({sub.obs[res_key].nunique()} clusters) ===")
    key_audit = ["FOXF1", "NKX6_1", "PITX2", "TBX5", "TBX4", "ISL1",
                 "PRRX1", "TWIST1", "IRX3", "TBX1", "BARX1", "HAND1", "NR2F2"]
    header = f"{'cl':>3s} {'n':>6s}"
    for m in key_audit:
        header += f"  {m:>7s}"
    print(header)
    print("-" * (10 + 9 * len(key_audit)))

    for cl in sorted(out[res_key].unique(), key=int):
        cl_data = out[out[res_key] == cl]
        n = len(cl_data)
        row = f"{cl:>3s} {n:>6d}"
        for m in key_audit:
            col = f"{m}_cp10k"
            if col in cl_data.columns:
                frac = (cl_data[col] > 0).mean() * 100
                row += f"  {frac:6.1f}%"
            else:
                row += f"      --"
        print(row)

    # Marker-discovery triad (Stein 2026 U-method + Wilcoxon + lineage
    # scoring).  Emits three TSVs consumed by
    # /tmp/d2g_audit/annotate_via_compendium.py for compendium-based
    # cluster annotation.  This framework replaces earlier threshold-rule
    # classifiers (priority artifacts, threshold-gap mis-fits, silent drops).
    compute_marker_discovery(sub, stage_label, outdir,
                             resolutions=(0.2, 0.3))


# ── Marker-discovery triad (used for compendium annotation) ──────────────

# Canonical lineage panels — each is (up_markers, down_markers).  Signed
# score = score(up) - score(down), so panels with clear negative markers
# can discriminate against co-expressing clusters.  Duck NKX6-1 is
# annotated as LOC113843492 in the reference; keep the LOC symbol here to
# match raw gene names.
#
# Down-markers were added after v1 (unsigned) panels produced systematic
# confusions: Splanchnic_vis vs Ventral_spl (both ZFPM2-high), BARX1+
# vs Splanchnic_vis (single-gene BARX1 panel lost to multi-gene
# Splanchnic).  Positive-direction markers augmented with U-method-
# discovered UEGs from the v1 run (e.g. WNT2B/PTCH1 for Dorsal, TBX20
# for cardiac-adj Splanchnic, IRX1 for Cardio_Syringeal,
# RUNX2/WNT5A/DPP10 for Differentiating LPM, LRIG1/SMOC1 for BARX1+).
# Stage gating (Early_LPM → HH13 only; Differentiating → HH25 only)
# is applied in the compendium scorer rather than at score computation.
LINEAGE_PANELS = {
    "Dorsal_subepi":    (["LOC113843492", "WNT2B", "PTCH1"],
                         ["HAND2", "GATA4"]),
    "Ventral_spl":      (["TBX5", "TBX4", "ZFPM2", "ISL1"],
                         ["GATA4", "BARX1", "HAND2"]),
    "BARX1_pos":        (["BARX1", "LRIG1", "SMOC1"],
                         ["TBX5", "GATA4"]),
    "Cardio_Syringeal": (["HAND2", "ISL1", "IRX1"],
                         ["TBX5", "GATA4", "BARX1", "LOC113843492"]),
    "Splanchnic_vis":   (["GATA4", "ZFPM2", "TBX20"],
                         ["TBX5", "TBX4", "BARX1", "HAND2"]),
    "Somatic":          (["PRRX1", "TWIST1", "EBF2", "EYA1", "EYA2"],
                         ["GATA4", "ZFPM2", "FOXF1", "HAND2"]),
    "Epicardial":       (["WT1", "TCF21"],
                         ["FOXF1", "TBX5", "HAND2", "GATA4", "PRRX1"]),
    "Differentiating":  (["RUNX2", "WNT5A", "DPP10"],
                         ["ZFPM2", "GATA4", "HAND2", "LOC113843492"]),
    "Early_LPM":        (["TBX5", "FOXF1", "ISL1"],
                         []),
    "Sclerotome":       (["PAX1", "PAX9", "FOXC2"],
                         []),
    "Dermomyotome":     (["PAX3", "PAX7"],
                         []),
}


def compute_marker_discovery(sub, stage_label, outdir, resolutions=(0.2, 0.3)):
    """U-method UEGs + Wilcoxon DE + lineage scores, per (resolution, cluster).

    sub.raw must hold CP10K log-normalized expression for ALL genes — it
    is set before HVG subsetting + scaling upstream.
    """
    from scipy.sparse import issparse

    raw = sub.raw
    assert raw is not None, "sub.raw must be populated"
    X_raw = raw.X
    var_names = list(raw.var_names)

    # Lineage scoring: per-cell signed score (up-panel minus down-panel) per
    # panel.  Handles single-gene up-panels by using raw CP10K; multi-gene
    # panels via sc.tl.score_genes.  Down-panel is subtracted when present.
    def _score_panel(genes: list[str], key: str) -> bool:
        """Write per-cell panel score to sub.obs[key]; return True if any."""
        present = [g for g in genes if g in var_names]
        if not present:
            sub.obs[key] = 0.0
            return False
        if len(present) >= 2:
            sc.tl.score_genes(sub, gene_list=present, score_name=key, use_raw=True)
        else:
            j = var_names.index(present[0])
            vals = X_raw[:, j]
            if issparse(vals):
                vals = vals.toarray().ravel()
            else:
                vals = np.asarray(vals).ravel()
            # Normalise single-gene panels to z-score scale so they compare
            # fairly against sc.tl.score_genes output (which centers on a
            # matched control set).  Without this, single-gene panels win
            # on magnitude alone for highly-expressed markers like BARX1.
            mu = float(np.mean(vals))
            sd = float(np.std(vals)) or 1.0
            sub.obs[key] = (vals - mu) / sd
        return True

    for lin, (up_genes, down_genes) in LINEAGE_PANELS.items():
        up_key = f"_score_{lin}_up"
        down_key = f"_score_{lin}_down"
        have_up = _score_panel(up_genes, up_key)
        if not have_up:
            print(f"  WARNING: lineage panel {lin} has no up-genes — skipping")
            sub.obs[f"score_{lin}"] = 0.0
            continue
        _score_panel(down_genes, down_key)
        # Signed panel = up - down; down contributes 0 if panel is empty
        sub.obs[f"score_{lin}"] = sub.obs[up_key] - sub.obs[down_key]

    ueg_rows = []
    wil_rows = []
    score_rows = []
    TOP_N = 50

    for res in resolutions:
        res_key = f"leiden_{res}"
        if res_key not in sub.obs.columns:
            continue
        clusters = sorted(sub.obs[res_key].astype(str).unique(), key=lambda x: int(x))

        # Detection-rate matrix → U-scores (Stein 2026)
        det = np.zeros((len(clusters), raw.n_vars), dtype=np.float32)
        cluster_sizes = {}
        for i, c in enumerate(clusters):
            mask = (sub.obs[res_key].astype(str) == c).values
            cluster_sizes[c] = int(mask.sum())
            sub_X = X_raw[mask]
            if issparse(sub_X):
                det[i] = np.asarray((sub_X > 0).mean(axis=0)).flatten()
            else:
                det[i] = (sub_X > 0).mean(axis=0)
        u = np.zeros_like(det)
        for i in range(len(clusters)):
            other_max = (np.delete(det, i, axis=0).max(axis=0)
                         if len(clusters) > 1 else np.zeros(det.shape[1]))
            u[i] = det[i] - other_max
        for i, c in enumerate(clusters):
            order = np.argsort(u[i])[::-1][:TOP_N]
            for rank, j in enumerate(order, start=1):
                ueg_rows.append({
                    "stage_label": stage_label,
                    "resolution": res,
                    "cluster": c,
                    "n_cells": cluster_sizes[c],
                    "rank": rank,
                    "gene": var_names[j],
                    "u_score": round(float(u[i, j]), 4),
                    "p_in": round(float(det[i, j]), 4),
                    "p_out_max": round(float(det[i, j] - u[i, j]), 4),
                })

        # Wilcoxon DE (scanpy native; operates on sub.raw via use_raw=True)
        try:
            sc.tl.rank_genes_groups(sub, groupby=res_key, method="wilcoxon",
                                    n_genes=TOP_N, use_raw=True)
            for c in clusters:
                df = sc.get.rank_genes_groups_df(sub, group=c).head(TOP_N)
                df = df.rename(columns={"names": "gene", "scores": "wilcoxon_score",
                                         "logfoldchanges": "lfc",
                                         "pvals": "pval", "pvals_adj": "pval_adj"})
                df.insert(0, "stage_label", stage_label)
                df.insert(1, "resolution", res)
                df.insert(2, "cluster", c)
                df.insert(3, "n_cells", cluster_sizes[c])
                df.insert(4, "rank", range(1, len(df) + 1))
                wil_rows.extend(df.to_dict(orient="records"))
        except Exception as exc:
            raise RuntimeError(f"Wilcoxon DE failed for res={res}") from exc

        # Per-cluster mean lineage scores
        for c in clusters:
            mask = (sub.obs[res_key].astype(str) == c).values
            rec = {"stage_label": stage_label, "resolution": res,
                   "cluster": c, "n_cells": cluster_sizes[c]}
            for lin in LINEAGE_PANELS:
                key = f"score_{lin}"
                if key in sub.obs.columns:
                    rec[lin] = round(float(sub.obs.loc[mask, key].mean()), 4)
            score_rows.append(rec)

    for out_rows, filename in [
        (ueg_rows,   f"d2g_uegs_{stage_label}.tsv"),
        (wil_rows,   f"d2g_wilcoxon_{stage_label}.tsv"),
        (score_rows, f"d2g_lineage_scores_{stage_label}.tsv"),
    ]:
        path = os.path.join(outdir, filename)
        pd.DataFrame(out_rows).to_csv(path, sep="\t", index=False)
        print(f"  Wrote {path} ({len(out_rows)} rows)")


if __name__ == "__main__":
    main()
