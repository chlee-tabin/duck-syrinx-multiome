#!/usr/bin/env python3
"""cAL — chick stage-restricted LPM subclustering (duck `d2g` analog), pooled Group A+B.

Faithful port of `pipeline/d2g_stage_restricted_subclustering.py`.  The whole-embryo
atlas (cAE) shows the LPM-lineage mesenchyme is ~99% STAGE-segregated, not D/V-
segregated ("Mesenchyme (early)" = HH13/HH15-16; "LPM mesenchyme" = HH25).  So the
dorsal/ventral split that feeds the PITX2 laterality readout (duck fig02 analog) only
emerges when we restrict to ONE stage and recompute HVG/PCA/clustering — exactly the
situation the duck `d2g` exists for.

Recipe (per the duck d2g, pooling A+B and Harmony-correcting on `library` to absorb the
ARC-vs-GEX platform batch + lane):
  1. Load cAE_atlas.h5ad (.X = CP10K log-norm all-genes; layers['counts'] = raw).
  2. Attach FIRM per-cell stage: donor -> stage from chick_stage_alignment.tsv, with the
     Group B resolved-12 relabel (cS donor11/donor2 splits) via _groupB_resolved12.
  3. Select LPM-lineage (res0.5 {"LPM mesenchyme","Mesenchyme (early)"}).
  4. Per FIRM stage HH13 / HH15-16 / HH25 SEPARATELY:
       re-HVG(2000, batch=library, mask W/Z sex genes) -> scale (regress pct_mito;
       +S/G2M for the ccreg variant) -> PCA(30) -> Harmony(library) -> UMAP ->
       Leiden sweep [0.1..1.0].
  5. Export the per-cell TSV (markers as *_cp10k) + the marker-discovery TRIAD
     (Stein U-method UEGs + Wilcoxon DE + signed lineage scores) consumed by
     cAM_compendium_annotate.py.

Both the default and `ccreg` (cell-cycle-regressed) variants are produced, matching the
duck — the wrapper picks `ccreg` as canonical (collider-bias-safe).  Chick adaptations vs
duck: (a) atlas h5ad source instead of step04; (b) NKX6-1 is a real symbol (duck = LOC…);
(c) W/Z sex genes masked from HVG (strong ZW confound; duck regresses sex at d1); (d) firm
stage from the settled 24-embryo staging, not e1b PC1.

Submit (ABSOLUTE scratch path; sbatch copies the wrapper to /var/spool/slurmd):
  sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAL_stage_subclustering.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os

# ── Pin BLAS/OpenMP threading BEFORE the numpy-pulling imports (matches d1/d2g) ──
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import traceback
import random

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import issparse

# ── Global random seeds + single-threaded scanpy (matches d1/d2g) ──
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
sc.settings.n_jobs = 1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _groupB_resolved12 import resolve_donor, apply_split_meta, SPLIT_EXCLUDED  # noqa: E402

# ── paths (scratch) ──────────────────────────────────────────────────────────
OUT = f"{SCRATCH_USER}/chick_demux/vireo_seeded/seed_qc"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"        # Group B cS split files
H5AD = f"{OUT}/cAE_atlas.h5ad"
STAGE_TSV = f"{OUT}/chick_stage_alignment.tsv"           # copied next to the h5ad
WZ_DIRS = [f"{SCRATCH_USER}/chick_demux",       # chick_W_genes.tsv / chick_Z_genes.tsv
           f"{SCRATCH_USER}/chick_demux_B"]

LPM_LINEAGE = {"LPM mesenchyme", "Mesenchyme (early)"}
CT_COL = "celltype_res0.5"
STAGES = ["HH13", "HH15-16", "HH25"]
N_HVG = 2000
N_PCS = 30
RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
DISCOVERY_RES = (0.2, 0.3)                # canonical-candidate resolutions for the triad

# Markers exported as {gene}_cp10k.  The override-required set (must be present for the
# compendium HH13 + WT1/BARX1 overrides) is the duck's: NKX6_1, FOXF1, GATA4, ZFPM2, TBX5,
# TBX4, BARX1, ISL1, HAND2, PRRX1, WT1, NR2F2.  NKX6-1 is exported under the column NKX6_1.
MARKERS = [
    "NKX6-1", "FOXF1", "PITX2", "FOXF2", "GATA4", "ZFPM2", "TBX5", "TBX4",
    "ISL1", "HAND1", "HAND2", "BARX1", "PRRX1", "TWIST1", "NR2F2", "WT1", "TCF21",
    "WNT2B", "PTCH1", "MSX1", "IRX1", "IRX3", "TBX1", "TBX20",
    "LRIG1", "SMOC1", "RUNX2", "WNT5A", "DPP10", "EBF2", "EYA1", "EYA2",
    "PAX1", "PAX9", "FOXC2", "PAX3", "PAX7",
    "COL1A1", "COL3A1", "POSTN", "LUM", "DCN", "ALDH1A2", "BMP4",
]


def _g(*xs):
    """Expand a gene list so NKX6-1 also tries the NKX6_1 underscore spelling."""
    out = []
    for x in xs:
        out.append(x)
        if x == "NKX6-1":
            out.append("NKX6_1")
    return out


# Canonical lineage panels (up, down) — duck d2g LINEAGE_PANELS with chick symbols
# (NKX6-1 in place of the duck LOC113843492).  Signed score = score(up) - score(down).
LINEAGE_PANELS = {
    "Dorsal_subepi":    (_g("NKX6-1", "WNT2B", "PTCH1"),
                         _g("HAND2", "GATA4")),
    "Ventral_spl":      (_g("TBX5", "TBX4", "ZFPM2", "ISL1"),
                         _g("GATA4", "BARX1", "HAND2")),
    "BARX1_pos":        (_g("BARX1", "LRIG1", "SMOC1"),
                         _g("TBX5", "GATA4")),
    "Cardio_Syringeal": (_g("HAND2", "ISL1", "IRX1"),
                         _g("TBX5", "GATA4", "BARX1", "NKX6-1")),
    "Splanchnic_vis":   (_g("GATA4", "ZFPM2", "TBX20"),
                         _g("TBX5", "TBX4", "BARX1", "HAND2")),
    "Somatic":          (_g("PRRX1", "TWIST1", "EBF2", "EYA1", "EYA2"),
                         _g("GATA4", "ZFPM2", "FOXF1", "HAND2")),
    "Epicardial":       (_g("WT1", "TCF21"),
                         _g("FOXF1", "TBX5", "HAND2", "GATA4", "PRRX1")),
    "Differentiating":  (_g("RUNX2", "WNT5A", "DPP10"),
                         _g("ZFPM2", "GATA4", "HAND2", "NKX6-1")),
    "Early_LPM":        (_g("TBX5", "FOXF1", "ISL1"),
                         []),
    "Sclerotome":       (_g("PAX1", "PAX9", "FOXC2"),
                         []),
    "Dermomyotome":     (_g("PAX3", "PAX7"),
                         []),
}


def log(m):
    print(m, flush=True)


def _stage_label(stage, regress_cc):
    lab = stage.replace("/", "").replace("-", "")          # HH13 / HH1516 / HH25
    return lab + ("_ccreg" if regress_cc else "")


def load_wz_genes():
    wz = set()
    for d in WZ_DIRS:
        for chrom in ("W", "Z"):
            p = f"{d}/chick_{chrom}_genes.tsv"
            if os.path.exists(p):
                wz |= set(pd.read_csv(p, header=None, sep="\t")[0].astype(str))
    log(f"[wz] {len(wz)} W/Z sex genes to mask from HVG "
        f"({'OK' if wz else 'NONE FOUND — sex genes will not be masked'})")
    return wz


def attach_firm_stage(adata):
    """Resolved-12 Group B relabel + firm donor->stage from chick_stage_alignment.tsv."""
    gbmask = (adata.obs["group"].astype(str) == "B").values
    donor = adata.obs["donor"].astype(str).copy()
    donor.loc[gbmask] = resolve_donor(donor[gbmask], WORK_B, log=log).values
    adata.obs["donor"] = donor.values
    # cS split doublet/unassigned cells -> SPLIT_EXCLUDED sentinel: drop them (never embryo 'a')
    _xs = (adata.obs["donor"].astype(str) == SPLIT_EXCLUDED).values
    log(f"  dropping {int(_xs.sum())} Group B split doublet/unassigned cells (SPLIT_EXCLUDED)")
    adata = adata[~_xs].copy()
    apply_split_meta(adata.obs, donor_col="donor", stage_col="stage", log=log)  # sex of 11a/11b/2a/2b

    align = pd.read_csv(STAGE_TSV, sep="\t")
    firm = dict(zip(align["group"].astype(str) + ":" + align["donor"].astype(str),
                    align["stage"].astype(str)))
    key = adata.obs["group"].astype(str) + ":" + adata.obs["donor"].astype(str)
    adata.obs["firm_stage"] = key.map(firm)
    n_drop = int(adata.obs["firm_stage"].isna().sum())
    if n_drop:
        bad = adata.obs.loc[adata.obs["firm_stage"].isna(), "donor"].value_counts()
        log(f"[stage] dropping {n_drop} phantom/unmapped cells "
            f"({dict(bad.head(8))})")
    adata = adata[adata.obs["firm_stage"].notna()].copy()
    log("[stage] firm per-cell stage attached:\n"
        + pd.crosstab(adata.obs["firm_stage"], adata.obs["group"]).to_string())
    return adata


def add_pct_ribo(adata):
    """Ribosomal UMI fraction (percent) from the counts layer — for the apoptotic check."""
    if "pct_ribo" in adata.obs.columns:
        return
    var = np.asarray(adata.var_names)
    ribo_idx = np.where([str(g).upper().startswith(("RPL", "RPS")) for g in var])[0]
    counts = adata.layers["counts"]
    counts = counts.tocsc() if hasattr(counts, "tocsc") else counts
    tot = np.asarray(counts.sum(axis=1)).ravel()
    rsum = (np.asarray(counts[:, ribo_idx].sum(axis=1)).ravel()
            if len(ribo_idx) else np.zeros_like(tot))
    adata.obs["pct_ribo"] = np.where(tot > 0, rsum / tot, 0.0) * 100.0
    log(f"[ribo] {len(ribo_idx)} ribo genes; median pct_ribo={np.median(adata.obs['pct_ribo']):.2f}%")


# ── marker-discovery triad (duck d2g compute_marker_discovery, chick panels) ──
def compute_marker_discovery(sub, stage_label, outdir, resolutions=DISCOVERY_RES):
    """U-method UEGs + Wilcoxon DE + signed lineage scores, per (resolution, cluster).

    sub.raw must hold CP10K log-normalized expression for ALL genes.
    """
    raw = sub.raw
    assert raw is not None, "sub.raw must be populated"
    X_raw = raw.X
    var_names = list(raw.var_names)

    def _score_panel(genes, key):
        present = [g for g in genes if g in var_names]
        if not present:
            sub.obs[key] = 0.0
            return False
        if len(present) >= 2:
            sc.tl.score_genes(sub, gene_list=present, score_name=key, use_raw=True)
        else:
            j = var_names.index(present[0])
            vals = X_raw[:, j]
            vals = vals.toarray().ravel() if issparse(vals) else np.asarray(vals).ravel()
            mu = float(np.mean(vals))
            sd = float(np.std(vals)) or 1.0
            sub.obs[key] = (vals - mu) / sd
        return True

    for lin, (up_genes, down_genes) in LINEAGE_PANELS.items():
        up_key, down_key = f"_score_{lin}_up", f"_score_{lin}_down"
        have_up = _score_panel(up_genes, up_key)
        if not have_up:
            log(f"  WARNING: lineage panel {lin} has no up-genes — skipping")
            sub.obs[f"score_{lin}"] = 0.0
            continue
        _score_panel(down_genes, down_key)
        sub.obs[f"score_{lin}"] = sub.obs[up_key] - sub.obs[down_key]

    ueg_rows, wil_rows, score_rows = [], [], []
    TOP_N = 50

    for res in resolutions:
        res_key = f"leiden_{res}"
        if res_key not in sub.obs.columns:
            continue
        clusters = sorted(sub.obs[res_key].astype(str).unique(), key=lambda x: int(x))

        det = np.zeros((len(clusters), raw.n_vars), dtype=np.float32)
        cluster_sizes = {}
        for i, c in enumerate(clusters):
            mask = (sub.obs[res_key].astype(str) == c).values
            cluster_sizes[c] = int(mask.sum())
            sub_X = X_raw[mask]
            det[i] = (np.asarray((sub_X > 0).mean(axis=0)).flatten()
                      if issparse(sub_X) else (sub_X > 0).mean(axis=0))
        u = np.zeros_like(det)
        for i in range(len(clusters)):
            other_max = (np.delete(det, i, axis=0).max(axis=0)
                         if len(clusters) > 1 else np.zeros(det.shape[1]))
            u[i] = det[i] - other_max
        for i, c in enumerate(clusters):
            order = np.argsort(u[i])[::-1][:TOP_N]
            for rank, j in enumerate(order, start=1):
                ueg_rows.append({
                    "stage_label": stage_label, "resolution": res, "cluster": c,
                    "n_cells": cluster_sizes[c], "rank": rank, "gene": var_names[j],
                    "u_score": round(float(u[i, j]), 4),
                    "p_in": round(float(det[i, j]), 4),
                    "p_out_max": round(float(det[i, j] - u[i, j]), 4),
                })

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
        (ueg_rows,   f"cAL_uegs_{stage_label}.tsv"),
        (wil_rows,   f"cAL_wilcoxon_{stage_label}.tsv"),
        (score_rows, f"cAL_lineage_scores_{stage_label}.tsv"),
    ]:
        path = os.path.join(outdir, filename)
        pd.DataFrame(out_rows).to_csv(path, sep="\t", index=False)
        log(f"  wrote {path} ({len(out_rows)} rows)")


# ── per-stage subclustering ──────────────────────────────────────────────────
def subcluster_stage(adata, stage, regress_cc, wz_genes):
    stage_label = _stage_label(stage, regress_cc)
    mask = (adata.obs["firm_stage"] == stage).values & adata.obs[CT_COL].isin(LPM_LINEAGE).values
    sub = adata[mask].copy()
    log(f"\n=== {stage_label}: {sub.n_obs:,} LPM-lineage cells "
        f"({sub.obs['group'].value_counts().to_dict()}, "
        f"{sub.obs['library'].nunique()} libraries) ===")
    if sub.n_obs < 100:
        log(f"  only {sub.n_obs} cells — skipped")
        return

    # drop libraries with <2 cells (HVG batch_key ZeroDivision guard)
    lc = sub.obs["library"].value_counts()
    small = lc[lc < 2].index.tolist()
    if small:
        log(f"  dropping {len(small)} libraries with <2 cells: {small}")
        sub = sub[~sub.obs["library"].isin(small)].copy()

    # .X is CP10K log-norm (all genes); preserve it as .raw for marker export + discovery
    sub.raw = sub.copy()

    sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, batch_key="library")
    if wz_genes:
        n_before = int(sub.var["highly_variable"].sum())
        sub.var.loc[sub.var_names.isin(wz_genes), "highly_variable"] = False
        log(f"  HVGs: {int(sub.var['highly_variable'].sum())} "
            f"(masked {n_before - int(sub.var['highly_variable'].sum())} W/Z sex genes)")
    sub = sub[:, sub.var["highly_variable"]].copy()

    regress_keys = ["pct_mito"] if "pct_mito" in sub.obs.columns else []
    if regress_cc:
        for k in ("S_score", "G2M_score"):
            if k in sub.obs.columns and sub.obs[k].abs().sum() > 0:
                regress_keys.append(k)
    if regress_keys:
        log(f"  regressing out: {regress_keys}")
        sc.pp.regress_out(sub, regress_keys)
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=N_PCS)

    log("  Harmony (batch=library)...")
    import harmonypy
    ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, "library",
                               max_iter_harmony=30, random_state=SEED)
    Z = ho.Z_corr
    sub.obsm["X_pca_harmony"] = Z if Z.shape[0] == sub.n_obs else Z.T

    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_pcs=N_PCS)
    sc.tl.umap(sub, random_state=42)
    for res in RESOLUTIONS:
        sc.tl.leiden(sub, resolution=res, key_added=f"leiden_{res}", random_state=42)
        log(f"  Leiden res={res}: {sub.obs[f'leiden_{res}'].nunique()} clusters")

    # ── per-cell export ──
    keep_obs = ["group", "library", "lane", "donor", "firm_stage", "sex",
                "pct_mito", "pct_ribo", "S_score", "G2M_score", "phase"]
    out = sub.obs[[c for c in keep_obs if c in sub.obs.columns]].copy()
    out = out.rename(columns={"firm_stage": "stage"})
    out["barcode"] = sub.obs_names
    out["UMAP1"] = sub.obsm["X_umap"][:, 0]
    out["UMAP2"] = sub.obsm["X_umap"][:, 1]
    for res in RESOLUTIONS:
        out[f"leiden_{res}"] = sub.obs[f"leiden_{res}"].values

    # NB: sub.raw.X is CP10K *log-norm* (log1p), so the exported `{gene}_cp10k` columns hold
    # log1p(CP10K), NOT linear cp10k. Detection (>0) is transform-invariant; any linear-space mean
    # must expm1 first (cAN/cb_floor_check do). The suffix is kept for continuity — see the §8
    # provenance note in docs/chick_control/laterality_fold_methods.md (not distributed).
    expr = sub.raw
    avail = [m for m in MARKERS if m in expr.var_names]
    miss = [m for m in MARKERS if m not in expr.var_names]
    log(f"  markers: {len(avail)}/{len(MARKERS)} available" + (f"; missing {miss}" if miss else ""))
    for m in avail:
        v = expr[:, m].X
        v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
        out[f"{m.replace('-', '_')}_cp10k"] = v   # NKX6-1 -> NKX6_1_cp10k (log1p(CP10K) value)

    fn = f"{OUT}/cAL_subclustering_{stage_label}.tsv.gz"
    out.to_csv(fn, sep="\t", index=False, compression="gzip")
    log(f"  saved {fn}  shape={out.shape}")

    # marker audit at res=0.3
    rk = "leiden_0.3"
    log(f"\n  --- marker audit {stage_label} at res=0.3 ({out[rk].nunique()} clusters) ---")
    audit = ["FOXF1", "NKX6_1", "PITX2", "TBX5", "TBX4", "ISL1", "ZFPM2",
             "BARX1", "HAND2", "WT1", "PRRX1", "PAX1", "PAX3"]
    hdr = f"  {'cl':>3s} {'n':>6s}" + "".join(f" {m:>7s}" for m in audit)
    log(hdr)
    for cl in sorted(out[rk].unique(), key=int):
        d = out[out[rk] == cl]
        row = f"  {cl:>3s} {len(d):>6d}"
        for m in audit:
            col = f"{m}_cp10k"
            row += f" {(d[col] > 0).mean() * 100:6.1f}%" if col in d.columns else "      --"
        log(row)

    compute_marker_discovery(sub, stage_label, OUT, resolutions=DISCOVERY_RES)


def main():
    log("=== cAL: chick stage-restricted LPM subclustering (duck d2g analog, pooled A+B) ===")
    adata = sc.read_h5ad(H5AD)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    assert "counts" in adata.layers, "cAE h5ad has no layers['counts']"
    assert CT_COL in adata.obs, f"{CT_COL} not in obs"
    if adata.raw is None:           # cAE sets .raw = CP10K log-norm; .X is the same
        adata.raw = adata
    adata = attach_firm_stage(adata)
    add_pct_ribo(adata)
    wz_genes = load_wz_genes()

    n_lpm = int(adata.obs[CT_COL].isin(LPM_LINEAGE).sum())
    log(f"[lpm] {n_lpm:,} LPM-lineage cells ({', '.join(sorted(LPM_LINEAGE))})")

    for regress_cc in (False, True):
        for stage in STAGES:
            subcluster_stage(adata, stage, regress_cc, wz_genes)   # a failed stage raises
    log("\n=== cAL done ===")


if __name__ == "__main__":
    main()
