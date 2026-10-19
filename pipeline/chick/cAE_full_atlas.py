#!/usr/bin/env python3
"""cAE — chick WHOLE-EMBRYO cell-type atlas (figS04 Panel B/D), Group A + Group B combined.

Faithful analog of the duck fig01 atlas pipeline (d1 normalize/PCA + d2 leiden + d2b annotation).
This v2 rewrite (2026-06-08, figS04 review) brings the chick atlas to METHODS-SECTION PARITY with
the duck d1->d2->d2b recipe, so a reader sees no analytic divergence except the one justified
chick-specific addition (Harmony-on-library, to correct the Group-A-vs-Group-B technical axis that
the single-cohort duck data does not have).

Duck-parity preprocessing (mirrors pipeline/d1_normalize_pca.py exactly):
  - SEED=42, BLAS/OpenMP pinned to 1 thread, sc.settings.n_jobs=1  (determinism)
  - per-cell sex_numeric (M=0 / F=1) from the established per-donor calls (Group A cT, Group B
    c5h; both REQUIRED) -> regressed out, just like duck. The 51 declared demux-remnant nuclei
    (_chick_design_rules.DEMUX_REMNANTS: Group A donor1/12/17, Group B donor5/13 survive QC) are
    KEPT in the atlas as "unassigned remnant" nuclei with the declared values stage/sex "?" and
    sex_numeric 0.5; any other donor missing from cT / c5h / the Group B stage map RAISES
  - HVG: seurat_v3_paper, n_top_genes=3000, batch_key="library", on RAW counts
  - normalize_total(1e4) + log1p ; .raw stored
  - cell-cycle scoring (duck S/G2M lists; chick uses the same avian symbols)
  - regress_out(["pct_mito", "sex_numeric", "S_score", "G2M_score"]) ; scale(max_value=10)
  - PCA n_comps=50, arpack, random_state=SEED
  - + Harmony on "library"  (THE chick-specific addition; duck has no A-vs-B batch)
  - neighbors(15, n_pcs=30, use_rep=X_pca_harmony) ; UMAP

Clustering + annotation (mirrors d2 + d2b):
  - leiden at res {0.2, 0.3, 0.5}, random_state=42  (0.2 = the duck fig01 atlas resolution)
  - annotation METHOD = duck d2b scoring (per-cluster score = mean fraction-expressing across a
    cell type's markers) but ARGMAX per cluster (NOT the duck greedy one-per-type, which mislabels
    the chick LPM blob "Splanchnic"); <0.01 floor -> "Unassigned_<cl>"; a "Mesenchyme (early)"
    set for the stage-split early mesenchyme; and a "Low-quality" override for clusters whose
    top-10 DE is ribo/mito-dominated. (Duck divergence: docs/chick_control/ANNOTATION_DIFF.md, not distributed.)
  - per-cluster rank_genes_groups DE (wilcoxon) + a curated paraxial/PSM probe score, exported so
    a suspect/Unassigned cluster can be adjudicated from markers (ground truth), not spatial guess.

Outputs (scratch chick_demux/vireo_seeded/seed_qc/):
  cAE_atlas.h5ad                        full annotated object (.raw + leiden_res* + celltype_res* + UMAP)
  cAE_atlas_umap.tsv.gz                 per-cell: barcode, group, library, lane, donor, stage, sex,
                                          leiden_res*, celltype_res*, UMAP1/2
  cAE_atlas_cluster_annotation_res{r}.tsv  cluster x celltype score matrix + d2b call + paraxial probe
  cAE_atlas_DE_res{r}.tsv               per-cluster top-25 rank_genes_groups DE
  cAE_atlas_celltype_markers.tsv        dotplot data (res=0.2 celltype x duck-panel marker)
  cAE_atlas_qc.png                      UMAP by celltype(res0.2)/group/leiden (quick review)

Submit (ABSOLUTE scratch paths; sbatch copies the wrapper to /var/spool/slurmd):
  sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAE_atlas.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os
# ── Pin BLAS/OpenMP threading for reproducibility (duck d1) ──────────────────
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import traceback, random
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Shared chick-atlas annotation (marker panel + argmax + Low-quality + Sclerotome promotion),
# imported so cAE (pooled A+B) and cAG (Group-B-only) label identically. Same-dir import works
# because the script's directory is on sys.path[0] when run directly (sbatch uses an absolute path).
from _chick_matrix_io import load_anndata          # shared cellranger-matrix loader
from _chick_atlas_annotation import (
    CELLTYPE_MARKERS, PARAXIAL_PROBE, MITO_EXACT,
    annotate_clusters, lowq_from_de, promote_sclerotome,
)
# Declared demux-remnant allowlist + their declared per-cell values (single source).
from _chick_design_rules import (
    DEMUX_REMNANTS, REMNANT_LABEL, REMNANT_STAGE, REMNANT_SEX, REMNANT_SEX_NUMERIC,
    is_remnant, lookup_or_remnant, c5h_table_path, read_donor_table, read_groupA_table,
)

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
sc.settings.n_jobs = 1
sc.settings.verbosity = 1

WORK_A = f"{SCRATCH_USER}/chick_demux"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"
OUT = f"{WORK_A}/vireo_seeded/seed_qc"
# matrix subdir to load: "matrix" (default, raw cellranger filtered) or "matrix_cb"
# (CellBender-corrected, re-staged by cb_restage.py). Set CHICK_MATRIX_SUBDIR=matrix_cb
# for the ambient-corrected / duck-parity atlas rerun.
MATRIX_SUBDIR = os.environ.get("CHICK_MATRIX_SUBDIR", "matrix")
CT = f"{OUT}/cT_groupA_stage_sex.tsv"        # Group A per-donor stage + sex (ARC0S-triangulated 4/4/4)
# Group B per-donor sex (N14/seed105): read from the c5h PRODUCER path ({WORK_B}/vireo_seeded/
# seed_qc), resolved (and checked against any stale copy in OUT) in load_donor_tables().
S_GENES_FILE = f"{OUT}/duck_s_genes.txt"
G2M_GENES_FILE = f"{OUT}/duck_g2m_genes.txt"

RESOLUTIONS = [0.2, 0.3, 0.5]
# res=0.5 is the chick's FEATURED resolution: at the chick's ~5x lower depth vs the duck,
# res=0.2 folds sclerotome + dorsal-neural into broad parents; res=0.5 separates them (and the
# Sclerotome margin-promotion fires there). Argmax homogenises the LPM into one label regardless,
# so res=0.5 does not over-fragment the displayed types — it recovers more of the duck roster.
DEFAULT_RES = 0.5

# (group, WORK, donor_ids, [(library, barcode_suffix), ...])
SOURCES = [
    ("A", WORK_A, f"{WORK_A}/vireo_seeded/N18/seed_117/donor_ids.tsv",
     [("ARC01", "1"), ("ARC02", "2")]),
    ("B", WORK_B, f"{WORK_B}/vireo_seeded/N14/seed_105/donor_ids.tsv",
     [("GEX09_GEX", "9"), ("GEX10_GEX", "10")]),
]
MIN_DEMUX_COVERAGE = 0.95    # min fraction of a library's cells that must join its donor_ids.tsv
LANE_NAME = {"ARC01": "ARC01", "ARC02": "ARC02", "GEX09_GEX": "GEX09", "GEX10_GEX": "GEX10"}

# CELLTYPE_MARKERS, PARAXIAL_PROBE, MITO_EXACT, annotate_clusters, lowq_from_de,
# promote_sclerotome are imported from _chick_atlas_annotation (shared with cAG).


def log(m): print(m, flush=True)


def load_library(work, library, suffix):
    """Load one library's GEX matrix -> AnnData (RAW counts), tagged with its lane name.
    `matrix` (raw cellranger) vs `matrix_cb` (CellBender) chosen via MATRIX_SUBDIR."""
    a = load_anndata(work, library, suffix, MATRIX_SUBDIR)
    a.obs["library"] = LANE_NAME[library]
    log(f"  [{library}] {a.n_obs} cells x {a.n_vars} GEX genes")
    return a


def group_b_stage_map():
    """Coarse (base-donor -> stage) for Group B (colouring attribute only)."""
    early = {"donor0", "donor1", "donor4", "donor6", "donor7", "donor9", "donor10"}
    late = {"donor2", "donor8"}
    out = {d: "early" for d in early}
    out.update({d: "HH25" for d in late})
    out["donor11"] = "early/late"
    return out


def load_donor_tables():
    """Per-donor stage + sex tables, keyed by group then base donor id. cT (Group A stage + sex)
    and c5h (Group B sex) are REQUIRED — a missing table raises instead of silently giving every
    cell sex '?' / sex_numeric 0.5. Group B stage = the coarse group_b_stage_map()."""
    a = read_groupA_table(CT)          # unique donors, stages in GROUPA_STAGES at 4/4/4, sex M/F
    b = read_donor_table(c5h_table_path(WORK_B, legacy_dir=OUT), ["donor", "sex"], "c5h Group B")
    stage = {"A": dict(zip(a["donor"], a["stage"])), "B": group_b_stage_map()}
    sex = {"A": dict(zip(a["donor"], a["sex"])),
           "B": dict(zip(b["donor"], b["sex"].str.upper()))}
    return stage, sex


# ════════════════════════════════════════════════════════════════════════════
# 1. load + concat (inner gene intersection across all 4 refs); keep RAW counts
# ════════════════════════════════════════════════════════════════════════════
parts = []
for grp, work, bdir, libs in SOURCES:
    dmap = pd.read_csv(bdir, sep="\t").set_index("cell").donor_id.to_dict()
    for library, suf in libs:
        a = load_library(work, library, suf)
        a.obs["group"] = grp
        a.obs["donor"] = a.obs_names.map(dmap)
        # base-demux join coverage (data integrity): a barcode-suffix mismatch must raise, not
        # silently drop the whole library at the donor.notna() filter below. Expected ~100% (per
        # lane the donor_ids row count equals the filtered-matrix cell count, cb_libqc_perlane.tsv).
        _cov = float(a.obs["donor"].notna().mean()) if a.n_obs else 0.0
        if _cov < MIN_DEMUX_COVERAGE:
            raise ValueError(
                f"{library}: only {int(a.obs['donor'].notna().sum())}/{a.n_obs} ({_cov:.1%}) cells join "
                f"{bdir} (< {MIN_DEMUX_COVERAGE:.0%}; barcode-suffix mismatch?); "
                f"e.g. atlas {list(a.obs_names[:3])} vs demux {list(dmap)[:3]}")
        parts.append(a)

common = sorted(set.intersection(*[set(a.var_names) for a in parts]))
log(f"[genes] {len(common)} common to all {len(parts)} libraries")
parts = [a[:, common].copy() for a in parts]
adata = ad.concat(parts, join="inner")
adata.obs_names_make_unique()
adata.obs["lane"] = adata.obs["library"]
log(f"[concat] {adata.n_obs} cells x {adata.n_vars} genes "
    f"(A={int((adata.obs.group=='A').sum())}, B={int((adata.obs.group=='B').sum())})")

# drop doublet / unassigned / no-donor (nothing stage/laterality-based: atlas = ALL demuxed A+B)
adata = adata[adata.obs.donor.notna() & ~adata.obs.donor.isin(["doublet", "unassigned"])].copy()
log(f"[demuxed] {adata.n_obs} cells")

# ── stage (Group A cT 'stage'; Group B coarse) + sex (cT + c5h -> sex_numeric) ──
# Declared rule: a real embryo MUST be in the per-donor tables (else raise); an allowlisted demux
# remnant gets the declared REMNANT_* values (kept in the atlas, excluded from per-embryo analyses).
stage_tab, sex_tab = load_donor_tables()
pair_stage, pair_sex = {}, {}
for g, d in set(zip(adata.obs.group.astype(str), adata.obs.donor.astype(str))):
    pair_stage[(g, d)] = lookup_or_remnant(stage_tab[g], g, d, "stage", REMNANT_STAGE)
    sx = lookup_or_remnant(sex_tab[g], g, d, "sex", REMNANT_SEX)
    if not is_remnant(g, d) and sx not in ("M", "F"):
        raise ValueError(f"sex: Group {g} {d} has non-M/F sex {sx!r} in its per-donor table")
    pair_sex[(g, d)] = sx
pairs = list(zip(adata.obs.group.astype(str), adata.obs.donor.astype(str)))
adata.obs["stage"] = [pair_stage[p] for p in pairs]
sex = [pair_sex[p] for p in pairs]
adata.obs["sex"] = sex
adata.obs["sex_numeric"] = [REMNANT_SEX_NUMERIC if is_remnant(*p) else {"M": 0.0, "F": 1.0}[s]
                            for p, s in zip(pairs, sex)]
rem = pd.Series([p for p in pairs if is_remnant(*p)], dtype=object).value_counts()
log(f"[remnant] {int(rem.sum())} pre-QC nuclei from declared demux remnants ({REMNANT_LABEL}; "
    f"stage/sex {REMNANT_STAGE!r}, sex_numeric {REMNANT_SEX_NUMERIC}): "
    + ", ".join(f"{g}:{d}={n}" for (g, d), n in rem.items()))
log(f"[sex] M={int((adata.obs.sex=='M').sum())}, F={int((adata.obs.sex=='F').sum())}, "
    f"?={int((adata.obs.sex=='?').sum())}")

# ════════════════════════════════════════════════════════════════════════════
# 2. QC + pct_mito (on raw counts) ; keep counts layer for seurat_v3_paper HVG
# ════════════════════════════════════════════════════════════════════════════
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
mito = [g for g in adata.var_names if g.upper() in MITO_EXACT or g.upper().startswith("MT-")]
if mito:
    tot = np.asarray(adata.X.sum(axis=1)).ravel()
    mt = np.asarray(adata[:, mito].X.sum(axis=1)).ravel()
    adata.obs["pct_mito"] = np.where(tot > 0, mt / tot, 0.0)
    log(f"[mito] {len(mito)} mito genes {sorted(mito)}; median pct_mito={np.median(adata.obs['pct_mito']):.4f}")
else:
    raise SystemExit("[mito] no mito genes matched (pct_mito regression would be a no-op)")
for g in ("EGFP", "eGFP", "GFP"):
    if g in adata.var_names:
        adata = adata[:, adata.var_names != g].copy()
adata.layers["counts"] = adata.X.copy()
log(f"[qc] {adata.n_obs} cells x {adata.n_vars} genes after filter")
_rq = adata.obs[[is_remnant(g, d) for g, d in zip(adata.obs.group, adata.obs.donor)]]
log(f"[remnant] {len(_rq)} post-QC {REMNANT_LABEL} nuclei kept in the atlas: "
    + ", ".join(f"{g}:{d}={n}" for (g, d), n in _rq.groupby(["group", "donor"], observed=True).size().items()))

# ── normalize + log1p ; store .raw (log-norm, all genes) ─────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ── cell-cycle scoring (duck S/G2M lists; chick = same avian symbols) ─────────
# no fallback: the published atlas (2026-06-11 h5ad) carries real S/G2M scores for all 37,362
# cells, so a missing gene list or too few matched genes must raise, never skip CC regression
def _read_genes(path):
    return [g.strip() for g in open(path) if g.strip() and g.strip() in adata.var_names]
s_genes, g2m_genes = _read_genes(S_GENES_FILE), _read_genes(G2M_GENES_FILE)
if len(s_genes) <= 3 or len(g2m_genes) <= 3:
    raise RuntimeError(f"[cc] too few CC genes matched (S={len(s_genes)} from {S_GENES_FILE}, "
                       f"G2M={len(g2m_genes)} from {G2M_GENES_FILE})")
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, random_state=SEED)
cc_scored = True
log(f"[cc] scored: {len(s_genes)} S + {len(g2m_genes)} G2M genes")

# ════════════════════════════════════════════════════════════════════════════
# 3. HVG (seurat_v3_paper on raw counts, per-library) — duck d1 parity
# ════════════════════════════════════════════════════════════════════════════
# no fallback: the published atlas recorded uns['hvg']['flavor'] = 'seurat_v3_paper' — let it raise
hvg_flavor = "seurat_v3_paper"
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3_paper",
                            layer="counts", batch_key="library")
n_hvg = int(adata.var.highly_variable.sum())
log(f"[hvg] {n_hvg} HVGs ({hvg_flavor}, batch_key=library)")

ad_h = adata[:, adata.var.highly_variable].copy()

# ── regress out (pct_mito + sex + cell cycle) + scale — duck d1 parity ───────
regress_keys = ["pct_mito", "sex_numeric"]
if cc_scored:
    regress_keys += ["S_score", "G2M_score"]
log(f"[regress] {regress_keys}")
sc.pp.regress_out(ad_h, regress_keys)
sc.pp.scale(ad_h, max_value=10)

# ── PCA (50, arpack, seeded) ─────────────────────────────────────────────────
sc.tl.pca(ad_h, n_comps=50, svd_solver="arpack", random_state=SEED)

# ════════════════════════════════════════════════════════════════════════════
# 4. Harmony on 'library' — THE chick-specific addition (A-vs-B + lane batch)
# ════════════════════════════════════════════════════════════════════════════
harmony_applied = False
rep = "X_pca"
try:
    import harmonypy
    ho = harmonypy.run_harmony(ad_h.obsm["X_pca"], ad_h.obs, ["library"], random_state=SEED)
    Z = np.asarray(ho.Z_corr).T
    if Z.shape[0] != ad_h.n_obs:
        Z = np.asarray(ho.Z_corr)
    assert Z.shape[0] == ad_h.n_obs, f"harmony Z {Z.shape} != n_obs {ad_h.n_obs}"
    ad_h.obsm["X_pca_harmony"] = Z
    rep = "X_pca_harmony"
    harmony_applied = True
    log(f"[harmony] integrated on 'library' (Z {Z.shape})")
except Exception as e:
    try:
        sc.external.pp.harmony_integrate(ad_h, "library")
        rep = "X_pca_harmony"
        harmony_applied = True
        log(f"[harmony] direct failed ({type(e).__name__}: {e}); scanpy wrapper ok")
    except Exception as e2:
        # no silent fallback to uncorrected X_pca: the atlas is published as Harmony-integrated
        raise RuntimeError(f"[harmony] both direct harmonypy ({type(e).__name__}: {e}) and the scanpy "
                           f"wrapper ({type(e2).__name__}: {e2}) failed") from e2

# ── neighbours (15, 30 PCs of the chosen rep) + UMAP ─────────────────────────
n_pcs = min(30, ad_h.obsm[rep].shape[1])
sc.pp.neighbors(ad_h, use_rep=rep, n_neighbors=15, n_pcs=n_pcs, random_state=SEED)
# no fallback to raw Harmony dims: the published S12B panel is a real UMAP — let sc.tl.umap raise
sc.tl.umap(ad_h, random_state=SEED)
adata.obsm["X_umap"] = ad_h.obsm["X_umap"]
umap_applied = True
adata.obsm["X_pca"] = ad_h.obsm["X_pca"]
if harmony_applied:
    adata.obsm["X_pca_harmony"] = ad_h.obsm["X_pca_harmony"]

# ════════════════════════════════════════════════════════════════════════════
# 5. leiden at {0.2,0.3,0.5} + d2b annotation + DE + paraxial probe (per res)
# ════════════════════════════════════════════════════════════════════════════
os.makedirs(OUT, exist_ok=True)
sc.tl.score_genes(adata, [g for g in PARAXIAL_PROBE if g in adata.raw.var_names],
                  score_name="paraxial_probe", use_raw=True)
celltype_cols = {}
for res in RESOLUTIONS:
    lk = f"leiden_res{res}"
    sc.tl.leiden(ad_h, resolution=res, key_added=lk, random_state=SEED)
    adata.obs[lk] = ad_h.obs[lk].values
    # DE first (wilcoxon, top 25): feeds BOTH the low-quality flag and the export below.
    # no fallback: a failed DE would silently skip the Low-quality flag — let it raise
    sc.tl.rank_genes_groups(adata, lk, method="wilcoxon", use_raw=True, n_genes=25)
    rgg = adata.uns["rank_genes_groups"]
    lowq = lowq_from_de(rgg)
    # ARGMAX annotation (+ Low-quality override) then Sclerotome margin-promotion (shared with cAG).
    # The promotion is resolution-dependent by construction: at res=0.2 the sclerotome is diluted
    # (score ~0.16, no promotion); at res=0.5 it concentrates (~0.37, PAX1 logFC ~6.7, CNMD/KERA
    # chondrogenic DE) and promotes — which is why res=0.5 is the chick's featured atlas.
    assignments, score_mat = annotate_clusters(adata, lk, lowq=lowq)
    paraxial_by_cl = adata.obs.groupby(lk, observed=True)["paraxial_probe"].mean()
    promote_sclerotome(assignments, score_mat, paraxial_by_cl, log=log)
    ck = f"celltype_res{res}"
    adata.obs[ck] = adata.obs[lk].map(assignments)
    celltype_cols[res] = ck
    # per-cluster annotation export (scores + paraxial probe + n + call)
    ann = score_mat.copy()
    ann.insert(0, "n_cells", adata.obs[lk].value_counts())
    ann["paraxial_probe"] = adata.obs.groupby(lk, observed=True)["paraxial_probe"].mean()
    ann["celltype_call"] = pd.Series(assignments)
    ann.index.name = "leiden"
    ann.round(4).to_csv(f"{OUT}/cAE_atlas_cluster_annotation_res{res}.tsv", sep="\t")
    n_un = sum(v.startswith("Unassigned") for v in assignments.values())
    n_lq = sum(v == "Low-quality" for v in assignments.values())
    log(f"[res {res}] {adata.obs[lk].nunique()} clusters -> "
        f"{adata.obs[ck].nunique()} celltypes ({n_un} Unassigned, {n_lq} Low-quality)")
    # DE export (reuse rgg computed above)
    if rgg is not None:
        de_rows = []
        for cl in rgg["names"].dtype.names:
            for rank in range(len(rgg["names"][cl])):
                de_rows.append(dict(cluster=cl, rank=rank + 1,
                                    gene=rgg["names"][cl][rank],
                                    logFC=round(float(rgg["logfoldchanges"][cl][rank]), 3),
                                    pval_adj=float(rgg["pvals_adj"][cl][rank]),
                                    celltype_call=assignments[cl]))
        pd.DataFrame(de_rows).to_csv(f"{OUT}/cAE_atlas_DE_res{res}.tsv", sep="\t", index=False)

# ════════════════════════════════════════════════════════════════════════════
# 6. exports
# ════════════════════════════════════════════════════════════════════════════
um = adata.obsm["X_umap"]
cells = pd.DataFrame({
    "barcode": adata.obs_names, "group": adata.obs.group.values,
    "library": adata.obs.library.values, "lane": adata.obs.lane.values,
    "donor": adata.obs.donor.values, "stage": adata.obs.stage.values,
    "sex": adata.obs.sex.values,
    "UMAP1": um[:, 0].round(3), "UMAP2": um[:, 1].round(3)})
for res in RESOLUTIONS:
    cells[f"leiden_res{res}"] = adata.obs[f"leiden_res{res}"].values
    cells[f"celltype_res{res}"] = adata.obs[f"celltype_res{res}"].values
cells.to_csv(f"{OUT}/cAE_atlas_umap.tsv.gz", sep="\t", index=False, compression="gzip")

# dotplot data (default-res celltype x duck-panel marker)
dck = celltype_cols[DEFAULT_RES]
dot_rows = []
raw = adata.raw
for ct in sorted(adata.obs[dck].dropna().unique()):
    sub_idx = np.where((adata.obs[dck] == ct).values)[0]
    if len(sub_idx) == 0:
        continue
    for mct, gs in CELLTYPE_MARKERS.items():
        for g in gs:
            if g not in raw.var_names:
                continue
            col = np.asarray(raw[:, g].X[sub_idx].todense()).ravel()
            dot_rows.append(dict(celltype=ct, marker=g, marker_celltype=mct,
                                 mean_expr=round(float(col.mean()), 4),
                                 pct_expr=round(float((col > 0).mean() * 100), 2),
                                 n_cells=int(len(sub_idx))))
pd.DataFrame(dot_rows).to_csv(f"{OUT}/cAE_atlas_celltype_markers.tsv", sep="\t", index=False)

# full annotated object — MANDATORY (cAK / cAJ / cAL read it): a write failure propagates, so a
# stale h5ad can never sit beside freshly written TSV exports under a success exit.
adata.write(f"{OUT}/cAE_atlas.h5ad")
log(f"[h5ad] wrote {OUT}/cAE_atlas.h5ad")

# QC figure
try:
    fig, axes = plt.subplots(1, 3, figsize=(21, 7))
    cmap = plt.get_cmap("tab20")
    cts = sorted(adata.obs[dck].dropna().unique())
    for i, ct in enumerate(cts):
        m = (adata.obs[dck] == ct).values
        axes[0].scatter(um[m, 0], um[m, 1], s=2, color=cmap(i % 20), label=ct, rasterized=True)
    axes[0].legend(markerscale=4, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
    axes[0].set_title(f"celltype res={DEFAULT_RES} (harmony={harmony_applied}, hvg={hvg_flavor})")
    for grp, col in [("A", "#1f77b4"), ("B", "#d62728")]:
        m = (adata.obs.group == grp).values
        axes[1].scatter(um[m, 0], um[m, 1], s=2, color=col, label=f"Group {grp}", rasterized=True)
    axes[1].legend(markerscale=4); axes[1].set_title("group (batch-mix check)")
    lk = f"leiden_res{DEFAULT_RES}"
    for cl in sorted(adata.obs[lk].unique(), key=int):
        m = (adata.obs[lk] == cl).values
        axes[2].scatter(um[m, 0], um[m, 1], s=2, color=cmap(int(cl) % 20), rasterized=True)
        axes[2].text(um[m, 0].mean(), um[m, 1].mean(), cl, fontsize=9, fontweight="bold")
    axes[2].set_title(f"leiden res={DEFAULT_RES}")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.suptitle("Chick whole-embryo atlas v2 (duck-d1 parity + d2b annotation; Harmony on library)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f"{OUT}/cAE_atlas_qc.png", dpi=150, bbox_inches="tight")
    log(f"[fig] wrote {OUT}/cAE_atlas_qc.png")
except Exception as e:
    raise SystemExit(f"[fig] FAILED: {type(e).__name__}: {e}\n{traceback.format_exc()}") from e

log(f"\n=== DONE. {adata.n_obs} cells. harmony={harmony_applied} umap={umap_applied} "
    f"hvg={hvg_flavor} cc={cc_scored} ===")
log("VALIDATE cAE_atlas_cluster_annotation_res*.tsv (d2b calls + paraxial probe) + "
    "cAE_atlas_DE_res*.tsv before finalising figures.")
