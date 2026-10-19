#!/usr/bin/env python3
"""cAD — chick PITX2-locus snRNA pileup, FAITHFUL duck d2i / fig02D port (STAGE x POPULATION).

Supersedes cAC (a from-scratch L/R version that did not mirror fig02D). This ports the duck
`pipeline/d2i_pitx2_rna_pileup.py` + `figures/fig02/generate_pitx2_pileup_multistage.R` template
to chick, with five additions over cAC:
  (a) marks internal polyA-priming sites (scans the chick PITX2-locus FASTA for A-rich tracts);
  (b) full gene model from the GTF: 3 isoform promoters (XM_025149516 / XM_025149515 / NM_205010,
      first exons 57,785,032 / 57,787,817 / 57,794,829 -> a / b / c-analog by promoter order) with
      UTR vs CDS vs gene-body (CDS + UTR + start/stop rows present in the galGal6a GTF);
  (c) per-cell UMI depth, emitted so the R plot can use a labeled shared/two-tier y-axis;
  (d) HH13 vs HH15-16 RESOLVED into separate stage tracks (cY maturation split within GFP class +
      design cardinality; TENTATIVE — Group B has no physical HH13/HH15 anchor);
  (e) per-SUBCLUSTER tracks at HH25 (Dorsal sub-epi. / Ventral spl. / Syringeal-Cardio), like duck.

Panel D = STAGE x POPULATION trajectory (NOT L/R; the designed L/R control is the panel-C analog
cP2/cR/cU). Sense-strand (-F 16, PITX2 + strand). Total PITX2 (isoform first-exons unresolved, cZ).
Uses ALL Group B mesenchyme cells per stage x population (coverage, not the WT-only laterality test).

Harmony (lane) — whole-B and per-stage — is a direct harmonypy.run_harmony(..., random_state=SEED)
call (cAE's port, with its Z_corr shape handling) and RAISES on failure. The previous
sc.external.pp.harmony_integrate call sat in a try/except that fell back to uncorrected X_pca
(lane axis NOT removed) with only a printed warning. SEED=0 equals harmonypy's own default
random_state, so the seed is UNCHANGED; the real differences are (i) no silent X_pca fallback and
(ii) BLAS/OpenMP pinned to one thread. If harmonypy imported fine in the published run, the output
should change only by thread-count float noise.

Outputs go to a SEPARATE directory so the published seed_qc/ copies are never overwritten:
  CAD_OUT_DIR (default {WORK}/vireo_seeded/seed_qc_cAD_harmonyfix/; refusing seed_qc/ itself).
Inputs (split donor_ids) are still read from seed_qc/; Group B stages from CAD_STAGE_TSV (chick_stage_alignment.tsv)..

Outputs (CAD_OUT_DIR):
  cAD_pitx2_pileup_depth.tsv.gz     per-base: position + <track>_raw + <track>_percell
  cAD_pitx2_gene_model.tsv          isoform/exon/CDS/UTR features (from GTF)
  cAD_pitx2_polya_priming.tsv       internal A-rich priming tract positions (from FASTA)
  cAD_pitx2_pileup_tracks.tsv       track -> stage, population, n_cells, total_depth

  sbatch pipeline/chick/cAD_groupB_pitx2_pileup_d2i.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip, subprocess, tempfile, traceback, re
# ── Pin BLAS/OpenMP threading for reproducibility (cAE / duck d1) — before numpy import ──
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _groupB_split_excluded import fold_split  # noqa: E402  (no {base}a default; fails loud)
import harmonypy

WORK = f"{SCRATCH_USER}/chick_demux_B"
# matrix subdir for the internal mesenchyme grouping clustering: "matrix" (raw cellranger) or
# "matrix_cb" (CellBender, duck-parity). Coverage itself is BAM-based (CB-invariant); this only
# affects which cells fall in each stage x population track. Set CHICK_MATRIX_SUBDIR=matrix_cb.
MATRIX_SUBDIR = os.environ.get("CHICK_MATRIX_SUBDIR", "matrix")
BDIR = f"{WORK}/vireo_seeded/N14/seed_105"
OUT = f"{WORK}/vireo_seeded/seed_qc"                       # INPUTS (splits, per-cell sex)
CAD_OUT = os.environ.get("CAD_OUT_DIR") or f"{WORK}/vireo_seeded/seed_qc_cAD_harmonyfix"
if os.path.realpath(CAD_OUT) == os.path.realpath(OUT):
    raise RuntimeError(f"CAD_OUT_DIR={CAD_OUT} is the published seed_qc/ dir — refusing to overwrite "
                       "the published cAD outputs; pick a new directory")
os.makedirs(CAD_OUT, exist_ok=True)
REF = f"{SCRATCH_USER}/chick_demux/ref"
GTF = f"{REF}/genes.gtf.gz"
FASTA = f"{REF}/genome.fa"
LIBS = [("GEX09_GEX", "9", f"{WORK}/bams/GEX09_GEX/possorted_genome_bam.bam"),
        ("GEX10_GEX", "10", f"{WORK}/bams/GEX10_GEX/possorted_genome_bam.bam")]
# Resolved-12 identity (same rule as _groupB_resolved12.py): donor11 and donor2 are split into
# a/b by their N=2 vireo runs (donor1 -> a, donor0 -> b); cells the split calls doublet or
# unassigned are DROPPED. donor8 is intentionally NOT split (weak split, r=0.51).
SPLITS = {"donor11": f"{OUT}/split_donor11/N2/donor_ids.tsv", "donor2": f"{OUT}/split_donor2/N2/donor_ids.tsv"}
SPLIT_SUBDONORS = ["donor11a", "donor11b", "donor2a", "donor2b"]
# split coverage: fold_split (pipeline/chick/_groupB_split_excluded.py) requires 100% (raises on any absent base cell)
# fraction of each lane's matrix cells that must join the base N14/seed_105 donor_ids.tsv (any
# call, incl. doublet/unassigned). Expected ~100%: per lane the donor_ids row count equals the
# filtered-matrix cell count (cb_libqc_perlane.tsv n_demux == n_cells);
# a barcode-suffix mismatch drops a whole lane to ~0% and must raise, not silently drop the lane.
MIN_DEMUX_COVERAGE = 0.95
STAGES = ("HH13", "HH15-16", "HH25")
# Group B per-embryo stages come from the published staging table (Fig S12C), not re-derived here.
STAGE_TSV = os.environ.get("CAD_STAGE_TSV") or ""
# Known demultiplexing remnants (tiny unassigned donor columns; see cAE DEMUX_REMNANTS): dropped.
GROUPB_REMNANTS = {"donor3", "donor5", "donor12", "donor13"}
MES = ["PDGFRA", "COL1A1", "COL3A1", "PRRX1", "TWIST1", "TWIST2", "POSTN", "LUM", "DCN",
       "TBX5", "TBX4", "FOXF1", "GATA4", "WT1", "HAND1", "HAND2", "FOXC2", "BMPER", "MEOX2"]
PITX2 = "PITX2"
CHROM_CANDIDATES = ["chr4", "NC_052535.1", "4"]
GENE_START, GENE_END = 57785032, 57798939   # GTF gene span; refined from GTF at runtime
FLANK = 1500
SUBCLUST_RES = 0.3
MIN_TRACK = 40        # min cells to render a track
SEED = 0
np.random.seed(SEED)
sc.settings.n_jobs = 1

# design cardinality for the HH13/HH15-16 split within each early GFP class (cY / design):
#   early WT  = 1 HH13 + 3 HH15-16 ; early GFP+ = 3 HH13 + 1 HH15-16  (numbering resets per stage)

# cW Yan/duck signed lineage panels (late/HH25 D/V + early pre-TES)
def _g(*xs):
    out = []
    for x in xs:
        out.append(x)
        if x == "NKX6-1": out.append("NKX6_1")
    return out
PANELS = {
    "Dorsal sub-epi.": (_g("NKX6-1", "WNT2B", "PTCH1"), _g("GATA4", "HAND2", "TBX4")),
    "Ventral spl.":    (_g("TBX5", "TBX4", "ZFPM2", "ISL1"), _g("NKX6-1", "BARX1", "WT1")),
    "Cardio/Syringeal":(_g("HAND2", "ISL1", "IRX1", "HAND1"), _g("NKX6-1", "WT1")),
    "Chondrogenic":    (_g("RUNX2", "WNT5A", "DPP10", "SOX9"), _g("GATA4")),
    "Early LPM":       (_g("TBX5", "FOXF1", "ISL1", "HAND1"), _g("PAX1", "PAX3", "WT1")),
    "BARX1+":          (_g("BARX1", "LRIG1", "SMOC1", "NR2F2"), _g("WT1")),
    "Mesothelium":     (_g("WT1", "TCF21", "UPK3B"), _g("FOXF1", "GATA4", "TBX5")),
    "Somatic":         (_g("PRRX1", "TWIST1", "TWIST2", "MSX1"), _g("GATA4", "FOXF1", "NKX6-1")),
    "Sclerotome":      (_g("PAX1", "PAX9", "FOXC2", "NKX3-2"), _g("GATA4", "FOXF1")),
    "Dermomyotome":    (_g("PAX3", "PAX7", "MYF5", "LBX1"), _g("GATA4", "FOXF1")),
}
GATE = {"early": {"Early LPM", "BARX1+", "Dorsal sub-epi.", "Cardio/Syringeal", "Mesothelium",
                  "Somatic", "Sclerotome", "Dermomyotome"},
        "late":  {"Ventral spl.", "Dorsal sub-epi.", "Cardio/Syringeal", "Chondrogenic", "BARX1+",
                  "Mesothelium", "Somatic", "Sclerotome", "Dermomyotome"}}
# which populations become rendered pileup tracks, per stage (mirror duck: pre-TES pooled-ish; HH25 D/V)
TRACK_POPS = {"HH13":    ["All LPM"],
              "HH15-16": ["Dorsal-biased", "Ventral (pooled)"],
              "HH25":    ["Dorsal sub-epi.", "Ventral spl.", "Syringeal-Cardio"]}


def log(m): print(m, flush=True)


def harmony_lane(a, what):
    """Direct seeded harmonypy on 'lane' (cAE port). Writes a.obsm['X_pca_harmony']; RAISES on any
    failure — never a silent fallback to uncorrected X_pca."""
    ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, ["lane"], random_state=SEED)
    Z = np.asarray(ho.Z_corr).T
    if Z.shape[0] != a.n_obs:          # harmonypy versions differ in Z_corr orientation
        Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != a.n_obs:
        raise RuntimeError(f"[harmony] {what}: Z {Z.shape} != n_obs {a.n_obs}")
    a.obsm["X_pca_harmony"] = Z
    log(f"[harmony] Harmony applied ({what}; Z shape {Z.shape}, random_state={SEED})")
    return "X_pca_harmony"


def safe(t):
    """Shell/column-safe key: alphanumeric + underscore only (track names have spaces/()/.)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", t).strip("_")


def sh(cmd, timeout=600):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)


# ── GTF gene model (3 isoforms, exon/CDS/UTR) ─────────────────────────────────
def parse_gene_model():
    rows = []
    gene = None
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or "PITX2" not in p[8]:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', p[8]))
            tid = attrs.get("transcript_id", "")
            rows.append(dict(chrom=p[0], feature=p[2], start=int(p[3]), end=int(p[4]),
                             strand=p[6], transcript_id=tid))
            if p[2] == "gene":
                gene = (p[0], int(p[3]), int(p[4]), p[6])
    return gene, pd.DataFrame(rows)


# ── polyA internal-priming scan on the sense strand of the locus ──────────────
def polya_priming(chrom, start, end, win=10, min_a=8, run=6):
    """A-rich tracts on the + (sense) strand within [start,end]: report tract midpoints.
    Criteria (internal-priming heuristic): a run of >=`run` consecutive A, OR a `win`-bp window
    with >=`min_a` A. Adjacent hits within 20 bp merged; midpoint reported."""
    r = sh(f"samtools faidx {FASTA} {chrom}:{start}-{end}")
    seq = "".join(r.stdout.split("\n")[1:]).upper()
    if r.returncode != 0 or not seq:
        raise RuntimeError(f"samtools faidx {chrom}:{start}-{end} failed (rc={r.returncode}): {r.stderr[:300]}")
    hits = set()
    # consecutive A runs
    for m in re.finditer(r"A{%d,}" % run, seq):
        for i in range(m.start(), m.end()):
            hits.add(i)
    # windowed A-richness
    a = np.array([1 if c == "A" else 0 for c in seq])
    csum = np.concatenate([[0], np.cumsum(a)])
    for i in range(0, len(seq) - win + 1):
        if csum[i + win] - csum[i] >= min_a:
            for j in range(i, i + win):
                hits.add(j)
    if not hits:
        return []
    sorted_hits = sorted(hits)
    tracts = []
    s = prev = sorted_hits[0]
    for h in sorted_hits[1:]:
        if h - prev <= 20:
            prev = h
        else:
            tracts.append((start + (s + prev) // 2)); s = prev = h
    tracts.append(start + (s + prev) // 2)
    return tracts


def detect_chrom_suffix(bam):
    chrom = None
    for c in CHROM_CANDIDATES:
        if sh(f"samtools view {bam} {c}:{GENE_START}-{GENE_END} 2>/dev/null | head -50").stdout.strip():
            chrom = c; break
    if chrom is None:
        raise RuntimeError(f"no reads at PITX2 locus in {bam}")
    sfx = ""
    for line in sh(f"samtools view {bam} {chrom}:{GENE_START}-{GENE_END} | head -3000").stdout.split("\n"):
        for tok in line.split("\t"):
            if tok.startswith("CB:Z:"):
                cb = tok[5:]; sfx = cb[16:] if len(cb) > 16 else ""; break
        if "CB:Z:" in line:
            break
    return chrom, sfx


def depth_for(bam, bcfile, region, flag="-F 16"):
    """flag='-F 16' = sense (PITX2 + strand, forward reads); '-f 16' = antisense control."""
    tmp = tempfile.mktemp(suffix=".bam", dir="/tmp", prefix="cad_")
    try:
        r1 = sh(f"samtools view -b {flag} -D CB:'{bcfile}' '{bam}' {region} > '{tmp}'")
        if r1.returncode != 0:
            raise RuntimeError(f"samtools view failed (rc={r1.returncode}): {r1.stderr[:300]}")
        ri = sh(f"samtools index '{tmp}'")
        if ri.returncode != 0:
            raise RuntimeError(f"samtools index failed (rc={ri.returncode}): {ri.stderr[:300]}")
        rd = sh(f"samtools depth -a -r {region} '{tmp}'")
        if rd.returncode != 0:
            raise RuntimeError(f"samtools depth failed (rc={rd.returncode}): {rd.stderr[:300]}")
        out = rd.stdout
        d = {}
        for line in out.strip().split("\n"):
            if line:
                pp = line.split("\t")
                if len(pp) >= 3:
                    d[int(pp[1])] = int(pp[2])
        return d
    finally:
        for fn in (tmp, tmp + ".bai"):
            try: os.unlink(fn)
            except OSError: pass


# ── load Group B GEX, mesenchyme, stage axis, GFP (cP2/cW logic) ──────────────
parts = []
for lib, suf, _ in LIBS:
    md = f"{WORK}/{MATRIX_SUBDIR}/{lib}"
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
    sym = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    typ = np.array([f[2] if len(f) > 2 else "Gene Expression" for f in feats])
    gex = typ == "Gene Expression"
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    egi = [i for i in range(len(sym)) if sym[i].strip().upper() == "EGFP"]
    erow = np.asarray(M[egi[0]].todense()).ravel() if egi else np.zeros(M.shape[1])
    a = ad.AnnData(sp.csr_matrix(M[gex].T), obs=pd.DataFrame(index=[f"{bc}-{suf}" for bc in bcs]),
                   var=pd.DataFrame(index=sym[gex]))
    a.var_names_make_unique(); a.obs["raw16"] = [b.split("-")[0] for b in bcs]
    a.obs["egfp"] = erow; a.obs["lane"] = "GEX09" if suf == "9" else "GEX10"; parts.append(a)
adata = ad.concat(parts, join="outer"); adata.obs_names_make_unique()
dmap = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t").set_index("cell").donor_id.to_dict()
adata.obs["donor"] = adata.obs_names.map(dmap)
# base-demux join coverage (data integrity), validated per lane BEFORE split coverage: split
# coverage only sees cells that already joined a base donor, so a lane-wide join failure would
# otherwise pass it at 100% while that lane's cells fall out as donor NaN below.
for _lane, _don in adata.obs.groupby("lane", observed=True)["donor"]:
    _cov = float(_don.notna().mean())
    if _cov < MIN_DEMUX_COVERAGE:
        raise ValueError(
            f"{_lane}: only {int(_don.notna().sum())}/{len(_don)} ({_cov:.1%}) matrix cells join "
            f"{BDIR}/donor_ids.tsv (< {MIN_DEMUX_COVERAGE:.0%}; barcode-suffix mismatch?); "
            f"e.g. atlas {list(_don.index[:3])} vs demux {list(dmap)[:3]}")
    log(f"[demux join] {_lane}: {int(_don.notna().sum())}/{len(_don)} cells joined base donor_ids")
# split fold via the shared validated helper: donor0 -> {base}b, donor1 -> {base}a;
# doublet/unassigned keep that status and are dropped below; a missing split file, a base cell
# absent from it (100% split coverage), a drifted clean assignment or exclusion set raises.
for base, f in SPLITS.items():
    dmap = fold_split(dmap, base, f, "B", log=log)
adata.obs["donor"] = adata.obs_names.map(dmap)
_empty = [d for d in SPLIT_SUBDONORS if not (adata.obs["donor"] == d).any()]
if _empty:
    raise ValueError(f"resolved-12 split sub-donors with zero cells: {_empty}")
n_rem = int(adata.obs.donor.isin(GROUPB_REMNANTS).sum())
adata = adata[adata.obs.donor.notna() & ~adata.obs.donor.isin(["doublet", "unassigned", "split_excluded"])
              & ~adata.obs.donor.isin(GROUPB_REMNANTS)].copy()
log(f"[resolved12] dropped {n_rem} demux-remnant cells; donors kept: {sorted(adata.obs.donor.unique())}")
if "EGFP" in adata.var_names: adata = adata[:, adata.var_names != "EGFP"].copy()
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata); adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
ad_h = adata[:, adata.var.highly_variable].copy(); sc.pp.scale(ad_h, max_value=10); sc.tl.pca(ad_h, n_comps=30)
rep = harmony_lane(ad_h, "whole-B")
sc.pp.neighbors(ad_h, use_rep=rep, n_neighbors=15); sc.tl.leiden(ad_h, resolution=0.5)
adata.obs["leiden"] = ad_h.obs["leiden"]
sc.tl.score_genes(adata, [g for g in MES if g in adata.raw.var_names], score_name="mes", use_raw=True)
clm = adata.obs.groupby("leiden")["mes"].mean().sort_values(ascending=False)
mes_cl = clm[clm > 0.05].index.tolist()
if not mes_cl:
    raise RuntimeError(f"no leiden cluster passes the mesenchyme-score floor 0.05: {clm.round(3).to_dict()}")
mes = adata[adata.obs.leiden.isin(mes_cl)].copy()
wz = set()
for p in [f"{WORK}/chick_W_genes.tsv", f"{WORK}/chick_Z_genes.tsv"]:
    if not os.path.exists(p):   # required: W/Z genes are excluded from the HVG set below
        raise FileNotFoundError(f"W/Z gene list missing: {p}")
    wz |= set(pd.read_csv(p, header=None, sep="\t")[0].astype(str))
# stage per donor: the published Group B staging (Fig S12C, chick_stage_alignment.tsv), required.
if not STAGE_TSV or not os.path.exists(STAGE_TSV):
    raise FileNotFoundError(f"CAD_STAGE_TSV must point to chick_stage_alignment.tsv (got {STAGE_TSV!r})")
_st = pd.read_csv(STAGE_TSV, sep="\t")
_need = {"group", "donor", "stage"} - set(_st.columns)
if _need:
    raise KeyError(f"{STAGE_TSV} lacks columns {sorted(_need)}")
_st = _st[_st.group == "B"]
if _st.donor.duplicated().any():
    raise ValueError(f"duplicate Group B donors in {STAGE_TSV}: {sorted(_st.donor[_st.donor.duplicated()])}")
_badst = _st[~_st.stage.isin(STAGES)]
if len(_badst):
    raise ValueError(f"Group B stages outside {STAGES} in {STAGE_TSV}: {_badst[['donor', 'stage']].to_dict('records')}")
_st = _st.set_index("donor").stage.to_dict()
_missing = sorted(set(mes.obs.donor) - set(_st))
if _missing:
    raise KeyError(f"Group B donors with mesenchyme cells but no published stage: {_missing}")
donor_stage = {d: _st[d] for d in sorted(set(mes.obs.donor))}
mes.obs["stage"] = mes.obs.donor.map(donor_stage)
if not mes.obs.stage.isin(STAGES).all():
    raise ValueError(f"mesenchyme cells without a valid stage: {mes.obs.stage.value_counts(dropna=False).to_dict()}")
mes.obs["PITX2_raw"] = np.asarray(mes.raw[:, PITX2].X.todense()).ravel()
log(f"[mesenchyme] {mes.n_obs} cells; donor_stage={donor_stage}")


# ── per-stage subcluster -> population (cW annotate) ──────────────────────────
def annotate(sub, stage_grp):
    admissible = [p for p in PANELS if p in GATE[stage_grp]]
    S = {}
    for pop in admissible:
        up, dn = PANELS[pop]
        u = [g for g in up if g in sub.raw.var_names]; d = [g for g in dn if g in sub.raw.var_names]
        if not u: continue
        sc.tl.score_genes(sub, u, score_name="_u", use_raw=True)
        sc.tl.score_genes(sub, d, score_name="_d", use_raw=True) if d else sub.obs.__setitem__("_d", 0.0)
        S[pop] = sub.obs.groupby("leiden", observed=True).apply(lambda x: (x["_u"] - x["_d"]).mean())
    S = pd.DataFrame(S)
    def pct(g):
        gg = [x for x in ([g, g.replace("-", "_")] if "-" in g else [g]) if x in sub.raw.var_names]
        if not gg: return pd.Series(0.0, index=S.index)
        v = np.asarray(sub.raw[:, gg[0]].X.todense()).ravel() > 0
        return pd.Series(v, index=sub.obs_names).groupby(sub.obs.leiden, observed=True).mean() * 100
    wt1, barx1 = pct("WT1"), pct("BARX1")
    call = {}
    for cl in S.index:
        c = S.loc[cl].idxmax()
        if wt1.get(cl, 0) > 35: c = "Mesothelium"
        elif barx1.get(cl, 0) > 40 and c != "Mesothelium": c = "BARX1+"
        call[cl] = c
    return call

mes.obs["population"] = pd.Series(index=mes.obs_names, dtype=object)
for stg in STAGES:
    sub = mes[mes.obs.stage == stg].copy()
    if sub.n_obs < 60:
        log(f"[{stg}] {sub.n_obs} cells -> all 'All LPM'"); mes.obs.loc[sub.obs_names, "population"] = "All LPM"; continue
    gate = "late" if stg == "HH25" else "early"
    sh2 = sub.copy(); sc.pp.highly_variable_genes(sh2, n_top_genes=2000)
    sh2.var.loc[sh2.var_names.isin(wz), "highly_variable"] = False
    sx = sh2[:, sh2.var.highly_variable].copy(); sc.pp.scale(sx, max_value=10); sc.tl.pca(sx, n_comps=30)
    srep = harmony_lane(sx, f"stage {stg}")
    sc.pp.neighbors(sx, use_rep=srep, n_neighbors=15); sc.tl.leiden(sx, resolution=SUBCLUST_RES)
    sub.obs["leiden"] = sx.obs["leiden"].values
    call = annotate(sub, gate)
    mes.obs.loc[sub.obs_names, "population"] = sub.obs.leiden.map(call).values
    log(f"[{stg}] " + ", ".join(f"{k}={v}" for k, v in pd.Series(sub.obs.leiden.map(call).values).value_counts().items()))

# map fine populations -> rendered track names per stage (mirror duck panel D)
def track_of(stage, pop):
    if pop in ("Sclerotome", "Dermomyotome"):
        return None   # paraxial mesoderm — excluded from all LPM tracks (user 2026-06-08)
    if stage == "HH13":
        return "HH13|LPM mesenchyme"
    if stage == "HH15-16":
        return f"HH15-16|{'Dorsal-biased' if pop == 'Dorsal sub-epi.' else 'Ventral (pooled)'}"
    if stage == "HH25":
        if pop == "Dorsal sub-epi.": return "HH25|Dorsal sub-epi."
        if pop == "Ventral spl.":    return "HH25|Ventral spl."
        if pop in ("Cardio/Syringeal", "Chondrogenic", "BARX1+"): return "HH25|Syringeal-Cardio"
        return None   # excluded (Mesothelium/Somatic) not rendered
    return None
mes.obs["track"] = [track_of(s, p) for s, p in zip(mes.obs.stage, mes.obs.population)]


# ── per-track barcodes (BAM CB form) + samtools sense pileup ──────────────────
bam_meta = {lib: (detect_chrom_suffix(bam)) for lib, suf, bam in LIBS}
for lib, (c, s) in bam_meta.items(): log(f"[{lib}] chrom={c} CB-suffix='{s}'")
BCDIR = f"{CAD_OUT}/cAD_bc"; os.makedirs(BCDIR, exist_ok=True)
gene, gm = parse_gene_model()
chrom0 = bam_meta[LIBS[0][0]][0]
g_start, g_end = (gene[1], gene[2]) if gene else (GENE_START, GENE_END)
region_pos = list(range(g_start - FLANK, g_end + FLANK + 1))
ncells = defaultdict(int); track_depth = {}
tracks_present = [t for t in mes.obs.track.dropna().unique()]
for track in tracks_present:
    cells = mes.obs[mes.obs.track == track]
    if len(cells) < MIN_TRACK:
        log(f"[skip] {track} n={len(cells)} < {MIN_TRACK}"); continue
    tot = defaultdict(int)
    for lib, suf, bam in LIBS:
        lane = "GEX09" if suf == "9" else "GEX10"
        lc = cells[cells.lane == lane]
        if not len(lc): continue
        c, sfx = bam_meta[lib]
        bcfile = f"{BCDIR}/{safe(track)}_{lib}.txt"
        with open(bcfile, "w") as fh: fh.write("\n".join(f"{b}{sfx}" for b in lc["raw16"]) + "\n")
        d = depth_for(bam, bcfile, f"{c}:{g_start - FLANK}-{g_end + FLANK}")
        for pos, v in d.items(): tot[pos] += v
    track_depth[track] = tot; ncells[track] = len(cells)
    log(f"[{track}] n={len(cells)} total_depth={sum(tot.values()):,}")

# ── antisense (-f 16) CONTROL: aggregate over ALL rendered-track cells ─────────
# Definitively rules out antisense contamination at the sense 5' peak: PITX2 is + strand,
# so genuine a/b-isoform transcription reads are sense (-F 16). If the 5' peak were an
# antisense ncRNA artifact it would appear here (-f 16). Expect ~0 at the sense 5' peak.
anti_tot = defaultdict(int); n_anti = 0
track_cells = mes.obs[mes.obs.track.notna()]
for lib, suf, bam in LIBS:
    lane = "GEX09" if suf == "9" else "GEX10"
    lc = track_cells[track_cells.lane == lane]
    if not len(lc):
        continue
    c, sfx = bam_meta[lib]
    bcfile = f"{BCDIR}/ANTISENSE_{lib}.txt"
    with open(bcfile, "w") as fh:
        fh.write("\n".join(f"{b}{sfx}" for b in lc["raw16"]) + "\n")
    d = depth_for(bam, bcfile, f"{c}:{g_start - FLANK}-{g_end + FLANK}", flag="-f 16")
    for pos, v in d.items():
        anti_tot[pos] += v
    n_anti += len(lc)
log(f"[antisense control] n={n_anti} cells, total_depth={sum(anti_tot.values()):,} "
    f"(-f 16; expect ~0 at the sense 5' peak)")

# ── outputs ───────────────────────────────────────────────────────────────────
order = [t for t in ["HH13|LPM mesenchyme", "HH15-16|Dorsal-biased", "HH15-16|Ventral (pooled)",
                     "HH25|Dorsal sub-epi.", "HH25|Ventral spl.", "HH25|Syringeal-Cardio"]
         if t in track_depth]
rows = []
for pos in region_pos:
    r = {"chrom": chrom0, "position": pos}
    for t in order:
        key = safe(t)
        r[f"{key}_raw"] = track_depth[t].get(pos, 0)
        r[f"{key}_percell"] = track_depth[t].get(pos, 0) / max(1, ncells[t])
    r["antisense_raw"] = anti_tot.get(pos, 0)
    r["antisense_percell"] = anti_tot.get(pos, 0) / max(1, n_anti)
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{CAD_OUT}/cAD_pitx2_pileup_depth.tsv.gz", sep="\t", index=False)
trk_rows = [{"track": t, "safe_key": safe(t), "stage": t.split("|")[0], "population": t.split("|")[1],
             "n_cells": ncells[t], "total_depth": sum(track_depth[t].values())} for t in order]
trk_rows.append({"track": "ANTISENSE control", "safe_key": "antisense", "stage": "control",
                 "population": "antisense (-f 16)", "n_cells": n_anti,
                 "total_depth": sum(anti_tot.values())})
pd.DataFrame(trk_rows).to_csv(f"{CAD_OUT}/cAD_pitx2_pileup_tracks.tsv", sep="\t", index=False)
gm.to_csv(f"{CAD_OUT}/cAD_pitx2_gene_model.tsv", sep="\t", index=False)
priming = polya_priming(chrom0, g_start, g_end)
pd.DataFrame({"position": priming}).to_csv(f"{CAD_OUT}/cAD_pitx2_polya_priming.tsv", sep="\t", index=False)
log(f"\ngene model: {len(gm)} features ({gm.transcript_id.nunique()} transcripts); "
    f"polyA priming tracts: {len(priming)} -> {priming}")
log(f"tracks rendered: {order}")
log(f"wrote {CAD_OUT}/cAD_pitx2_pileup_depth.tsv.gz + _tracks.tsv + _gene_model.tsv + _polya_priming.tsv")
