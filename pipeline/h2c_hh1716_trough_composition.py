#!/usr/bin/env python3
"""h2c: left/right cell-composition decomposition of the HH17/16 inference pool.

At HH17/16 (the PITX2 low between the two waves) the snRNA inference pool is two
subclusters (Cardiopharyngeal LPM + BARX1+; read from ``POOLS`` via
``get_pool_map``, never hardcoded).  AHR and CHRM2 are right-biased in the pooled
L/R contrast, and both are expressed far more in Cardiopharyngeal LPM than in
BARX1+ (per-embryo median 7.1x and 13.2x on the 2026-10-14 build).
This step asks how much of each pooled right-minus-left difference is the
Cardiopharyngeal/BARX1+ balance shifting between the two halves (composition)
and how much is a real per-subcluster expression change (within-subcluster).

Cohort = the h2 POOL paired cohort, reproduced exactly:
  * pseudobulk unit = (genotype, library_m, side) summed over the inference-pool
    subclusters of stage_anchor == "HH17/16";
  * floor = max(20, round(0.05 x median pooled n_cells)) (h2 "principled" mode);
  * keep genotypes (= embryos) with both an L and an R pseudobulk at or above the
    floor (at or above 20 pool nuclei on the current build; the test is n >= floor).
Counts are clipped at 0 and rounded, as h2 does before DESeq2.

Decomposition (per paired embryo, per gene; k = subcluster, s = side):
  w_k,s = UMI share of subcluster k in side s's pooled pseudobulk
  c_k,s = CPM of the gene in subcluster k, side s
  P_s   = sum_k w_k,s * c_k,s          (= pooled CPM; asserted against sum G / sum U)
  composition = sum_k (w_k,R - w_k,L) * (c_k,R + c_k,L) / 2
  within      = sum_k (w_k,R + w_k,L) / 2 * (c_k,R - c_k,L)
  composition + within == P_R - P_L    (exact; asserted on every embryo x gene)
A subcluster absent on one side has weight 0 there and borrows the other side's
rate (exact: a weight-0 term contributes nothing to that side's P).  A subcluster
absent on both sides contributes exactly 0 to both parts (rates left NaN).

Summary statistics: paired tests are exact two-sided Wilcoxon signed-rank tests by
exhaustive enumeration (figures/_shared/signed_rank.py: tie-correct mid-ranks, zero
differences dropped and counted).  A statistic that is genuinely undefined on valid
data (no non-zero difference; a 0/0 enrichment fold) is written as NaN/inf next to
its count, never an abort; panel F refuses to draw an undefined value.

Inputs (all explicit; no fallback chain; --h1-dir has NO default and must be given;
h1 metadata sample_ids and count columns must match one-to-one):
  <h1-dir>/h1_pseudobulk_snrna_counts.tsv.gz    genes x pseudobulks (h1)
  <h1-dir>/h1_pseudobulk_snrna_metadata.tsv     one row per pseudobulk (h1)
  --e1b  figures/data/e1b_stage_assignments.tsv  per-genotype sex (context column only)
  figures/_shared/subcluster_map.py POOLS / get_pool_map   inference-pool membership

Arguments: --h1-dir is required (no default).  --e1b defaults to
figures/data/e1b_stage_assignments.tsv and --out-dir to figures/data; both are
printed back by the build so the run log records exactly what was read/written.

Outputs (<out-dir>, default figures/data):
  h2c_hh1716_trough_composition.tsv          one row per (embryo, gene, subcluster)
  h2c_hh1716_trough_composition_summary.tsv  headline statistics (long: metric/gene/subcluster/value)

Consumer: figures/fig02_supp/generate_trough_composition.py (figS02 / Fig S9 panel F).

Run (local; conda env duckfig):
  conda run -n duckfig python3 pipeline/h2c_hh1716_trough_composition.py --h1-dir figures/data
  conda run -n duckfig python3 pipeline/h2c_hh1716_trough_composition.py --self-check
(--h1-dir is required for a build: a bare run would silently pick whichever h1
sits in figures/data, which can be older than the h1 run being summarised.)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _paths import REPO_ROOT  # portable roots: docs/audit_hardcoded_paths.md

sys.path.insert(0, str(REPO_ROOT / "figures"))
from _shared.signed_rank import self_test as signed_rank_self_test  # noqa: E402
from _shared.signed_rank import signed_rank_exact  # noqa: E402
from _shared.subcluster_map import POOLS, get_pool_map  # noqa: E402

STAGE = "HH17/16"
GENES = ("AHR", "CHRM2", "PITX2")      # PITX2 = reference (left-biased) gene
# The panel contrasts the two HH17/16 inference-pool subclusters by name; build()
# asserts the pool (from POOLS / get_pool_map) is exactly these two.
FOCAL = "Cardiopharyngeal LPM"
OTHER = "BARX1+"
SIDES = ("L", "R")
H1_COUNTS = "h1_pseudobulk_snrna_counts.tsv.gz"
H1_META = "h1_pseudobulk_snrna_metadata.tsv"
OUT_MAIN = "h2c_hh1716_trough_composition.tsv"
OUT_SUMMARY = "h2c_hh1716_trough_composition_summary.tsv"
META_COLS = ("sample_id", "genotype", "library_m", "side", "subcluster",
             "stage_anchor", "n_cells")
RTOL, ATOL = 1e-9, 1e-6


# ── Decomposition (pure function; exercised by --self-check) ──────────────
def decompose(U: dict, G: dict) -> dict:
    """Exact two-part decomposition of a pooled CPM difference.

    U[s], G[s]: float arrays over subclusters (total UMIs, gene UMIs) for side s.
    Returns per-subcluster weights, rates (after weight-0 borrowing), the
    per-subcluster composition/within terms, and the pooled P_L/P_R.
    """
    w, c, borrowed = {}, {}, {}
    for s in SIDES:
        tot = U[s].sum()
        if not tot > 0:
            raise ValueError(f"side {s}: pooled pseudobulk has zero UMIs")
        w[s] = U[s] / tot
        with np.errstate(divide="ignore", invalid="ignore"):
            c[s] = np.where(U[s] > 0, G[s] / U[s] * 1e6, np.nan)
    for a, b in (("L", "R"), ("R", "L")):
        borrowed[a] = np.isnan(c[a]) & ~np.isnan(c[b])
    for a, b in (("L", "R"), ("R", "L")):
        c[a] = np.where(borrowed[a], c[b], c[a])
    both_absent = np.isnan(c["L"]) & np.isnan(c["R"])
    cL = np.where(both_absent, 0.0, c["L"])
    cR = np.where(both_absent, 0.0, c["R"])
    P = {s: G[s].sum() / U[s].sum() * 1e6 for s in SIDES}
    for s, cs in (("L", cL), ("R", cR)):
        if not np.isclose((w[s] * cs).sum(), P[s], rtol=RTOL, atol=ATOL):
            raise AssertionError(f"side {s}: sum_k w*c = {(w[s]*cs).sum()} != pooled CPM {P[s]}")
    comp_k = (w["R"] - w["L"]) * (cR + cL) / 2
    within_k = (w["R"] + w["L"]) / 2 * (cR - cL)
    dP = P["R"] - P["L"]
    if not np.isclose(comp_k.sum() + within_k.sum(), dP, rtol=RTOL, atol=ATOL):
        raise AssertionError(f"composition {comp_k.sum()} + within {within_k.sum()} != P_R - P_L {dP}")
    return dict(w=w, c=c, borrowed=borrowed, comp_k=comp_k, within_k=within_k,
                P_L=P["L"], P_R=P["R"], delta=dP,
                composition=float(comp_k.sum()), within=float(within_k.sum()))


def self_check() -> None:
    """Synthetic cases with known answers; raises on any failure."""
    rng = np.random.default_rng(0)
    # 1. random weights/rates: exact sum
    for _ in range(200):
        K = int(rng.integers(2, 5))
        U = {s: rng.integers(1_000, 1_000_000, K).astype(float) for s in SIDES}
        G = {s: np.floor(U[s] * rng.uniform(0, 1e-3, K)) for s in SIDES}
        decompose(U, G)                       # asserts internally
    # 2. identical rates on both sides -> within == 0, all of dP is composition
    U = {"L": np.array([7e5, 3e5]), "R": np.array([3e5, 7e5])}
    G = {"L": np.array([7.0, 300.0]), "R": np.array([3.0, 700.0])}   # rates 10 / 1000 CPM
    r = decompose(U, G)
    assert abs(r["within"]) < ATOL and np.isclose(r["composition"], r["delta"]), r
    # 3. identical weights -> composition == 0, all of dP is within
    U = {"L": np.array([5e5, 5e5]), "R": np.array([5e5, 5e5])}
    G = {"L": np.array([5.0, 50.0]), "R": np.array([10.0, 200.0])}
    r = decompose(U, G)
    assert abs(r["composition"]) < ATOL and np.isclose(r["within"], r["delta"]), r
    # 4. subcluster absent on one side: weight 0, rate borrowed, still exact
    U = {"L": np.array([1e6, 0.0]), "R": np.array([4e5, 6e5])}
    G = {"L": np.array([20.0, 0.0]), "R": np.array([8.0, 300.0])}
    r = decompose(U, G)
    assert r["borrowed"]["L"].tolist() == [False, True] and r["c"]["L"][1] == r["c"]["R"][1], r
    # 5. absent on both sides: contributes exactly 0
    U = {"L": np.array([1e6, 0.0]), "R": np.array([2e6, 0.0])}
    G = {"L": np.array([20.0, 0.0]), "R": np.array([80.0, 0.0])}
    r = decompose(U, G)
    assert r["comp_k"][1] == 0 and r["within_k"][1] == 0, r
    print("h2c self-check: 5/5 synthetic decomposition cases passed "
          "(exact sum, pure-composition, pure-within, one-side-absent borrow, both-absent)")
    signed_rank_self_test()


# ── Data loading / cohort (mirrors h2 POOL paired cohort) ─────────────────
def _require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"missing input: {path}")
    return path


def load_inputs(h1_dir: Path):
    counts_path = _require_file(h1_dir / H1_COUNTS)
    meta_path = _require_file(h1_dir / H1_META)
    meta = pd.read_csv(meta_path, sep="\t")
    missing_cols = [c for c in META_COLS if c not in meta.columns]
    if missing_cols:
        raise KeyError(f"{meta_path}: missing columns {missing_cols}")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    if counts.index.duplicated().any():
        raise ValueError(f"{counts_path}: duplicated gene ids "
                         f"{counts.index[counts.index.duplicated()][:5].tolist()}")
    missing_genes = [g for g in GENES if g not in counts.index]
    if missing_genes:
        raise KeyError(f"{counts_path}: genes absent {missing_genes}")
    # Metadata rows and count columns must match one-to-one BEFORE any pseudobulk is
    # read: build() treats a (side, subcluster) with no metadata row as structurally
    # absent (0 nuclei, rate borrowed), so a dropped metadata row would otherwise turn
    # a real pseudobulk into a silent zero.
    dup_meta = sorted(meta.loc[meta["sample_id"].duplicated(), "sample_id"].unique())
    if dup_meta:
        raise ValueError(f"{meta_path}: duplicated sample_id {dup_meta}")
    dup_cols = sorted(counts.columns[counts.columns.duplicated()].unique())
    if dup_cols:
        raise ValueError(f"{counts_path}: duplicated count column {dup_cols}")
    meta_only = sorted(set(meta["sample_id"]) - set(counts.columns))
    counts_only = sorted(set(counts.columns) - set(meta["sample_id"]))
    if meta_only or counts_only:
        raise KeyError(f"h1 metadata sample_id and count columns are not one-to-one: "
                       f"{len(meta_only)} metadata-only {meta_only}; "
                       f"{len(counts_only)} count-column-only {counts_only}")
    return counts, meta


def paired_cohort(meta: pd.DataFrame, pool: list[str]):
    """Return (stage metadata restricted to paired embryos, floor, paired genotypes)."""
    m = meta[meta["stage_anchor"] == STAGE].copy()
    if m.empty:
        raise ValueError(f"no h1 pseudobulks with stage_anchor == {STAGE!r}")
    present = set(m["subcluster"])
    if present != set(pool):
        raise ValueError(f"{STAGE} h1 subclusters {sorted(present)} != inference pool "
                         f"{sorted(pool)} (POOLS drift or a non-pool subcluster in h1)")
    bad_side = sorted(set(m["side"]) - set(SIDES))
    if bad_side:
        raise ValueError(f"unexpected side values {bad_side}")
    # one library per (genotype, side): same-side duplicates must already be pooled by h1
    nlib = m.groupby(["genotype", "side"])["library_m"].nunique()
    if (nlib > 1).any():
        raise ValueError(f"(genotype, side) spread over >1 library_m: {nlib[nlib > 1].to_dict()}")
    if m.duplicated(["genotype", "side", "subcluster"]).any():
        raise ValueError("duplicate (genotype, side, subcluster) pseudobulks at HH17/16")
    pooled = m.groupby(["genotype", "library_m", "side"])["n_cells"].sum()
    floor = max(20, int(round(float(pooled.median()) * 0.05)))
    keep = pooled[pooled >= floor].reset_index()
    sides = keep.groupby("genotype")["side"].agg(set)
    paired = sorted(g for g, s in sides.items() if s == set(SIDES))
    if len(paired) < 2:
        raise ValueError(f"only {len(paired)} paired embryos at or above the floor ({floor} pool nuclei)")
    m = m[m["genotype"].isin(paired)].copy()
    return m, floor, paired, int(pooled.median())


def load_sex(e1b_path: Path, genotypes: list[str]) -> dict:
    e1b = pd.read_csv(_require_file(e1b_path), sep="\t")
    for col in ("genotype_id", "sex"):
        if col not in e1b.columns:
            raise KeyError(f"{e1b_path}: missing column {col!r}")
    if e1b["genotype_id"].duplicated().any():
        raise ValueError(f"{e1b_path}: duplicated genotype_id")
    sex = dict(zip(e1b["genotype_id"].astype(str), e1b["sex"]))
    missing = [g for g in genotypes if g not in sex]
    if missing:
        raise KeyError(f"{e1b_path}: no sex for paired genotypes {missing}")
    bad = {g: sex[g] for g in genotypes if sex[g] not in ("M", "F")}
    if bad:
        raise ValueError(f"{e1b_path}: sex outside {{M, F}} for paired genotypes {bad}")
    return {g: sex[g] for g in genotypes}


# ── Main ──────────────────────────────────────────────────────────────────
def build(h1_dir: Path, e1b_path: Path):
    pool = list(POOLS[STAGE]["inference"])
    pool_from_tsv = set(get_pool_map(STAGE, "inference").values())
    if pool_from_tsv != set(pool):
        raise ValueError(f"get_pool_map({STAGE}) labels {sorted(pool_from_tsv)} != POOLS {pool}")
    if set(pool) != {FOCAL, OTHER}:
        raise ValueError(f"{STAGE} inference pool {pool} != {{{FOCAL!r}, {OTHER!r}}}: "
                         f"the two-subcluster contrast no longer applies")
    counts, meta = load_inputs(h1_dir)
    m, floor, paired, med = paired_cohort(meta, pool)
    sex = load_sex(e1b_path, paired)
    print(f"[h2c] {STAGE} pooled-pseudobulk floor = {floor} nuclei (median pooled n = {med}); "
          f"paired embryos = {len(paired)}")

    ids = m["sample_id"].tolist()
    sub = counts[ids].clip(lower=0).round()
    umi_total = sub.sum(axis=0)
    gene_umi = sub.loc[list(GENES)]

    rows = []
    for gt in paired:
        g = m[m["genotype"] == gt].set_index(["side", "subcluster"])
        lib = {s: g.xs(s, level="side")["library_m"].iloc[0] for s in SIDES}
        # a subcluster with no nuclei in one half is a real zero (weight 0), not a default
        n_nuc = {s: np.array([float(g.loc[(s, k), "n_cells"]) if (s, k) in g.index else 0.0
                              for k in pool]) for s in SIDES}
        U = {s: np.array([float(umi_total[g.loc[(s, k), "sample_id"]]) if (s, k) in g.index else 0.0
                          for k in pool]) for s in SIDES}
        for gene in GENES:
            G = {s: np.array([float(gene_umi.loc[gene, g.loc[(s, k), "sample_id"]]) if (s, k) in g.index
                              else 0.0 for k in pool]) for s in SIDES}
            r = decompose(U, G)
            for i, k in enumerate(pool):
                Ub = U["L"][i] + U["R"][i]
                rows.append(dict(
                    genotype=gt, sex=sex[gt], stage_anchor=STAGE,
                    library_L=lib["L"], library_R=lib["R"], gene=gene, subcluster=k,
                    n_nuclei_L=int(n_nuc["L"][i]), n_nuclei_R=int(n_nuc["R"][i]),
                    pool_nuclei_L=int(n_nuc["L"].sum()), pool_nuclei_R=int(n_nuc["R"].sum()),
                    share_nuclei_L=n_nuc["L"][i] / n_nuc["L"].sum(),
                    share_nuclei_R=n_nuc["R"][i] / n_nuc["R"].sum(),
                    umi_L=U["L"][i], umi_R=U["R"][i],
                    weight_L=r["w"]["L"][i], weight_R=r["w"]["R"][i],
                    gene_umi_L=G["L"][i], gene_umi_R=G["R"][i],
                    cpm_L=r["c"]["L"][i], cpm_R=r["c"]["R"][i],
                    rate_borrowed_L=bool(r["borrowed"]["L"][i]),
                    rate_borrowed_R=bool(r["borrowed"]["R"][i]),
                    cpm_both_halves=((G["L"][i] + G["R"][i]) / Ub * 1e6) if Ub > 0 else np.nan,
                    composition_k=r["comp_k"][i], within_k=r["within_k"][i],
                    P_L=r["P_L"], P_R=r["P_R"], delta=r["delta"],
                    composition=r["composition"], within=r["within"],
                ))
    out = pd.DataFrame(rows)
    # Re-assert the exact sum on the assembled table (per embryo x gene), in memory,
    # before serialization; main() re-checks the TSV as written (6 significant digits).
    eg = out.groupby(["genotype", "gene"]).agg(
        comp=("composition_k", "sum"), within=("within_k", "sum"), delta=("delta", "first"))
    if not np.allclose(eg["comp"] + eg["within"], eg["delta"], rtol=RTOL, atol=ATOL):
        raise AssertionError("assembled table: composition + within != delta")
    return out, floor, med


def verify_written_table(path: Path) -> float:
    """Re-read the main TSV as written and check sum_k(composition_k + within_k) = delta
    per embryo x gene, to within the 6-significant-digit rounding of the written terms
    (each value carries <= 5e-6 relative error; the bound is 2x the summed worst case).
    Returns the largest residual in CPM."""
    t = pd.read_csv(path, sep="\t")
    t["abs_terms"] = t["composition_k"].abs() + t["within_k"].abs()
    eg = t.groupby(["genotype", "gene"]).agg(
        comp=("composition_k", "sum"), within=("within_k", "sum"), delta=("delta", "first"),
        mag=("abs_terms", "sum"))
    resid = (eg["comp"] + eg["within"] - eg["delta"]).abs()
    bound = 2 * 5e-6 * (eg["mag"] + eg["delta"].abs())
    if (resid > bound).any():
        raise AssertionError(f"{path}: composition + within != delta beyond 6-digit rounding: "
                             f"{eg.loc[resid > bound].index.tolist()}")
    return float(resid.max())


def summarise(out: pd.DataFrame, floor: int, med: int) -> pd.DataFrame:
    recs = []

    def add(metric, value, gene="", subcluster=""):
        recs.append(dict(metric=metric, gene=gene, subcluster=subcluster, value=value))

    n = out["genotype"].nunique()
    add("n_paired_embryos", n)
    add("pool_floor_nuclei", floor)
    add("pool_median_nuclei", med)
    emb = out[out["gene"] == GENES[0]]
    for k, d in emb.groupby("subcluster", sort=False):
        diff = d["share_nuclei_R"].to_numpy() - d["share_nuclei_L"].to_numpy()
        # Exact signed-rank p by enumeration (tie-correct; zero differences dropped and
        # counted -- convention in figures/_shared/signed_rank.py).  An undefined p (no
        # non-zero difference) is written as NaN; panel F refuses to print it.
        p, _, n_zero = signed_rank_exact(diff)
        add("share_nuclei_median_L", d["share_nuclei_L"].median(), subcluster=k)
        add("share_nuclei_median_R", d["share_nuclei_R"].median(), subcluster=k)
        add("share_nuclei_n_R_gt_L", int((diff > 0).sum()), subcluster=k)
        add("share_nuclei_wilcoxon_p_exact", p, subcluster=k)
        add("share_nuclei_n_zero_diff", n_zero, subcluster=k)
        add("share_umi_median_L", d["weight_L"].median(), subcluster=k)
        add("share_umi_median_R", d["weight_R"].median(), subcluster=k)
    per = out.groupby(["gene", "genotype"], sort=False).first().reset_index()
    for gene, d in per.groupby("gene", sort=False):
        add("pooled_cpm_mean_L", d["P_L"].mean(), gene)
        add("pooled_cpm_mean_R", d["P_R"].mean(), gene)
        add("delta_mean", d["delta"].mean(), gene)
        add("composition_mean", d["composition"].mean(), gene)
        add("within_mean", d["within"].mean(), gene)
        add("composition_pct_of_delta", 100 * d["composition"].mean() / d["delta"].mean(), gene)
        add("n_embryos_R_gt_L", int((d["delta"] > 0).sum()), gene)
        # panel F prints these two counts above its composition / within columns ("k/n > 0")
        add("within_n_R_gt_L", int((d["within"] > 0).sum()), gene)
        add("composition_n_R_gt_L", int((d["composition"] > 0).sum()), gene)
    both = out.pivot_table(index=["gene", "genotype"], columns="subcluster",
                           values="cpm_both_halves", dropna=False)
    for gene, d in both.groupby(level="gene", sort=False):
        if set(d.columns) != {FOCAL, OTHER}:
            raise KeyError(f"{gene}: both-halves CPM columns {sorted(d.columns)} != {FOCAL}, {OTHER}")
        focal, other = d[FOCAL].to_numpy(), d[OTHER].to_numpy()
        # Per-embryo fold FOCAL / OTHER on CPM over both halves (gene UMIs summed over
        # L + R / total UMIs summed over L + R).  Undefined -- excluded from the median
        # and counted -- when it is 0/0 (gene absent from both subclusters) or a
        # subcluster has no nuclei in either half (CPM NaN).  Zero CPM in OTHER with
        # > 0 in FOCAL is an infinitely higher fold: kept as +inf (the median is
        # rank-based) and counted, not hidden.  The median itself is written as is:
        # +inf when most embryos have a zero denominator, NaN when every fold is
        # undefined; both are explicit in the summary next to their counts.
        undefined = np.isnan(focal) | np.isnan(other) | ((focal == 0) & (other == 0))
        zero_den = ~undefined & (other == 0)
        f, o = focal[~undefined], other[~undefined]
        fold = np.where(o == 0, np.inf, f / np.where(o == 0, 1.0, o))
        add("n_embryos_focal_gt_other", int((f > o).sum()), gene,
            f"{FOCAL} > {OTHER}")
        add("fold_focal_over_other_median", float(np.median(fold)) if len(fold) else np.nan,
            gene, f"{FOCAL} / {OTHER}")
        add("fold_focal_over_other_n_zero_denominator", int(zero_den.sum()), gene,
            f"{FOCAL} / {OTHER}")
        add("fold_focal_over_other_n_undefined", int(undefined.sum()), gene,
            f"{FOCAL} / {OTHER}")
    for (gene, k), d in out.groupby(["gene", "subcluster"], sort=False):
        add("cpm_both_halves_mean", d["cpm_both_halves"].mean(), gene, k)
        add("cpm_both_halves_median", d["cpm_both_halves"].median(), gene, k)
        add("cpm_median_L", d["cpm_L"].median(), gene, k)
        add("cpm_median_R", d["cpm_R"].median(), gene, k)
        add("within_k_mean", d["within_k"].mean(), gene, k)
        add("composition_k_mean", d["composition_k"].mean(), gene, k)
        add("gene_umi_sum_L", float(d["gene_umi_L"].sum()), gene, k)
        add("gene_umi_sum_R", float(d["gene_umi_R"].sum()), gene, k)
        # Per-subcluster L/R expression test (context only, not plotted) on embryos
        # with an own rate on both sides: not borrowed (a borrowed rate is identical on
        # both sides by construction) and not absent from both halves (rate NaN).
        nb = d[~d["rate_borrowed_L"] & ~d["rate_borrowed_R"]
               & d["cpm_L"].notna() & d["cpm_R"].notna()]
        diff = nb["cpm_R"].to_numpy() - nb["cpm_L"].to_numpy()
        # Zero differences (e.g. the gene has 0 UMIs in this subcluster in both halves)
        # are valid data: dropped from the test and counted; with no non-zero
        # difference the p is undefined and written as NaN, not an abort.
        p, _, n_zero = signed_rank_exact(diff)
        add("cpm_lr_n_embryos_unborrowed", len(nb), gene, k)
        add("cpm_lr_median_L_unborrowed", nb["cpm_L"].median(), gene, k)
        add("cpm_lr_median_R_unborrowed", nb["cpm_R"].median(), gene, k)
        add("cpm_lr_n_R_gt_L_unborrowed", int((diff > 0).sum()), gene, k)
        add("cpm_lr_wilcoxon_p_exact_unborrowed", p, gene, k)
        add("cpm_lr_n_zero_diff_unborrowed", n_zero, gene, k)
    return pd.DataFrame(recs)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--h1-dir", type=Path, default=None,
                    help="REQUIRED for a build (no default): directory holding the h1 snRNA "
                         "inference-pool pseudobulk counts + metadata, e.g. figures/data")
    ap.add_argument("--e1b", type=Path, default=REPO_ROOT / "figures" / "data" / "e1b_stage_assignments.tsv",
                    help="per-genotype sex table, context column only "
                         "(default figures/data/e1b_stage_assignments.tsv)")
    ap.add_argument("--out-dir", type=Path, default=REPO_ROOT / "figures" / "data",
                    help="output directory (default figures/data, where the panel F script reads)")
    ap.add_argument("--self-check", action="store_true",
                    help="run the synthetic decomposition checks and exit")
    a = ap.parse_args()
    if not a.self_check and a.h1_dir is None:
        ap.error("--h1-dir is required (no default: name the h1 you mean to decompose)")
    self_check()
    if a.self_check:
        return
    out, floor, med = build(a.h1_dir.resolve(), a.e1b.resolve())
    summ = summarise(out, floor, med)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out_dir / OUT_MAIN, sep="\t", index=False, float_format="%.6g")
    summ.to_csv(a.out_dir / OUT_SUMMARY, sep="\t", index=False, float_format="%.6g",
                na_rep="NaN")   # an undefined statistic is written explicitly, never blank
    max_resid = verify_written_table(a.out_dir / OUT_MAIN)
    print(f"[h2c] written table: composition + within = delta to within 6-digit rounding "
          f"(max residual {max_resid:.2g} CPM)")
    print(f"[h2c] h1 input: {a.h1_dir.resolve()}")
    print(f"[h2c] e1b input: {a.e1b.resolve()}")
    print(f"[h2c] wrote {a.out_dir / OUT_MAIN} ({len(out)} rows) and {a.out_dir / OUT_SUMMARY}")
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
