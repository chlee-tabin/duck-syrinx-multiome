#!/usr/bin/env python3
"""hh21_trans_motif.py — does PITX2-motif content explain the HH21 DEG->promoter-DA residual?

BACKGROUND.  HH21 L/R RNA DEGs show only
partial concordance with L/R accessibility at their own promoters (July substrate:
gene-level 62.3%, per-peak 58.0%, 2/215 promoters clearing per-feature FDR; on the
2026-08 seeded relabel with PITX2's own locus elements held out the HH21 rows are
36/64 genes and 50/90 promoters, and 2 of 94 unique promoters clear the per-feature
HH21 FDR — PITX2c and the +5.4 kb element; at the TROUGH only PITX2c clears the
response-family FDR. Quote figures/data/hh21_deg_da_concordance*.tsv, not these
background figures).  The trans-skew reading: the
promoter is bilaterally OPEN and the asymmetry is supplied in trans by left-restricted
PITX2 — so a chromatin-flat promoter of an asymmetric gene is the EXPECTED case, not a
detection failure.

TEST (pre-registered).  Split the response promoters three ways on the SIGNED DESeq2
z-statistic (stat = log2FC/lfcSE, HH21 unpaired h7) rather than on significance (only 2
would survive) or on raw log2FC (DESeq2 inflates |log2FC| at low baseMean, so a raw-lfc
tercile sorts on coverage — which covaries with promoter GC class):

    concordant     top tercile |z|, sign agrees with the RNA direction
    flat           bottom tercile |z|                       <- the "residual"
    anti-concordant top tercile |z|, sign opposes the RNA direction  <- internal noise null

Each is compared against BILATERAL promoters (HH21 Ventral-spl. genes with padj>=0.1),
STRATIFIED BY DIRECTION (L-up / R-up).  PITX2 is left-restricted, so only the L-up arm can
be a PITX2 story; an R-up arm that enriches equally means a generic DEG-promoter-class
(AT-richness) effect, not PITX2.

    PRIMARY COEFFICIENT: L-up_flat vs bilateral, PITX2 motif, one-sided (> 0).

Everything else in the model is context or control.  The full g3f 12-motif JASPAR panel is
scanned in the same job so the claim can be a SPECIFICITY claim ("PITX2 is the only motif
in the panel that separates them"), BH-adjusted across the 12.  A raw TAATCC/GGATTA count
rides along as the GC-immune sanity check (same guardrail idiom as g3f's raw AGAC SBE).

Unit of analysis = the unique promoter PEAK.  A peak serving both a response and a
bilateral gene is assigned to the response set; a peak serving two response genes takes the
direction of the smaller RNA padj.  PITX2's own locus (its c/a/b promoters and the +5.4 kb
element) is EXCLUDED — those are PITX2 cis-elements, not downstream response promoters.

Two stages:
  prep  (cluster) GTF -> transcript TSS -> containing ArchR union peak -> response/bilateral
                  labels; writes the promoter map + a BED for the R motif scan.
  test  (local)   joins the motif scan, applies the z split, fits the GLM, writes stats.

Run:  prep  sbatch pipeline/run_hh21_trans_motif.sh          (does prep -> R scan)
      test  conda run -n duckfig python3 figures/figS03/hh21_trans_motif.py test
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

import grn_memory_node_scan as gm  # TSS parsing / peak mapping / BH — reused verbatim

HERE = os.path.dirname(os.path.abspath(__file__))
NW = os.path.dirname(os.path.dirname(HERE))
DATA = f"{NW}/figures/data"

DEG_TSV = f"{DATA}/h2_deg_SUB_Ventral_spl__HH21.tsv"   # same response universe as grn_memory_node_scan
H4_TSV = f"{DATA}/h4_peak_empirical_p_HH21.tsv.gz"     # HH21 union peaks (coords live in peak_id)
H7_TSV = f"{DATA}/h7_unpaired_deseq2_HH21.tsv.gz"      # HH21 promoter L/R DESeq2 (stat, lfcSE, baseMean)

OUT_MAP = f"{DATA}/hh21_trans_promoter_map.tsv"
OUT_BED = f"{DATA}/hh21_trans_promoters.bed"
SCAN_TSV = f"{DATA}/hh21_trans_promoter_motifs.tsv"    # written by hh21_trans_motif_scan.R
OUT_STATS = f"{DATA}/hh21_trans_motif_stats.tsv"
OUT_GROUPS = f"{DATA}/hh21_trans_motif_groups.tsv"

PADJ = 0.10          # response cut, matches grn_memory_node_scan
PRIMARY_MOTIF = "PITX2"
PRIMARY_GROUP = "L-up_flat"


# ── prep (cluster: needs the cdw GTF) ────────────────────────────────────────────────
def prep(gtf_tx: str) -> None:
    d = pd.read_csv(DEG_TSV, sep="\t")
    d = d[d["padj"].notna()].copy()
    resp = d["padj"] < PADJ
    d["set"] = np.where(resp, "response", "bilateral")
    d["direction"] = np.where(resp, np.where(d["log2FoldChange"] > 0, "L-up", "R-up"), "bilateral")
    print(f"[prep] genes: response={int(resp.sum())} "
          f"(L-up={int((d.direction=='L-up').sum())} R-up={int((d.direction=='R-up').sum())}) "
          f"bilateral={int((~resp).sum())}")

    # reuse gm's parser/mapper, re-pointed: the full transcript GTF (not gm's response-only grep)
    # and the HH21 h4 table (gm defaults to HH17_16 — same union peak set, but be explicit).
    gm.GTF_TX = gtf_tx
    gm.H4_TSV = H4_TSV
    tss = gm.parse_transcript_tss(set(d["gene_id"]))
    peaks = gm.load_union_peaks()
    tmap = gm.map_tss_to_containing_peak(tss, peaks)
    has = tmap[tmap["peak_id"].notna()].merge(
        d[["gene_id", "set", "direction", "log2FoldChange", "padj", "baseMean"]], on="gene_id")
    has = has.rename(columns={"log2FoldChange": "rna_log2FC", "padj": "rna_padj",
                              "baseMean": "rna_baseMean"})
    print(f"[prep] TSS inside a called peak: {len(has)} rows, "
          f"{has.gene_id.nunique()} genes, {has.peak_id.nunique()} peaks")

    # one row per (gene, peak); then collapse to one row per PEAK
    per_gp = (has.sort_values("rna_padj")
                 .drop_duplicates(["gene_id", "peak_id"])
                 .reset_index(drop=True))
    # response wins any peak it shares with a bilateral gene; within response, smallest padj wins
    per_gp["_rank"] = np.where(per_gp["set"] == "response", 0, 1)
    per_peak = (per_gp.sort_values(["_rank", "rna_padj"])
                      .drop_duplicates("peak_id")
                      .drop(columns="_rank")
                      .reset_index(drop=True))
    shared = int((per_gp.groupby("peak_id")["set"].nunique() > 1).sum())
    print(f"[prep] unique promoter peaks: {len(per_peak)} "
          f"(response={int((per_peak['set'] == 'response').sum())}, "
          f"bilateral={int((per_peak['set'] == 'bilateral').sum())}; {shared} peaks shared across sets)")

    coords = per_peak["peak_id"].str.extract(gm.PEAKID_RE)   # same parser as grn_memory_node_scan
    if coords.isna().any().any():
        raise RuntimeError("unparseable peak_id(s): " +
                           str(per_peak.loc[coords[0].isna(), "peak_id"].head().tolist()))
    per_peak["chrom"] = coords[0]
    per_peak["start"] = coords[1].astype(int)
    per_peak["end"] = coords[2].astype(int)

    # drop PITX2's own cis-elements — they are the reference locus, not response promoters
    c, lo, hi = gm.load_locus_window()
    in_locus = (per_peak["chrom"] == c) & (per_peak["end"] >= lo) & (per_peak["start"] <= hi)
    print(f"[prep] dropping {int(in_locus.sum())} peaks inside the PITX2 locus {c}:{lo}-{hi}")
    # Fail loud: PITX2's TSS is inside this window by construction, so a zero match means the
    # contig was renamed or the assembly moved. Silently admitting PITX2's own a/b/c promoters
    # into the L-up arm would manufacture the very result this analysis tests for.
    # Assert against the DEG UNIVERSE, not against per_peak — per_peak is the table this guard
    # protects, and PITX2 can be dedup'd out of it by a co-located gene with a smaller padj, which
    # would make the guard silently unfirable exactly under the drift it exists to catch.
    if "PITX2" not in set(d["gene_id"]):
        raise RuntimeError("PITX2 is not in the DESeq2-tested universe — this analysis is defined "
                           "relative to it; check DEG_TSV before trusting any exclusion")
    pitx2_mapped = has.loc[has["gene_id"] == "PITX2", "peak_id"].unique()
    if len(pitx2_mapped) == 0:
        raise RuntimeError("PITX2 has no TSS inside a called peak — cannot verify its own "
                           "promoters were excluded from the response set")
    # Check what the filter LEAVES BEHIND, not what it is about to remove: `in_locus` peaks are
    # dropped two lines down, so testing membership alone fires on the normal case every time.
    surviving = per_peak["peak_id"].isin(pitx2_mapped) & ~in_locus
    if surviving.any():
        raise RuntimeError(f"{int(surviving.sum())} PITX2 promoter peak(s) would survive the locus "
                           f"filter {c}:{lo}-{hi} — PITX2's own promoters would enter the L-up arm")
    per_peak = per_peak[~in_locus].reset_index(drop=True)

    per_peak.to_csv(OUT_MAP, sep="\t", index=False)
    per_peak[["chrom", "start", "end", "peak_id"]].to_csv(
        OUT_BED, sep="\t", index=False, header=False)
    print(f"[prep] wrote {OUT_MAP} and {OUT_BED} ({len(per_peak)} peaks)")


# ── test (local) ─────────────────────────────────────────────────────────────────────
def _assign_groups(m: pd.DataFrame) -> pd.DataFrame:
    """Three-way signed-z split of the response peaks; bilateral peaks stay one group."""
    resp = m["set"] == "response"
    z = m["atac_z"]
    hi = np.nanquantile(z[resp].abs(), 2 / 3)
    lo = np.nanquantile(z[resp].abs(), 1 / 3)
    want = np.where(m["direction"] == "L-up", 1.0, -1.0)   # RNA direction the chromatin should follow
    same = np.sign(z) == want

    grp = pd.Series("bilateral", index=m.index, dtype=object)
    grp[resp & (z.abs() <= lo)] = "flat"
    grp[resp & (z.abs() >= hi) & same] = "concordant"
    grp[resp & (z.abs() >= hi) & ~same] = "anti"
    grp[resp & (z.abs() > lo) & (z.abs() < hi)] = "mid"    # kept in the model, never interpreted
    m = m.copy()
    m["group"] = np.where(resp, m["direction"].astype(str) + "_" + grp, "bilateral")
    m["z_tercile_lo"], m["z_tercile_hi"] = lo, hi
    return m


def bh_strict(pvals: np.ndarray, where: str) -> np.ndarray:
    """BH, but refuse to silently shrink the family on a non-finite p.

    gm._bh sizes the family by the NON-NaN count, so one failed IRLS fit would quietly drop a
    motif from a pre-registered 12-motif family and inflate every surviving q. Verified:
    gm._bh([0.01, 0.02, nan, nan]) -> [0.02, 0.02, nan, nan], i.e. n=2, not 4.
    """
    pv = np.asarray(pvals, float)
    if not np.isfinite(pv).all():
        raise RuntimeError(f"{int((~np.isfinite(pv)).sum())} non-finite p in the BH family at "
                           f"{where} — a GLM failed to converge; the correction would be "
                           "anti-conservative. Fix the fit or mark the row excluded explicitly.")
    return gm._bh(pv)


def fit_panel(m: pd.DataFrame, motifs: list[str], rhs: str,
              offset: np.ndarray | None = None,
              family: str = "poisson") -> pd.DataFrame:
    """GLM (HC0) per motif -> tidy rows, one per C(group) contrast vs the reference.

    Shared by the ATAC-peak analysis (test(), below) and the TSS-anchored sequence analysis
    (hh21_trans_motif_seq.py), so both report identical statistics. `m` must carry a
    Categorical `group` whose FIRST category is the reference level.

    RETURNED SCHEMA (fixed — does NOT vary with the reference level):
      motif, group, n, mean_sites, mean_sites_ref, ref_group, log_rate_ratio, se, z,
      p_two_sided, p_one_sided_gt, p_one_sided_lt, degenerate
    `ref_group` names the comparator on every row, so a caller that re-levels the model (see
    hh21_trans_motif_seq._contrast_l_vs_r) cannot leave a column named for the wrong reference.
    An earlier version keyed the column as f"mean_sites_{ref}", which shipped the R-up mean under
    a `mean_sites_bilateral` header in two committed tables.

    `family` = "poisson" (counts) or "binomial" (0/1 presence). For a saturating binary outcome a
    count model with a length offset is misspecified, so presence callers must pass "binomial"
    and no offset. `degenerate` flags a fit where the outcome is ~constant in the group or the
    reference (min(rate, 1-rate) * n < 5) — for an abundant motif the presence readout saturates
    and manufactures tiny p-values off a few percentage points.
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from scipy import stats as sps

    if family == "binomial" and offset is not None:
        raise ValueError("binomial + offset is misspecified: a log(bp) offset on a saturating 0/1 "
                         "outcome models a 'presence rate per base pair', which is not a quantity. "
                         "Pass log(width) as an ordinary covariate if scan opportunity needs it.")
    fam = {"poisson": sm.families.Poisson, "binomial": sm.families.Binomial}[family]()
    ref = m["group"].cat.categories[0]
    ref_mask = m["group"] == ref
    rows = []
    for mo in motifs:
        y = f"n_{mo}"
        # Poisson + HC0 robust SE: sandwich covers the mild overdispersion of site counts
        # without committing to an NB dispersion per motif.
        fit = smf.glm(f"{y} ~ {rhs}", data=m, family=fam, offset=offset).fit(cov_type="HC0")
        for term in fit.params.index:
            if not term.startswith("C(group)"):
                continue
            grp = term.split("T.")[1].rstrip("]")
            g_mask = m["group"] == grp
            zval = fit.tvalues[term]
            degen = False
            if family == "binomial":
                degen = any(min(r, 1 - r) * n < 5 for r, n in
                            ((m.loc[k, y].mean(), int(k.sum())) for k in (g_mask, ref_mask)))
            rows.append({
                "motif": mo, "group": grp, "n": int(g_mask.sum()),
                "mean_sites": float(m.loc[g_mask, y].mean()),
                "mean_sites_ref": float(m.loc[ref_mask, y].mean()), "ref_group": str(ref),
                "log_rate_ratio": float(fit.params[term]), "se": float(fit.bse[term]),
                "z": float(zval), "p_two_sided": float(fit.pvalues[term]),
                "p_one_sided_gt": float(sps.norm.sf(zval)),
                "p_one_sided_lt": float(sps.norm.cdf(zval)),
                "degenerate": bool(degen),
            })
    return pd.DataFrame(rows)


def test() -> None:

    m = pd.read_csv(OUT_MAP, sep="\t")
    scan = pd.read_csv(SCAN_TSV, sep="\t")
    # peak mode derives width from the map; drop the scanner's per-interval copies so the
    # merge cannot silently produce width_x/width_y (the scanner grew these for sequence mode).
    scan = scan.drop(columns=[c for c in ("width", "gc_bases") if c in scan.columns])
    h7 = pd.read_csv(H7_TSV, sep="\t", usecols=["peak_id", "baseMean", "stat", "log2FoldChange"])
    h7 = h7.rename(columns={"baseMean": "atac_baseMean", "stat": "atac_z",
                            "log2FoldChange": "atac_log2FC"})

    m = m.merge(scan, on="peak_id", how="inner").merge(h7, on="peak_id", how="inner")
    m = m[m["atac_z"].notna() & (m["atac_baseMean"] > 0) & (m["rna_baseMean"] > 0)].reset_index(drop=True)
    m = _assign_groups(m)
    m["width"] = m["end"] - m["start"]
    m["log_atac"] = np.log2(m["atac_baseMean"])
    m["log_rna"] = np.log2(m["rna_baseMean"])
    # ArchR peaks are fixed-width (501 bp), so width carries no information and would make the
    # design singular. Assert rather than silently omit — a variable-width peak set NEEDS it.
    if m["width"].nunique() != 1:
        raise RuntimeError(f"peak widths are not constant ({sorted(m['width'].unique())[:5]}...) "
                           "— add log2(width) back as a model covariate")

    counts = m["group"].value_counts().to_dict()
    print(f"[test] peaks in model: {len(m)}")
    print(f"[test] z terciles |z| <= {m.z_tercile_lo.iloc[0]:.3f} (flat) / "
          f">= {m.z_tercile_hi.iloc[0]:.3f} (concordant|anti)")
    for g in sorted(counts):
        print(f"         {g:22s} {counts[g]}")
    if counts.get(PRIMARY_GROUP, 0) < 10:
        raise RuntimeError(f"primary group {PRIMARY_GROUP} has n={counts.get(PRIMARY_GROUP,0)} — too few")

    motifs = [c[2:] for c in scan.columns if c.startswith("n_")]
    panel = [x for x in motifs if x != "TAATCC_raw"]   # BH over the 12 JASPAR motifs only;
    #                                                    TAATCC_raw is the GC-immune guardrail
    m["group"] = pd.Categorical(m["group"],
                                categories=["bilateral"] + sorted(set(m["group"]) - {"bilateral"}))
    st = fit_panel(m, motifs, "C(group) + gc + log_atac + log_rna")

    # BH across the JASPAR panel, within the PRIMARY group only (the specificity claim)
    prim = st["group"] == PRIMARY_GROUP
    inpanel = prim & st["motif"].isin(panel)
    st["padj_motif_panel"] = np.nan
    st.loc[inpanel, "padj_motif_panel"] = bh_strict(
        st.loc[inpanel, "p_one_sided_gt"].to_numpy(), "hh21_trans_motif.test")
    st = st.sort_values(["group", "p_one_sided_gt"]).reset_index(drop=True)
    prim = st["group"] == PRIMARY_GROUP          # re-derive: reset_index invalidated the old mask
    st.to_csv(OUT_STATS, sep="\t", index=False, float_format="%.6g")
    m[["peak_id", "gene_id", "set", "direction", "group", "rna_log2FC", "rna_padj",
       "atac_log2FC", "atac_z", "atac_baseMean", "rna_baseMean", "gc", "width"]
      + [f"n_{x}" for x in motifs]].to_csv(OUT_GROUPS, sep="\t", index=False, float_format="%.6g")

    print(f"\n[test] wrote {OUT_STATS} / {OUT_GROUPS}")
    p = st[prim & (st["motif"] == PRIMARY_MOTIF)].iloc[0]
    print(f"\n=== PRIMARY: {PRIMARY_MOTIF} in {PRIMARY_GROUP} vs bilateral ===")
    print(f"  n={p.n}  mean sites {p.mean_sites:.3f} vs {p.ref_group} {p.mean_sites_ref:.3f}")
    print(f"  log rate ratio {p.log_rate_ratio:+.4f} (se {p.se:.4f})  "
          f"one-sided p={p.p_one_sided_gt:.4g}  BH over {len(panel)} motifs = {p.padj_motif_panel:.4g}")
    print(f"\n=== {PRIMARY_GROUP} panel (BH-adjusted, one-sided; TAATCC_raw = raw-hexamer comparator, NOT GC-immune) ===")
    print(st[prim][["motif", "mean_sites", "mean_sites_ref", "ref_group", "log_rate_ratio",
                    "p_one_sided_gt", "padj_motif_panel"]].to_string(index=False))
    print(f"\n=== {PRIMARY_MOTIF} across every group (direction control: R-up should be null) ===")
    print(st[st["motif"] == PRIMARY_MOTIF][["group", "n", "mean_sites", "log_rate_ratio",
                                            "p_two_sided"]].to_string(index=False))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "test"
    if mode == "prep":
        prep(sys.argv[2])
    elif mode == "test":
        test()
    else:
        sys.exit(f"unknown mode {mode!r} (prep <tx.gtf> | test)")
