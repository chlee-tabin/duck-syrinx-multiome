#!/usr/bin/env python3
"""g4h_response_locus_trough_test.py — is the HH21 LEFT-response PROGRAM epigenetically
pre-lateralized in the HH17/16 TROUGH, or is PITX2 the sole cis-memory node?

Hypothesis (the figure this powers): in the inter-wave trough the PITX2c promoter is
ALREADY L-open (chromatin memory), but the genes that only become L-biased at HH21 (the
"response" program, downstream of PITX2 reactivation) are NOT yet lateralized in their
own chromatin.  => PITX2 is the single upstream memory locus; the response cascade merely
follows its reactivation.  This is simultaneously the specificity control for the fig03
keystone AND a GRN-circuitry result (one upstream node -> instructed cascade).

DESIGN = ACCEPT-THE-NULL for the response set, so we never lean on a bare p>0.05:
  - response set       : TOST equivalence  (|mean L/R log2FC| bounded within +/-DELTA)
  - PITX2c promoter    : difference test    (mean log2FC > 0) — the matched-power positive
  - baseMean-matched peak BOOTSTRAP : response/PITX2c aggregate vs a null of size- AND
                         depth-matched random peaks.  THIS is what turns "symmetric"
                         into "genuinely null", not "underpowered".
  - per-locus FDR scan : surfaces any single response locus that IS pre-lateralized
                         (the "PITX2 + N nodes" discovery case — a feature, not a bug)
  - dorsal control     : same response peaks in the DORSAL trough -> expect symmetric

Reads ONLY existing pipeline outputs (h2 snRNA DEGs + a genome-wide per-peak TROUGH L/R
table); NO cluster compute.  The trough L/R source defaults to the paired-balanced h6
table; point TROUGH_LR_TSV at a VENTRAL-compartment-restricted genome-wide trough L/R if
you have one (ideal — the response program is ventral; see NOTE at bottom).

peak_id convention = "chrom:start-end" (matches g4e/g4f/archr_union_peaks.bed.gz).
SIDE convention      = log2FoldChange > 0 == higher on LEFT (h2/g4f use levels R<L).
                       Set SIDE_POS_IS_L=0 if your trough table is signed the other way.

Run (local):  conda run -n duckfig python3 pipeline/g4h_response_locus_trough_test.py
Depends on the HH17/16 trough at the low (20-cell) floor — confirm the trough, not just
HH13, gets the floor relaxation before trusting the equivalence call.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "figures", "data")

# ── inputs (all overridable by env) ──────────────────────────────────────────
H2_DEG_TSV   = os.environ.get("H2_DEG_TSV",   f"{DATA}/h2_deg_POOL_HH21.tsv")          # response source
# Bare-script DEFAULT = the genome-wide h6 table; the committed run OVERRIDES this via run_g4h.sh to the
# VENTRAL g4h_prep trough tables (g4h_lr_ventral_*.tsv) — the response program is ventral (see NOTE at EOF).
TROUGH_LR_TSV= os.environ.get("TROUGH_LR_TSV",f"{DATA}/h6_paired_deseq2_HH17_16_balanced.tsv")
DORSAL_LR_TSV= os.environ.get("DORSAL_LR_TSV","")    # optional dorsal-trough per-peak L/R (compartment ctrl)
GENE_BED     = os.environ.get("GENE_BED",     f"{DATA}/duck_gene_tss.tsv")             # gene_id,chrom,tss,strand
RESPONSE_PEAKS_TSV = os.environ.get("RESPONSE_PEAKS_TSV", "")  # precomputed gene_id,peak_id (else map here)
PITX2_LOCUS_TSV    = os.environ.get("PITX2_LOCUS_TSV", f"{DATA}/h0_lvr_peaks_pitx2_locus.tsv")
PITX2C_PEAKS_TXT   = os.environ.get("PITX2C_PEAKS_TXT", "")    # peak_ids of the c-promoter (positive ctrl)
ABPROM_PEAKS_TXT   = os.environ.get("ABPROM_PEAKS_TXT", "")    # peak_ids of the a/b promoter (within-locus neg)
OUT_DIR      = os.environ.get("OUT_DIR", DATA)

# ── parameters ───────────────────────────────────────────────────────────────
DEG_PADJ   = float(os.environ.get("DEG_PADJ", "0.1"))     # h2 DEG significance for the response set
WINDOW     = int(os.environ.get("WINDOW", "50000"))       # +/- bp of TSS for cCRE assignment (fallback mapper)
PROM_BP    = int(os.environ.get("PROM_BP", "2000"))       # +/- bp of TSS counted as 'promoter'
DELTA      = float(os.environ.get("DELTA", "0.20"))       # TOST equivalence margin (log2FC units)
N_BOOT     = int(os.environ.get("N_BOOT", "2000"))
SEED       = int(os.environ.get("SEED", "0"))
SIDE_POS_IS_L = os.environ.get("SIDE_POS_IS_L", "1") == "1"


def ts(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _col(df: pd.DataFrame, *cands: str) -> str:
    """First matching column (case-insensitive) — tolerate h6/h2 schema drift."""
    low = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    sys.exit(f"ERROR: none of {cands} in columns {list(df.columns)}")


def load_trough_lr(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    pid = _col(df, "peak_id", "peak", "name")
    l2f = _col(df, "log2FoldChange", "log2fc", "lfc")
    bm  = _col(df, "baseMean", "basemean", "mean_cpm")
    padj= _col(df, "padj", "padj_locus")   # ADJUSTED only — never silently treat a raw pvalue as padj
    out = df[[pid, l2f, bm, padj]].copy()
    out.columns = ["peak_id", "log2FC", "baseMean", "padj"]
    out["peak_id"] = out["peak_id"].astype(str)
    if not SIDE_POS_IS_L:
        out["log2FC"] = -out["log2FC"]
    out = out.dropna(subset=["log2FC"])
    ts(f"  trough L/R: {len(out):,} peaks  (log2FC>0 == LEFT)")
    return out


def response_genes() -> list[str]:
    d = pd.read_csv(H2_DEG_TSV, sep="\t")
    g = _col(d, "gene_id", "gene"); lf = _col(d, "log2FoldChange", "log2fc"); pj = _col(d, "padj")
    sig = d[(d[lf] > 0) & (d[pj] < DEG_PADJ)]          # LEFT-biased response genes at HH21
    genes = [x for x in sig[g].astype(str) if x.upper() != "PITX2"]
    ts(f"  response genes (HH21 L-up, padj<{DEG_PADJ}, PITX2 dropped): {len(genes)}")
    return genes


def parse_peakid(pid: str):
    chrom, rng = pid.rsplit(":", 1); s, e = rng.split("-")
    return chrom, (int(s) + int(e)) // 2


def map_genes_to_peaks(genes: list[str], trough_peaks: set[str]) -> pd.DataFrame:
    """gene_id -> peak_id.  Prefer a precomputed map (e.g. ArchR peak2gene); else fall
    back to a TSS-window overlap on the ArchR union (self-contained, coord-based)."""
    if RESPONSE_PEAKS_TSV and os.path.isfile(RESPONSE_PEAKS_TSV):
        m = pd.read_csv(RESPONSE_PEAKS_TSV, sep="\t")
        m = m.rename(columns={_col(m, "gene_id", "gene"): "gene_id",
                              _col(m, "peak_id", "peak"): "peak_id"})
        m = m[m["gene_id"].astype(str).isin(genes)].copy()
        ts(f"  peak map: {len(m)} (gene,peak) pairs from {os.path.basename(RESPONSE_PEAKS_TSV)}")
    else:
        ts("  no RESPONSE_PEAKS_TSV — TSS-window mapping over the trough peak set")
        gb = pd.read_csv(GENE_BED, sep="\t")
        gb = gb.rename(columns={_col(gb, "gene_id", "gene"): "gene_id",
                                _col(gb, "chrom", "chr", "seqname"): "chrom",
                                _col(gb, "tss", "start"): "tss"})
        gb = gb[gb["gene_id"].astype(str).isin(genes)]
        # Build (chrom, mid, peak_id) straight from the trough peak_ids (parse_peakid) — this
        # guarantees the mapped peak_ids are exactly the ones tested in the trough table; reading
        # a separate union bed risks a coord/format/off-by-one mismatch (and the bed column count
        # is not fixed). No UNION_BED dependency.
        op = pd.DataFrame([(*parse_peakid(p), p) for p in trough_peaks],
                          columns=["chrom", "mid", "peak_id"])
        rows = []
        for _, r in gb.iterrows():
            near = op[(op["chrom"] == r["chrom"]) & ((op["mid"] - r["tss"]).abs() <= WINDOW)]
            for _, p in near.iterrows():
                role = "promoter" if abs(p["mid"] - r["tss"]) <= PROM_BP else "ccre"
                rows.append({"gene_id": r["gene_id"], "peak_id": p["peak_id"], "role": role})
        m = pd.DataFrame(rows)
        ts(f"  peak map: {len(m)} (gene,peak) pairs within +/-{WINDOW//1000}kb of TSS "
           f"(from {len(gb)}/{len(genes)} response genes that have a TSS in GENE_BED)")
    m = m[m["peak_id"].astype(str).isin(trough_peaks)].drop_duplicates(["gene_id", "peak_id"])
    ts(f"  peak map after restricting to trough-tested peaks: {len(m)} pairs / "
       f"{m['peak_id'].nunique()} unique peaks / {m['gene_id'].nunique()} genes")
    return m


def _read_peak_list(path: str) -> set[str]:
    if path and os.path.isfile(path):
        return set(pd.read_csv(path, header=None)[0].astype(str))
    return set()


def aggregate(name: str, peaks: set[str], lr: pd.DataFrame) -> dict:
    sub = lr[lr["peak_id"].isin(peaks)]
    if not len(sub):
        return {"set": name, "n_peaks": 0}
    return {"set": name, "n_peaks": len(sub),
            "mean_log2FC": float(sub["log2FC"].mean()),
            "median_log2FC": float(sub["log2FC"].median()),
            "frac_L": float((sub["log2FC"] > 0).mean()),
            "n_L_sig": int(((sub["log2FC"] > 0) & (sub["padj"] < 0.1)).sum()),
            "n_R_sig": int(((sub["log2FC"] < 0) & (sub["padj"] < 0.1)).sum())}


def tost_zero(vals: np.ndarray, delta: float) -> dict:
    """Two one-sided t-tests: is mean(vals) within +/-delta of 0 (equivalence)?"""
    from scipy import stats
    n = len(vals); m = float(np.mean(vals))
    if n < 2:
        return {"n": n, "mean": m, "se": np.nan, "p_tost": np.nan, "equivalent": False}
    se = float(np.std(vals, ddof=1) / np.sqrt(n))
    if se == 0:
        return {"n": n, "mean": m, "se": 0.0,
                "equivalent": abs(m) < delta, "p_tost": 0.0 if abs(m) < delta else 1.0}
    # TOST: reject the union of two one-sided nulls to conclude |mean| < delta.
    #   lower null  H0: mean <= -delta  -> reject when t large POSITIVE -> UPPER tail (sf)
    #   upper null  H0: mean >= +delta  -> reject when t large NEGATIVE -> LOWER tail (cdf)
    # (the lower-bound test uses the upper tail; the prior code had both tails swapped,
    #  which inverted p_tost — a genuinely symmetric set would have read as NOT equivalent.)
    p_lo = stats.t.sf((m - (-delta)) / se, df=n - 1)        # H0: mean <= -delta
    p_hi = stats.t.cdf((m - delta) / se, df=n - 1)          # H0: mean >= +delta
    p = max(p_lo, p_hi)
    return {"n": n, "mean": m, "se": se, "p_tost": float(p), "equivalent": bool(p < 0.05)}


def matched_bootstrap(target_peaks: set[str], lr: pd.DataFrame, rng) -> np.ndarray:
    """Null = mean log2FC of random peak sets matched in SIZE and baseMean-decile to the
    target.  Tells us how extreme the target's L/R aggregate is vs depth-matched chance."""
    lr = lr.copy()
    lr["bm_bin"] = pd.qcut(lr["baseMean"].rank(method="first"), 10, labels=False)
    tgt = lr[lr["peak_id"].isin(target_peaks)]
    bin_counts = tgt["bm_bin"].value_counts().to_dict()
    pools = {b: lr[lr["bm_bin"] == b]["log2FC"].to_numpy() for b in bin_counts}
    null = np.empty(N_BOOT)
    for i in range(N_BOOT):
        draws = [rng.choice(pools[b], size=k, replace=True) for b, k in bin_counts.items()]
        null[i] = float(np.concatenate(draws).mean())
    return null


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ts("=== g4h_response_locus_trough_test ===")
    for p in (H2_DEG_TSV, TROUGH_LR_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing input {p}")

    lr = load_trough_lr(TROUGH_LR_TSV)
    trough_peaks = set(lr["peak_id"])
    genes = response_genes()
    gmap = map_genes_to_peaks(genes, trough_peaks)
    resp_peaks = set(gmap["peak_id"])
    # resp_peaks is the response set AS MAPPED (±50 kb TSS window). It can include the PITX2c peak
    # when a PITX2-locus neighbour gene (LOC106017192/LOC113843539) is itself an HH21 DEG — see
    # ISSUE 4. That overlap is immaterial to the response aggregate (~0.0009 in mean log2FC over the
    # ~1.6k-peak set) and the per-locus scan reports it explicitly, so it is retained rather than
    # special-cased; the random background below DOES exclude resp+locus so the null is never drawn
    # from the response or locus peaks.

    # control peak sets
    locus = pd.read_csv(PITX2_LOCUS_TSV, sep="\t")
    locus_peaks = set(locus[_col(locus, "peak_id", "peak")].astype(str)) & trough_peaks
    c_peaks = _read_peak_list(PITX2C_PEAKS_TXT) & trough_peaks
    if not c_peaks:
        sys.exit("ERROR: PITX2C_PEAKS_TXT must name the c-promoter peak(s) present in the trough "
                 f"table (the matched-power positive control). PITX2C_PEAKS_TXT={PITX2C_PEAKS_TXT!r} "
                 "-> empty intersection with trough peaks. Refusing to silently fall back to the "
                 "whole locus (a silent fallback would dilute the positive control).")
    ab_peaks = _read_peak_list(ABPROM_PEAKS_TXT) & trough_peaks
    # clamp the random-background size against the ACTUAL sampling pool (not len(trough_peaks)), so
    # replace=False can never request more than the pool holds on a smaller/restricted trough table.
    rand_pool = sorted(trough_peaks - resp_peaks - locus_peaks)
    rand_peaks = set(rng.choice(rand_pool, size=min(len(resp_peaks), len(rand_pool)), replace=False))

    # ── aggregates ───────────────────────────────────────────────────────────
    sets = {"PITX2c_promoter": c_peaks, "response_program": resp_peaks,
            "ab_promoter": ab_peaks, "random_matched": rand_peaks}
    agg = pd.DataFrame([aggregate(n, p, lr) for n, p in sets.items() if p])
    ts(f"\n{agg.to_string(index=False)}\n")

    # ── equivalence (response) vs difference (PITX2c) ─────────────────────────
    resp_vals = lr[lr["peak_id"].isin(resp_peaks)]["log2FC"].to_numpy()
    c_vals    = lr[lr["peak_id"].isin(c_peaks)]["log2FC"].to_numpy()
    tost = tost_zero(resp_vals, DELTA)
    from scipy import stats
    c_t, c_p = stats.ttest_1samp(c_vals, 0.0) if len(c_vals) > 1 else (np.nan, np.nan)
    ts(f"  response TOST (delta={DELTA}): mean={tost['mean']:+.3f} p_tost={tost.get('p_tost'):.3g} "
       f"-> {'EQUIVALENT to 0 (no memory)' if tost['equivalent'] else 'NOT equivalent'}")
    c_mean = float(np.mean(c_vals))
    c_verdict = ("single peak (n=1): t-test n/a — read the bootstrap pct" if np.isnan(c_p)
                 else "L-biased (memory)" if (c_p < 0.05 and c_mean > 0) else "n.s.")
    ts(f"  PITX2c diff-test: mean={c_mean:+.3f} t={c_t:.2f} p={c_p:.3g} -> {c_verdict}")

    # ── matched-power bootstrap ───────────────────────────────────────────────
    # NB: `null` is matched to the RESPONSE set's size + baseMean deciles, so it is a distribution
    # of ~|resp_peaks|-peak MEANS (tight). resp_pct (response vs its own null) is a calibrated
    # self-comparison. c_pct scores PITX2c against this SAME response-sized null — for the single-peak
    # c-set that is a MAGNITUDE statement ("exceeds any size-matched random aggregate"), NOT a
    # per-locus significance level (a 1-peak value has far wider spread than this null). See LEGEND caveat 4.
    null = matched_bootstrap(resp_peaks, lr, rng)
    resp_pct = float((null >= np.mean(resp_vals)).mean())
    c_pct    = float((null >= c_mean).mean())
    ts(f"  bootstrap null mean={null.mean():+.3f} sd={null.std():.3f} | "
       f"response pct={resp_pct:.3f} (expect ~0.5) | PITX2c pct={c_pct:.4f} (expect ~0/extreme)")

    # ── per-locus FDR scan (the PITX2 + N nodes case) ─────────────────────────
    per_locus = (gmap.merge(lr, on="peak_id")
                 .sort_values("padj").groupby("gene_id")
                 .first().reset_index()[["gene_id", "peak_id", "log2FC", "padj"]])
    hits = per_locus[(per_locus["padj"] < 0.1)]
    ts(f"  per-locus scan: {len(hits)}/{per_locus['gene_id'].nunique()} response genes "
       f"with a trough L/R-sig peak (padj<0.1) — inspect these as candidate extra memory nodes")

    # ── dorsal control (optional) ─────────────────────────────────────────────
    dorsal_row = None
    if DORSAL_LR_TSV and os.path.isfile(DORSAL_LR_TSV):
        dlr = load_trough_lr(DORSAL_LR_TSV)
        dorsal_row = aggregate("response_program__DORSAL", resp_peaks & set(dlr["peak_id"]), dlr)
        ts(f"  dorsal control: {dorsal_row}")

    # ── write results ─────────────────────────────────────────────────────────
    agg.to_csv(f"{OUT_DIR}/g4h_aggregate.tsv", sep="\t", index=False, float_format="%.5g")
    per_locus.to_csv(f"{OUT_DIR}/g4h_per_locus_scan.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame({
        "metric": ["response_mean_log2FC", "response_tost_p", "response_equivalent",
                   "pitx2c_mean_log2FC", "pitx2c_diff_p", "response_boot_pct", "pitx2c_boot_pct",
                   "delta", "n_response_peaks", "n_response_genes", "n_per_locus_hits"],
        "value": [tost["mean"], tost.get("p_tost"), int(tost["equivalent"]),
                  float(np.mean(c_vals)), c_p, resp_pct, c_pct,
                  DELTA, len(resp_peaks), gmap["gene_id"].nunique(), len(hits)],
    }).to_csv(f"{OUT_DIR}/g4h_summary.tsv", sep="\t", index=False, float_format="%.5g")
    np.savetxt(f"{OUT_DIR}/g4h_bootstrap_null.txt", null, fmt="%.5f")
    # the numeric results are already written above; a panel failure still fails the run
    _panel(lr, sets, null, np.mean(resp_vals), c_mean)
    ts("=== g4h DONE ===")


def _panel(lr, sets, null, resp_mean, c_mean):
    """B&W per house style: per-peak L/R effect strip + the matched-power null."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))
    order = [k for k in ("PITX2c_promoter", "response_program", "ab_promoter", "random_matched")
             if sets.get(k)]
    for i, k in enumerate(order):
        v = lr[lr["peak_id"].isin(sets[k])]["log2FC"].to_numpy()
        x = np.full(len(v), i) + (np.random.default_rng(i).random(len(v)) - 0.5) * 0.3
        ax[0].scatter(x, v, s=6, c="0.25" if "PITX2c" in k else "0.6", alpha=0.6, lw=0)
        ax[0].plot([i - 0.25, i + 0.25], [v.mean()] * 2, c="k", lw=2)
    ax[0].axhline(0, c="k", lw=0.8); ax[0].set_xticks(range(len(order)))
    ax[0].set_xticklabels(order, rotation=20, ha="right", fontsize=8)
    ax[0].set_ylabel("trough L/R log2FC  (>0 = Left)")
    ax[0].set_title("Per-peak laterality in the HH17/16 trough")
    ax[1].hist(null, bins=40, color="0.8", edgecolor="0.5")
    ax[1].axvline(resp_mean, c="k", ls="--", lw=1.5, label="response program")
    ax[1].axvline(c_mean, c="k", lw=1.5, label="PITX2c promoter")
    ax[1].set_xlabel("mean L/R log2FC (baseMean-matched null)")
    ax[1].set_title("Matched-power bootstrap"); ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/panel_response_locus_trough.png", dpi=200)
    ts(f"  -> {OUT_DIR}/panel_response_locus_trough.png")


if __name__ == "__main__":
    main()

# ─────────────────────────────────────────────────────────────────────────────
# NOTE on the trough L/R source.  h6_paired_deseq2_HH17_16_balanced is genome-wide and
# paired but may be POOLED-LPM, not ventral-restricted.  The response program is ventral,
# so the cleaner input is a ventral-compartment trough L/R run — i.e. g4f's lr_dynamics()
# (which already filters to the VENTRAL pool by stage) generalized OFF the 359-locus slice
# to genome-wide.  Two ways to get it:
#   (a) add a "genome_wide=True" path to g4f that skips the locus restriction at HH17/16
#       and writes a g4f_lr_genomewide_HH17_16.tsv (peak_id, log2FoldChange, padj, baseMean);
#       point TROUGH_LR_TSV at it.  Light cluster step, reuses g4f machinery verbatim.
#   (b) build RESPONSE_PEAKS_TSV cluster-side from ArchR peak2gene (mirror g4e), so the
#       gene->peak mapping uses real linkages instead of the +/-50kb TSS window here.
# GENE_BED (gene_id,chrom,tss,strand) builds from the canonical GTF in `cdw`:
#   zcat $cdw/crarc_zju1_PITX2/genes/genes.gtf.gz | awk -F'\t' '$3=="transcript"{...}'  > duck_gene_tss.tsv
