#!/usr/bin/env python3
"""Locus-restricted BH-FDR for h6 paired + h7 unpaired DESeq2 outputs.

Restricts BH multiple-testing correction to the PITX2 focal-locus window.
The default window is the ±500 kb convention window extended modestly
upstream to the FAM241A gene-locus border (359 peaks,
`pitx2_locus_tad_peaks.tsv`, from derive_pitx2_locus_tad_peaks.py — the
filename keeps a legacy `tad` tag, but we do NOT claim a duck-measured TAD).
The strict ±500 kb set (329 peaks, `h0_lvr_peaks_pitx2_locus.tsv`, a subset)
stays selectable via `LOCUS_INDEX=h0`.  Window rationale:
docs/pitx2_locus_window_rationale.md (not distributed).
Parallels scripts/compute_pitx2_locus_fdr.py (g2c/h0 binom screening).

Outer-joins against the full locus index so peaks dropped upstream
(e.g. h7 MIN_PEAK_COUNT=10 filter at low-depth HH13) appear as explicit
NaN rows in the output rather than silently absent. BH-padj denominator
is still the count of TESTED hypotheses (`n_locus_tested` column in
locus_fdr_summary_h6_h7.tsv), not the locus-index size — peaks that were
never tested don't contribute to the multiple-testing burden.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "figures" / "data"

# Focal-locus window for the BH denominator.  Default = the ±500 kb convention
# window extended upstream to the FAM241A gene-locus border (359 peaks; NOT a
# duck-measured TAD — the filename keeps a legacy `tad` tag only); the strict
# ±500 kb set (329, a subset) via LOCUS_INDEX=h0; or an explicit path via
# LOCUS_INDEX=/abs/path.tsv.  The resolved window label is written to the
# `locus_window` column of locus_fdr_summary_h6_h7.tsv.
_LOCUS_CHOICES = {
    "default": ("pitx2_500kb_fam241a_359", DATA / "pitx2_locus_tad_peaks.tsv"),
    "h0":      ("pitx2_500kb_329", DATA / "h0_lvr_peaks_pitx2_locus.tsv"),
}
# Lowercase ONLY for matching the named choices; a custom path must keep its
# original case (scratch paths like /n/scratch/.../FC_08656/... are case-sensitive,
# so lowercasing the whole value would break load_locus_index()).
_locus_raw = os.environ.get("LOCUS_INDEX", "default")
_locus_key = _locus_raw.lower()
if _locus_key in _LOCUS_CHOICES:
    LOCUS_WINDOW, LOCUS_INDEX = _LOCUS_CHOICES[_locus_key]
else:
    LOCUS_WINDOW, LOCUS_INDEX = (f"custom:{Path(_locus_raw).name}", Path(_locus_raw))
CANDIDATES_BED = REPO / "references" / "pitx2_candidate_peaks.bed"

# Headline candidate set — the ONLY named peaks reported in the locus-FDR
# candidate table.  Reduced 2026-06-19 from the 8-peak grab-bag to the peaks
# with a standing basis: the PITX2c promoter (a-priori — it IS the promoter),
# the +5.4 kb peak (Peak_3; independent motif basis, see fig04), and the
# PITX2a promoter (bilateral within-locus control).  ASE_region / SMAD5plus /
# Peak_2 / Peak_4 / Peak_5 are retired from the headline (significance-driven
# or since-reframed selection); the BED keeps them for the exploratory motif /
# co-accessibility scans that document why they are weak.  Map = BED name ->
# clearer display label.
HEADLINE_CANDIDATES = {
    "Peak_1": "PITX2c_promoter",
    "Peak_3": "PITX2_+5.4kb",
    "PITX2a_promoter": "PITX2a_promoter",
}

# h6 inputs are EXPECTED to be locus-pre-filtered (the producer is
# scripts/extract_h6_locus_subset.py, which restricts to the focal-locus
# window — the FAM241A-extended ±500 kb set by default — before writing).  We re-apply isin(locus_ids)
# defensively so a future replacement of these files with the full
# genome-wide h6 output would be caught (BH on 200K peaks mislabeled
# as locus padj would silently invalidate the focal claim).
H6_INPUTS = {
    "HH17_16": DATA / "h6_paired_deseq2_HH17_16_balanced_locus.tsv",
    "HH25":    DATA / "h6_paired_deseq2_HH25_balanced_locus.tsv",
}
H7_INPUTS = {
    "HH13":    DATA / "h7_unpaired_deseq2_HH13.tsv.gz",
    "HH17_16": DATA / "h7_unpaired_deseq2_HH17_16.tsv.gz",
    "HH21":    DATA / "h7_unpaired_deseq2_HH21.tsv.gz",
    "HH25":    DATA / "h7_unpaired_deseq2_HH25.tsv.gz",
}


def load_locus_index() -> pd.DataFrame:
    df = pd.read_csv(LOCUS_INDEX, sep="\t")
    return df[["peak_id", "chrom", "start", "end"]].copy()


def load_candidates() -> pd.DataFrame:
    rows = []
    seen = set()
    # encoding pinned: the BED header carries non-ASCII (§, —, ⚠) so a C/POSIX
    # locale would otherwise UnicodeDecodeError on a comment line we only skip.
    with open(CANDIDATES_BED, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"{CANDIDATES_BED.name}: malformed BED line (need >=4 cols): {line!r}")
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if name not in HEADLINE_CANDIDATES:
                continue  # retired-from-headline candidate (kept in BED for exploratory scans)
            seen.add(name)
            rows.append({"chrom": chrom, "start": start, "end": end,
                         "name": HEADLINE_CANDIDATES[name],
                         "peak_id": f"{chrom}:{start}-{end}"})
    missing = set(HEADLINE_CANDIDATES) - seen
    if missing:
        # Fail loud (repo rule: a missing key is a real error, never a silent
        # default): a renamed/typo'd/retired BED row would otherwise silently
        # shrink the headline candidate table feeding a PUBLISHED figure.
        raise KeyError(
            f"{CANDIDATES_BED.name}: HEADLINE_CANDIDATES not found in BED: "
            f"{sorted(missing)} — fix the BED peak name or HEADLINE_CANDIDATES.")
    return pd.DataFrame(rows)


def compute_locus_fdr(per_peak_tsv: Path, locus_index: pd.DataFrame,
                       pval_col: str = "pvalue",
                       require_locus_subset: bool = False,
                       producer_ids: "set[str] | None" = None,
                       expected_extra_ids: "set[str] | None" = None) -> tuple[pd.DataFrame, int, int]:
    """Restrict to locus peaks, BH-recompute padj, outer-join missing rows as NaN.

    Returns (locus_df, n_tested, n_dropped) where n_tested is rows with
    valid pvalue (used as BH denominator) and n_dropped is locus peaks
    absent from the input.

    `require_locus_subset=True` asserts the input is already restricted to the
    PRODUCER locus window (`producer_ids` — the 359-peak set
    extract_h6_locus_subset.py pre-filters h6 to), used for h6.  The check is
    against the producer window, NOT `locus_index`, so a narrower FDR index
    (`LOCUS_INDEX=h0`, a 329-peak strict subset) can still be applied to a
    359-peak h6 input as a sensitivity tier — the in-producer-not-in-h0 peaks
    are simply dropped from the BH denominator below, not treated as a fatal
    genome-wide input.  A row OUTSIDE the producer window means the file is
    mislabeled (e.g. the full genome-wide h6 output), which would silently
    invalidate the locus padj denominator — so that fails loud.  If
    `producer_ids` is None the check falls back to `locus_index` (legacy).
    """
    if str(per_peak_tsv).endswith(".gz"):
        df = pd.read_csv(per_peak_tsv, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(per_peak_tsv, sep="\t")
    n_input = len(df)
    if require_locus_subset:
        # h6 path: the input must already be restricted to the producer locus
        # window.  Check against the PRODUCER window (not the chosen FDR
        # locus_index), so LOCUS_INDEX=h0 can select a strict 329-peak subset of
        # a 359-peak h6 input without the in-producer-not-in-h0 peaks tripping a
        # false genome-wide abort.  A row outside the producer window IS
        # genome-wide contamination and fails loud (BH on 200K peaks mislabeled
        # as locus padj would silently invalidate the focal claim).
        check_ids = producer_ids if producer_ids is not None else set(locus_index["peak_id"])
        n_outside = int((~df["peak_id"].isin(check_ids)).sum())
        if n_outside > 0:
            raise RuntimeError(
                f"{per_peak_tsv.name}: expected a locus-pre-filtered h6 input "
                f"(subset of the {len(check_ids)}-peak producer window) but found "
                f"{n_outside} rows outside it out of {n_input} — looks genome-wide; "
                f"extract_h6_locus_subset.py may have regressed.")
        # Guard the OPPOSITE error: a STALE NARROWER extract (e.g. a pre-migration
        # h0/329 _locus.tsv fed under the 359 default).  It is a valid SUBSET of
        # the producer window (so the check above passes), but it LACKS the peaks
        # the chosen window adds beyond the strict h0 set, so the BH denominator
        # would silently revert to the narrower window while the output is labeled
        # with the broader one.  A correct extract contains those peaks (26-30 of
        # 30 at the real h6 stages); a stale narrower one contains none.
        if expected_extra_ids and not (set(df["peak_id"]) & expected_extra_ids):
            raise RuntimeError(
                f"{per_peak_tsv.name}: the chosen window adds {len(expected_extra_ids)} "
                f"peaks beyond the strict h0 set, but this input contains NONE of them "
                f"— it looks like a STALE narrower extract (extracted against a smaller "
                f"window, then fed under the broader default), which would silently "
                f"revert the BH denominator. Re-extract h6 against the chosen window.")
    locus_ids = set(locus_index["peak_id"])
    sub = df[df["peak_id"].isin(locus_ids)].copy()
    n_kept = len(sub)

    # BH on the tested peaks
    valid = sub[pval_col].notna()
    n_tested = int(valid.sum())
    sub["padj_locus"] = np.nan
    if n_tested > 0:
        _, padj_locus, _, _ = multipletests(sub.loc[valid, pval_col].values, method="fdr_bh")
        sub.loc[valid, "padj_locus"] = padj_locus

    # Outer-join against the full locus index so dropped peaks appear
    # as explicit NaN rows
    full = locus_index.merge(sub.drop(columns=[c for c in ("chrom", "start", "end") if c in sub.columns]),
                              on="peak_id", how="left")
    n_dropped = len(full) - len(sub)

    print(f"  {per_peak_tsv.name}: input n={n_input}, locus kept n={n_kept}, "
          f"tested (valid {pval_col}) n={n_tested}, "
          f"locus peaks absent from input n={n_dropped}")
    if n_dropped > 0:
        # Surface the drop loudly so callers know n_tested < n_locus
        # (HH13 h7 commonly drops most locus peaks — 322/359 under the default
        # window, 294/329 under LOCUS_INDEX=h0 — because MIN_PEAK_COUNT=10
        # filters most peaks at the 3-pseudobulk depth)
        print(f"    NOTE: BH-FDR denominator = {n_tested} tested peaks "
              f"(NOT {len(locus_index)} locus peaks). {n_dropped} peaks "
              f"failed upstream count filter.")
    return full, n_tested, n_dropped


def main() -> None:
    locus_index = load_locus_index()
    n_locus = len(locus_index)
    print(f"[i] {n_locus} locus peaks from {LOCUS_INDEX.name} (window={LOCUS_WINDOW})")
    candidates = load_candidates()
    print(f"[i] {len(candidates)} named candidates from {CANDIDATES_BED.name}")

    # load_candidates() guards the BED-name path; this guards the COORDINATE path.
    # The per-stage candidate table is built by an inner-join on peak_id, so a
    # candidate whose BED coord isn't in the chosen locus index (e.g. an index
    # regen shifts it out of the window, or LOCUS_INDEX=h0 excludes it) would be
    # silently dropped from the published table.  Fail loud instead.
    missing_in_index = set(candidates["peak_id"]) - set(locus_index["peak_id"])
    if missing_in_index:
        raise KeyError(
            f"headline candidate(s) absent from locus index {LOCUS_INDEX.name} "
            f"(window={LOCUS_WINDOW}): {sorted(missing_in_index)} — the candidate-"
            f"table inner-join would silently drop them.")

    # Producer window for the h6 containment check: the 359-peak default set
    # extract_h6_locus_subset.py pre-filters h6 to, unioned with the chosen FDR
    # index.  Lets LOCUS_INDEX=h0 run the sensitivity tier directly on the
    # committed 359-peak h6 inputs (restrict 359->329 for BH) instead of aborting.
    producer_ids = set(locus_index["peak_id"]) | set(
        pd.read_csv(_LOCUS_CHOICES["default"][1], sep="\t")["peak_id"])
    # Peaks the chosen window adds beyond the strict h0 (329) set — used to detect
    # a stale NARROWER h6 extract (empty when LOCUS_INDEX=h0, so nothing to check).
    h0_ids = set(pd.read_csv(_LOCUS_CHOICES["h0"][1], sep="\t")["peak_id"])
    chosen_extra = set(locus_index["peak_id"]) - h0_ids

    summary_rows = []
    candidate_rows = []

    print("\n=== h6 paired DESeq2 (locus-pre-filtered inputs) ===")
    for stage, src in H6_INPUTS.items():
        if not src.exists():
            print(f"  [{stage}] missing {src.name} — sync from scratch via "
                  f"scripts/extract_h6_locus_subset.py, then rerun")
            continue
        full, n_tested, n_dropped = compute_locus_fdr(
            src, locus_index, pval_col="pvalue", require_locus_subset=True,
            producer_ids=producer_ids, expected_extra_ids=chosen_extra)
        out_path = DATA / f"h6_paired_deseq2_{stage}_balanced_locus_fdr.tsv"
        full.to_csv(out_path, sep="\t", index=False)
        print(f"  [{stage}] wrote {out_path.name}: n={len(full)} locus peaks "
              f"({n_tested} tested, {n_dropped} not tested)")

        n_l_010 = ((full["log2FoldChange"] > 0) & (full["padj_locus"] < 0.10)).sum()
        n_r_010 = ((full["log2FoldChange"] < 0) & (full["padj_locus"] < 0.10)).sum()
        n_l_005 = ((full["log2FoldChange"] > 0) & (full["padj_locus"] < 0.05)).sum()
        n_r_005 = ((full["log2FoldChange"] < 0) & (full["padj_locus"] < 0.05)).sum()
        summary_rows.append({
            "stage": stage, "classifier": "h6_paired",
            "locus_window": LOCUS_WINDOW,
            "n_locus_in_index": n_locus, "n_locus_tested": n_tested,
            "n_locus_dropped": n_dropped,
            "n_L_locus_padj_010": int(n_l_010), "n_R_locus_padj_010": int(n_r_010),
            "n_L_locus_padj_005": int(n_l_005), "n_R_locus_padj_005": int(n_r_005),
        })

        named = full.merge(candidates[["peak_id", "name"]], on="peak_id", how="inner")
        for _, r in named.iterrows():
            candidate_rows.append({
                "stage": stage, "classifier": "h6_paired",
                "name": r["name"], "peak_id": r["peak_id"],
                "log2FC": r.get("log2FoldChange"), "pvalue": r.get("pvalue"),
                "padj_genomewide": r.get("padj"), "padj_locus": r["padj_locus"],
            })

    print("\n=== h7 unpaired DESeq2 (genome-wide inputs) ===")
    for stage, src in H7_INPUTS.items():
        if not src.exists():
            print(f"  [{stage}] missing {src.name}")
            continue
        full, n_tested, n_dropped = compute_locus_fdr(
            src, locus_index, pval_col="pvalue", require_locus_subset=False)
        out_path = DATA / f"h7_unpaired_deseq2_{stage}_locus_fdr.tsv"
        full.to_csv(out_path, sep="\t", index=False)
        print(f"  [{stage}] wrote {out_path.name}: n={len(full)} locus peaks "
              f"({n_tested} tested, {n_dropped} not tested)")

        n_l_010 = ((full["log2FoldChange"] > 0) & (full["padj_locus"] < 0.10)).sum()
        n_r_010 = ((full["log2FoldChange"] < 0) & (full["padj_locus"] < 0.10)).sum()
        n_l_005 = ((full["log2FoldChange"] > 0) & (full["padj_locus"] < 0.05)).sum()
        n_r_005 = ((full["log2FoldChange"] < 0) & (full["padj_locus"] < 0.05)).sum()
        summary_rows.append({
            "stage": stage, "classifier": "h7_unpaired",
            "locus_window": LOCUS_WINDOW,
            "n_locus_in_index": n_locus, "n_locus_tested": n_tested,
            "n_locus_dropped": n_dropped,
            "n_L_locus_padj_010": int(n_l_010), "n_R_locus_padj_010": int(n_r_010),
            "n_L_locus_padj_005": int(n_l_005), "n_R_locus_padj_005": int(n_r_005),
        })

        named = full.merge(candidates[["peak_id", "name"]], on="peak_id", how="inner")
        for _, r in named.iterrows():
            candidate_rows.append({
                "stage": stage, "classifier": "h7_unpaired",
                "name": r["name"], "peak_id": r["peak_id"],
                "log2FC": r.get("log2FoldChange"), "pvalue": r.get("pvalue"),
                "padj_genomewide": r.get("padj"), "padj_locus": r["padj_locus"],
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = DATA / "locus_fdr_summary_h6_h7.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[i] wrote {summary_path.name}")

    cand_df = pd.DataFrame(candidate_rows)
    cand_path = DATA / "locus_fdr_candidates_h6_h7.tsv"
    cand_df.to_csv(cand_path, sep="\t", index=False)
    print(f"[i] wrote {cand_path.name}")

    print("\n=== Per-stage / per-classifier locus-FDR counts ===")
    print(summary_df.to_string(index=False))

    print("\n=== Per-candidate locus-padj across stages/classifiers ===")
    pivot = cand_df.pivot_table(index="name", columns=["classifier", "stage"],
                                 values="padj_locus", aggfunc="first")
    print(pivot.to_string(float_format=lambda x: f"{x:.4g}" if pd.notna(x) else "-"))


if __name__ == "__main__":
    main()
