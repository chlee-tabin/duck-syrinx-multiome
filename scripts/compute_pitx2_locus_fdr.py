#!/usr/bin/env python3
"""
Locus-restricted BH-FDR on the 329 PITX2 ±500 kb peaks.

The HH17/16 INFERENCE_POOL ArchR aggregate (`g2c_aggregated_HH1716_
INFERENCE_POOL_archr.tsv`) carries a binom_padj computed under the
GENOME-WIDE multiple-testing burden (n=283,818 peaks).  At the PITX2
locus, only 3 L>R peaks survive genome-wide padj<0.1 — including
PITX2c promoter (Peak #1), PITX2c CDS (Peak #2, gene-body, likely not
a regulator), and a 3' peak (Peak #3).  A 4th candidate, the duck
ortholog of the chick ASE region (`NC_051775.1:3974421-3974922`,
log2FC=+1.05, binom_p=0.003), is buried at genome-wide padj=0.74.

This script asks: when the multiple-testing burden is restricted to the
329 peaks within ±500 kb of the PITX2 gene body (a defensible
locus-level analysis, not a genome-wide screen), which additional
L-biased peaks clear padj<0.1?  Specifically, does the ASE-region peak
clear?

Workflow
--------
1. Read figures/data/h0_lvr_peaks_pitx2_locus.tsv (329 peaks).
2. Apply Benjamini-Hochberg correction over only the 329 raw binom_p
   values (statsmodels.stats.multitest.multipletests, method='fdr_bh').
   This reproduces the same BH variant used genome-wide in
   pipeline/g2c_aggregated_descriptive.py.
3. Write the same TSV with a new `binom_padj_locus_fdr` column,
   sorted by binom_padj_locus_fdr ascending.
4. Stdout: list L-biased (log2FC > 0) peaks that newly clear locus-FDR
   padj<0.1 but did NOT clear genome-wide padj<0.1, and annotate each
   with its position relative to the PITX2 gene model.

The PITX2 gene model is taken from references/pitx2_promoter_regions.bed
and references/pitx2_locus_50kb_features.tsv: PITX2 gene body
3,961,387 - 3,976,254 (+ strand) with three TSS isoforms; the locus
also contains LOC106016804 (3,891,849 - 3,954,422, - strand, upstream
lncRNA), LOC113843538 (3,963,664 - 3,964,114, - strand, inside PITX2
intron), LOC113843539 (3,973,586 - 3,976,858, - strand, overlapping
PITX2 3' end), LOC106017192 (3,976,875 - 4,033,262, + strand, 3' of
PITX2), and LOC101795212 (4,003,281 - 4,037,227, - strand).  See the
BED + locus_50kb_features files for the GTF provenance and the
PITX2a/b/c TSS-assignment correction history.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
IN_TSV = REPO / "figures/data/h0_lvr_peaks_pitx2_locus.tsv"
OUT_TSV = REPO / "figures/data/h0_lvr_peaks_pitx2_locus_locus_fdr.tsv"

# PITX2 gene structure (+ strand, NC_051775.1, ZJU1.0 / GCF_015476345.1).
PITX2_CHROM = "NC_051775.1"
PITX2_GENE_START = 3_961_387
PITX2_GENE_END = 3_976_254
PITX2_TSS = {
    "PITX2a": 3_961_387,
    "PITX2b": 3_964_253,
    "PITX2c": 3_971_385,
}
# All PITX2 exons across the 3 isoforms (start/end, 1-based inclusive
# from the GTF; converted to 0-based half-open at use).
PITX2_EXONS = [
    (3_961_387, 3_961_458),  # PITX2a exon 1 / gene-body 5' start
    (3_964_253, 3_964_353),  # PITX2b exon 1
    (3_966_591, 3_966_646),  # shared exon
    (3_971_385, 3_972_064),  # PITX2c exon 1 (long, contains "Peak #2" range)
    (3_973_108, 3_973_313),  # shared exon
    (3_974_789, 3_976_254),  # last exon (3'-UTR / coding tail)
]

# Other genes/lncRNAs in ±500 kb that we care about for context.  Each
# row: (gene, start, end, strand, role).
NEIGHBOR_GENES = [
    ("LOC106016804", 3_891_849, 3_954_422, "-", "upstream lncRNA"),
    ("LOC113843538", 3_963_664, 3_964_114, "-",
     "antisense lncRNA inside PITX2 intron 1 (between PITX2a-TSS and PITX2b-TSS)"),
    ("LOC113843539", 3_973_586, 3_976_858, "-",
     "antisense lncRNA overlapping PITX2 3' end / 3'-UTR"),
    ("LOC106017192", 3_976_875, 4_033_262, "+",
     "downstream lncRNA, starts ~620 bp 3' of PITX2"),
    ("LOC101795212", 4_003_281, 4_037_227, "-", "downstream protein-coding (-)"),
]


def annotate_gene_structure(peak_start: int, peak_end: int) -> str:
    """One-line structural context for a peak interval.

    Reports whether the peak overlaps PITX2 gene body (and where —
    promoter, exonic, intronic), and any neighbor lncRNA / coding-gene
    overlaps.  Coordinates from the locus features file.
    """
    tags = []
    overlaps_pitx2 = peak_start < PITX2_GENE_END and peak_end > PITX2_GENE_START
    if overlaps_pitx2:
        # Promoter calls: TSS ± 300 bp matches the promoter|... entries
        # in references/pitx2_promoter_regions.bed.
        for iso, tss in PITX2_TSS.items():
            if peak_start <= tss + 300 and peak_end >= tss - 300:
                tags.append(f"overlaps PITX2 {iso} promoter (TSS±300 bp)")
                break
        else:
            # Check if any PITX2 exon overlaps; otherwise intronic.
            exon_hit = False
            for ex_s, ex_e in PITX2_EXONS:
                if peak_start < ex_e and peak_end > ex_s:
                    tags.append(f"overlaps PITX2 exon {ex_s}-{ex_e}")
                    exon_hit = True
                    break
            if not exon_hit:
                tags.append("intronic to PITX2")
    else:
        # Upstream or downstream of PITX2 gene body.
        if peak_end <= PITX2_GENE_START:
            distance = PITX2_GENE_START - peak_end
            tags.append(f"5' upstream of PITX2 gene body ({distance:,} bp)")
        else:
            distance = peak_start - PITX2_GENE_END
            tags.append(f"3' downstream of PITX2 gene body ({distance:,} bp)")

    # Always report any neighbor overlap (LOC genes), even if the peak
    # also overlaps PITX2 (since some lncRNAs are antisense within
    # PITX2's exons/introns).
    for name, s, e, strand, role in NEIGHBOR_GENES:
        if peak_start < e and peak_end > s:
            tags.append(f"overlaps {name} ({strand}, {role})")

    return "; ".join(tags) if tags else "intergenic (no PITX2 / LOC overlap)"


def main() -> None:
    df = pd.read_csv(IN_TSV, sep="\t")
    n = len(df)
    print(f"[i] loaded {n} peaks from {IN_TSV.relative_to(REPO)}")

    # ---- data-integrity assertions --------------------------------------
    assert n == 329, f"expected 329 peaks at PITX2 ±500 kb, got {n}"
    assert df["binom_p"].isna().sum() == 0, "binom_p has NaN"
    assert df["binom_padj"].isna().sum() == 0, "binom_padj has NaN"
    assert df["log2FC"].isna().sum() == 0, "log2FC has NaN"
    assert ((df["binom_p"] >= 0) & (df["binom_p"] <= 1)).all(), "binom_p out of [0,1]"
    assert ((df["binom_padj"] >= 0) & (df["binom_padj"] <= 1)).all(), \
        "binom_padj out of [0,1]"
    assert (df["chrom"] == PITX2_CHROM).all(), \
        f"expected all peaks on {PITX2_CHROM}"

    # ---- locus-restricted BH --------------------------------------------
    _, padj_locus, _, _ = multipletests(df["binom_p"].values, method="fdr_bh")
    df["binom_padj_locus_fdr"] = padj_locus

    # Sanity: locus padj <= genome-wide padj on every peak (since the
    # multiple-testing burden is smaller).  This must hold for BH on a
    # nested set.  Allow a tiny float tolerance.  (Equality occurs when
    # binom_p is identical and both BH paths happen to pick the same
    # max-rank step-up.)
    delta = df["binom_padj"] - df["binom_padj_locus_fdr"]
    violations = (delta < -1e-9).sum()
    if violations:
        print(f"[!] WARNING: {violations} peaks have locus_fdr > gw_fdr "
              "(unexpected — BH on a subset of tests cannot increase padj). "
              "Check that the locus_subset binom_p values exactly match the "
              "GW-table binom_p column.")
    else:
        print(f"[i] all {n} locus-FDR padj <= genome-wide padj (expected; "
              f"max reduction = {delta.max():.3g}, mean = {delta.mean():.3g})")

    # ---- emerging peaks: locus padj<0.1 but NOT GW padj<0.1 -------------
    locus_sig_l = df.query("log2FC > 0 and binom_padj_locus_fdr < 0.1")
    gw_sig_l = df.query("log2FC > 0 and binom_padj < 0.1")
    n_locus_l = len(locus_sig_l)
    n_gw_l = len(gw_sig_l)
    print()
    print(f"[i] L-biased peaks clearing padj<0.1 at LOCUS-FDR: {n_locus_l}")
    print(f"[i] L-biased peaks clearing padj<0.1 at GENOME-WIDE FDR: {n_gw_l}")
    print(f"[i] newly emerging (locus<0.1 & gw>=0.1): {n_locus_l - n_gw_l}")

    # Pretty-print every L-biased peak that clears locus padj<0.1, ordered
    # by locus padj ascending.
    print()
    print("=" * 76)
    print(f"All L-biased peaks at PITX2 ±500 kb with binom_padj_locus_fdr < 0.1")
    print("=" * 76)
    cols_view = [
        "peak_id", "start", "end", "dist_genebody",
        "log2FC", "binom_p", "binom_padj", "binom_padj_locus_fdr",
        "n_pitx_motifs", "n_smad23_motifs",
    ]
    for _, r in locus_sig_l.sort_values("binom_padj_locus_fdr").iterrows():
        emerged = " <-- NEW (gw padj>=0.1)" if r["binom_padj"] >= 0.1 else ""
        print(f"  {r['peak_id']:38s}")
        print(f"    log2FC={r['log2FC']:+.3f}  binom_p={r['binom_p']:.4g}  "
              f"gw_padj={r['binom_padj']:.4g}  "
              f"LOCUS_padj={r['binom_padj_locus_fdr']:.4g}{emerged}")
        print(f"    dist_to_PITX2_gene_body = {int(r['dist_genebody']):+,} bp  "
              f"PITX-motifs={int(r['n_pitx_motifs'])}  "
              f"SMAD2/3-motifs={int(r['n_smad23_motifs'])}")
        print(f"    structure: {annotate_gene_structure(int(r['start']), int(r['end']))}")
        print()

    # ASE-region peak (3,974,421 - 3,974,922) — explicit verdict.
    print("=" * 76)
    print("ASE-region peak verdict (duck syntenic to chick PITX2 ASE region)")
    print("=" * 76)
    ase = df[(df["start"] == 3_974_421) & (df["end"] == 3_974_922)]
    if ase.empty:
        print("[!] ASE peak NC_051775.1:3,974,421-3,974,922 NOT found in input.")
    else:
        r = ase.iloc[0]
        verdict_locus = "CLEARS" if r["binom_padj_locus_fdr"] < 0.1 else "DOES NOT CLEAR"
        verdict_gw = "CLEARS" if r["binom_padj"] < 0.1 else "DOES NOT CLEAR"
        print(f"  peak_id: {r['peak_id']}")
        print(f"  log2FC = {r['log2FC']:+.4f} (L-biased)")
        print(f"  binom_p (raw) = {r['binom_p']:.4g}")
        print(f"  genome-wide BH padj = {r['binom_padj']:.4g} -> {verdict_gw} padj<0.1")
        print(f"  locus-FDR BH padj   = {r['binom_padj_locus_fdr']:.4g} "
              f"-> {verdict_locus} padj<0.1")
        print(f"  dist to PITX2 gene body: {int(r['dist_genebody']):+,} bp")
        print(f"  structure: {annotate_gene_structure(int(r['start']), int(r['end']))}")
        print(f"  PITX motifs = {int(r['n_pitx_motifs'])}, "
              f"SMAD2/3 motifs = {int(r['n_smad23_motifs'])}")

    # ---- write output ----------------------------------------------------
    out = df.sort_values("binom_padj_locus_fdr").reset_index(drop=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print()
    print(f"[i] wrote {OUT_TSV.relative_to(REPO)} ({len(out)} rows, "
          f"{len(out.columns)} cols incl. new binom_padj_locus_fdr)")


if __name__ == "__main__":
    main()
