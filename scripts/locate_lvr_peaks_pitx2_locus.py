#!/usr/bin/env python3
"""
Localize the HH17/16 INFERENCE_POOL ArchR aggregated L>R peaks to the duck
PITX2 locus, and report PITX-family + SMAD2/3 motif content for each one.

Input
-----
- figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv
  (g2c aggregated descriptive — peak_id, log2FC, binom_padj, etc.)
- figures/data/g4_hh1716_peak_motif_matrix_archr.tsv.gz
  (binary peak × 21-motif matrix from motifmatchr, p<5e-4)
- references/pitx2_promoter_regions.bed
  (PITX2 gene-body + per-isoform TSS coordinates, NC_051775.1, + strand)

Output
------
- figures/data/h0_lvr_peaks_pitx2_locus.tsv
  per-peak distance + motif counts within PITX2 ±500 kb
- stdout summary with bin counts at ±50/±100/±500 kb under two padj cutoffs.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
AGG_TSV = REPO / "figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv"
MOTIF_TSV = REPO / "figures/data/g4_hh1716_peak_motif_matrix_archr.tsv.gz"
OUT_TSV = REPO / "figures/data/h0_lvr_peaks_pitx2_locus.tsv"

PITX2_CHROM = "NC_051775.1"
PITX2_GENE_START = 3_961_386  # BED 0-based start of gene|PITX2
PITX2_GENE_END = 3_976_254    # BED end of gene|PITX2 (exclusive)
TSS = {
    "PITX2a": 3_961_387,
    "PITX2b": 3_964_253,
    "PITX2c": 3_971_385,
}

PITX_MOTIF_COLS = ["PITX1", "PITX2", "PITX3"]
SMAD23_MOTIF_COLS = ["SMAD2_3_4_human", "SMAD3", "Smad2_3_mouse"]
# SMAD2_3_4_human and Smad2_3_mouse are NODAL-pathway transducer composites;
# SMAD3 alone is the lone-monomer scan.


def parse_peak_id(peak_id: str) -> tuple[str, int, int]:
    chrom, coords = peak_id.split(":")
    start_s, end_s = coords.split("-")
    return chrom, int(start_s), int(end_s)


def signed_distance(start: int, end: int) -> int:
    """Signed distance from peak interval to PITX2 gene body.

    0  = peak overlaps gene body
    +  = peak is downstream of gene end (3' of gene on + strand)
    -  = peak is upstream of gene start (5' of gene)
    """
    if end <= PITX2_GENE_START:
        return start - PITX2_GENE_START  # negative (upstream)
    if start >= PITX2_GENE_END:
        return start - PITX2_GENE_END + 1  # positive (downstream)
    return 0  # overlap


def main() -> None:
    agg = pd.read_csv(AGG_TSV, sep="\t")
    print(f"[i] aggregated peaks total: {len(agg):,}")

    # L>R sets at two cutoffs
    headline = agg.query("log2FC > 0 and binom_padj < 0.1")
    strict = agg.query("log2FC > 0 and binom_padj < 0.05")
    print(f"[i] L>R peaks padj<0.1: {len(headline)}  | padj<0.05: {len(strict)}")

    # Parse coords once on the full aggregate; reused for both sets.
    coords = agg["peak_id"].str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    agg["chrom"] = coords["chrom"]
    agg["start"] = coords["start"].astype(int)
    agg["end"] = coords["end"].astype(int)

    # Restrict to PITX2 chromosome and ±500 kb of gene body for the
    # locus-level report.
    on_chrom = agg[agg["chrom"] == PITX2_CHROM].copy()
    on_chrom["dist_genebody"] = on_chrom.apply(
        lambda r: signed_distance(r["start"], r["end"]), axis=1
    )
    on_chrom["abs_dist_genebody"] = on_chrom["dist_genebody"].abs()
    # Distance from each isoform TSS (signed: peak midpoint - TSS).
    on_chrom["peak_mid"] = (on_chrom["start"] + on_chrom["end"]) // 2
    for iso, tss in TSS.items():
        on_chrom[f"dist_{iso}_TSS"] = on_chrom["peak_mid"] - tss
    locus500 = on_chrom[on_chrom["abs_dist_genebody"] <= 500_000].copy()
    print(
        f"[i] peaks on {PITX2_CHROM} within PITX2 ±500 kb (gene body): "
        f"{len(locus500):,} (any direction, any sig)"
    )

    # Motif matrix — load only rows we care about to save memory.
    keep_ids = set(locus500["peak_id"])
    motif_rows = []
    with gzip.open(MOTIF_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        peak_idx = header.index("peak_id")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[peak_idx] in keep_ids:
                motif_rows.append(parts)
    motif_df = pd.DataFrame(motif_rows, columns=header)
    for c in header[1:]:
        motif_df[c] = motif_df[c].astype(int)
    print(f"[i] motif matrix rows joined: {len(motif_df)}")

    locus500 = locus500.merge(motif_df, on="peak_id", how="left")
    locus500["n_pitx_motifs"] = locus500[PITX_MOTIF_COLS].sum(axis=1)
    locus500["n_smad23_motifs"] = locus500[SMAD23_MOTIF_COLS].sum(axis=1)
    locus500["any_pitx"] = (locus500["n_pitx_motifs"] > 0).astype(int)
    locus500["any_smad23"] = (locus500["n_smad23_motifs"] > 0).astype(int)
    locus500["any_pitx_or_smad23"] = (
        (locus500["any_pitx"] + locus500["any_smad23"]) > 0
    ).astype(int)

    # Flag which L>R sig set each peak falls into.
    locus500["sig_padj_0p1"] = (
        (locus500["log2FC"] > 0) & (locus500["binom_padj"] < 0.10)
    ).astype(int)
    locus500["sig_padj_0p05"] = (
        (locus500["log2FC"] > 0) & (locus500["binom_padj"] < 0.05)
    ).astype(int)

    # Sort by distance.
    locus500 = locus500.sort_values(["abs_dist_genebody", "binom_padj"])

    cols = (
        ["peak_id", "chrom", "start", "end",
         "dist_genebody", "abs_dist_genebody",
         "dist_PITX2a_TSS", "dist_PITX2b_TSS", "dist_PITX2c_TSS",
         "L_count", "R_count", "L_cpm", "R_cpm",
         "log2FC", "binom_p", "binom_padj",
         "sig_padj_0p1", "sig_padj_0p05",
         "n_pitx_motifs", "n_smad23_motifs",
         "any_pitx", "any_smad23", "any_pitx_or_smad23"]
        + PITX_MOTIF_COLS + SMAD23_MOTIF_COLS
    )
    locus500[cols].to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[i] wrote {OUT_TSV.relative_to(REPO)}")

    # ---- bin counts -----------------------------------------------------
    def bin_summary(df: pd.DataFrame, label: str) -> None:
        print(f"\n== {label} ==")
        print(f"  N peaks within PITX2 ±500 kb on {PITX2_CHROM}: {len(df)}")
        for window in (50_000, 100_000, 500_000):
            sub = df[df["abs_dist_genebody"] <= window]
            n = len(sub)
            n_pitx = int(sub["any_pitx"].sum()) if n else 0
            n_smad = int(sub["any_smad23"].sum()) if n else 0
            n_any = int(sub["any_pitx_or_smad23"].sum()) if n else 0
            print(
                f"  ±{window//1000:>3d} kb: n={n:>3d}  "
                f"PITX-motif+={n_pitx:>2d}  SMAD2/3-motif+={n_smad:>2d}  "
                f"(any PITX or SMAD2/3)={n_any:>2d}"
            )
        if len(df):
            print("\n  --- per-peak (±500 kb, sorted by distance) ---")
            for _, r in df.iterrows():
                tag = ""
                if r["sig_padj_0p05"]:
                    tag = "**"
                elif r["sig_padj_0p1"]:
                    tag = "*"
                print(
                    f"    {r['peak_id']:38s} "
                    f"dist={int(r['dist_genebody']):+9d}  "
                    f"log2FC={r['log2FC']:+.2f}  "
                    f"padj={r['binom_padj']:.3g}  "
                    f"PITX={int(r['n_pitx_motifs'])} SMAD23={int(r['n_smad23_motifs'])}"
                    f"  {tag}"
                )

    # Headline = the 99 sig L>R peaks at padj<0.1
    bin_summary(
        locus500[locus500["sig_padj_0p1"] == 1],
        "HEADLINE: L>R peaks padj<0.1 (n=99 total genome-wide)"
    )
    bin_summary(
        locus500[locus500["sig_padj_0p05"] == 1],
        "STRICT: L>R peaks padj<0.05 (n=54 total genome-wide)"
    )

    # Also report ALL peaks within ±500 kb (regardless of sig) so we can see
    # the baseline locus accessibility — this is the denominator the locus
    # asymmetry has to clear.
    bin_summary(
        locus500,
        "REFERENCE (all directions, any sig): every peak within ±500 kb"
    )


if __name__ == "__main__":
    main()
