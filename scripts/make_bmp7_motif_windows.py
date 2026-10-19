#!/usr/bin/env python3
"""Build the duck BMP7-locus window BED for the g3h chick lift + PITX-motif scan.

Answers the question Table S3 raises but does not test: the duck BMP7 locus
carries three left-open distal elements with PITX2-family motifs — do their
CHICK orthologs carry the motifs too, and does the chick BMP7 promoter?

Window set (name = <label>|<class>, 0-based, matching g4f_domain_motif_windows.bed):
  Lopen_motif  the 3 genome-wide-significant Left-biased, PITX2-motif+ peaks
               (= the Table S3 claim; padj 0.0068 / 0.013 / 0.017)
  bg           background peaks from the same +-100 kb window, sampled evenly
               over the locus and stratified by motif status, so the chick
               scan gets its own background rate rather than borrowing duck's
               (duck background in this window is 38/119 = 31.9% PITX2-motif+,
               so single-window presence/absence is uninformative without it)
The BMP7 promoter window is appended by the runner, which reads strand + TSS
from the duck GTF (not available locally).

Run: /opt/miniconda3/envs/duckfig/bin/python scripts/make_bmp7_motif_windows.py
Writes: figures/data/g3h_bmp7_motif_windows.bed
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
SRC = os.path.join(REPO, "figures", "data", "m5_mitogen_loci_pitx2_motif_peaks.tsv")
OUT = os.path.join(REPO, "figures", "data", "g3h_bmp7_motif_windows.bed")
N_BG = 16


def main():
    d = pd.read_csv(SRC, sep="\t")
    b = d[d.gene == "BMP7"].copy()
    if b.empty:
        raise SystemExit("no BMP7 rows in the mitogen motif table")
    pk = b.peak_id.str.extract(r"^(\S+):(\d+)-(\d+)$")
    b["chrom"], b["lo"], b["hi"] = pk[0], pk[1].astype(int), pk[2].astype(int)

    claim = b[(b.padj < 0.05) & (b.log2FC_LvsR > 0)].sort_values("padj")
    if len(claim) != 3 or int(claim.PITX2_motif.sum()) != 3:
        raise SystemExit(f"expected the 3 motif+ gw-L peaks, got {len(claim)} "
                         f"({int(claim.PITX2_motif.sum())} motif+) — Table S3 drifted")

    rest = b[~b.peak_id.isin(claim.peak_id)].sort_values("lo")
    # stratify by motif status so the chick background is not all-positive or
    # all-negative, then spread each stratum evenly across the locus
    bg = []
    for motif_status, share in ((1, N_BG // 2), (0, N_BG - N_BG // 2)):
        pool = rest[rest.PITX2_motif == motif_status]
        if len(pool) <= share:
            bg.append(pool)
        else:
            step = len(pool) / share
            bg.append(pool.iloc[[int(i * step) for i in range(share)]])
    bg = pd.concat(bg).sort_values("lo")

    rows = []
    for _, r in claim.iterrows():
        rows.append((r.chrom, r.lo, r.hi,
                     f"BMP7_L{int(r.dist_bp)//1000}kb|Lopen_motif"))
    # Labels MUST be unique: g3g keys its scan matrices by window label, so a
    # collision silently double-counts one window and drops its twin. dist_bp is
    # 0 for every peak inside the gene body, so it cannot disambiguate alone —
    # carry an index.
    for i, (_, r) in enumerate(bg.iterrows()):
        rows.append((r.chrom, r.lo, r.hi,
                     f"BMP7_bg{i:02d}_{int(r.dist_bp)//1000}kb_m{int(r.PITX2_motif)}|bg"))
    labels = [n.split("|")[0] for _, _, _, n in rows]
    if len(set(labels)) != len(labels):
        raise SystemExit(f"duplicate window labels: {sorted(labels)}")

    with open(OUT, "w") as fh:
        fh.write("# g3h duck BMP7-locus windows for the chick lift + PITX-motif scan.\n")
        fh.write("# name=<label>|<class>. 0-based. See scripts/make_bmp7_motif_windows.py.\n")
        fh.write(f"# duck background in this window: {int(b.PITX2_motif.sum())}/{len(b)} "
                 f"= {b.PITX2_motif.mean():.3f} PITX2-motif+ per 501-bp peak.\n")
        fh.write("# The promoter window is appended by pipeline/run_g3h_bmp7_chick.sh (needs the GTF).\n")
        for chrom, lo, hi, name in rows:
            fh.write(f"{chrom}\t{lo}\t{hi}\t{name}\t0\t+\n")
    print(f"wrote {OUT} ({len(rows)} windows: {len(claim)} claim + {len(bg)} background)")
    print(f"duck background PITX2-motif rate in the BMP7 window: "
          f"{b.PITX2_motif.mean():.3f} ({int(b.PITX2_motif.sum())}/{len(b)})")


if __name__ == "__main__":
    main()
