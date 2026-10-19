"""Descriptive check: do peaks near the asymmetric mitogen ligands carry PITX2
motifs, and are any of those peaks Left-biased at HH25? Includes a
discriminatory panel: same-class expressed-but-null controls (WNT11, FRZB, DKK3,
RSPO3), a same-class new positive (BMP5), and the genome-wide Right-biased
control (added 2026-08-23).

Substrate (all committed, local):
  h6_paired_deseq2_HH25_balanced.tsv          genome-wide paired L/R DESeq2 (+log2FC = Left)
  g4_hh1716_peak_motif_matrix_archr.tsv.gz    ArchR motif matches, full 283,818-peak union set
  m5_mitogen_gene_coords_zju1.tsv             NCBI gene locationhist for GCF_015476345.1
                                              (fetched 2026-08-23; role column = panel arm)

Sign-convention anchor: the PITX2c TSS peak reproduces ~+1.33 at HH25
(Peak #1; re-verified at code level 2026-08-23).
`dist_bp` in the output = distance from peak midpoint to the NEAREST GENE-BODY
EDGE (not the TSS — strand-aware TSS distances are not computed here).
Full per-peak export: EVERY union peak in each window (no preview truncation).
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
WIN = 100_000

h6 = pd.read_csv(DATA / "h6_paired_deseq2_HH25_balanced.tsv", sep="\t")
mm = pd.read_csv(DATA / "g4_hh1716_peak_motif_matrix_archr.tsv.gz", sep="\t",
                 usecols=["peak_id", "PITX1", "PITX2", "PITX3"])
if len(mm) < 280_000:
    raise SystemExit(f"motif matrix unexpectedly small: {len(mm)}")

ANCHOR_ID = "NC_051775.1:3971182-3971683"
best = h6[h6.peak_id == ANCHOR_ID]
if len(best) != 1:
    raise SystemExit(f"anchor peak {ANCHOR_ID} not found in h6 table")
best = best.iloc[0]
print(f"sign anchor (PITX2c TSS peak): {best.peak_id} "
      f"log2FC={best.log2FoldChange:+.2f} padj={best.padj:.2g}")
if not (0.8 < best.log2FoldChange < 2.2):
    raise SystemExit("sign-convention anchor failed — do not trust directions")

d = h6.merge(mm, on="peak_id", how="left")
if d.PITX2.isna().any():
    raise SystemExit(f"{int(d.PITX2.isna().sum())} h6 peaks missing from motif matrix")

pk = d.peak_id.str.extract(r"^(\S+):(\d+)-(\d+)$")
d["chrom"], d["p_lo"], d["p_hi"] = pk[0], pk[1].astype(int), pk[2].astype(int)

bg = d.PITX2.mean()
gw_L = d[(d.padj < 0.05) & (d.log2FoldChange > 0)]
gw_R = d[(d.padj < 0.05) & (d.log2FoldChange < 0)]
print(f"\nbackground: PITX2-motif fraction over all {len(d):,} tested peaks = {bg:.3f}")
print(f"Left-biased  (padj<0.05): {len(gw_L):>3} | PITX2-motif+: {int(gw_L.PITX2.sum())} "
      f"({gw_L.PITX2.mean():.3f}, {gw_L.PITX2.mean()/bg:.2f}x)")
print(f"Right-biased (padj<0.05): {len(gw_R):>3} | PITX2-motif+: {int(gw_R.PITX2.sum())} "
      f"({gw_R.PITX2.mean():.3f}, {gw_R.PITX2.mean()/bg:.2f}x)  <- specificity control")

genes = pd.read_csv(DATA / "m5_mitogen_gene_coords_zju1.tsv", sep="\t")
if "role" not in genes.columns:
    raise SystemExit("coords TSV lacks the role column — refetch")

print(f"\nwindows: gene body +-{WIN//1000} kb; dist = midpoint to nearest gene-body edge\n")
rows = []
summary = []
for _, g in genes.iterrows():
    w = d[(d.chrom == g.chracc) & (d.p_hi > g.start - WIN) & (d.p_lo < g.end + WIN)].copy()
    nomL = w[(w.pvalue < 0.01) & (w.log2FoldChange > 0)]
    nomR = w[(w.pvalue < 0.01) & (w.log2FoldChange < 0)]
    gwLw = w[(w.padj < 0.05) & (w.log2FoldChange > 0)]
    summary.append(dict(gene=g.gene, role=g.role, n_peaks=len(w),
                        motif_rate=round(float(w.PITX2.mean()), 3),
                        gw_L=len(gwLw), gw_L_motif=int(gwLw.PITX2.sum()),
                        nom_L=len(nomL), nom_L_motif=int(nomL.PITX2.sum()),
                        nom_R=len(nomR)))
    print(f"== {g.gene} [{g.role}]  ({g.chracc}:{g.start:,}-{g.end:,})")
    print(f"   peaks {len(w)} | motif rate {w.PITX2.mean():.2f} | gw-L {len(gwLw)} "
          f"(motif+ {int(gwLw.PITX2.sum())}) | nom-L {len(nomL)} | nom-R {len(nomR)}")
    for _, r in w.sort_values("pvalue").iterrows():
        mid = (r.p_lo + r.p_hi) // 2
        dist = 0 if g.start <= mid <= g.end else min(abs(mid - g.start), abs(mid - g.end))
        rows.append(dict(gene=g.gene, role=g.role, peak_id=r.peak_id, dist_bp=dist,
                         log2FC_LvsR=round(r.log2FoldChange, 3),
                         pvalue=r.pvalue, padj=r.padj,
                         PITX2_motif=int(r.PITX2), PITX1=int(r.PITX1),
                         PITX3=int(r.PITX3)))
    for _, r in w.sort_values("pvalue").head(5).iterrows():
        mid = (r.p_lo + r.p_hi) // 2
        dist = 0 if g.start <= mid <= g.end else min(abs(mid - g.start), abs(mid - g.end))
        flag = " *L(gw)*" if (r.padj < 0.05 and r.log2FoldChange > 0) else (
               " *L(nom)*" if (r.pvalue < 0.01 and r.log2FoldChange > 0) else "")
        print(f"     {r.peak_id}  d={dist/1000:.1f}kb  lfc={r.log2FoldChange:+.2f} "
              f"p={r.pvalue:.3g} padj={r.padj:.2g} PITX2motif={int(r.PITX2)}{flag}")
    print()

pd.DataFrame(rows).to_csv(DATA / "m5_mitogen_loci_pitx2_motif_peaks.tsv",
                          sep="\t", index=False)
print(pd.DataFrame(summary).to_string(index=False))
print(f"\nwrote {DATA/'m5_mitogen_loci_pitx2_motif_peaks.tsv'} "
      f"({len(rows)} rows — ALL window peaks, no truncation)")
