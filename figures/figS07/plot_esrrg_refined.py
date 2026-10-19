#!/usr/bin/env python3
"""plot_esrrg_refined.py — ESRRG intron-inflation mitigation (refined; replaces the log bars).

Two stacked panels:
  TOP  — per-embryo intron fraction for the 6 Fig-6a estrogen receptors. One dot per embryo
         per gene (pseudobulk, summed across library), x = % intronic UMIs (linear, no log),
         dot area ∝ total UMIs, ESRRG highlighted. Pooled %-intronic shown as a black tick;
         the pooled exonic anchor (UMIs) annotated at right. ESR1 = big exonic anchor at a
         moderate intron fraction; ESRRG = the intron-inflated outlier (few exonic UMIs).
  BOTTOM — ESRRG gene model (16 exons, − strand, ~404 kb) with the major scRNA coverage peaks
         dotted as priming sites, coloured intronic vs exonic. The peaks sit in introns, not at
         the 3′ terminal exon → internal priming on pre-mRNA explains the broad intronic signal.

B&W + a single accent for ESRRG (figure-style: minimise colour); y/x start at 0 / locus start.

Inputs (shared old-duck data home figures/data/scrna_e6e9e11/):
  er_perembryo_ie.tsv      gene, library, embryo, umi_exon, umi_intron
  er_intron_exon_umi.tsv   per-gene pooled exon/intron UMIs, pct_intronic, span_kb
  esrrg_coverage_binned.tsv  pos, depth   (pooled cr6-BAM coverage, 200 bp bins)
  esrrg_exons.bed          chrom, start, end  (16 exons)
Output: panel_esrrg_refined.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, os.pardir, "data", "scrna_e6e9e11")   # shared old-duck data home
OUT = os.path.join(HERE, "panel_esrrg_refined.png")

ACC = "#B2182B"        # ESRRG accent (the flagged gene, top panel)
GREY = "#7A7A7A"       # other genes
# ESRRG locus (NC_051774.1, − strand): 3′ terminal exon at the LOW coordinate
GS, GE = 20177932, 20581985
MB = 1e6


def per_embryo():
    df = pd.read_csv(os.path.join(D, "er_perembryo_ie.tsv"), sep="\t")
    g = (df.groupby(["gene", "embryo"], as_index=False)[["umi_exon", "umi_intron"]].sum())
    g["total"] = g["umi_exon"] + g["umi_intron"]
    g = g[g["total"] > 0].copy()
    g["pct_intronic"] = 100.0 * g["umi_intron"] / g["total"]
    return g


def pooled():
    p = pd.read_csv(os.path.join(D, "er_intron_exon_umi.tsv"), sep="\t")
    return p.set_index("gene")


def find_peaks(pos, depth, pct=85):
    """All coverage-peak (priming) sites: bins above the pct-th percentile, grouped within 2 kb,
    one peak (group-max) per group. Returns genomic-position-sorted bin indices."""
    thr = np.percentile(depth, pct)
    idx = np.where(depth >= thr)[0]
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 10:          # 10 bins * 200 bp = 2 kb
            cur.append(i)
        else:
            groups.append(cur); cur = [i]
    groups.append(cur)
    return sorted(g[int(np.argmax(depth[g]))] for g in groups)


def main():
    pe = per_embryo()
    pl = pooled()
    genes = pe.groupby("gene")["pct_intronic"].median().sort_values(ascending=False).index.tolist()
    n = len(genes)
    maxtot = pe["total"].max()

    def dot_size(tot):
        return 12 + 300.0 * (tot / maxtot)        # area ∝ total UMIs

    fig = plt.figure(figsize=(11.0, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.42)
    axT = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    # ── TOP: per-embryo % intronic dots ──
    for rank, gene in enumerate(genes):
        y = n - 1 - rank                              # highest-pct gene at top
        is_e = (gene == "ESRRG")
        col = ACC if is_e else GREY
        if is_e:
            axT.axhspan(y - 0.45, y + 0.45, color=ACC, alpha=0.07, zorder=0)
        sub = pe[pe["gene"] == gene]
        rng = np.random.default_rng(rank)            # deterministic jitter
        jit = (rng.random(len(sub)) - 0.5) * 0.5
        axT.scatter(sub["pct_intronic"], y + jit, s=[dot_size(t) for t in sub["total"]],
                    facecolor=col, edgecolor="white", linewidth=0.4, alpha=0.85, zorder=3)
        # pooled %-intronic reference tick (black)
        pooled_pct = float(pl.loc[gene, "pct_intronic"])
        axT.plot([pooled_pct, pooled_pct], [y - 0.42, y + 0.42], color="black", lw=1.4, zorder=4)
    axT.set_yticks(range(n))
    axT.set_yticklabels([genes[n - 1 - y] for y in range(n)], fontsize=9)
    for tl, gene in zip(axT.get_yticklabels(), [genes[n - 1 - y] for y in range(n)]):
        if gene == "ESRRG":
            tl.set_color(ACC); tl.set_fontweight("bold")
    axT.set_xlim(-7, 107); axT.set_ylim(-0.6, n - 0.4)   # margin so 0%/100% embryo dots aren't clipped
    axT.set_xticks([0, 25, 50, 75, 100])
    axT.set_xlabel("% intronic UMIs")
    for sp in ("top", "right"):
        axT.spines[sp].set_visible(False)
    # size legend
    leg_tot = [t for t in (100, 1000, 4000) if t <= maxtot]
    handles = [axT.scatter([], [], s=dot_size(t), facecolor=GREY, edgecolor="white",
                           linewidth=0.4, label=f"{t:,}") for t in leg_tot]   # axT (not plt) — keep off axB
    axT.legend(handles=handles, title="total UMIs", loc="upper left", frameon=False,
               fontsize=7, title_fontsize=7.5, labelspacing=1.3, handletextpad=1.0)

    # ── BOTTOM: ESRRG gene model + priming-site dots ──
    cov = pd.read_csv(os.path.join(D, "esrrg_coverage_binned.tsv"), sep="\t",
                      header=None, names=["pos", "depth"])
    ex = pd.read_csv(os.path.join(D, "esrrg_exons.bed"), sep="\t",
                     header=None, names=["chrom", "start", "end"])
    pos = cov["pos"].to_numpy(); depth = cov["depth"].to_numpy(); dmax = depth.max()
    pk = find_peaks(pos, depth, pct=85)        # coverage-peak (priming) sites; density across the gene body

    GENE_Y = -0.12
    ARROW_Y = 0.98                             # row of internal-priming-site arrowheads (fig02D style)
    # the UMI pileup (coverage) — light backdrop, scaled below the arrowhead row
    axB.fill_between(pos / MB, depth / dmax * 0.82, step="mid", color="0.88", linewidth=0, zorder=0)
    # internal-priming-site markers: downward arrowheads (▼) at each peak — cf. fig02 panel D.
    # Their horizontal density shows where internal priming concentrates (across the introns).
    axB.scatter([pos[i] / MB for i in pk], [ARROW_Y] * len(pk), marker="v",
                s=10, facecolor="0.45", edgecolor="none", zorder=4, clip_on=False)
    # gene model (16 exons, − strand; transcription high→low, 3′ polyA at the LOW coordinate)
    axB.hlines(GENE_Y, GS / MB, GE / MB, color="black", lw=1.1, zorder=2)
    for _, r in ex.iterrows():
        axB.add_patch(Rectangle((r["start"] / MB, GENE_Y - 0.045),
                                max((r["end"] - r["start"]) / MB, 0.0015), 0.09,
                                facecolor="black", edgecolor="none", zorder=3))
    axB.annotate("", xy=(GE / MB - 0.012, GENE_Y), xytext=(GE / MB, GENE_Y),
                 arrowprops=dict(arrowstyle="-|>", color="black", lw=1.0))
    axB.text(GS / MB, GENE_Y - 0.12, "3′", fontsize=8, ha="left", va="top", color="0.3")
    axB.text(GE / MB, GENE_Y - 0.12, "5′", fontsize=8, ha="right", va="top", color="0.3")
    axB.text((GS + GE) / 2 / MB, 1.14, "ESRRG", fontsize=9.5, ha="center", va="bottom", fontstyle="italic")
    axB.set_xlabel("NC_051774.1 position (Mb)")
    axB.set_xlim((GS - 8000) / MB, (GE + 8000) / MB); axB.set_ylim(GENE_Y - 0.26, 1.22)
    axB.set_yticks([])
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT, dpi=180, bbox_inches="tight")
    print(f"[esrrg_refined] wrote {OUT}")
    print(f"  genes (by median %intronic): {genes}")
    print(f"  embryos/gene: {pe.groupby('gene')['embryo'].nunique().to_dict()}")
    print(f"  priming-site peaks dotted: {len(pk)}")


if __name__ == "__main__":
    main()
