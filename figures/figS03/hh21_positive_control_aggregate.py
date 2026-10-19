"""Aggregate (set-level) positive-control statistics — stronger than per-promoter.
Uses only committed grn_memory_node_promoter_scan.tsv + h7_unpaired_deseq2_HH21.tsv.gz."""
import numpy as np, pandas as pd
from scipy import stats

import os
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
scan = pd.read_csv(f"{DATA}/grn_memory_node_promoter_scan.tsv", sep="\t")
h7 = pd.read_csv(f"{DATA}/h7_unpaired_deseq2_HH21.tsv.gz", sep="\t")[
    ["peak_id", "log2FoldChange", "padj", "mean_L", "mean_R", "baseMean"]]
d = scan.merge(h7, on="peak_id", how="left")
# Every PITX2-locus element out (is_ref covers only PITX2's own promoters; the +5.4 kb
# element is PITX2's cis-element too), matching panel H's independent-responder set.
d = d[d.log2FoldChange.notna() & ~d.pitx2_locus.astype(bool)].copy()
L = d[d.group == "L-up"]; R = d[d.group == "R-up"]
print(f"L-up promoter peaks: {len(L)}   R-up promoter peaks: {len(R)}")

# 1) COMPETITIVE aggregate: do L-up promoters have more-Left chromatin than R-up promoters?
u, up = stats.mannwhitneyu(L.log2FoldChange, R.log2FoldChange, alternative="greater")
t, tp = stats.ttest_ind(L.log2FoldChange, R.log2FoldChange, equal_var=False)
# rank-biserial effect size
rbc = 2 * u / (len(L) * len(R)) - 1
print(f"\n[1] Competitive L-up vs R-up chromatin log2FC:")
print(f"    L-up mean {L.log2FoldChange.mean():+.3f} vs R-up {R.log2FoldChange.mean():+.3f} "
      f"(diff {L.log2FoldChange.mean()-R.log2FoldChange.mean():+.3f})")
print(f"    Mann-Whitney one-sided p={up:.2e}, rank-biserial={rbc:+.3f}; Welch t p={tp:.2e}")

# 2) PEAK-PAIRED module accessibility: within each module, is L-side > R-side (paired over peaks)?
for name, s, side_hi in (("L-up", L, "L"), ("R-up", R, "R")):
    hi, lo = (s.mean_L, s.mean_R) if side_hi == "L" else (s.mean_R, s.mean_L)
    w, wp = stats.wilcoxon(hi, lo, alternative="greater")
    frac = (hi.values > lo.values).mean()
    print(f"\n[2] {name} module: mean access on its OWN side vs opposite (paired over {len(s)} peaks):")
    print(f"    {side_hi}-side>{('R' if side_hi=='L' else 'L')}-side in {frac*100:.0f}% of peaks; "
          f"Wilcoxon signed-rank p={wp:.2e}")

# 3) PERMUTATION null on the competitive mean-difference (shuffle L-up/R-up labels)
obs = L.log2FoldChange.mean() - R.log2FoldChange.mean()
vals = d.log2FoldChange.values; nL = len(L)
rng = np.random.default_rng(0)
perm = np.array([(lambda idx: vals[idx[:nL]].mean() - vals[idx[nL:]].mean())(rng.permutation(len(vals)))
                 for _ in range(10000)])
pp = (perm >= obs).mean()
print(f"\n[3] Permutation (10k) on L-up minus R-up mean chromatin log2FC:")
print(f"    observed {obs:+.3f}, null mean {perm.mean():+.3f} sd {perm.std():.3f}, empirical p={pp:.4f}")
