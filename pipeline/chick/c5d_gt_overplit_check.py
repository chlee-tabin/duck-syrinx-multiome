#!/usr/bin/env python3
"""c5d — over-split diagnostic: GT-distance at higher N.

At N=12 the 12 vireo donors are all GT-distinct (max pairwise r=0.47), matching the
12-embryo design. But the best seed at N=14 fills 13 donors and at N=16 fills 15, with
the substantial-donor count creeping up (9 -> 10 -> 11). Two explanations:
  (a) over-split: forcing more donors fragments a real embryo -> the fragment is
      GT-CLOSE (r>0.90) to its parent;
  (b) genuine: there are really >12 genotypes (contamination / design undercount) ->
      the extra donors are GT-distinct from all others.
This prints, per requested seed dir, the donors with >= MIN_CELLS cells and any pair with
r>0.90 (over-split signature). A close pair at N=14/16 confirms N=12 is the right cut.

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5d_gt_overplit_check.py \
         <seed_dir1> [<seed_dir2> ...]
"""
import os, sys, gzip
import numpy as np
import pandas as pd

MIN_CELLS = 300


def _open(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def parse_gt(vcf):
    donors, rows = [], []
    with _open(vcf) as f:
        for ln in f:
            if ln.startswith("##"): continue
            if ln.startswith("#CHROM"):
                donors = ln.rstrip("\n").split("\t")[9:]; continue
            x = ln.rstrip("\n").split("\t"); fmt = x[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            ds = []
            for s in x[9:]:
                g = s.split(":")[gi].replace("|", "/").split("/")
                try: ds.append(sum(int(a) for a in g))
                except ValueError: ds.append(np.nan)
            rows.append(ds)
    return donors, np.array(rows, float).T


def analyze(bdir):
    gt = next((f"{bdir}/{e}" for e in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf")
               if os.path.exists(f"{bdir}/{e}")), None)
    if not gt:
        print(f"[{bdir}] no GT_donors VCF — skipped"); return
    did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t")
    counts = did[~did.donor_id.isin(["doublet", "unassigned"])].donor_id.value_counts()
    donors, G = parse_gt(gt); n = len(donors)
    C = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            m = ~(np.isnan(G[i]) | np.isnan(G[j]))
            if m.sum() >= 100:
                C[i, j] = C[j, i] = np.corrcoef(G[i][m], G[j][m])[0, 1]
    sub = sorted([d for d in counts.index if counts[d] >= MIN_CELLS],
                 key=lambda d: -counts[d])
    print(f"\n=== {os.path.basename(os.path.dirname(bdir))}/{os.path.basename(bdir)} ===")
    print(f"  {len(sub)} substantial donors (>= {MIN_CELLS} cells): "
          + ", ".join(f"{d}({counts[d]})" for d in sub))
    di = {d: k for k, d in enumerate(donors)}
    pairs = []
    for a in range(len(sub)):
        for b in range(a + 1, len(sub)):
            r = C[di[sub[a]], di[sub[b]]]
            if not np.isnan(r):
                pairs.append((r, sub[a], sub[b]))
    pairs.sort(reverse=True)
    close = [p for p in pairs if p[0] > 0.90]
    print(f"  max pairwise r among substantial donors = {pairs[0][0]:.3f} "
          f"({pairs[0][1]}~{pairs[0][2]})" if pairs else "  (no pairs)")
    if close:
        print(f"  *** {len(close)} OVER-SPLIT pair(s) r>0.90 (confirms N too high): "
              + "; ".join(f"{a}~{b}:{r:.3f}" for r, a, b in close))
    else:
        print("  no r>0.90 pair — substantial donors all GT-distinct at this N")
    print("  top pairs: " + " | ".join(f"{a}~{b}:{r:.2f}" for r, a, b in pairs[:4]))


for d in sys.argv[1:]:
    analyze(d.rstrip("/"))
