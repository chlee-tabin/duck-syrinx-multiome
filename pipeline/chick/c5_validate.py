#!/usr/bin/env python3
"""c5 — deep per-donor validation (both groups), the duck-grade gate for clutch siblings.

Runs after c3/c4 (multi-seed vireo + seed-QC). Computes:
  1. GT-DISTANCE between donors (donor x donor genotype-dosage correlation from
     GT_donors.vireo.vcf.gz). Genetically-close clutch-mates show high correlation —
     flag pairs r>0.90 as potential over-split (should merge) or barely-separable
     siblings needing orthogonal corroboration. This is the signal that caught the
     duck's hard-to-discriminate pairs.
  2. ZW SEX per donor (chrW vs chrZ pseudobulk depth). Bird sex: female=ZW (W present),
     male=ZZ (W absent). A donor with intermediate/ambiguous W is a candidate M+F merger.
Sex + GT-distance + (Group B) GFP + cross-lane are the orthogonal tie-breakers that
resolve what genotype-only clustering cannot.

  conda run -n vireo python3 pipeline/chick/c5_validate.py <WORK_DIR> [expected_N]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux"
EXPECTED = int(sys.argv[2]) if len(sys.argv) > 2 else 12
SEEDED = f"{WORK}/vireo_seeded"; OUT = f"{SEEDED}/seed_qc"; CELL = f"{WORK}/cellsnp/combined"
os.makedirs(OUT, exist_ok=True)

rep = pd.read_csv(f"{OUT}/seed_qc_report.tsv", sep="\t")
pick = rep[(rep.N == EXPECTED) & (rep.n_donors == EXPECTED)]
if pick.empty: pick = rep[rep.n_donors == EXPECTED]
if pick.empty: pick = rep
best = pick.sort_values("logLik", ascending=False).iloc[0]
bdir = best["dir"]
print(f"Best seed for validation: N={best['N']} seed={best['seed']}\n  {bdir}")


def _open(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def parse_gt(vcf):
    donors, rows = [], []
    with _open(vcf) as f:
        for ln in f:
            if ln.startswith("##"): continue
            if ln.startswith("#CHROM"): donors = ln.rstrip("\n").split("\t")[9:]; continue
            x = ln.rstrip("\n").split("\t"); fmt = x[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            ds = []
            for s in x[9:]:
                g = s.split(":")[gi].replace("|", "/").split("/")
                try: ds.append(sum(int(a) for a in g))
                except ValueError: ds.append(np.nan)
            rows.append(ds)
    return donors, np.array(rows, float).T  # donors x variants

# ── 1. GT-distance / close-pair detection ──────────────────────────────────
gtvcf = next((f"{bdir}/{e}" for e in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf")
              if os.path.exists(f"{bdir}/{e}")), None)
if gtvcf:
    donors, G = parse_gt(gtvcf); n = len(donors)
    C = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i, n):
            m = ~(np.isnan(G[i]) | np.isnan(G[j]))
            if m.sum() >= 100:
                C[i, j] = C[j, i] = np.corrcoef(G[i][m], G[j][m])[0, 1]
    pd.DataFrame(C, index=donors, columns=donors).to_csv(f"{OUT}/gt_pairwise_corr.tsv", sep="\t")
    close = [(donors[i], donors[j], C[i, j]) for i in range(n) for j in range(i + 1, n)
             if not np.isnan(C[i, j]) and C[i, j] > 0.90]
    print(f"\n[GT-distance] {n} donors. Close pairs r>0.90 (clutch/over-split candidates):")
    print("\n".join(f"  {a} ~ {b}: r={r:.3f}" for a, b, r in sorted(close, key=lambda t: -t[2]))
          or "  none — donors well-separated")
else:
    print("\n[GT-distance] no GT_donors VCF found; skipped")

# ── 2. ZW sex per donor ─────────────────────────────────────────────────────
base = next((f"{CELL}/{e}" for e in ("cellSNP.base.vcf.gz", "cellSNP.base.vcf")
             if os.path.exists(f"{CELL}/{e}")), None)
chrom = np.array([ln.split("\t", 1)[0] for ln in _open(base) if not ln.startswith("#")])
DP = sio.mmread(f"{CELL}/cellSNP.tag.DP.mtx").tocsr()              # variants x cells
samples = [l.strip() for l in open(f"{CELL}/cellSNP.samples.tsv") if l.strip()]
did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t")
Wdp = np.asarray(DP[chrom == "chrW"].sum(0)).ravel()
Zdp = np.asarray(DP[chrom == "chrZ"].sum(0)).ravel()
col = {bc: i for i, bc in enumerate(samples)}
rows = []
for d in sorted(set(did.donor_id) - {"doublet", "unassigned"}):
    idx = [col[c] for c in did.loc[did.donor_id == d, "cell"] if c in col]
    w, z = Wdp[idx].sum(), Zdp[idx].sum()
    rows.append({"donor": d, "n_cells": len(idx), "chrW_DP": int(w), "chrZ_DP": int(z),
                 "W_over_Z": round(w / (z + 1e-9), 4)})
sdf = pd.DataFrame(rows).sort_values("W_over_Z")
sdf.to_csv(f"{OUT}/per_donor_sex.tsv", sep="\t", index=False)
print("\n[ZW sex] chrW/chrZ pseudobulk depth (high=female ZW, ~0=male ZZ); "
      "look for a clean bimodal gap — straddlers = candidate M+F mergers:")
print(sdf.to_string(index=False))
print(f"\nWrote {OUT}/gt_pairwise_corr.tsv, {OUT}/per_donor_sex.tsv")
