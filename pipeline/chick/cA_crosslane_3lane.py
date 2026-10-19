#!/usr/bin/env python3
"""cA 3-lane cross-lane concordance — Group A ARC01 / ARC02 / ARC0S donor genotype correlations.

Group A has THREE lanes that sample the SAME 12 embryos:
  - ARC01, ARC02 : complementary L/R halves, jointly demuxed (cellsnp/combined + N18/seed_117).
  - ARC0S        : the L-only HH15/HH25 "symmetry-breaking" lane (8 left half-embryos), genotyped at
                   the IDENTICAL 283,435-variant base and PROJECTED onto the Group A donors via a
                   reference-VCF vireo run (--donorFile).  arc0s/cellsnp + arc0s/vireo_donorfile/seed_117.

For a correctly-assigned donor, its pseudobulk allele frequencies agree across every lane it appears in
(same embryo genotype): ARC01~ARC02 (all 12) and ARC0S~ARC01, ARC0S~ARC02 (the 8 HH15/HH25 left-halves).
The 4 HH13 donors are absent from ARC0S by design -> a built-in negative control (no ARC0S left-half).

This supersedes the 2-lane cA_crosslane_matrix.py (ARC01 x ARC02 only; removed) by adding the ARC0S
lane, and writes a tidy long-form concordance table + per-pair wide matrices + a compact per-donor
diagonal table for the figS04_supp demux-concordance panel. The ARC01 x ARC02 block reproduces the
prior 2-lane result exactly (cross-check).

  sbatch -p priority pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cA_crosslane_3lane.py [WORK] [BDIR]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux"
# production Group-A seed (matches cT_groupA_stage_sex / cAE_full_atlas).
BDIR = sys.argv[2] if len(sys.argv) > 2 else f"{WORK}/vireo_seeded/N18/seed_117"
SEED = os.path.basename(BDIR.rstrip("/"))             # e.g. "seed_117" — the Group A donor run in use
COMB = f"{WORK}/cellsnp/combined"                      # ARC01+ARC02 merged cellData (suffix -1 / -2)
ARC0S_CELL = f"{WORK}/arc0s/cellsnp"                   # ARC0S genotyped at the SAME 283,435 base
# ARC0S --donorFile projection MUST be onto the SAME Group A donor run as BDIR: vireo donor labels are
# run-specific, so a mismatched run would correlate unrelated donor IDs. Derive the path from SEED (not
# a hardcoded seed_117) so overriding BDIR either finds the matching projection or fails loud (no file).
ARC0S_BDIR = f"{WORK}/arc0s/vireo_donorfile/{SEED}"
OUT = f"{WORK}/vireo_seeded/seed_qc"
MIN_DP = 15            # min depth (per lane) for a variant to enter a correlation
MIN_PAIR_VAR = 100     # min shared informative variants for a correlation
MIN_CELLS = 20         # min cells to form a per-donor pseudobulk AF vector


def load_mtx(cell_dir):
    AD = sio.mmread(f"{cell_dir}/cellSNP.tag.AD.mtx").tocsc()
    DP = sio.mmread(f"{cell_dir}/cellSNP.tag.DP.mtx").tocsc()
    samples = [l.strip() for l in open(f"{cell_dir}/cellSNP.samples.tsv") if l.strip()]
    return AD, DP, {bc: i for i, bc in enumerate(samples)}


def base_positions(cell_dir):
    """(CHROM, POS) per variant = the row identity of this dir's AD/DP matrices (cellSNP base order)."""
    import gzip
    out = []
    with gzip.open(f"{cell_dir}/cellSNP.base.vcf.gz", "rt") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                c = ln.split("\t", 2)
                out.append((c[0], c[1]))
    return out


AD_c, DP_c, col_c = load_mtx(COMB)
AD_s, DP_s, col_s = load_mtx(ARC0S_CELL)
assert AD_c.shape[0] == AD_s.shape[0], (
    f"variant-row mismatch combined({AD_c.shape[0]}) vs arc0s({AD_s.shape[0]}) — base VCFs differ, "
    f"AF vectors would be misaligned. Genotype ARC0S with -R cellsnp/combined/cellSNP.base.vcf.gz.")
# Equal row COUNT is necessary but not sufficient — the AD/DP rows must also be the SAME variants in the
# SAME order, else every cross-lane correlation is silently computed over mismatched positions. Enforce
# position identity in-code (replaces the prior out-of-band md5 check) so any future ARC0S re-genotyping
# against a different base fails loud here rather than producing wrong-but-plausible r values.
assert base_positions(COMB) == base_positions(ARC0S_CELL), (
    "combined and arc0s cellSNP base VCFs are NOT position-identical (same length, different variants/"
    "order) — AF rows misaligned. Genotype ARC0S with -R cellsnp/combined/cellSNP.base.vcf.gz.")
print(f"[load] {AD_c.shape[0]} variants (position-verified)  ·  combined {AD_c.shape[1]} cells  ·  "
      f"arc0s {AD_s.shape[1]} cells  ·  Group A run {SEED}", flush=True)

# combined (ARC01+ARC02): lane from the final barcode suffix, donor from the joint demux
did = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["lane"] = did.cell.str.rsplit("-", n=1).str[-1].map({"1": "ARC01", "2": "ARC02"})

# ARC0S: donor from the reference-VCF projection (single lane; drop doublet/unassigned = duck carry-over)
dis = pd.read_csv(f"{ARC0S_BDIR}/donor_ids.tsv", sep="\t")
n_proj_total = len(dis)
dis = dis[~dis.donor_id.isin(["doublet", "unassigned"])].copy()
print(f"[arc0s] {n_proj_total} projected cells -> {len(dis)} assigned singlets "
      f"(median prob_max {dis.prob_max.astype(float).median():.3f})", flush=True)


def af(cells, col, AD, DP):
    idx = [col[c] for c in cells if c in col]
    if len(idx) < MIN_CELLS:
        return None
    a = np.asarray(AD[:, idx].sum(1)).ravel()
    d = np.asarray(DP[:, idx].sum(1)).ravel()
    return a, d


def corr(p1, p2):
    if p1 is None or p2 is None:
        return np.nan, 0
    a1, d1 = p1
    a2, d2 = p2
    m = (d1 >= MIN_DP) & (d2 >= MIN_DP)
    if m.sum() < MIN_PAIR_VAR:
        return np.nan, int(m.sum())
    return float(np.corrcoef(a1[m] / d1[m], a2[m] / d2[m])[0, 1]), int(m.sum())


def dkey(x):
    s = "".join(ch for ch in str(x) if ch.isdigit())
    return (int(s) if s else 9999, str(x))


donors = sorted(set(did.donor_id.unique()) | set(dis.donor_id.unique()), key=dkey)
af01 = {d: af(did.loc[(did.donor_id == d) & (did.lane == "ARC01"), "cell"].tolist(), col_c, AD_c, DP_c) for d in donors}
af02 = {d: af(did.loc[(did.donor_id == d) & (did.lane == "ARC02"), "cell"].tolist(), col_c, AD_c, DP_c) for d in donors}
af0S = {d: af(dis.loc[dis.donor_id == d, "cell"].tolist(), col_s, AD_s, DP_s) for d in donors}
afmap = {"ARC01": af01, "ARC02": af02, "ARC0S": af0S}
present = {L: [d for d in donors if afmap[L][d] is not None] for L in afmap}
for L in ("ARC01", "ARC02", "ARC0S"):
    print(f"[present] {L}: {len(present[L])} donors {present[L]}", flush=True)

# long-form concordance for the 3 lane pairs
PAIRS = [("ARC01", "ARC02"), ("ARC0S", "ARC01"), ("ARC0S", "ARC02")]
rows = []
for la, lb in PAIRS:
    for i in present[la]:
        for j in present[lb]:
            r, nv = corr(afmap[la][i], afmap[lb][j])
            rows.append({"lane_a": la, "donor_a": i, "lane_b": lb, "donor_b": j,
                         "r": round(r, 4) if r == r else np.nan, "n_var": nv,
                         "same_donor": int(i == j)})
long = pd.DataFrame(rows)
os.makedirs(OUT, exist_ok=True)
long.to_csv(f"{OUT}/cA_crosslane_3lane_long.tsv", sep="\t", index=False)

print("\n=== cross-lane concordance summary ===", flush=True)
for la, lb in PAIRS:
    sub = long[(long.lane_a == la) & (long.lane_b == lb)]
    w = sub.pivot(index="donor_a", columns="donor_b", values="r")
    w = w.reindex(index=sorted(w.index, key=dkey), columns=sorted(w.columns, key=dkey))
    w.round(4).to_csv(f"{OUT}/cA_crosslane_3lane_{la}_vs_{lb}.tsv", sep="\t")
    diag = sub[sub.same_donor == 1].r.values
    off = sub[sub.same_donor == 0].r.values
    print(f"{la} vs {lb}: same-embryo r median {np.nanmedian(diag):.3f} "
          f"(min {np.nanmin(diag):.3f}, n={int(np.sum(~np.isnan(diag)))})  ·  "
          f"different-embryo r median {np.nanmedian(off):.3f} (max {np.nanmax(off):.3f})", flush=True)

# full 3-lane block matrix: every (lane, donor) pseudobulk-AF vector vs every other.
# Axis order = ARC01 donors | ARC02 donors | ARC0S donors (grouped by lane). The ARC01xARC02
# sub-block reproduces cA_crosslane_3lane_ARC01_vs_ARC02; the within-lane blocks carry the
# same-lane different-embryo null (diagonal self-correlation = 1).
blk_order = ([("ARC01", d) for d in present["ARC01"]]
             + [("ARC02", d) for d in present["ARC02"]]
             + [("ARC0S", d) for d in present["ARC0S"]])
blk_labels = [f"{L}:{d}" for L, d in blk_order]
blk = pd.DataFrame(index=blk_labels, columns=blk_labels, dtype=float)
for La, da in blk_order:
    for Lb, db in blk_order:
        r, _ = corr(afmap[La][da], afmap[Lb][db])
        blk.loc[f"{La}:{da}", f"{Lb}:{db}"] = r
blk.index.name = "lane_donor"
blk.round(4).to_csv(f"{OUT}/cA_crosslane_3lane_block.tsv", sep="\t")
print(f"[block] wrote cA_crosslane_3lane_block.tsv ({len(blk_order)}x{len(blk_order)})", flush=True)

# compact per-donor diagonal table (easy to read back; carries ARC0S cell counts).
# Restrict to donors substantial in at least one lane — drops the vireo over-split phantom columns
# (e.g. donor1/12/17, 2-11 cells each) that would otherwise appear as all-NaN rows indistinguishable
# from a real HH13 embryo that legitimately has no ARC0S left-half.
present_any = [d for d in donors if any(d in present[L] for L in present)]
diagrows = []
for d in present_any:
    rec = {"donor": d}
    for la, lb in PAIRS:
        s = long[(long.lane_a == la) & (long.lane_b == lb) & (long.donor_a == d) & (long.donor_b == d)]
        rec[f"{la}_vs_{lb}"] = float(s.r.iloc[0]) if len(s) and s.r.iloc[0] == s.r.iloc[0] else np.nan
    rec["n_arc0s_cells"] = int((dis.donor_id == d).sum())
    diagrows.append(rec)
diagtab = pd.DataFrame(diagrows)
diagtab.to_csv(f"{OUT}/cA_crosslane_3lane_diag.tsv", sep="\t", index=False)
print("\n=== per-donor same-embryo cross-lane r (ARC0S cols NaN = no left-half in ARC0S = HH13) ===", flush=True)
print(diagtab.round(3).to_string(index=False), flush=True)
print(f"\nwrote {OUT}/cA_crosslane_3lane_long.tsv  (+ 3 per-pair wide matrices + _diag.tsv)", flush=True)
