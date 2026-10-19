#!/usr/bin/env python3
"""c4 — chick vireo seed QC + cross-lane concordance.

Adapts a4_vireo_seed_qc.py. For each swept N, rank the 20 seeds by logLik, count
donors, flag against expected (12). For the best seed of each N, build the
per-donor x per-lane (ARC01/-1, ARC02/-2, ARC0S/-S) cell table: the 12 true embryos
appear in BOTH ARC01 and ARC02 (same embryos in both lanes) — the strongest
internal check that the demux is real (replaces the duck's external baseline).

Run:  conda run -n vireo python3 pipeline/chick/c4_seed_qc.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, glob
import pandas as pd

WORK = f"{SCRATCH_USER}/chick_demux"
SEEDED = f"{WORK}/vireo_seeded"
OUT = f"{SEEDED}/seed_qc"
EXPECTED = 12
LANE = {"1": "ARC01", "2": "ARC02", "S": "ARC0S-chick"}
os.makedirs(OUT, exist_ok=True)


def loglik(p):
    if not os.path.isfile(p):
        return None
    for ln in open(p):
        if ln.startswith("logLik:"):
            return float(ln.split(":")[1])
    return None


def n_donors(summary_path):
    if not os.path.isfile(summary_path):
        return None, None
    df = pd.read_csv(summary_path, sep="\t")
    # cellSNP/vireo summary: Var1 = donor/doublet/unassigned, Freq = count
    sing = df[~df["Var1"].isin(["doublet", "unassigned"])]
    return len(sing), int(sing["Freq"].sum())


rows = []
for ndir in sorted(glob.glob(f"{SEEDED}/N*")):
    N = os.path.basename(ndir).replace("N", "")
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        seed = os.path.basename(sd).replace("seed_", "")
        ll = loglik(f"{sd}/_log.txt")
        nd, nsing = n_donors(f"{sd}/summary.tsv")
        if ll is None:
            continue
        rows.append({"N": int(N), "seed": int(seed), "logLik": ll,
                     "n_donors": nd, "n_singlets": nsing, "dir": sd})

if not rows:
    raise SystemExit("No completed vireo seeds found under " + SEEDED)

df = pd.DataFrame(rows).sort_values(["N", "logLik"], ascending=[True, False])
df.to_csv(f"{OUT}/seed_qc_report.tsv", sep="\t", index=False)

print("=== Seed ranking per N (top by logLik) ===")
for N in sorted(df["N"].unique()):
    g = df[df["N"] == N]
    stable = (g["n_donors"] == EXPECTED).mean()
    best = g.iloc[0]
    flag = "" if best["n_donors"] == EXPECTED else f"  <-- best seed has {best['n_donors']} donors!"
    print(f"  N={N:>2}: best logLik={best['logLik']:.4e} seed={best['seed']} "
          f"n_donors={best['n_donors']} singlets={best['n_singlets']} | "
          f"{stable*100:.0f}% of seeds gave {EXPECTED} donors{flag}")

# Cross-lane concordance for the best seed at the expected N (fallback: best overall)
pick = df[(df["N"] == EXPECTED) & (df["n_donors"] == EXPECTED)]
if pick.empty:
    pick = df[df["n_donors"] == EXPECTED]
if pick.empty:
    print("\nNo seed recovered the expected donor count; skipping cross-lane table.")
    raise SystemExit(0)

best = pick.iloc[0]
did = pd.read_csv(f"{best['dir']}/donor_ids.tsv", sep="\t")
did["lane"] = did["cell"].str.rsplit("-", n=1).str[-1].map(LANE).fillna("?")
sing = did[~did["donor_id"].isin(["doublet", "unassigned"])]
tab = pd.crosstab(sing["donor_id"], sing["lane"])
tab.to_csv(f"{OUT}/best_seed_donor_by_lane.tsv", sep="\t")

print(f"\n=== Cross-lane composition (best seed: N={best['N']} seed={best['seed']}, "
      f"{best['dir']}) ===")
print(tab.to_string())
both = ((tab.get("ARC01", 0) > 0) & (tab.get("ARC02", 0) > 0)).sum()
print(f"\nDonors present in BOTH ARC01 and ARC02: {both}/{len(tab)} "
      f"(expect {EXPECTED} — same embryos in both lanes)")
print(f"\nWrote: {OUT}/seed_qc_report.tsv, {OUT}/best_seed_donor_by_lane.tsv")
