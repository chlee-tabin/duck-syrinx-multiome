#!/usr/bin/env python3
"""c6b — Group B donor -> embryo cross-check (sex + GFP + lane vs design).

Joins the per-donor orthogonal signals into one identity table and checks them against the
12-embryo design (design/chick_groupB_lr_design.tsv):
  - GFP-status  (c5b: per_donor_gfp.tsv)
  - ZW sex      (c5c: per_donor_sex_wgene.tsv)
  - lane split  (c4b: best_seed_donor_by_lane.tsv)  — every embryo is L/R-split across lanes
Full naming to a specific embryo_id also needs STAGE (HH13/15/25 via mesenchyme PCA, a later
c6 step). Until then this maps at the (GFP x sex) level, confirms the 6 GFP / 6 WT balance,
and flags under-captured donors (the ones needing lenient barcode recovery to characterize).

  python3 pipeline/chick/c6b_assign_embryo.py [IN_DIR] [DESIGN_TSV]
    IN_DIR holds the three per-donor TSVs (pulled from scratch seed_qc/).
"""
import os, sys
import pandas as pd

# repo root from this script's location (pipeline/chick/<this>) — portable across checkouts/users,
# unlike a hardcoded worktree path.
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = sys.argv[1] if len(sys.argv) > 1 else f"{REPO}/docs/chick_control/sex_wgene"
DESIGN = sys.argv[2] if len(sys.argv) > 2 else f"{REPO}/design/chick_groupB_lr_design.tsv"
MIN_CELLS = 100

gfp = pd.read_csv(f"{IN}/per_donor_gfp.tsv", sep="\t")[["donor", "n_cells", "frac_EGFP_pos", "call"]]
gfp = gfp.rename(columns={"call": "gfp_call"})
sex = pd.read_csv(f"{IN}/per_donor_sex_wgene.tsv", sep="\t")[
    ["donor", "median_pct_w", "median_pct_z", "sex_W", "sex_Z", "sex_call"]]
lane = pd.read_csv(f"{IN}/best_seed_donor_by_lane.tsv", sep="\t")
lane.columns = ["donor"] + list(lane.columns[1:])
lc = [c for c in lane.columns if c != "donor"]
lane["gex09_n"] = lane[[c for c in lc if "09" in c or "A" in c]].sum(axis=1)
lane["gex10_n"] = lane[[c for c in lc if "10" in c or "B" in c]].sum(axis=1)

d = gfp.merge(sex, on="donor", how="outer").merge(
    lane[["donor", "gex09_n", "gex10_n"]], on="donor", how="left")
d["well_covered"] = d["n_cells"] >= MIN_CELLS
d["gfp_clean"] = d["gfp_call"].str.replace(r"\s*\(low-cell\)", "", regex=True)
d["class"] = d["gfp_clean"].astype(str) + " / " + d["sex_W"].astype(str)
d = d.sort_values(["well_covered", "gfp_clean", "sex_W", "n_cells"],
                  ascending=[False, True, True, False]).reset_index(drop=True)

design = pd.read_csv(DESIGN, sep="\t")
exp_gfp = (design.gfp_status == "GFP").sum()
exp_wt = (design.gfp_status == "WT").sum()

print("=== Group B per-donor identity (sex + GFP + lane) ===")
cols = ["donor", "n_cells", "gex09_n", "gex10_n", "gfp_clean", "frac_EGFP_pos",
        "sex_W", "sex_Z", "median_pct_w", "median_pct_z", "well_covered"]
print(d[cols].to_string(index=False))

wc = d[d.well_covered]
print(f"\n=== Design check (design = {exp_gfp} GFP + {exp_wt} WT, 4/stage HH13/15-16/25) ===")
print(f"well-covered donors (>= {MIN_CELLS} cells): {len(wc)}")
print("  (GFP x sex) classes among well-covered:")
print(wc.groupby(["gfp_clean", "sex_W"]).size().to_string())
under = d[~d.well_covered]
print(f"\nunder-captured donors (< {MIN_CELLS} cells, GT-distinct = real but low-yield): "
      f"{len(under)}")
for _, r in under.iterrows():
    print(f"  {r['donor']}: {int(r['n_cells'])} cells, gfp~{r['gfp_clean']}, "
          f"W%={r['median_pct_w']:.4f} (sex unreliable at this depth)")

print(f"\nTotal distinct genotypes recovered: {len(d)} (= 12 design embryos if GT-distinct). "
      f"Well-characterized: {len(wc)}; need recovery: {len(under)}.")
print("Stage (HH13/15/25) naming -> mesenchyme-PCA staging (c6), not run here.")
out = f"{IN}/groupB_donor_identity.tsv"
d.to_csv(out, sep="\t", index=False)
print(f"\nWrote {out}")
