#!/usr/bin/env python3
"""stats_cycling_lr_hh25.py — P-S1 primary + sensitivity arms (supp S16).

Implements Session A of `docs/m5_esr1_prolif_stats_plan.md` (not distributed) for the
proliferation-timing supplement.  Everything here runs off committed tables.

PRIMARY (plan §2.3)
  Composition-standardised cycling fraction (`phase != "G1"`) per embryo x side
  over the HH25 inference pool, exact Wilcoxon signed-rank on the 8 paired
  differences, TWO-SIDED — the claim is an ABSENCE, so there is no a-priori
  direction to spend.

  THE CLAIM VEHICLE IS THE INTERVAL, NOT THE p-VALUE.  The licensed sentence
  "no Left excess at HH25" requires the exact CI upper bound to fall below
  delta = +0.137 (the E09 Left cycling effect the later-stage claim rests on;
  pre-specified).  If it does not, the supplement must say the
  multiome cannot yet exclude an E09-sized Left excess, and may NOT fall back
  on the null p-value.  "Symmetric" is prohibited as a positive claim.  The
  Right-side tilt is reported two-sided with its CI as a trend, never asserted
  as a finding.

COLLIDER SAFETY
  This runs on the ccreg (cell-cycle-regressed) partition, whose clusters are
  therefore NOT defined by cycling; `phase`/`S_score`/`G2M_score` are computed
  upstream of that regression.  Running the contrast on the non-ccreg partition
  would be a collider artefact and is prohibited.

ALSO REPORTED
  The §2.1 falsification check; sensitivity arms S1, S2, S5, S6; the
  HH13/HH17-16/HH21 secondary family (Holm-adjusted, supporting the temporal
  emergence framing only); and the post-hoc HH21 Ventral spl. contrast,
  descriptive.

Inputs  (all under figures/data/):
  d2g_lpm_subclustering_<stage>_ccreg.tsv.gz   per-cell leiden + phase + S/G2M
  lr_inference_multiome.tsv                    canonical L/R sidecar (ONLY source)
  step03_sex_per_genotype.tsv                  per-genotype sex
Output: figures/data/m5_cycling_lr_hh25_stats.tsv  (tidy, long)

Run:  python figures/figS16/stats_cycling_lr_hh25.py
      --selfcheck  runs the shared numeric self-test first
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, FIG)
from _shared.paired_standardised import (  # noqa: E402
    DATA, MIN_POOL_PER_SIDE, MIN_STRATUM, crude, eligible_strata, fixed_weights,
    load_pool_cells, paired_arm, paired_cohort, paired_diffs, standardise,
    stats_row, write_stats, _selfcheck)

OUT = DATA / "m5_cycling_lr_hh25_stats.tsv"
FRAC = "fraction of pool cells"
DELTA = 0.137          # pre-specified: the E09 Left cycling effect the CI must exclude
SECONDARY = ["HH13", "HH17/16", "HH21"]

EXPECT = dict(n_pairs=8, pool_L=6383, pool_R=10137, w_syringeal=0.146,
              w_ventral=0.854, min_pool_per_side=122, n_male=6, n_female=2,
              stratum_rule_fires_for=["g00"],
              # scoping values the plan quotes in §2.3 (crude, unstandardised pool)
              crude_k_L_gt_R=2, crude_p=0.0546875,
              hh21_ventral_p=0.0234375)


def check(log: list[dict], name: str, got, want, tol: float = 0.0) -> None:
    ok = (abs(float(got) - float(want)) <= tol) if isinstance(want, (int, float)) else (got == want)
    log.append(dict(name=name, got=got, want=want, ok=bool(ok)))
    print(f"  [{'ok ' if ok else 'FAIL'}] {name}: got {got!r}, plan says {want!r}")


def cycling_cells(stage: str) -> pd.DataFrame:
    """Inference-pool singlets at `stage` with side, sex and cell-cycle readouts."""
    cells = load_pool_cells(stage, "inference",
                            extra_cols=("phase", "S_score", "G2M_score"))
    bad = sorted(set(cells["phase"].astype(str)) - {"G1", "S", "G2M"})
    if bad:
        raise RuntimeError(f"{stage}: unexpected cell-cycle phase value(s) {bad} — the "
                           "cycling definition (phase != 'G1') would silently mis-split")
    cells["cycling"] = (cells["phase"] != "G1").astype(float)
    return cells


def stage_arm(stage: str) -> tuple[np.ndarray, dict, list[str]]:
    """Standardised cycling-fraction L-R differences for one stage.

    Returns (differences, cohort-rule drops, standardisation-infeasible embryos).

    The plan's two §2.1 rules can conflict away from HH25: the cohort rule keeps
    an embryo with >= 30 pool cells per side, but with k strata the stratum rule
    needs up to k x 20 cells on the thinner side, so an embryo can clear the
    cohort rule and still have NO eligible stratum.  Per the 2026-08-22 §2.1
    tie-break amendment, such embryos are EXCLUDED and the
    exclusion disclosed by name — never silently dropped.
    """
    cells = cycling_cells(stage)
    cohort, dropped = paired_cohort(cells)
    if len(cohort) < 4:
        raise RuntimeError(f"{stage}: only {len(cohort)} complete pairs clear "
                           f"{MIN_POOL_PER_SIDE} cells/side — too few for a paired test")
    cells = cells[cells["genotype"].isin(cohort)].copy()
    n = (cells.groupby(["genotype", "population", "inferred_side"]).size()
              .unstack("inferred_side", fill_value=0))
    for s in ("L", "R"):
        if s not in n.columns:
            n[s] = 0
    ok = n[(n["L"] >= MIN_STRATUM) & (n["R"] >= MIN_STRATUM)]
    infeasible = sorted(set(cohort) - set(ok.index.get_level_values("genotype")))
    keep = [g for g in cohort if g not in infeasible]
    cells = cells[cells["genotype"].isin(keep)].copy()
    std, log = standardise(cells, "cycling", fixed_weights(cells))
    for s in log:
        print(f"    {stage} substitution: {s}")
    return paired_diffs(std, keep).values, dropped, infeasible


def main() -> int:
    if "--selfcheck" in sys.argv:
        _selfcheck()

    xchecks: list[dict] = []
    rows: list[dict] = []

    print("== substrate (HH25) ==")
    cells = cycling_cells("HH25")
    cohort, dropped = paired_cohort(cells)
    if dropped:
        print(f"  cohort drops: {dropped}")
    cells = cells[cells["genotype"].isin(cohort)].copy()
    sex = cells.groupby("genotype")["sex"].first()
    males = sorted(sex[sex == "M"].index)
    weights = fixed_weights(cells)
    per_side = cells.groupby(["genotype", "inferred_side"]).size()

    print("== cross-checks against the pre-specified plan ==")
    check(xchecks, "complete L/R pairs", len(cohort), EXPECT["n_pairs"])
    check(xchecks, "pool cells Left", int((cells["inferred_side"] == "L").sum()), EXPECT["pool_L"])
    check(xchecks, "pool cells Right", int((cells["inferred_side"] == "R").sum()), EXPECT["pool_R"])
    check(xchecks, "fixed weight Syringeal", round(float(weights["Syringeal"]), 3),
          EXPECT["w_syringeal"], tol=0.0005)
    check(xchecks, "fixed weight Ventral spl.", round(float(weights["Ventral spl."]), 3),
          EXPECT["w_ventral"], tol=0.0005)
    check(xchecks, "min pool cells per side", int(per_side.min()), EXPECT["min_pool_per_side"])
    check(xchecks, "sex M", int((sex == "M").sum()), EXPECT["n_male"])
    check(xchecks, "sex F", int((sex == "F").sum()), EXPECT["n_female"])
    elig, subst = eligible_strata(cells, weights, MIN_STRATUM)
    fires = sorted(g for g, e in elig.items() if len(e) < len(weights))
    check(xchecks, f"embryos where the >={MIN_STRATUM}-cells-both-sides rule fires",
          fires, EXPECT["stratum_rule_fires_for"])
    for s in subst:
        print(f"    substitution: {s}")

    print("== falsification check (plan §2.1, reported, never gating) ==")
    tot = cells.groupby(["genotype", "inferred_side"]).size()
    syr = (cells[cells["population"] == "Syringeal"]
           .groupby(["genotype", "inferred_side"]).size().reindex(tot.index, fill_value=0))
    d_share = paired_diffs((syr / tot).rename("value").reset_index(), cohort)
    rows.append(paired_arm(
        "falsification_check", "Syringeal share of the inference pool (L - R)",
        d_share.values, "two-sided", "percentage points", scale=100,
        notes="pre-specified §2.1; a non-null result is what the standardisation absorbs "
              "and does NOT invalidate the primary"))

    print("== primary (plan §2.3) ==")
    std, _ = standardise(cells, "cycling", weights)
    d = paired_diffs(std, cohort)
    piv = std.pivot(index="genotype", columns="inferred_side", values="value")
    print(piv.round(4).assign(L_minus_R=d.round(4)).to_string())

    prim = paired_arm(
        "primary", "cycling fraction (phase != G1), composition-standardised (L - R)",
        d.values, "two-sided", FRAC,
        notes="PRIMARY, claim-bearing; exact Wilcoxon signed-rank; n = 8 embryos "
              "(16 half-embryo pseudobulks); ccreg partition so the clusters are NOT "
              "defined by cycling (collider-safe)")
    excludes = bool(prim["ci_high"] < DELTA)
    licensed = ('LICENSED: "no Left excess at HH25" — the exact CI upper bound '
                f'{prim["ci_high"]:+.4f} lies below delta = {DELTA:+.3f}'
                if excludes else
                'NOT LICENSED: the CI upper bound {:+.4f} does not clear delta = {:+.3f}; '
                'the supplement must say the multiome cannot yet exclude an E09-sized Left '
                'excess, and may NOT fall back on the null p-value'.format(prim["ci_high"], DELTA))
    prim["notes"] += (f"; delta={DELTA:+.3f} (E09 Left cycling effect, pre-specified); "
                      f"excludes_delta={excludes}; {licensed}; "
                      f'"symmetric" is PROHIBITED as a positive claim; '
                      f"median(L-R)={float(np.median(d)):+.4f}")
    rows.append(prim)
    print(f"  HL DeltaL-R = {prim['estimate']:+.4f} "
          f"[{prim['ci_low']:+.4f}, {prim['ci_high']:+.4f}], p = {prim['p']:.4g}, "
          f"{prim['direction_count']}")
    print(f"  {licensed}")
    rows.append(stats_row(
        "benchmark_delta", "E09 Left cycling effect used as the exclusion benchmark",
        unit=FRAC, estimate=DELTA,
        notes="figS07 panel F, Syringeal chondrocyte (renamed 2026-08-22, ex-'Chondrogenic "
              "LPM'), p = 0.031 (old-scRNA, dissected sides); "
              "pre-specified so 'no Left excess yet at HH25' is a bounded, "
              "falsifiable statement rather than an acceptance of the null"))

    print("== sensitivity arms ==")
    dc = paired_diffs(crude(cells, "cycling"), cohort)
    crude_row = paired_arm(
        "S1_crude", "cycling fraction, UNstandardised (L - R)", dc.values, "two-sided", FRAC,
        notes=f"shows the size of the composition correction; confounded by stratum mix; "
              f"median(L-R)={float(np.median(dc)):+.4f}")
    rows.append(crude_row)
    check(xchecks, "crude pool direction count (L>R)", int((dc > 0).sum()),
          EXPECT["crude_k_L_gt_R"])
    check(xchecks, "crude pool exact two-sided p", round(float(crude_row["p"]), 6),
          EXPECT["crude_p"], tol=1e-6)

    for pop in list(weights.index):
        p_cells = cells[cells["population"] == pop]
        pcoh, pdrop = paired_cohort(p_cells)
        if len(pcoh) < 4:
            rows.append(stats_row(f"S2_{pop}", f"cycling fraction within {pop}",
                                  n_embryos=len(pcoh),
                                  notes=f"NOT RUN — only {len(pcoh)} embryos clear "
                                        f"{MIN_POOL_PER_SIDE} cells/side"))
            continue
        rows.append(paired_arm(
            f"S2_{pop}", f"cycling fraction within {pop} (L - R)",
            paired_diffs(crude(p_cells[p_cells["genotype"].isin(pcoh)], "cycling"),
                         pcoh).values, "two-sided", FRAC, require_ci=False,
            notes=f"conditions on the boundary the 2026-08 seeded relabel moved; drops {sorted(pdrop)}"))

    rows.append(paired_arm(
        "S5_males_only", "cycling fraction, standardised, males only (L - R)",
        d.loc[males].values, "two-sided", FRAC,
        notes=f"n = {len(males)} male embryos ({sorted(males)}); underpowered; confirms the "
              "R-tilt is not driven by the two females. Male-specificity is NOT a multiome "
              "claim — with 6 M / 2 F no sex-interaction test is performed; that arm is "
              "carried entirely by the E09/E11 scRNA analysis"))

    for col, lab in [("S_score", "S"), ("G2M_score", "G2M")]:
        sstd, _ = standardise(cells, col, weights)
        rows.append(paired_arm(
            f"S6_mean_{lab}_score", f"mean {lab} score, standardised (L - R)",
            paired_diffs(sstd, cohort).values, "two-sided", f"{lab} score (arbitrary units)",
            notes=f"continuous analogue of the phase fraction; computed upstream of the "
                  "ccreg regression, so unaffected by it"))

    print("== secondary family (plan §2.3, Holm over three stages) ==")
    clean = []
    for stage in SECONDARY:
        ds, drops, infeasible = stage_arm(stage)   # a failing secondary arm raises
        note = (f"secondary family; supports the TEMPORAL EMERGENCE framing only; "
                f"cohort drops {sorted(drops)}")
        if infeasible:
            # 2026-08-22 plan amendment (§2.1 tie-break): an embryo that
            # clears the cohort rule but has NO stratum with >= MIN_STRATUM cells on both
            # sides has undefined standardisation and is EXCLUDED, disclosed.
            note += (f"; additionally excludes {infeasible} under the 2026-08-22 §2.1 "
                     f"tie-break amendment (clears >={MIN_POOL_PER_SIDE} pool cells/side "
                     f"but no stratum holds >={MIN_STRATUM} cells on both sides, so the "
                     "standardisation is undefined — excluded, disclosed)")
            print(f"  {stage}: excluding {infeasible} per the 2026-08-22 tie-break amendment")
        clean.append((stage, ds, note))
    if clean:
        raw = [paired_arm(f"secondary_{stage}", "cycling fraction, standardised (L - R)",
                          ds, "two-sided", FRAC, require_ci=(len(ds) >= 6), notes=note)
               for stage, ds, note in clean]
        fam = ", ".join(s for s, _, _ in clean)
        holm = multipletests([r["p"] for r in raw], method="holm")[1]
        for r, q in zip(raw, holm):
            r["notes"] += f"; holm_p={q:.4g} (family of {len(raw)}: {fam})"
            print(f"  {r['arm']}: {r['estimate']:+.4f} "
                  f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] p={r['p']:.4g} holm={q:.4g} "
                  f"{r['direction_count']} n={int(r['n_embryos'])}")
            rows.append(r)

    print("== post-hoc (descriptive) ==")
    h21 = cycling_cells("HH21")
    h21 = h21[h21["population"] == "Ventral spl."]
    vcoh, vdrop = paired_cohort(h21)
    vrow = paired_arm(
        "posthoc_HH21_ventral", "cycling fraction within HH21 Ventral spl. (L - R)",
        paired_diffs(crude(h21[h21["genotype"].isin(vcoh)], "cycling"), vcoh).values,
        "two-sided", FRAC,
        notes="POST-HOC across populations — descriptive only, no multiplicity frame; "
              f"cohort drops {sorted(vdrop)}")
    rows.append(vrow)
    print(f"  {vrow['arm']}: {vrow['estimate']:+.4f} p={vrow['p']:.4g} "
          f"{vrow['direction_count']} n={int(vrow['n_embryos'])}")
    check(xchecks, "HH21 Ventral spl. exact two-sided p", round(float(vrow["p"]), 6),
          EXPECT["hh21_ventral_p"], tol=1e-6)

    for c in xchecks:
        rows.append(stats_row(
            "crosscheck", f"plan cross-check: {c['name']}",
            notes=f"{'PASS' if c['ok'] else 'DISCREPANCY'} — recomputed {c['got']!r}, "
                  f"plan pre-registers {c['want']!r}"))

    df = write_stats(rows, OUT)
    pd.set_option("display.width", 250)
    print(f"\n== wrote {OUT} ==")
    print(df[["arm", "estimate", "ci_low", "ci_high", "p", "sidedness",
              "n_embryos", "direction_count"]].to_string(index=False))

    bad = [c["name"] for c in xchecks if not c["ok"]]
    if bad:
        print(f"\n!! PLAN CROSS-CHECK DISCREPANCIES: {bad} — recorded in {OUT.name}; "
              "the affected arm must be re-read against the plan, not forced")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
