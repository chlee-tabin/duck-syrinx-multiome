#!/usr/bin/env python3
"""stats_esr1_lr_hh25.py — E-S1 primary + sensitivity arms (supp S17).

Implements Session A of `docs/m5_esr1_prolif_stats_plan.md` (not distributed) for the ESR1
laterality supplement.  Everything here runs off committed tables.

PRIMARY (plan §2.2)
  Composition-standardised ESR1+ detection fraction per embryo x side over the
  HH25 inference pool (Ventral spl. + Syringeal), exact Wilcoxon signed-rank on
  the 8 paired differences, ONE-SIDED (Left > Right — direction fixed a priori
  by the independently dissected E09/E11 dataset).  Reported as a
  Hodges-Lehmann effect in percentage points with an exact distribution-free CI.
  Multiplicity frame: one test, one gene, one stage.  ESR1 is a pre-specified
  target of this manuscript, not a discovery.

ALSO REPORTED
  §2.1 falsification check (per-embryo Syringeal share difference) — reported,
  never gating; the pre-specified Dorsal sub-epi. bilateral negative control;
  the committed pool-level paired DESeq2 row (a table read, no refit); and
  sensitivity arms S1-S8 (all non-claim-bearing).

Inputs  (all under figures/data/):
  d2g_lpm_subclustering_HH25_ccreg.tsv.gz   per-cell leiden + QC (ccreg = the
                                            partition the annotation SSOT keys on)
  d2g_signaling_HH25.tsv.gz                 per-cell ESR1_cp10k
  lr_inference_multiome.tsv                 canonical L/R sidecar (ONLY source)
  step03_sex_per_genotype.tsv               per-genotype sex
  h2_deg_POOL_HH25.tsv                      committed paired DESeq2 (supporting; S4 parity)
  h1_pseudobulk_snrna_{counts.tsv.gz,metadata.tsv}  inference-pool pseudobulks (S4 refit)
  hh25_percell_counts_esr1.tsv.gz           per-cell raw ESR1 UMI + depth (S8)
Output: figures/data/m5_esr1_lr_hh25_stats.tsv  (tidy, long)

Run:  conda run -n duckfig python figures/figS17/stats_esr1_lr_hh25.py
      --selfcheck  runs the shared numeric self-test first
      (env from envs/duckfig_portable.yml: S4 needs DESeq2 and S8 glmmTMB in the R
      that DUCK_RSCRIPT or PATH resolves)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, FIG)
from _shared.paired_standardised import (  # noqa: E402
    DATA, MIN_POOL_PER_SIDE, MIN_STRATUM, crude, eligible_strata, fixed_weights,
    load_pool_cells, paired_arm, paired_cohort, paired_diffs, standardise,
    stats_row, write_stats, _selfcheck)
from _shared.subcluster_map import POOLS  # noqa: E402

OUT = DATA / "m5_esr1_lr_hh25_stats.tsv"
PP = "percentage points"
# R for the S8 glmmTMB fit: DUCK_RSCRIPT if set, else the Rscript on PATH (the figure env).
RSCRIPT = os.environ.get("DUCK_RSCRIPT") or shutil.which("Rscript")

# Values the plan pre-registers (§2.1/§2.2).  Recomputed and compared; a
# mismatch is recorded as a `crosscheck` row and exits non-zero rather than
# being forced through.
EXPECT = dict(n_pairs=8, pool_L=6383, pool_R=10137, w_syringeal=0.146,
              w_ventral=0.854, min_pool_per_side=122, n_male=6, n_female=2,
              stratum_rule_fires_for=["g00"],
              deseq2=dict(log2FoldChange=0.343, lfcSE=0.142, pvalue=0.0154,
                          padj=0.254, n_paired_genotypes=8))


def check(log: list[dict], name: str, got, want, tol: float = 0.0) -> None:
    ok = (abs(float(got) - float(want)) <= tol) if isinstance(want, (int, float)) else (got == want)
    log.append(dict(name=name, got=got, want=want, ok=bool(ok)))
    print(f"  [{'ok ' if ok else 'FAIL'}] {name}: got {got!r}, plan says {want!r}")


def esr1_cells() -> pd.DataFrame:
    """HH25 inference-pool singlets with side, sex, ESR1 detection and level."""
    cells = load_pool_cells("HH25", "inference")
    sig = pd.read_csv(DATA / "d2g_signaling_HH25.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", "ESR1_cp10k"])
    n_before = len(cells)
    cells = cells.merge(sig, on=["genotype", "library", "barcode"], how="left")
    if len(cells) != n_before:
        raise RuntimeError(f"ESR1 join changed the cell count {n_before} -> {len(cells)} "
                           "— duplicate (genotype, library, barcode) keys in the signaling export")
    if cells["ESR1_cp10k"].isna().any():
        bad = cells.loc[cells["ESR1_cp10k"].isna(), ["genotype", "library"]].drop_duplicates()
        raise RuntimeError(f"pool cells with no ESR1_cp10k row: {bad.to_dict('records')}")
    cells["esr1_pos"] = (cells["ESR1_cp10k"] > 0).astype(float)
    cells["esr1_lin"] = np.expm1(cells["ESR1_cp10k"])   # linear CP10K (S3)
    return cells


def population_share(cells: pd.DataFrame, population: str) -> pd.DataFrame:
    """Per embryo x side share of the pool held by `population` (falsification check)."""
    if population not in set(cells["population"]):
        raise RuntimeError(f"population {population!r} absent from the pool "
                           f"({sorted(set(cells['population']))}) — d2g relabel drift")
    tot = cells.groupby(["genotype", "inferred_side"]).size()
    k = (cells[cells["population"] == population]
         .groupby(["genotype", "inferred_side"]).size().reindex(tot.index, fill_value=0))
    return (k / tot).rename("value").reset_index()


def deseq2_support_row() -> dict:
    """The committed pool-level paired DESeq2 ESR1 row — a table read, no refit."""
    d = pd.read_csv(DATA / "h2_deg_POOL_HH25.tsv", sep="\t")
    r = d[d["gene_id"] == "ESR1"]
    if len(r) != 1:
        raise RuntimeError(f"h2_deg_POOL_HH25.tsv has {len(r)} ESR1 rows, expected 1")
    r = r.iloc[0]
    n_tested = int(d["padj"].notna().sum())
    return stats_row(
        "support_deseq2_h2", "ESR1 pool-level paired DESeq2 (committed, no refit)",
        unit="log2 fold change (L/R)", estimate=float(r["log2FoldChange"]),
        ci_low=float(r["log2FoldChange"] - 1.96 * r["lfcSE"]),
        ci_high=float(r["log2FoldChange"] + 1.96 * r["lfcSE"]),
        p=float(r["pvalue"]), sidedness="two-sided (Wald)",
        n_embryos=int(r["n_paired_genotypes"]),
        notes=(f"lfcSE={r['lfcSE']:.4g}; baseMean={r['baseMean']:.4g}; Wald CI (not exact); "
               f"genome-wide padj={r['padj']:.4g} over {n_tested} genes clearing DESeq2 "
               "independent filtering — reported here under the targeted single-gene frame, "
               "with the genome-wide padj disclosed in the same sentence. Different estimator "
               "on the SAME cells: agreement is a robustness statement, not independent evidence"))


def s4_deseq2_refit() -> dict:
    """S4 — local DESeq2 refit ~ genotype + syringeal_fraction + side on the committed
    h1 inference-pool counts, pooled exactly as h2's POOL scope.  The same pooling fitted
    ~ genotype + side must first reproduce the committed h2 ESR1 row (parity gate);
    a mismatch means the pooling differs from h2 and the arm is not reported."""
    if not RSCRIPT or not os.path.exists(RSCRIPT):
        raise FileNotFoundError(f"S4: Rscript not found ({RSCRIPT!r}); set DUCK_RSCRIPT or put Rscript on PATH")
    pops = POOLS["HH25"]["inference"]
    with tempfile.TemporaryDirectory() as td:
        dst = os.path.join(td, "s4.tsv")
        p = subprocess.run([RSCRIPT, os.path.join(HERE, "s4_deseq2_refit.R"),
                            str(DATA / "h1_pseudobulk_snrna_counts.tsv.gz"),
                            str(DATA / "h1_pseudobulk_snrna_metadata.tsv"),
                            "HH25", "|".join(pops), "ESR1", dst],
                           capture_output=True, text=True)
        print(p.stdout, end="")
        if p.returncode != 0 or not os.path.exists(dst):
            tail = "\n".join((p.stderr or "").strip().splitlines()[-6:])
            raise RuntimeError(f"S4: DESeq2 refit failed (rc={p.returncode}):\n{tail}")
        fit = pd.read_csv(dst, sep="\t").set_index("model")
    if set(fit.index) != {"parity", "composition"}:
        raise RuntimeError(f"s4_deseq2_refit.R returned models {sorted(fit.index)}")

    # parity gate against the committed h2 POOL row (same cells, same pooling, same design)
    h2 = pd.read_csv(DATA / "h2_deg_POOL_HH25.tsv", sep="\t")
    ref = h2[h2["gene_id"] == "ESR1"].iloc[0]
    par = fit.loc["parity"]
    for col in ("log2FoldChange", "lfcSE", "pvalue"):
        if not np.isclose(par[col], ref[col], rtol=1e-6, atol=0):
            raise RuntimeError(f"S4 parity: ~ genotype + side refit {col}={par[col]!r} != "
                               f"committed h2 {ref[col]!r} — pooling differs from h2")
    if int(par["n_genes_padj_nonNA"]) != int(h2["padj"].notna().sum()):
        raise RuntimeError("S4 parity: genes clearing independent filtering differ from h2")
    if int(par["n_paired_genotypes"]) != int(ref["n_paired_genotypes"]):
        raise RuntimeError("S4 parity: paired-genotype count differs from h2")

    r = fit.loc["composition"]
    print(f"  S4 parity OK (log2FC {par['log2FoldChange']:.6g} vs h2 {ref['log2FoldChange']:.6g}); "
          f"composition fit log2FC={r['log2FoldChange']:.4g} p={r['pvalue']:.4g}")
    return stats_row(
        "S4_deseq2_refit",
        "ESR1, local DESeq2 refit ~ genotype + syringeal_fraction + side on h1 pool counts",
        unit="log2 fold change (L/R)", estimate=float(r["log2FoldChange"]),
        ci_low=float(r["log2FoldChange"] - 1.96 * r["lfcSE"]),
        ci_high=float(r["log2FoldChange"] + 1.96 * r["lfcSE"]),
        p=float(r["pvalue"]), sidedness="two-sided (Wald)",
        n_embryos=int(r["n_paired_genotypes"]),
        notes=(f"lfcSE={r['lfcSE']:.4g}; baseMean={r['baseMean']:.4g}; Wald CI (not exact); "
               f"genome-wide padj={r['padj']:.4g} over {int(r['n_genes_padj_nonNA'])} genes "
               "clearing DESeq2 independent filtering; composition-covariate variant of the "
               "committed h2 pool arm: the committed h1 inference-pool counts summed per embryo "
               f"x side as h2 does (cell floor {int(r['floor_cells'])}, poscounts size factors, "
               "parametric dispersion), syringeal_fraction = Syringeal / pool nuclei per half; "
               "parity: the same pooling fitted ~ genotype + side reproduces the committed h2 "
               f"ESR1 row (log2FC {par['log2FoldChange']:.3f}, p {par['pvalue']:.4f}, "
               f"{int(par['n_genes_tested'])} genes tested); SENSITIVITY only, "
               "non-claim-bearing; NB/Wald assumptions at n = 8 embryos are heavier than the "
               "exact paired test"))


def s8_depth_adjusted(cells: pd.DataFrame, weights: pd.Series) -> list[dict]:
    """S8 — depth-adjusted per-cell detection: binomial GLMM, log-depth offset,
    embryo random intercept, on exactly the cells the standardised estimate uses."""
    pcf = DATA / "hh25_percell_counts_esr1.tsv.gz"
    # S8 is required (the figS17 panel C renderer reads its rows): fail loud, never a NOT RUN row
    if not pcf.exists():
        raise FileNotFoundError(f"S8: {pcf} absent")
    if not RSCRIPT or not os.path.exists(RSCRIPT):
        raise FileNotFoundError(f"S8: Rscript not found ({RSCRIPT!r}); set DUCK_RSCRIPT or put Rscript on PATH")

    elig, _ = eligible_strata(cells, weights, MIN_STRATUM)
    keep = [p in elig[g] for g, p in zip(cells["genotype"], cells["population"])]
    sub = cells[keep].copy()

    pc = pd.read_csv(pcf, sep="\t")
    n_before = len(sub)
    sub = sub.merge(pc[["genotype", "library", "barcode", "ESR1_count", "n_counts"]],
                    on=["genotype", "library", "barcode"], how="left")
    if len(sub) != n_before:
        raise RuntimeError(f"per-cell depth join changed the row count {n_before} -> {len(sub)}")
    if sub[["ESR1_count", "n_counts"]].isna().any().any():
        bad = sub.loc[sub["n_counts"].isna(), ["genotype", "library"]].drop_duplicates()
        raise RuntimeError(f"pool cells absent from hh25_percell_counts_esr1: {bad.to_dict('records')}")
    if (sub["n_counts"] <= 0).any():
        raise RuntimeError("non-positive n_counts — log-depth offset undefined")
    # the raw count and the committed log1p-CP10K must agree on DETECTION
    disagree = int(((sub["ESR1_count"] > 0) != (sub["esr1_pos"] > 0)).sum())
    if disagree:
        raise RuntimeError(f"{disagree} cells where ESR1_count>0 disagrees with ESR1_cp10k>0 "
                           "— the two exports are from different cell sets or vintages")
    sub["detected"] = (sub["ESR1_count"] > 0).astype(int)

    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "cells.tsv")
        dst = os.path.join(td, "fit.tsv")
        sub[["genotype", "inferred_side", "population", "detected", "n_counts"]].to_csv(
            src, sep="\t", index=False)
        p = subprocess.run([RSCRIPT, os.path.join(HERE, "s8_glmm.R"), src, dst],
                           capture_output=True, text=True)
        print(p.stdout, end="")
        if p.returncode != 0 or not os.path.exists(dst):
            tail = "\n".join((p.stderr or "").strip().splitlines()[-6:])
            raise RuntimeError(f"S8: glmmTMB fit failed (rc={p.returncode}):\n{tail}")
        fit = pd.read_csv(dst, sep="\t")

    n_emb = sub["genotype"].nunique()
    suffix = {"side": "", "side_pop": "_pop_adj", "side_slope": "_embryo_slope"}
    extra = {
        "side": "the plan's S8 spec; side is a CELL-level fixed effect, so the Wald p is "
                "anticonservative by orders of magnitude — the replication unit is the "
                "embryo (n=8), not the cell. Read the estimate, not this p",
        "side_pop": "+ population fixed effect = the GLMM analogue of the composition "
                    "standardisation; same cell-level pseudoreplication caveat",
        "side_slope": "random side slope by embryo — the side effect is tested at the EMBRYO "
                      "level, so this is the row comparable to the primary",
    }
    rows = []
    for _, r in fit.iterrows():
        model = str(r["model"])
        if model not in suffix:
            raise RuntimeError(f"s8_glmm.R returned an unknown model label {model!r}")
        rows.append(stats_row(
            f"S8_depth_adjusted{suffix[model]}",
            f"ESR1 detected (per cell), binomial GLMM [{model}] + offset(log depth)",
            unit="log odds (L vs R)", estimate=float(r["estimate"]),
            ci_low=float(r["ci_low"]), ci_high=float(r["ci_high"]), p=float(r["p"]),
            sidedness="two-sided (Wald)", n_embryos=n_emb,
            notes=(f"binomial GLMM (glmmTMB), {int(r['n_cells'])} pool cells over {n_emb} "
                   f"embryos; SENSITIVITY only, non-claim-bearing; {extra[model]}")))
    return rows


def main() -> int:
    if "--selfcheck" in sys.argv:
        _selfcheck()

    xchecks: list[dict] = []
    rows: list[dict] = []

    print("== substrate ==")
    cells = esr1_cells()
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
    share = population_share(cells, "Syringeal")
    d_share = paired_diffs(share, cohort)
    print(f"  per-embryo Syringeal share L-R: "
          f"{dict(zip(d_share.index, np.round(d_share.values, 4)))}")
    rows.append(paired_arm(
        "falsification_check", "Syringeal share of the inference pool (L - R)",
        d_share.values, "two-sided", PP, scale=100,
        notes="pre-specified §2.1; a non-null result is what the standardisation absorbs "
              "and does NOT invalidate the primary"))

    print("== primary (plan §2.2) ==")
    std, _ = standardise(cells, "esr1_pos", weights)
    d = paired_diffs(std, cohort)
    piv = std.pivot(index="genotype", columns="inferred_side", values="value")
    print((piv * 100).round(2).assign(L_minus_R=(d * 100).round(2)).to_string())
    rows.append(paired_arm(
        "primary", "ESR1+ detection fraction, composition-standardised (L - R)",
        d.values, "greater", PP, scale=100,
        notes="PRIMARY, claim-bearing; exact Wilcoxon signed-rank; direction fixed a priori "
              "by the independently dissected E09/E11 dataset; n = 8 embryos "
              "(16 half-embryo pseudobulks); one test on one gene at one stage — no "
              f"genome-wide correction applies; exact one-sided p floor at n=8 = {1/256:.4f}"))

    print("== negative control (plan §2.2) ==")
    dor = load_pool_cells("HH25", "bilateral")
    sig = pd.read_csv(DATA / "d2g_signaling_HH25.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", "ESR1_cp10k"])
    dor = dor.merge(sig, on=["genotype", "library", "barcode"], how="left")
    if dor["ESR1_cp10k"].isna().any():
        raise RuntimeError("Dorsal sub-epi. cells with no ESR1_cp10k row")
    dor["esr1_pos"] = (dor["ESR1_cp10k"] > 0).astype(float)
    dcoh, ddrop = paired_cohort(dor)
    print(f"  dorsal cohort n={len(dcoh)}; dropped {ddrop}")
    dstd, _ = standardise(dor[dor["genotype"].isin(dcoh)], "esr1_pos",
                          fixed_weights(dor[dor["genotype"].isin(dcoh)]))
    rows.append(paired_arm(
        "negcontrol_dorsal", "ESR1+ detection fraction, Dorsal sub-epi. control (L - R)",
        paired_diffs(dstd, dcoh).values, "two-sided", PP, scale=100,
        notes=f"pre-specified bilateral control, expected null; single-stratum pool so "
              f"standardisation is vacuous; same >={MIN_POOL_PER_SIDE} cells/side cohort rule "
              f"as the primary drops {sorted(ddrop)}"))
    rows.append(paired_arm(
        "negcontrol_dorsal_all8", "ESR1+ detection fraction, Dorsal sub-epi. control, all 8 embryos",
        paired_diffs(crude(dor, "esr1_pos"), sorted(dor['genotype'].unique())).values,
        "two-sided", PP, scale=100,
        notes="cohort rule waived (crude fraction) so the control is also reported on the "
              "full 8-embryo set; two embryos contribute <30 dorsal cells on one side"))

    print("== supporting DESeq2 (committed table read) ==")
    dr = deseq2_support_row()
    print(f"  ESR1 log2FC={dr['estimate']:.4g} p={dr['p']:.4g} | {dr['notes']}")
    check(xchecks, "h2 DESeq2 ESR1 log2FoldChange", round(dr["estimate"], 3),
          EXPECT["deseq2"]["log2FoldChange"], tol=0.0005)
    check(xchecks, "h2 DESeq2 ESR1 pvalue", round(dr["p"], 4),
          EXPECT["deseq2"]["pvalue"], tol=5e-5)
    check(xchecks, "h2 DESeq2 n_paired_genotypes", int(dr["n_embryos"]),
          EXPECT["deseq2"]["n_paired_genotypes"])
    rows.append(dr)

    print("== sensitivity arms ==")
    # S1 crude
    dc = paired_diffs(crude(cells, "esr1_pos"), cohort)
    rows.append(paired_arm("S1_crude", "ESR1+ detection fraction, UNstandardised (L - R)",
                           dc.values, "greater", PP, scale=100,
                           notes="shows the size of the composition correction; confounded by "
                                 "stratum mix"))
    # S2 per population
    for pop in list(weights.index):
        p_cells = cells[cells["population"] == pop]
        pcoh, pdrop = paired_cohort(p_cells)
        if len(pcoh) < 4:
            rows.append(stats_row(f"S2_{pop}", f"ESR1+ detection fraction within {pop}",
                                  n_embryos=len(pcoh),
                                  notes=f"NOT RUN — only {len(pcoh)} embryos clear "
                                        f"{MIN_POOL_PER_SIDE} cells/side"))
            continue
        rows.append(paired_arm(
            f"S2_{pop}", f"ESR1+ detection fraction within {pop} (L - R)",
            paired_diffs(crude(p_cells[p_cells["genotype"].isin(pcoh)], "esr1_pos"),
                         pcoh).values, "greater", PP, scale=100, require_ci=False,
            notes=f"conditions on the boundary the 2026-08 seeded relabel moved; drops {sorted(pdrop)}"))
    # S3 level
    lstd, _ = standardise(cells, "esr1_lin", weights)
    rows.append(paired_arm("S3_level_expm1",
                           "ESR1 mean linear CP10K (expm1), standardised (L - R)",
                           paired_diffs(lstd, cohort).values, "greater", "CP10K",
                           notes="estimator-dependent; the committed export is log1p-CP10K so the "
                                 "back-transform is forced"))
    # S4 DESeq2 refit (composition covariate) on the committed h1 pool counts
    rows.append(s4_deseq2_refit())
    # S5 males only
    rows.append(paired_arm("S5_males_only",
                           "ESR1+ detection fraction, standardised, males only (L - R)",
                           d.loc[males].values, "greater", PP, scale=100,
                           notes=f"n = {len(males)} male embryos ({sorted(males)}); underpowered; "
                                 "guards against a two-female artefact"))
    # S7 verdict
    rows.append(stats_row(
        "S7_markerscore_gate", "Cluster-free marker-score gate within the pool",
        notes="STILL DEFERRED — the d2g lineage score is sc.tl.score_genes, whose control bins "
              "need the FULL gene matrix; only 29 per-cell CP10K columns are committed, so any "
              "local gate is a DIFFERENT estimator with no source of truth. Per-cell depth is "
              "now available (S8), which removes the depth objection but not this one: the "
              "gate's source-of-truth score is cluster-side. A PITX2+ gate is circular outright; "
              "the identity panel (TBX5/TBX4/ZFPM2/ISL1 up, GATA4/BARX1/HAND2 down) is not "
              "verified L/R-neutral. No substitute gate was invented"))
    # S8 depth-adjusted GLMM
    rows.extend(s8_depth_adjusted(cells, weights))

    for c in xchecks:
        rows.append(stats_row(
            "crosscheck", f"plan cross-check: {c['name']}",
            notes=f"{'PASS' if c['ok'] else 'DISCREPANCY'} — recomputed {c['got']!r}, "
                  f"plan pre-registers {c['want']!r}"))

    df = write_stats(rows, OUT)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
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
