#!/usr/bin/env python3
"""cAN — chick PITX2 laterality: orientation-free Panel F + the directional control (cU unified).

Two readouts from ONE source of truth (the cAM per-cell populations), folded by the duck-faithful
chick POOLS (`_chick_subcluster_pools.py`; rationale in docs/chick_control/laterality_fold_methods.md, not distributed):

(1) **Orientation-free half-vs-half (Panel F, duck fig02C analog).**  Per embryo × spatial DOMAIN,
    PITX2 in one half vs the other; distance from the y=x diagonal = asymmetry REGARDLESS of which
    half is Left or Right.  Three facets, folded C-strict (a literal duck port):
        "LPM mes. (HH13)"          = HH13 inference (D/V unresolved pre-TES; no dorsal)
        "Ventral splanchnic pool"  = HH15-16 + HH25 inference
        "Dorsal sub-epi."          = HH15-16 + HH25 bilateral (the within-figure bilateral control)
    Skeletogenic LPM (was Somatic) / Mesothelium / Venous-pole spl. / paraxial are dropped at EVERY
    stage, as the duck does (the original cAN kept Skeletogenic LPM+Mesothelium at HH13 — fixed
    here).  EVERY embryo contributes
    (Group A + B, WT + GFP+ identically); the only per-embryo drop is data-quality (donor8).

(2) **Directional control (Group B only) — the unified, de-staled cU.**  Orient each Group B embryo
    by the wetlab-DESIGNED LEFT lane (chick_groupB_lr_design.tsv) using firm 3-stage × GFP design
    classes (uniform / singleton / elimination — NOT the old early/late binary), then test whether
    PITX2 is higher in that designed-LEFT half.  Runs on the INFERENCE pool only (where laterality
    lives), with a per-embryo Fisher test.  GFP+ counted identically to WT; a class is left out of the
    oriented count only when un-orientable (non-uniform design class) or not significantly lateralised
    — never as a quality tier (project convention: "Chick Control Conventions").

Reads : cAM_percell_populations.tsv.gz (group, lane/library, donor, stage, population, PITX2_cp10k)
          minus groupB_split_excluded_cells.tsv (163 Group B cS-split doublet/unassigned cells;
          92 of them in cAM — see _groupB_split_excluded.py)
        chick_stage_alignment.tsv ; design/chick_groupB_lr_design.tsv ; cAL_subcluster_annotations.tsv
Writes: cAN_halfvshalf.tsv          (Panel F input: one row per group×donor×domain)
        cAN_domain_counts.tsv       (kept/dropped audit — the dot accounting)
        cAN_directional_perembryo.tsv (Group B per-embryo: PITX2-high lane, designed-LEFT, p, concordant)
        cAN_directional_summary.tsv (per stage×gfp: n_orientable / n_sig / n_sig_concordant)

Run locally (cAM per-cell in repo) or on cluster (CHICK_D2G_DATA=<seed_qc>):
  conda run -n duckfig python3 pipeline/chick/cAN_stage_laterality.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chick_subcluster_pools import (  # noqa: E402
    facet_of, is_inference, FACET_ORDER, validate_against_annotations)
from _groupB_split_excluded import drop_split_excluded  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
# every cAN input AND output is figure-consumed, so the whole step reads/writes figures/data.
DATA = os.environ.get("CHICK_D2G_DATA",
                      os.path.join(REPO, "figures", "data", "chick_control"))
PERCELL = os.path.join(DATA, "cAM_percell_populations.tsv.gz")
STAGE_TSV = os.path.join(DATA, "chick_stage_alignment.tsv")
ANNOT = os.path.join(REPO, "figures", "data", "chick_control",
                     "cAL_subcluster_annotations.tsv")
DESIGN = os.path.join(REPO, "design", "chick_groupB_lr_design.tsv")
OUTDIR = DATA

MIN_HALF = 5                          # min cells per half for a per-embryo point (chick is sparse)
STAGES = ["HH13", "HH15-16", "HH25"]
EXCLUDE_GROUPB_DONORS = {"donor8"}    # data-quality (one half ~97% under-captured); NOT GFP-based
LANES_A = ("ARC01", "ARC02")
LANES_B = ("GEX09", "GEX10")


def log(m):
    print(m, flush=True)


def lane_of(s):
    """Map a library/lane string to one of the four known lanes; fail closed on anything else.

    Matches the full lane token (exact / prefix / suffix), NOT a bare 2-digit substring — the old
    `tag[-2:] in s` test mis-mapped e.g. "donor10" -> GEX10 and "donor01" -> ARC01. Raises rather
    than silently returning the raw string (which the `pc["half"]` assignment would otherwise lump
    into half "h2").
    """
    s = str(s)
    for tag in (*LANES_A, *LANES_B):
        if s == tag or s.startswith(tag) or s.endswith(tag):
            return tag
    raise ValueError(f"unrecognized library/lane string: {s!r} "
                     f"(expected one of {(*LANES_A, *LANES_B)})")


# ── designed-LEFT orientation (firm 3-stage × GFP classes; uniform / singleton / elimination) ──
def build_designed_left(design, percell_B, log=log):
    """Return {donor: designed_LEFT_lane or 'ambiguous'} for Group B, by firm-stage × GFP class.

    DESIGN-ALONE orientation — genuinely non-circular.  A donor is orientable only when its
    (firm stage × GFP) design class is UNIFORM (all candidate embryos share an L lane) or a
    SINGLETON; non-uniform classes are 'ambiguous'.  We deliberately do NOT resolve non-uniform
    classes by PITX2-elimination (cU did): assigning the minority embryo by its PITX2-high lane
    and then testing PITX2 concordance is CIRCULAR for that embryo.  So the directional headline
    rests only on embryos oriented by the wetlab design, with zero PITX2 input.  (See
    docs/chick_control/laterality_fold_methods.md §6, not distributed.)
    """
    des = design.copy()
    des["gfp"] = des["gfp_status"].map({"GFP": "GFP+", "WT": "WT"}).fillna(des["gfp_status"])
    des["L_lane"] = np.where(des["gex09_side"] == "L", "GEX09", "GEX10")
    donor_class = {donor: (g["stage"].iloc[0], g["gfp"].iloc[0])
                   for donor, g in percell_B.groupby("donor", observed=True)}
    designed_left = {}
    for (st, gf), grp in des.groupby(["stage", "gfp"]):
        lanes = set(grp["L_lane"])
        dl = next(iter(lanes)) if len(lanes) == 1 else "ambiguous"   # uniform/singleton -> orient
        for d in (d for d, c in donor_class.items() if c == (st, gf)):
            designed_left[d] = dl
    return designed_left


def directional_control(percell_B, designed_left, log=log):
    """Per Group B embryo: inference-pool PITX2 L/R vs designed-LEFT + Fisher significance."""
    rows = []
    for donor, g in percell_B.groupby("donor", observed=True):
        if donor in EXCLUDE_GROUPB_DONORS:
            continue
        st, gf = g["stage"].iloc[0], g["gfp"].iloc[0]
        inf = g[[is_inference(s, p) for s, p in zip(g["stage"], g["population"])]]
        g9, g10 = inf[inf.lane == "GEX09"], inf[inf.lane == "GEX10"]
        n9, n10 = len(g9), len(g10)
        pos9, pos10 = int((g9.PITX2.values > 0).sum()), int((g10.PITX2.values > 0).sum())
        f9 = pos9 / n9 if n9 else np.nan
        f10 = pos10 / n10 if n10 else np.nan
        testable = n9 >= MIN_HALF and n10 >= MIN_HALF
        p = (fisher_exact([[pos9, n9 - pos9], [pos10, n10 - pos10]])[1]
             if testable else np.nan)
        hi = ("GEX09" if (f9 > f10) else "GEX10") if (n9 and n10 and f9 != f10) else \
             ("GEX09" if (n9 and not n10) else "GEX10" if (n10 and not n9) else "")
        dl = designed_left.get(donor, "ambiguous")
        conc = ("yes" if hi == dl else "no") if dl in ("GEX09", "GEX10") and hi else ""
        rows.append(dict(donor=donor, stage=st, gfp=gf, n_GEX09=n9, n_GEX10=n10,
                         posfrac_GEX09=round(f9, 4) if n9 else np.nan,
                         posfrac_GEX10=round(f10, 4) if n10 else np.nan,
                         pitx2_high_lane=hi, designed_LEFT=dl,
                         testable=testable, fisher_p=p,
                         significant=bool(testable and p < 0.05),
                         concordant=conc))
    return pd.DataFrame(rows)


def main():
    pc = pd.read_csv(PERCELL, sep="\t")
    # Group B cS-split doublet/unassigned cells are still in the committed pre-fix cAM (as
    # donor11a/donor2a); drop them via the committed groupB_split_excluded_cells.tsv (0 hits once
    # cAM is regenerated with the fixed _groupB_resolved12 overlay).
    pc = drop_split_excluded(pc, "cAM_percell_populations", log=log)
    assert "PITX2_cp10k" in pc.columns, "cAM per-cell TSV lacks PITX2_cp10k"
    assert "population" in pc.columns, "cAM per-cell TSV lacks population"
    pc["group"] = pc["group"].astype(str)
    pc["donor"] = pc["donor"].astype(str)
    pc["lane"] = (pc["library"] if "library" in pc.columns else pc["lane"]).map(lane_of)
    # NB: PITX2_cp10k holds log1p(CP10K) (cAL exports the .X log-norm value despite the name);
    # detection (>0) is transform-invariant, but the linear means below must expm1 first (the
    # np.expm1(...).mean() in the Panel F loop).
    pc["PITX2"] = pc["PITX2_cp10k"].astype(float)

    # transparency: which observed populations fall outside CHICK_POOLS (would silently drop)
    validate_against_annotations(ANNOT, log=log)

    align = pd.read_csv(STAGE_TSV, sep="\t")
    align["group"] = align["group"].astype(str)
    align["donor"] = align["donor"].astype(str)
    stage_map = {(r.group, r.donor): r.stage for r in align.itertuples()}
    gfp_map = {(r.group, r.donor): r.gfp for r in align.itertuples()}

    # every (group, donor) must be in the stage-alignment SSOT: no per-cell stage or "WT"
    # default (the Group B directional control uses gfp to pick the designed-LEFT lane)
    unmapped = sorted({(g, d) for g, d in zip(pc["group"], pc["donor"]) if (g, d) not in stage_map})
    if unmapped:
        raise SystemExit(f"[align] (group, donor) absent from {os.path.basename(STAGE_TSV)}: {unmapped}")
    pc["stage"] = pc.apply(lambda r: stage_map[(r["group"], r["donor"])], axis=1)
    pc["gfp"] = pc.apply(lambda r: gfp_map[(r["group"], r["donor"])], axis=1)
    bad_stage = pc.loc[~pc["stage"].isin(STAGES), "stage"]
    pc = pc[pc["stage"].isin(STAGES)].copy()
    if len(bad_stage):                       # signature silent-drop guard — never drop unaudited
        log(f"[stage-filter] dropped {len(bad_stage)} cells with stage not in {STAGES}: "
            f"{bad_stage.value_counts().to_dict()}")
    pc["domain"] = [facet_of(s, p) for s, p in zip(pc["stage"], pc["population"])]

    # arbitrary orientation: lower lane (ARC01 / GEX09) = half1 (the duck uses arbitrary Library A/B)
    pc["half"] = np.where(pc["lane"].isin(("ARC01", "GEX09")), "h1", "h2")

    # ── directional control (Group B, inference pool) ─────────────────────────
    design = pd.read_csv(DESIGN, sep="\t")
    pcB = pc[pc.group == "B"].copy()
    designed_left = build_designed_left(design, pcB)
    direc = directional_control(pcB, designed_left)
    direc.to_csv(os.path.join(OUTDIR, "cAN_directional_perembryo.tsv"), sep="\t", index=False)

    # significance + orientability summary per (stage × gfp); headline = orientable & significant
    drows = []
    for (st, gf), g in direc.groupby(["stage", "gfp"]):
        orient = g[g.designed_LEFT.isin(("GEX09", "GEX10"))]
        sig_or = orient[orient.significant]
        drows.append(dict(stage=st, gfp=gf, n_donors=len(g),
                          n_orientable=len(orient),
                          n_sig=int(g.significant.sum()),
                          n_sig_orientable=len(sig_or),
                          n_sig_concordant=int((sig_or.concordant == "yes").sum())))
    dsum = pd.DataFrame(drows)
    dsum.to_csv(os.path.join(OUTDIR, "cAN_directional_summary.tsv"), sep="\t", index=False)
    tot_so = int(dsum.n_sig_orientable.sum()); tot_c = int(dsum.n_sig_concordant.sum())
    log("\n=== directional control (Group B, inference pool, firm 3-stage) ===")
    log(dsum.to_string(index=False))
    log(f"  HEADLINE: {tot_c}/{tot_so} orientable + significantly-lateralised embryos are "
        f"PITX2-high in the designed-LEFT half (GFP+ counted identically to WT)")

    # ── orientation-free Panel F (per embryo × facet half-vs-half) ────────────
    pcf = pc[pc["domain"].notna()].copy()
    rows, diag = [], []
    for (group, donor, domain), g in pcf.groupby(["group", "donor", "domain"], observed=True):
        g1, g2 = g[g.half == "h1"], g[g.half == "h2"]
        excluded = group == "B" and donor in EXCLUDE_GROUPB_DONORS
        kept = (not excluded) and len(g1) >= MIN_HALF and len(g2) >= MIN_HALF
        reason = ("data-quality (donor8)" if excluded else "" if kept
                  else f"n<{MIN_HALF} in a half (h1={len(g1)},h2={len(g2)})")
        diag.append(dict(group=group, donor=donor, stage=stage_map[(group, donor)],
                         gfp=gfp_map[(group, donor)], domain=domain,
                         n1=len(g1), n2=len(g2), kept=kept, drop_reason=reason))
        if not kept:
            continue
        m1, m2 = float(np.expm1(g1.PITX2).mean()), float(np.expm1(g2.PITX2).mean())
        f1, f2 = float((g1.PITX2 > 0).mean()), float((g2.PITX2 > 0).mean())
        # design concordance for the Panel F label: in THIS facet, is the higher-frac half the
        # designed-LEFT half?  (designed_LEFT is the embryo's orientation from the inference pool.)
        dl = designed_left.get(donor, "ambiguous") if group == "B" else "groupA-na"
        hi_half_lane = ("ARC01" if group == "A" else "GEX09") if f1 > f2 else \
                       ("ARC02" if group == "A" else "GEX10")
        conc = ("yes" if hi_half_lane == dl else "no") if dl in ("GEX09", "GEX10") and f1 != f2 else ""
        rows.append(dict(
            group=group, donor=donor, stage=stage_map[(group, donor)],
            gfp=gfp_map[(group, donor)], domain=domain,
            half1_lane=("ARC01" if group == "A" else "GEX09"),
            half2_lane=("ARC02" if group == "A" else "GEX10"),
            designed_LEFT=dl, concordant=conc,
            mean_h1=round(m1, 4), mean_h2=round(m2, 4),
            frac_h1=round(f1, 4), frac_h2=round(f2, 4),
            n1=len(g1), n2=len(g2), min_cells=int(min(len(g1), len(g2))),
            abs_log2_mean=round(float(abs(np.log2((m1 + .03) / (m2 + .03)))), 3)))
    out = pd.DataFrame(rows)
    out["domain"] = pd.Categorical(out["domain"], FACET_ORDER, ordered=True)
    out = out.sort_values(["domain", "group", "stage", "donor"])
    out.to_csv(os.path.join(OUTDIR, "cAN_halfvshalf.tsv"), sep="\t", index=False)
    log(f"\n[cAN] {len(out)} embryo×domain points "
        f"(A={int((out.group=='A').sum())}, B={int((out.group=='B').sum())}; "
        f"GFP+={int((out.gfp=='GFP+').sum())}, WT={int((out.gfp=='WT').sum())})")

    dg = pd.DataFrame(diag)
    dg["domain"] = pd.Categorical(dg["domain"], FACET_ORDER, ordered=True)
    dg = dg.sort_values(["domain", "group", "stage", "donor"])
    dg.to_csv(os.path.join(OUTDIR, "cAN_domain_counts.tsv"), sep="\t", index=False)
    log("\n=== dot accounting: kept / total per domain×group (dropped = sparse half) ===")
    acc = (dg.groupby(["domain", "group"], observed=True)
           .agg(total=("donor", "nunique"), kept=("kept", "sum")).reset_index())
    log(acc.to_string(index=False))
    drops = dg[~dg["kept"] & (dg["drop_reason"] != "data-quality (donor8)")]
    if len(drops):
        log("\ndropped (sparse) embryo×domain:")
        log(drops[["group", "donor", "stage", "domain", "n1", "n2"]].to_string(index=False))


if __name__ == "__main__":
    main()
