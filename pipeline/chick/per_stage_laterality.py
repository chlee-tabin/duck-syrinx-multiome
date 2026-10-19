#!/usr/bin/env python3
"""Per-stage PITX2 lateral divergence — TWO metrics × TWO modalities × TWO species.

Per the discussion: the chick has no external L/R for most embryos, so the ABSOLUTE
(orientation-free) |log2(L/R)| is the chick's BEST CASE (counts divergence either way → can
only over-state laterality). Where the GFP/WT design fixes orientation we ALSO report the
SIGNED log2(L/R) (real direction) WITHOUT subsetting the absolute analysis. Per stage we ask:
  - ABSOLUTE  : per-embryo |log2 half-ratio| vs the within-embryo shuffle null (ATAC) / magnitude (RNA)
  - SIGNED    : among design-orientable embryos, do they point the SAME way (sign test vs 50/50)?

Two model checks:
  (1) CHICK — is the ventral pool non-lateralized at EVERY stage (incl. HH13), in both modalities?
  (2) DUCK  — does RNA laterality DAMPEN at HH17/16 while ATAC laterality PERSISTS (chromatin memory)?

ATAC: figures/data/{duck,chick}_locus per-cell (side L/R = design/d2l-anchored; duck HH17/16
      re-mapped onto d2n via apply_d2n_hh1716_override; lanes = halves).
RNA : duck = canonical d2g pools (figures/data/d2g_lpm_subclustering_*_ccreg.tsv.gz, fig02's
      grouping — the res=0.5 'Splanchnic core/transitional' tiers are DROPPED); chick =
      figures/data/chick_control/cAN_halfvshalf.tsv (chick d2g facets, h1/h2 + designed_LEFT).

  conda run -n duckfig python3 pipeline/chick/per_stage_laterality.py
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binomtest, wilcoxon

L_COLOR, R_COLOR = "#2c5fa8", "#c0392b"  # project rule: Left = blue, Right = red


def make_panel(dk_sgn, ck_sgn, rdk, rck, out_path):
    """Per-stage SIGNED PITX2 L/R divergence (d2g grouping): ATAC (top) + RNA (bottom),
    duck (filled) vs chick (open), each point an embryo colored by direction (L>R blue, R>L red)."""
    cols = ["HH13", "mid", "HH21", "HH25"]
    xlab = ["HH13", "HH17/16\n(chick HH15-16)", "HH21\n(duck only)", "HH25"]
    smap = {"HH13": "HH13", "HH17/16": "mid", "HH15-16": "mid", "HH21": "HH21", "HH25": "HH25"}
    rng = np.random.default_rng(0)

    def strip(ax, df, col, scale, filled):
        for _, r in df.iterrows():
            c = smap.get(str(r["stage"]))
            if c not in cols or not np.isfinite(r[col]):
                continue
            x = cols.index(c) + (-0.17 if filled else 0.17) + rng.uniform(-0.04, 0.04)
            v = r[col] * scale
            col_dir = L_COLOR if v >= 0 else R_COLOR
            ax.scatter(x, v, marker=("o" if filled else "s"), s=36, zorder=3, linewidth=1.0,
                       facecolor=(col_dir if filled else "none"), edgecolor=col_dir)

    fig, axes = plt.subplots(2, 1, figsize=(7, 7.2), sharex=True)
    strip(axes[0], dk_sgn, "signed_diff_L_R", 1e5, True)
    strip(axes[0], ck_sgn, "signed_diff_L_R", 1e5, False)
    axes[0].set_ylabel("ATAC  signed ΔFRiP (×1e-5)\n(+ = L>R,  − = R>L)")
    strip(axes[1], rdk, "signed", 1.0, True)
    strip(axes[1], rck, "signed", 1.0, False)
    axes[1].set_ylabel("RNA  signed Δcp10k\n(+ = L>R,  − = R>L)")
    for ax in axes:
        ax.axhline(0, color="0.5", ls="--", lw=0.8, zorder=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(xlab, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    hh = [plt.Line2D([], [], marker="o", color="k", ls="", mfc="k", label="duck"),
          plt.Line2D([], [], marker="s", color="k", ls="", mfc="none", label="chick"),
          plt.Line2D([], [], marker="o", color=L_COLOR, ls="", mfc=L_COLOR, label="L>R"),
          plt.Line2D([], [], marker="o", color=R_COLOR, ls="", mfc=R_COLOR, label="R>L")]
    axes[0].legend(handles=hh, fontsize=6.5, frameon=False, loc="upper right", ncol=2)
    fig.suptitle("Per-stage signed PITX2 L/R (d2g; HH17/16 on chromatin/d2n sides) — duck ATAC "
                 "L>R every stage;\nRNA flips at the HH17/16 trough (snRNA unreliable there); chick flat/sparse",
                 fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
HERE = os.path.dirname(os.path.abspath(__file__))            # pipeline/chick
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)   # cross_species_halfvshalf / _chick_subcluster_pools are siblings
from cross_species_halfvshalf import load_species, halfvshalf, facet_of, lane_of, FACET_VENTRAL  # noqa
from _chick_subcluster_pools import CHICK_POOLS  # noqa
sys.path.insert(0, os.path.join(REPO, "figures"))
from _shared.subcluster_map import POOLS as DUCK_POOLS  # noqa
from _shared.library_naming import collapse_lib  # noqa
from _shared.lr_firmness import resolve_sides, same_side_duplicates  # noqa
EPS = 1e-9
DATA = os.path.join(REPO, "figures", "data")
DOC = os.path.join(DATA, "chick_control")
# diagnostic-panel home: docs/chick_control/ (not distributed; no shipped figure reads the
# diagnostic panel). The directory is created on demand.
XS = os.path.join(REPO, "docs", "chick_control")
os.makedirs(XS, exist_ok=True)


def duck_d2n_hh1716_sides():
    """HH17/16 chromatin (d2n) L/R sides — AUTHORITATIVE at the trough.

    At HH17/16 the snRNA (d2l) side-call is unreliable (PITX2c at the noise floor;
    total-PITX2 is a/b-dominated and bilateral), and it fails the experimental-design
    pair-constraint (figS02-D = 5 reds).  The chromatin call (d2n) satisfies the design
    (fig03-E = 0 reds) and is adopted as authoritative — fig03 panel E overrides snRNA
    with it.  We mirror that here so the cross-species ATAC/RNA panel is self-consistent
    with fig03 at the trough.  Returns {(genotype, library_m): 'L'/'R'}.
    """
    d = pd.read_csv(os.path.join(DATA, "d2n_chromatin_lr_inference_hh1716.tsv"), sep="\t")
    # design-fit (pair_constraint_ok) is the INDEPENDENT anchor that justifies adopting d2n here;
    # drop design-conflicting calls so we don't smuggle back the d2l sides we switched away from.
    # The firmness gate mirrors the d2l convention (direct/placement only, no singleton_rank),
    # so a weak singleton chromatin guess can never overwrite a firm d2l side.
    # G3_unfiltered_donor0/3/9 leaked into HH17/16 RNA.
    # DELIBERATELY reads the raw d2n TSV, not the canonical d2o merge (lr_inference_multiome.tsv):
    # the d2o merge at the canonical filter carries only 26 HH17/16 keys and OMITS the
    # g01/g08/g20 pairs present in the locus extract. Routing registry:
    # docs/laterality_sidecar_routing.md §3 (not distributed).
    # Shared firmness rule (figures/_shared/lr_firmness.py): chromatin_ prefix stripped, unknown
    # values raise; a same-side duplicate (lane-loading-record split tube, design side, not gated by
    # pair_constraint_ok) is keyed on its POOLED lane label -- the lane the locus extract gives it.
    d = resolve_sides(d, {"direct", "placement"}, source="d2n_chromatin_lr_inference_hh1716.tsv",
                      require_pair_constraint=True)
    return {(g, k): s for g, k, s in zip(d.genotype, d.library_key, d.inferred_side)}


def apply_d2n_hh1716_override(pc, d2n=None):
    """Re-map a DUCK per-cell frame's `side` at the HH17/16 trough onto the chromatin
    (d2n) call — the raw-d2l side baked into locus_percell_counts is circular/stale
    there (the laterality-sidecar routing rule).  Cells whose (genotype, lane) has no
    design-consistent d2n call become "?" and drop out of any L/R split (every
    consumer filters side.isin(["L","R"]); the print below makes the drop visible).

    On the current committed extract this is a NO-OP: the d2l pair-constraint gate in
    duck_pitx2_locus_percell_extract.py already excludes every genotype whose d2l and
    d2n sides disagree, so all baked HH17/16 sides equal d2n.  The remap guards against
    future divergence; it did not produce the current numbers.

    Duck only — chick sides are design-derived and must not pass through here (raises
    on the chick stage label).  Fails loud when the remap resolves ZERO HH17/16 cells:
    that is a key-alphabet/plumbing failure, never a legitimate design exclusion.
    """
    if (pc["stage"] == "HH15-16").any():
        raise ValueError("apply_d2n_hh1716_override is duck-only; got chick stage 'HH15-16'")
    d2n = d2n if d2n is not None else duck_d2n_hh1716_sides()
    pc = pc.copy()
    m = (pc["stage"] == "HH17/16").to_numpy()
    old = list(pc.loc[m, "side"])  # plain-list snapshot — .to_numpy() can alias the
    # frame's own buffer and mutate under the assignment below, zeroing the flip count
    new = [d2n.get((g, collapse_lib(lane_of(c))), "?")
           for g, c in zip(pc.loc[m, "genotype"], pc.loc[m, "cell"])]
    pc.loc[m, "side"] = new
    n_hh = int(m.sum())
    if n_hh:
        n_miss = sum(1 for s in new if s == "?")
        if n_miss == n_hh:
            miss = sorted({(g, collapse_lib(lane_of(c)))
                           for g, c in zip(pc.loc[m, "genotype"], pc.loc[m, "cell"])})
            raise RuntimeError(
                f"d2n override resolved 0 of {n_hh} HH17/16 cells — key-alphabet mismatch "
                f"(stage label / lane / genotype drift), not a design exclusion. Keys: {miss}")
        n_flip = int((old != np.array(new, dtype=object)).sum())
        miss = sorted({(g, collapse_lib(lane_of(c)))
                       for g, c, s in zip(pc.loc[m, "genotype"], pc.loc[m, "cell"], new)
                       if s == "?"})
        print(f"[d2n override] HH17/16 {n_hh:,} cells: {n_flip:,} changed, {n_miss:,} -> '?'"
              + (f" (dropped keys: {miss})" if miss else ""))
    return pc


def sign_test(signs):
    """signs: array of +1/-1 (drop 0). Returns (n_pos, n, two-sided p vs 0.5)."""
    s = signs[signs != 0]
    if len(s) == 0:
        return 0, 0, np.nan
    npos = int((s > 0).sum())
    return npos, len(s), binomtest(npos, len(s), 0.5).pvalue


# ── ATAC: per-embryo signed ΔFRiP (L−R) + absolute |ΔFRiP| via halfvshalf(metric="diff") ──
def atac_per_stage(pc, pools, peak="PITX2c_prom_pooled", min_cells=10):
    # absolute magnitude = DIFFERENCE (self-gates small counts; no zero-inflation), vs shuffle null
    hv = halfvshalf(pc, peak, FACET_VENTRAL, pools, min_cells=5, metric="diff")
    # signed: per (genotype, stage) with side L/R — sign(ΔFRiP)==sign(log2 ratio), so DIRECTION
    # is metric-invariant; we keep ΔFRiP for unit-consistency with the absolute metric.
    # Duck HH17/16 sides are d2n-remapped UPSTREAM (main() applies apply_d2n_hh1716_override
    # to the whole frame); chick frames arrive on design/d2l sides untouched.
    keep = [facet_of(s, p, pools) == FACET_VENTRAL for s, p in zip(pc.stage, pc.population)]
    sub = pc[keep].copy()
    sub = sub[sub.side.isin(["L", "R"])]
    sgn = []
    for (g, st), d in sub.groupby(["genotype", "stage"], observed=True):
        L, R = d[d.side == "L"], d[d.side == "R"]
        if len(L) < min_cells or len(R) < min_cells:
            continue
        fL = L[peak].sum() / max(L.nFrags.sum(), EPS)
        fR = R[peak].sum() / max(R.nFrags.sum(), EPS)
        sgn.append({"genotype": g, "stage": st, "signed_diff_L_R": fL - fR, "fL": fL, "fR": fR,
                    "nL": len(L), "nR": len(R)})
    sgn = pd.DataFrame(sgn)
    return hv, sgn


def rna_duck_per_stage(d2n_side=None):
    """CANONICAL d2g grouping (fig02 pools) — the res=0.5 'Splanchnic core/+transitional' tiers
    are DROPPED (they fed the stale xs1 path; see duck_rna_d2g_per_stage.py for the standalone).
    Per genotype × stage in the d2g inference (Ventral splanchnic) pool: PITX2_cp10k DIFFERENCE
    between the two library-halves (abs) + signed (mean_L − mean_R via the side lookup).

    Side source: d2l (snRNA) for HH13/HH21/HH25; chromatin (d2n) for HH17/16 when provided,
    so the RNA is oriented on the SAME authoritative sides as the ATAC.  On d2n sides the RNA
    flips for the trough-misled (a/b-dominated) embryos — the intended chromatin-value exhibit."""
    from _shared.subcluster_map import get_leiden_col, get_pool_map
    stage_file = {"HH13": "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
                  "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
                  "HH21": "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
                  "HH25": "d2g_lpm_subclustering_HH25_ccreg.tsv.gz"}
    lc = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1",
          "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}
    d2l_raw = pd.read_csv(os.path.join(DATA, "d2l_lr_inference.tsv"), sep="\t")
    d2l = resolve_sides(d2l_raw, {"direct", "placement"}, source="d2l_lr_inference.tsv",
                        paired=True, require_pair_constraint=True)
    side = {(g, lm): s for g, lm, s in zip(d2l.genotype, d2l.library_m, d2l.inferred_side)}
    # a same-side duplicate's two libraries hold ONE half: never a library-pair difference
    dups = set(same_side_duplicates(d2l_raw, "d2l_lr_inference.tsv"))
    print(f"[rna_duck_per_stage] same-side duplicates excluded from the library pairs: "
          f"{sorted(dups)}")
    rows = []
    for st, fn in stage_file.items():
        sidemap = d2n_side if (d2n_side is not None and st == "HH17/16") else side
        df = pd.read_csv(os.path.join(DATA, fn), sep="\t")
        df["cl"] = df[get_leiden_col(st)].astype(str)
        df = df[df.cl.isin(set(get_pool_map(st, "inference").keys()))]
        df["lm"] = df.library.map(lambda x: lc.get(x, x))
        for g, d in df.groupby("genotype", observed=True):
            if g in dups:
                continue
            libs = sorted(d.lm.unique())
            if len(libs) != 2:
                continue
            a, b = libs
            da, db = d[d.lm == a], d[d.lm == b]
            if len(da) < 10 or len(db) < 10:
                continue
            ma, mb = da.PITX2_cp10k.mean(), db.PITX2_cp10k.mean()
            sa, sb = sidemap.get((g, a)), sidemap.get((g, b))
            signed = ((ma - mb) if sa == "L" else (mb - ma)) if {sa, sb} == {"L", "R"} else np.nan
            rows.append({"stage": st, "abs": abs(ma - mb), "signed": signed})
    return pd.DataFrame(rows)


def rna_chick_per_stage():
    df = pd.read_csv(os.path.join(DOC, "cAN_halfvshalf.tsv"), sep="\t")
    df = df[df.domain.isin(["Ventral splanchnic pool", "LPM mes. (HH13)"])].copy()
    df["abs"] = (df["mean_h1"] - df["mean_h2"]).abs()   # cp10k DIFFERENCE

    def signed(r):  # orient mean_h1/mean_h2 by designed_LEFT (Group B firm only)
        dl = str(r["designed_LEFT"])
        if dl == str(r["half1_lane"]):
            return r["mean_h1"] - r["mean_h2"]
        if dl == str(r["half2_lane"]):
            return r["mean_h2"] - r["mean_h1"]
        return np.nan
    df["signed"] = df.apply(signed, axis=1)
    return df


def report(label, abs_vals, signed_vals, null_vals=None, scale=1.0, unit=""):
    n = len(abs_vals)
    med = np.median(abs_vals) * scale if n else np.nan
    line = f"  {label:30s} n={n:2d} | |Δ|{unit} med={med:6.3f}"
    if null_vals is not None and n:
        gap = np.median(np.asarray(abs_vals) - np.asarray(null_vals)) * scale
        try:
            wp = wilcoxon(abs_vals, null_vals, alternative="greater").pvalue
        except ValueError:
            wp = np.nan
        line += f" | obs-null gap={gap:+6.3f} (p={wp:.3f})"
    sv = np.asarray([s for s in signed_vals if np.isfinite(s)])
    if len(sv):
        npos, ns, sp = sign_test(np.sign(sv))
        line += f" | SIGNED L>R {npos}/{ns} (p={sp:.3f})"
    else:
        line += " | SIGNED (none orientable)"
    print(line)


def main():
    (duck_pc, _), _ = (load_species("duck", DUCK_POOLS), None)
    (chick_pc, _), _ = (load_species("chick", CHICK_POOLS), None)
    duck_d2n = duck_d2n_hh1716_sides()  # HH17/16 chromatin-authoritative sides (duck only)
    duck_pc = apply_d2n_hh1716_override(duck_pc, d2n=duck_d2n)
    dk_hv, dk_sgn = atac_per_stage(duck_pc, DUCK_POOLS)
    ck_hv, ck_sgn = atac_per_stage(chick_pc, CHICK_POOLS)  # chick: no d2n (design/d2l sides)
    rdk, rck = rna_duck_per_stage(d2n_side=duck_d2n), rna_chick_per_stage()
    panel = os.path.join(XS, "per_stage_signed_two_metric.png")
    make_panel(dk_sgn, ck_sgn, rdk, rck, panel)
    print(f"wrote {os.path.relpath(panel, REPO)}")

    print("\n############ ATAC — Ventral splanchnic pool, per stage (metric = ΔFRiP, ×1e-5) ############")
    print("(absolute = orientation-free best case, DIFFERENCE not ratio; SIGNED = design-anchored direction)")
    for sp_name, hv, sgn, stages in [("DUCK", dk_hv, dk_sgn, ["HH13", "HH17/16", "HH21", "HH25"]),
                                     ("CHICK", ck_hv, ck_sgn, ["HH13", "HH15-16", "HH25"])]:
        print(f"\n{sp_name} ATAC:")
        for st in stages:
            h = hv[hv.stage == st]
            s = sgn[sgn.stage == st]["signed_diff_L_R"].to_numpy() if len(sgn) else np.array([])
            if not len(h) and not len(s):
                print(f"  {st:30s} (no embryos ≥5/half)"); continue
            report(st, h.abs_log2_div.to_numpy(), s, h.null_median_div.to_numpy(),
                   scale=1e5, unit="(×1e-5)")

    # Q3: name any duck HH17/16 minority-sign embryo (historically 11/12 pre-relabel;
    # 7/7 since the 2026-08 duck_locus relabel — the block is a regression tripwire now)
    o = dk_sgn[dk_sgn.stage == "HH17/16"]
    if len(o):
        rev = o[o.signed_diff_L_R < 0]
        print(f"\n  [Q3] DUCK HH17/16 signed: {(o.signed_diff_L_R>0).sum()}/{len(o)} L>R. "
              f"Minority (R>L): {list(rev.genotype)} "
              + (f"→ {rev.iloc[0].genotype}: nL={int(rev.iloc[0].nL)} nR={int(rev.iloc[0].nR)}, "
                 f"FRiP_L={rev.iloc[0].fL*1e5:.2f} vs FRiP_R={rev.iloc[0].fR*1e5:.2f} (×1e-5)"
                 if len(rev) else ""))

    print("\n############ RNA (snRNA PITX2) — splanchnic, per stage (metric = Δcp10k) ############")
    for sp_name, rdf, stages in [("DUCK", rdk, ["HH13", "HH17/16", "HH21", "HH25"]),
                                 ("CHICK", rck, ["HH13", "HH15-16", "HH25"])]:
        print(f"\n{sp_name} RNA:")
        for st in stages:
            d = rdf[rdf.stage == st]
            if not len(d):
                print(f"  {st:30s} (no embryos)"); continue
            report(st, d["abs"].to_numpy(), d["signed"].to_numpy(), unit="(cp10k)")

    print("\nMODEL CHECKS:")
    print("  (1) CHICK lateralized at any stage? → SIGNED L>R sign-test p<0.05 AND obs-null gap>0.")
    print("  (2) DUCK chromatin memory? → ATAC obs-null gap + signed PERSIST at HH17/16 where RNA is noisier.")
    print("  NOTE: duck HH13 RNA is cp10k (depth-erases the Wave-1 peak; absolute counts needed to see it).")


if __name__ == "__main__":
    main()
