#!/usr/bin/env python3
"""build_orientfree_isoform_data — per-embryo ORIENTATION-FREE half-vs-half divergence across the
four fig04 readouts (ATAC PITX2c / PITX2a / PITX2b promoters + total PITX2 RNA), duck & chick.

Data builder for fig04 panel B (the orientation-free main panel; renderer =
figures/fig04/panel_B_xspecies_orientfree.py reads the TSV this writes).  Mirrors the isoform-panel
split: the per-cell-dependent COMPUTE lives here in pipeline/chick (next to the cross_species_halfvshalf engine);
the figure script is a thin TSV consumer.

For each readout, per embryo, we ask whether its two physical halves diverge MORE than its own
within-embryo cell-shuffle null — the L/R label DISCARDED, so the magnitude is non-circular for
BOTH species and EVERY embryo contributes (incl. the orientation-ambiguous chick Group-B donors):

    gap = obs( half-A vs half-B divergence )  −  median( within-embryo shuffle null )

ATAC divergence = null-corrected |log2(FRiP_a/FRiP_b)| at the pooled promoter peaks (ratio metric);
RNA  divergence = null-corrected |mean(PITX2_cp10k)_a − _b| (a level → a mean difference).  emp_p =
fraction of the per-embryo null ≥ obs (scale-free rank test; the figure-relevant DUCK-vs-CHICK
contrast is a one-sided Fisher on the per-embryo divergent fraction emp_p<0.05, computed at render).

Output: figures/data/chick_control/xspecies_orientfree_isoform.tsv
  cols: modality, isoform, species, genotype, stage, n_min, obs, null, gap, emp_p

  conda run -n duckfig python3 pipeline/chick/build_orientfree_isoform_data.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))            # pipeline/chick
sys.path.insert(0, HERE)   # cross_species_halfvshalf / _chick_subcluster_pools are siblings
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "figures"))
from cross_species_halfvshalf import load_species, halfvshalf, facet_of, FACET_VENTRAL  # noqa: E402
from _chick_subcluster_pools import CHICK_POOLS  # noqa: E402
from _groupB_split_excluded import drop_split_excluded  # noqa: E402
from _shared.subcluster_map import POOLS as DUCK_POOLS, get_leiden_col, get_pool_map  # noqa: E402
from _shared.lr_firmness import normalize_firmness, same_side_duplicates  # noqa: E402

DATA = os.path.join(REPO, "figures", "data")
DOC = os.path.join(DATA, "chick_control")

# isoform promoter windows (1-based incl.) — MUST match references/{,chick_}pitx2_promoter_regions.bed
# and build_isoform_data.py ISO_PROM (a known duplication; keep in sync).
ISO_PROM = {
    "duck": {"PITX2a": ("NC_051775.1", 3_961_086, 3_961_687),
             "PITX2b": ("NC_051775.1", 3_963_952, 3_964_553),
             "PITX2c": ("NC_051775.1", 3_971_084, 3_971_685)},
    "chick": {"PITX2a": ("chr4", 57_784_731, 57_785_332),
              "PITX2b": ("chr4", 57_787_516, 57_788_117),
              "PITX2c": ("chr4", 57_794_528, 57_795_129)},
}
LC = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1", "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}
D2G_FILE = {"HH13": "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
            "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
            "HH21": "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
            "HH25": "d2g_lpm_subclustering_HH25_ccreg.tsv.gz"}
READOUTS = [("ATAC", "PITX2c"), ("ATAC", "PITX2a"), ("ATAC", "PITX2b"), ("RNA", "total")]
# Per-half cell floor = 5 (DELIBERATELY below the floor-20 census used by the oriented panel G /
# panel A): the orientation-free metric is null-CORRECTED per embryo (gap = obs − own cell-shuffle
# null), which self-calibrates for small n — a 5-cell half has a wide null so it rarely clears
# emp_p<0.05 — and a low floor maximises the "every embryo" inclusion that is this panel's whole
# point (esp. the small chick Group-B halves).  Point area ∝ min-half cells encodes the reliability.
MIN_CELLS = 5


def invalid_duck_genotypes():
    """Genotypes flagged d2l `pair_constraint_ok == False` = the G3_unfiltered demux duplicates that
    leak into HH17/16 RNA (the oriented sibling drops them via its d2l/d2n side gate).
    The pair-constraint is a library-PAIRING validity property, NOT an L/R-side property, so excluding
    it is orientation-INDEPENDENT and leaves the 'every (valid) embryo' premise intact — valid
    side-ambiguous embryos (weak firmness) are KEPT; only the pair-invalid duplicates are dropped.
    Verified: this set is exactly the 5 G3_unfiltered_donor* genotypes, disjoint from the valid
    side-ambiguous embryos (G2_genotype_09, g19) the orientation-free panel must keep."""
    d2l = pd.read_csv(os.path.join(DATA, "d2l_lr_inference.tsv"), sep="\t")
    normalize_firmness(d2l["firmness"], "d2l_lr_inference.tsv")   # unknown values raise
    return set(d2l[~d2l.pair_constraint_ok].genotype)


def same_side_duplicate_genotypes():
    """Same-side duplicates (lane-loading-record split tube: ONE half sequenced in two libraries).
    Their two lanes hold the same tissue, so they are never an orientation-free half-vs-half
    pair -- excluded by this explicit rule (figures/_shared/lr_firmness.py), independent of
    whether pair_constraint_ok also happens to flag them."""
    d2l = pd.read_csv(os.path.join(DATA, "d2l_lr_inference.tsv"), sep="\t")
    return set(same_side_duplicates(d2l, "d2l_lr_inference.tsv"))


def pool_iso(pc, coords, species, iso):
    ch, ps, pe = ISO_PROM[species][iso]
    overlap = coords[(coords.chrom == ch) & (coords.start < pe) & (coords.end > ps)].peak_id.tolist()
    hit = [p for p in overlap if p in pc.columns]
    if not hit:
        raise SystemExit(f"{species} {iso}: no locus peak overlaps {ch}:{ps}-{pe}")
    if len(hit) < len(overlap):   # partial drop: coords vs per-cell-matrix peak universes disagree
        print(f"WARN {species} {iso}: {len(overlap) - len(hit)} overlapping peak(s) absent from the "
              f"per-cell matrix — pooling only the {len(hit)} present")
    pc = pc.copy()
    pc[f"{iso}_pool"] = pc[hit].sum(axis=1).to_numpy(float)
    return pc


def mean_halfvshalf(df, val, half, key_cols, min_cells=MIN_CELLS, k_null=2000, seed=0):
    """Orientation-free |mean_a − mean_b| with a within-embryo cell-shuffle null (RNA level).

    Shares the null convention of cross_species_halfvshalf.halfvshalf (the ATAC engine) — same
    default_rng(0) seed, k_null=2000, pp[:na]/pp[na:] split, emp_p=(#null≥obs + 1)/(k_null+1) — so
    the RNA and ATAC emp_p are directly comparable; only the per-half statistic differs (column-mean
    here vs depth-normalised FRiP there).  Keep the two in lockstep if either changes (a shared `shuffle_null` helper is the cleaner long-term shape — deferred to a follow-up)."""
    rng = np.random.default_rng(seed)
    rows = []
    for kv, d in df.groupby(key_cols, observed=True):
        halves = sorted(d[half].unique())
        if len(halves) != 2:
            continue
        a, b = halves
        da, db = d[d[half] == a], d[d[half] == b]
        if len(da) < min_cells or len(db) < min_cells:
            continue
        obs = abs(da[val].mean() - db[val].mean())
        allv = d[val].to_numpy(float)
        na = len(da)
        idx = np.arange(len(d))
        nulls = np.array([abs(allv[(pp := rng.permutation(idx))[:na]].mean() - allv[pp[na:]].mean())
                          for _ in range(k_null)])
        emp = (np.sum(nulls >= obs) + 1) / (k_null + 1)
        # join the full group key so chick "group:donor" is unique (donor IDs reused across A/B);
        # duck keys on genotype alone (already globally unique) → single-element join is a no-op.
        geno = ":".join(str(x) for x in kv) if isinstance(kv, tuple) else str(kv)
        rows.append(dict(genotype=geno, obs=obs, null=float(np.median(nulls)),
                         gap=obs - float(np.median(nulls)), emp_p=emp, n_min=min(len(da), len(db))))
    return pd.DataFrame(rows)


def atac_readout(species, pools, pc, coords, iso):
    pc2 = pool_iso(pc, coords, species, iso)
    s = halfvshalf(pc2, f"{iso}_pool", FACET_VENTRAL, pools, min_cells=MIN_CELLS, metric="ratio")
    if not len(s):
        return s
    s = s.rename(columns={"abs_log2_div": "obs", "null_median_div": "null"}).copy()
    s["gap"] = s.obs - s.null
    s["n_min"] = np.minimum(s.n_a, s.n_b)
    if species == "chick":   # chick donor IDs are reused across Group A/B → key by "group:donor"
        # Strip any existing prefix before re-applying, so this is IDEMPOTENT: since
        # the chick locus artifact already carries the group-qualified
        # id (A:donorN), and blindly prefixing produced "A:A:donorN".  Same bare-donor
        # idiom as build_isoform_data.py and isoform_probe/chick_atac_isoform_promoters.py.
        s["genotype"] = (s["group"].astype(str) + ":"
                         + s["genotype"].astype(str).str.split(":").str[-1])
    return s[["genotype", "stage", "obs", "null", "gap", "emp_p", "n_min"]]


def duck_rna(invalid):
    out = []
    for st, fn in D2G_FILE.items():
        df = pd.read_csv(os.path.join(DATA, fn), sep="\t")
        df = df[df[get_leiden_col(st)].astype(str).isin(set(get_pool_map(st, "inference")))].copy()
        df = df[~df.genotype.isin(invalid)]   # drop the d2l pair-invalid G3_unfiltered duplicates
        df["lm"] = df.library.map(lambda x: LC.get(x, x))
        s = mean_halfvshalf(df, "PITX2_cp10k", "lm", ["genotype"])
        if len(s):
            s["stage"] = st
            out.append(s)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def chick_rna():
    cam = pd.read_csv(os.path.join(DOC, "cAM_percell_populations.tsv.gz"), sep="\t")
    cam = drop_split_excluded(cam, "cAM_percell_populations")   # Group B split doublet/unassigned
    cam = cam[[facet_of(s, p, CHICK_POOLS) == FACET_VENTRAL
               for s, p in zip(cam.stage, cam.population)]].copy()
    if cam.empty:
        raise SystemExit("chick_rna: 0 cells folded to the ventral pool — CHICK_POOLS vs the cAM "
                         "population labels disagree (a chick subcluster rerun may have renamed them)")
    out = []
    for st, d in cam.groupby("stage", observed=True):
        s = mean_halfvshalf(d, "PITX2_cp10k", "lane", ["group", "donor"])
        if len(s):
            s["stage"] = st
            out.append(s)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    duck_pc, duck_coords = load_species("duck", DUCK_POOLS)
    chick_pc, chick_coords = load_species("chick", CHICK_POOLS)

    # drop the d2l pair-constraint-invalid duck genotypes (the G3_unfiltered leak) from BOTH
    # modalities so duck ATAC and duck RNA see the same valid embryo set (ATAC was only incidentally
    # clean via the 2-lane gate; make it explicit).
    pc_invalid = invalid_duck_genotypes()
    dups = same_side_duplicate_genotypes()
    invalid = pc_invalid | dups
    excluded = sorted(invalid & set(duck_pc.genotype))
    duck_pc = duck_pc[~duck_pc.genotype.isin(invalid)].copy()
    print(f"[duck] excluded {len(excluded)} genotype(s) present in the locus extract "
          f"(d2l pair_constraint_ok==False, or same-side duplicate): {excluded}")
    print(f"[duck] same-side duplicates excluded from every half-vs-half pair: {sorted(dups)}")

    data = {}
    for iso in ("PITX2c", "PITX2a", "PITX2b"):
        data[("ATAC", iso)] = {
            "duck": atac_readout("duck", DUCK_POOLS, duck_pc, duck_coords, iso),
            "chick": atac_readout("chick", CHICK_POOLS, chick_pc, chick_coords, iso)}
    data[("RNA", "total")] = {"duck": duck_rna(invalid), "chick": chick_rna()}

    rows = []
    print("\n=== divergent-embryo fraction (emp_p<0.05) per readout ===")
    for key in READOUTS:
        for sp in ("duck", "chick"):
            for _, r in data[key][sp].iterrows():
                rows.append(dict(modality=key[0], isoform=key[1], species=sp, genotype=r.genotype,
                                 stage=r.stage, n_min=int(r.n_min), obs=round(r.obs, 4),
                                 null=round(r.null, 4), gap=round(r.gap, 4), emp_p=round(r.emp_p, 4)))
        dfr, dfc = data[key]["duck"], data[key]["chick"]
        print(f"  {key[0]}·{key[1]:6s}  duck {float((dfr.emp_p < 0.05).mean()):.2f} (n={len(dfr)}) | "
              f"chick {float((dfc.emp_p < 0.05).mean()):.2f} (n={len(dfc)})")
    emb = pd.DataFrame(rows)
    os.makedirs(DOC, exist_ok=True)
    tsv = os.path.join(DOC, "xspecies_orientfree_isoform.tsv")
    emb.to_csv(tsv, sep="\t", index=False)
    print(f"\nwrote {os.path.relpath(tsv, REPO)}: {emb.shape}")


if __name__ == "__main__":
    main()
