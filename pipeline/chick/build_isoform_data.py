#!/usr/bin/env python3
"""Consolidated cross-species PITX2 ISOFORM x MODALITY figure data.

Per (modality, isoform, stage, group, embryo): signed L-R divergence + n_cells (min half).
  ATAC  PITX2c   : in_promoter_c pooled FRiP                 -- asymmetric isoform (DIRECTION / laterality)
  ATAC  PITX2a   : PITX2a promoter FRiP                      -- BMP/SMAD5 cis-anchor; expected BILATERAL
  ATAC  PITX2b   : PITX2b promoter FRiP                      -- bystander; expected BILATERAL (within-locus ctrl)
  RNA   total    : total PITX2 cp10k                         -- FISH-comparable LEVEL (two-wave/HH17 trough)
  RNA   PITX2a/b : chick a/b-body (cZ2) / duck d2j pitx2ab   -- DROPPED from the figure (see main())

Duck sides: d2l (snRNA) at HH13/HH21/HH25, chromatin d2n at the HH17/16 trough (authoritative,
fig03/POC).  Chick: design-legal orientation envelope (best/worst), NEVER a cherry-picked side.
Same d2g grouping (get_pool_map / CHICK_POOLS facet); DIFFERENCE not ratio; no silent fallbacks
(a missing peak / side / pool key raises).

  conda run -n duckfig python3 pipeline/chick/build_isoform_data.py
Output: figures/data/chick_control/xspecies_pitx2_isoform_laterality.tsv
        (rendered by figures/figS04/panel_G_xspecies_isoform.R = the ORIENTED figS04 panel G; the
         orientation-free fig04 panel B reads the separate xspecies_orientfree_isoform.tsv; also
         poc/isoform_beeswarm.R)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))            # pipeline/chick
REPO = os.path.dirname(os.path.dirname(HERE))
# cross_species_halfvshalf / per_stage_laterality / _chick_subcluster_pools are siblings
for p in (HERE, os.path.join(REPO, "figures")):
    sys.path.insert(0, p)
from cross_species_halfvshalf import load_species, facet_of, lane_of, FACET_VENTRAL  # noqa: E402
from _chick_subcluster_pools import CHICK_POOLS  # noqa: E402
from _shared.subcluster_map import get_leiden_col, get_pool_map  # noqa: E402
from _shared.subcluster_map import POOLS as DUCK_POOLS  # noqa: E402
from _shared.lr_firmness import resolve_sides, same_side_duplicates  # noqa: E402
import per_stage_laterality as P  # noqa: E402

DATA = os.path.join(REPO, "figures", "data")
DOC = os.path.join(DATA, "chick_control")
DES = os.path.join(REPO, "design")
CAM = os.path.join(DOC, "cAM_percell_populations.tsv.gz")                       # chick per-cell meta
EPS = 1e-9
MINC = int(os.environ.get("MINC", "20"))                    # >=N cells/half/embryo. DEFAULT 20 =
#   floor-20, matching the committed figS04 panel G TSV + LEGEND (a plain regen reproduces the
#   oriented panel this builder feeds). NB the orientation-free fig04 panel B is a SEPARATE builder
#   (build_orientfree_isoform_data.py) at floor 5. Override MINC=10 only for the relaxed power sweep.
SD = {"HH13": "HH13", "HH15-16": "HH17/16", "HH25": "HH25"}  # chick internal -> figure stage
DUCK_STAGES = ["HH13", "HH17/16", "HH21", "HH25"]
# a and b shown SEPARATELY (not pooled).  The chick galGal6a GTF annotates all three isoforms
# (XM_025149516.1=PITX2a, XM_025149515.1=PITX2b, NM_205010.1=PITX2c) and the ATAC peaks resolve
# a vs b, so the earlier "PITX2a/b" pooling was a presentation choice, not an annotation limit.
# a = the BMP/SMAD5 cis-anchor (mechanistically distinct, ties to the a-vs-c promoter model);
# b = bystander.  Both expected BILATERAL (L≈R) vs the asymmetric PITX2c.
ISO_COLS = [("PITX2c", "PITX2c_prom"), ("PITX2a", "PITX2a_prom"), ("PITX2b", "PITX2b_prom")]

# isoform promoter windows (1-based inclusive); the c-overlap is asserted == in_promoter_c.
# CANONICAL SOURCE = references/pitx2_promoter_regions.bed (duck) + references/chick_pitx2_promoter_regions.bed
# (chick), promoter|PITX2{a,b,c} rows. These literals MUST match those BEDs — duplication is a known
# drift risk; a follow-up should parse the BEDs here and in isoform_probe/.
ISO_PROM = {
    "duck": {"PITX2a": ("NC_051775.1", 3_961_086, 3_961_687),
             "PITX2b": ("NC_051775.1", 3_963_952, 3_964_553),
             "PITX2c": ("NC_051775.1", 3_971_084, 3_971_685)},
    "chick": {"PITX2a": ("chr4", 57_784_731, 57_785_332),
              "PITX2b": ("chr4", 57_787_516, 57_788_117),
              "PITX2c": ("chr4", 57_794_528, 57_795_129)},
}


def map_iso_peaks(coords, promoters):
    """isoform -> overlapping union-peak ids; assert c == in_promoter_c (fail loud on drift)."""
    out = {}
    for iso, (ch, ps, pe) in promoters.items():
        hit = coords[(coords.chrom == ch) & (coords.start < pe) & (coords.end > ps)]
        if hit.empty:
            raise RuntimeError(f"{iso}: no union peak overlaps {ch}:{ps}-{pe}")
        out[iso] = hit.peak_id.tolist()
    c_flag = set(coords.loc[coords.in_promoter_c, "peak_id"])
    if set(out["PITX2c"]) != c_flag:
        raise RuntimeError(f"PITX2c overlap {out['PITX2c']} != in_promoter_c {sorted(c_flag)}")
    return out


def add_iso_cols(pc, coords, species):
    peaks = map_iso_peaks(coords, ISO_PROM[species])
    pc = pc.copy()
    for iso, pl in peaks.items():
        miss = [p for p in pl if p not in pc.columns]
        if miss:
            raise RuntimeError(f"{species} {iso}: peak cols absent from per-cell matrix: {miss}")
        pc[f"{iso}_prom"] = pc[pl].sum(axis=1).to_numpy(float)
    return pc  # PITX2a_prom / PITX2b_prom / PITX2c_prom — a & b emitted separately (not pooled)


# ── duck ATAC: per-embryo signed ΔFRiP (sides already on pc, d2n @ HH17/16) ──
def duck_atac_signed(pc, iso_col, stage):
    keep = np.array([facet_of(s, p, DUCK_POOLS) == FACET_VENTRAL
                     for s, p in zip(pc.stage, pc.population)])
    sub = pc[keep & (pc.stage.to_numpy() == stage) & pc.side.isin(["L", "R"]).to_numpy()]
    out = []
    for g, d in sub.groupby("genotype", observed=True):
        L, R = d[d.side == "L"], d[d.side == "R"]
        if len(L) < MINC or len(R) < MINC:
            continue
        fL = L[iso_col].sum() / max(L.nFrags.sum(), EPS)
        fR = R[iso_col].sum() / max(R.nFrags.sum(), EPS)
        out.append(dict(embryo=g, signed=(fL - fR) * 1e5, n_cells=int(min(len(L), len(R)))))
    return out


# ── chick design-legal orientation envelope (design-table best/worst over legal lane->side maps) ──
def _ng(g):
    return "GFP" if str(g).strip() in ("GFP", "GFP+") else "WT"


ga = pd.read_csv(f"{DES}/chick_groupA_lr_design.tsv", sep="\t")
gb = pd.read_csv(f"{DES}/chick_groupB_lr_design.tsv", sep="\t")
gb["gfp_n"] = gb.gfp_status.map(_ng)
gaA = {st: list(s.arc01_side) for st, s in ga.groupby("stage")}
gbB = {(st, g): list(s.gex09_side) for (st, g), s in gb.groupby(["stage", "gfp_n"])}
can = pd.read_csv(f"{DOC}/cAN_halfvshalf.tsv", sep="\t")
can["gfp_n"] = can.gfp.map(_ng)
ven = can[can.domain.isin(["Ventral splanchnic pool", "LPM mes. (HH13)"])]


def chick_atac_embs(chick_pc, cst, iso_col):
    keep = [facet_of(s, p, CHICK_POOLS) == FACET_VENTRAL
            for s, p in zip(chick_pc.stage, chick_pc.population)]
    sub = chick_pc[[k and s == cst for k, s in zip(keep, chick_pc.stage)]]
    rost = {(r.group, r.donor): r.gfp_n for r in ven[ven.stage == cst].itertuples()}  # (group,donor) key: chick donor IDs are reused across groups (A:donorN != B:donorN)
    out = []
    for (grp, geno), e in sub.groupby(["group", "genotype"], observed=True):
        lanes = sorted(e.lane.unique())
        if len(lanes) != 2:
            continue
        a, b = lanes
        ea, eb = e[e.lane == a], e[e.lane == b]
        if len(ea) < MINC or len(eb) < MINC:
            continue
        fa = ea[iso_col].sum() / max(ea.nFrags.sum(), EPS)
        fb = eb[iso_col].sum() / max(eb.nFrags.sum(), EPS)
        out.append(dict(e=geno, grp=grp, gfp=rost[(grp, geno.split(":")[-1])],  # (group, bare-donor) key; raises if absent (no silent default)
                        delta=(fa - fb) * 1e5, n=int(min(len(ea), len(eb)))))
    return out


def chick_rna_total_embs(cst):
    out = []
    d = ven[(ven.stage == cst) & (ven.n1 >= MINC) & (ven.n2 >= MINC)]
    for r in d.itertuples():
        # group-qualified id, same alphabet as chick_atac_embs / the locus artifact
        out.append(dict(e=f"{r.group}:{r.donor}", grp=r.group, gfp=r.gfp_n,
                        delta=r.mean_h1 - r.mean_h2, n=int(min(r.n1, r.n2))))
    return out


def orientations(embs, cst):
    classes = []
    A = [e for e in embs if e["grp"] == "A"]
    if A:
        asg = sorted(set(itertools.permutations(gaA[cst], len(A))))
        classes.append((A, [[s == "L" for s in a] for a in asg]))
    for g in ["WT", "GFP"]:
        B = [e for e in embs if e["grp"] == "B" and e["gfp"] == g]
        if not B:
            continue
        asg = sorted(set(itertools.permutations(gbB[(cst, g)], len(B))))
        classes.append((B, [[s == "L" for s in a] for a in asg]))
    if not classes:
        return []
    oris = []
    for combo in itertools.product(*[a for _, a in classes]):
        sgn = {}
        for (el, _), ch in zip(classes, combo):
            for e, h1L in zip(el, ch):
                sgn[e["e"]] = e["delta"] if h1L else -e["delta"]
        oris.append(sgn)
    return oris


def chick_cases(embs, cst):
    """{'chick best','chick worst'} signed-dicts = the design-legal mean-signed extremes."""
    oris = orientations(embs, cst)
    if not oris:
        return {}
    means = [np.mean(list(o.values())) for o in oris]
    order = np.argsort(means)
    return {"chick worst": oris[order[0]], "chick best": oris[order[-1]]}


def main():
    rows = []

    # ── DUCK ATAC (c + a/b), d2n sides at HH17/16 ──
    duck_coords = pd.read_csv(f"{DATA}/duck_locus/locus_peak_coords.tsv", sep="\t")
    duck_pc, _ = load_species("duck", DUCK_POOLS)
    duck_pc = add_iso_cols(duck_pc, duck_coords, "duck")
    duck_d2n = P.duck_d2n_hh1716_sides()
    duck_pc = P.apply_d2n_hh1716_override(duck_pc, d2n=duck_d2n)
    for iso, col in ISO_COLS:
        for st in DUCK_STAGES:
            for r in duck_atac_signed(duck_pc, col, st):
                rows.append(dict(modality="ATAC", isoform=iso, stage=st, group="duck", **r))

    # ── CHICK ATAC (c + a/b), design-legal envelope ──
    chick_coords = pd.read_csv(f"{DATA}/chick_locus/locus_peak_coords.tsv", sep="\t")
    chick_pc, _ = load_species("chick", CHICK_POOLS)
    chick_pc = add_iso_cols(chick_pc, chick_coords, "chick")
    chick_pc["lane"] = chick_pc["cell"].map(lane_of)
    for iso, col in ISO_COLS:
        for cst, fst in SD.items():
            embs = chick_atac_embs(chick_pc, cst, col)
            nmap = {e["e"]: e["n"] for e in embs}
            for case, o in chick_cases(embs, cst).items():
                for e, v in o.items():
                    rows.append(dict(modality="ATAC", isoform=iso, stage=fst, group=case,
                                     embryo=e, signed=v, n_cells=nmap[e]))

    # ── DUCK RNA total (d2g PITX2_cp10k; d2n @ HH17/16) ──
    lc = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1",
          "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}
    # shared firmness rule (figures/_shared/lr_firmness.py): unknown values raise; this is
    # an embryo-PAIRED signed L-R readout, so a same-side duplicate (one half in two
    # libraries, no opposite half) is excluded -- from the side map AND the library loop.
    d2l_raw = pd.read_csv(f"{DATA}/d2l_lr_inference.tsv", sep="\t")
    d2l = resolve_sides(d2l_raw, {"direct", "placement"}, source="d2l_lr_inference.tsv",
                        paired=True, require_pair_constraint=True)
    side = {(g, lm): s for g, lm, s in zip(d2l.genotype, d2l.library_m, d2l.inferred_side)}
    dups = set(same_side_duplicates(d2l_raw, "d2l_lr_inference.tsv"))
    sf = {"HH13": "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
          "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
          "HH21": "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
          "HH25": "d2g_lpm_subclustering_HH25_ccreg.tsv.gz"}
    for st, fn in sf.items():
        sm = duck_d2n if st == "HH17/16" else side
        df = pd.read_csv(f"{DATA}/{fn}", sep="\t")
        df = df[df[get_leiden_col(st)].astype(str).isin(set(get_pool_map(st, "inference")))]
        df["lm"] = df.library.map(lambda x: lc.get(x, x))
        for g, d in df.groupby("genotype", observed=True):
            if g in dups:
                continue
            libs = sorted(d.lm.unique())
            if len(libs) != 2:
                continue
            a, b = libs
            da, db = d[d.lm == a], d[d.lm == b]
            if len(da) < MINC or len(db) < MINC:
                continue
            ma, mb = da.PITX2_cp10k.mean(), db.PITX2_cp10k.mean()
            sa, sb = sm.get((g, a)), sm.get((g, b))
            if {sa, sb} != {"L", "R"}:
                continue
            rows.append(dict(modality="RNA", isoform="total", stage=st, group="duck", embryo=g,
                             signed=(ma - mb) if sa == "L" else (mb - ma),
                             n_cells=int(min(len(da), len(db)))))

    # ── CHICK RNA total (cAN ventral; design-legal envelope) ──
    for cst, fst in SD.items():
        embs = chick_rna_total_embs(cst)
        nmap = {e["e"]: e["n"] for e in embs}
        for case, o in chick_cases(embs, cst).items():
            for e, v in o.items():
                rows.append(dict(modality="RNA", isoform="total", stage=fst, group=case,
                                 embryo=e, signed=v, n_cells=nmap[e]))

    # ── CHICK RNA PITX2a/b: DROPPED from the figure (build_pooled_rna_ab.py showed the
    #    depth-robust a/b-body FRACTION is c-mirrored, not an independent control: duck pooled
    #    Δ(ab/shared) is R>L at HH13/HH21/HH25 [= the arithmetic shadow of c], bilateral only at
    #    the HH17/16 trough — matching the prior AB_FRAC analysis).  The clean a/b bilateral
    #    control is the chromatin ATAC a/b panel.  The per-embryo a/b-body computation, the
    #    pooled c-mirror evidence and the LEGEND a/b-capture stat (ab/shared 0.55-0.82) come
    #    from an internal script (not distributed). ──

    out = pd.DataFrame(rows)
    # Canonical output lives with the other figure-consumed chick data in
    # figures/data/chick_control/ (the same home figS04 panel_F reads cAN_halfvshalf.tsv
    # from) — consumed by the figS04 panel_G renderer.
    out_tsv = os.path.join(DOC, "xspecies_pitx2_isoform_laterality.tsv")
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {os.path.relpath(out_tsv, REPO)}: {out.shape}")
    g = (out.groupby(["modality", "isoform", "stage", "group"])
         .agg(n=("signed", "size"), n_LgtR=("signed", lambda s: int((s > 0).sum())),
              med=("signed", "median")).round(3))
    pd.set_option("display.width", 200)
    print(g.to_string())


if __name__ == "__main__":
    main()
