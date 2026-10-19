"""Same-side duplicate libraries: one embryo-half tube split across two lanes.

The lane-loading record (``design/lane_loading_record_layout.tsv``) records that the
first multiome batch split ONE left-half tube -- lane-loading-record group ``17-B(2)`` (not repo embryo 17-B), repo
embryos 17-5 and 17-6 -- across ARC04 (GEX04) and ARC0S (GEX0S-duck).  Those two
embryos therefore carry the SAME side (L) in two distinct libraries: a technical
duplicate of one tissue, never an L/R pair.  Their right halves were never
sequenced (lane ARC03 was struck).

In ``design/embryo_lr_design.tsv`` this appears as an embryo with rows in two
DISTINCT ``library_m`` that share a side.  That is a different case from the
GEX03-1A/1B and PILOTGEX4A/B technical replicates, which ``collapse_lib`` folds
into ONE ``library_m`` before any inference.

Which GENOTYPES are those embryos is a data question (genotype -> embryo is never
assumed a priori).  It is answered here WITHOUT any PITX2 / laterality signal, by
the tube-split identity test: two aliquots of one nuclei suspension must have the
same cell-type composition up to multinomial noise, whereas a left half and a
right half (different tubes, dissected and dissociated separately) need not.
For each design duplicate library pair holding ``k`` embryos, among the genotypes
with >= ``MIN_CELLS`` nuclei in both libraries, EXACTLY ``k`` must pass the
composition-homogeneity chi-square (p > ``P_SAME``) and every other candidate must
fail it decisively (p < ``P_DIFF``).  When all ``k`` duplicate embryos came from
one lane-loading-record tube, their mutual nuclei ratio must also be preserved across the two
lanes (p > ``P_SAME``).  Anything else raises.

Where the answer lives
----------------------
The test reads ``figures/data/e3_composition.tsv``, an output of the downstream
step e3.  d2l / d2n must not depend on e3 at run time, so the answer is
COMMITTED as ``references/same_side_duplicate_genotypes.tsv`` (every candidate
genotype with its evidence and the ``is_duplicate`` call).  d2l / d2n read that
reference (``load_duplicate_reference``: no e3 read; the reference must still
reconcile with the design and the lane-loading record).  ``pipeline/validate_same_side_duplicates.py``
recomputes the evidence from e3 (``identify_duplicate_genotypes``) and raises on
any disagreement with the reference; it runs after e3 in the cluster runner
``FC_08656/run_regenerate_figure_data.sh`` (not distributed) and in ``make check``.

Consumer contract (sidecar semantics)
-------------------------------------
* d2l / d2n write BOTH library rows of a duplicate genotype with
  ``inferred_side`` = the design side, ``firmness`` = ``SAME_SIDE_DUP``
  (d2n: ``chromatin_`` + it), ``votes`` = ``DUP_VOTES``, ``is_singleton`` = False,
  ``partner_lib`` = the other library (display only -- never used to derive an
  opposite side; placement raises if it ever reaches a duplicate genotype).
* A duplicate genotype is never an L/R pair: paired contrasts must drop it
  (it has no R), and unpaired contrasts must count it ONCE (pool its two
  libraries -- see ``merged_library_label``), as the tech-rep collapse does.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

REPO = Path(__file__).resolve().parent.parent
FIGDATA = REPO / "figures" / "data"
DESIGN_PATH = REPO / "design" / "embryo_lr_design.tsv"
LANE_LAYOUT_PATH = REPO / "design" / "lane_loading_record_layout.tsv"
COMPOSITION_PATH = FIGDATA / "e3_composition.tsv"
DUPLICATE_REFERENCE = REPO / "references" / "same_side_duplicate_genotypes.tsv"
EVIDENCE_COLS = ["genotype", "lib_a", "lib_b", "dup_side", "n_a", "n_b", "chi2", "dof",
                 "p_composition", "tv_distance", "is_duplicate", "dup_embryos",
                 "p_tube_ratio"]
_EXACT_COLS = ["genotype", "lib_a", "lib_b", "dup_side", "n_a", "n_b", "dof",
               "is_duplicate", "dup_embryos"]
_FLOAT_COLS = ["chi2", "p_composition", "tv_distance", "p_tube_ratio"]
REFERENCE_RTOL = 1e-6   # the reference is written at 10 significant digits

sys.path.insert(0, str(REPO / "figures"))
from _shared.library_naming import ARC_TO_D2L_LIBRARY, collapse_lib  # noqa: E402
# the firmness alphabet + pooled-library label are shared with every consumer
from _shared.lr_firmness import SAME_SIDE_DUP, merged_library_label  # noqa: E402,F401

DUP_VOTES = "dup"

MIN_CELLS = 100      # nuclei per library for a genotype to be a duplicate candidate
P_SAME = 0.05        # duplicate: composition homogeneous across the two lanes
P_DIFF = 1e-6        # distinct halves: composition decisively heterogeneous


def design_duplicate_pairs(design: pd.DataFrame) -> list[dict]:
    """Embryos whose design rows put the SAME side in two distinct library_m.

    ``design`` must already carry ``lib_m`` (tech-rep collapsed) and be deduped
    on (embryo, lib_m, side).  Returns one record per (library pair, side):
    {"libs": frozenset, "side": "L"/"R", "embryos": [...], "stages": {...}}.
    """
    for col in ("embryo", "stage", "lib_m", "side"):
        if col not in design.columns:
            raise KeyError(f"design frame lacks column {col!r}")
    pairs: dict[tuple[frozenset, str], dict] = {}
    for emb, grp in design.groupby("embryo"):
        by_lib = grp.groupby("lib_m")["side"].agg(lambda s: set(s))
        multi = [lib for lib, sides in by_lib.items() if len(sides) > 1]
        if multi:
            raise ValueError(f"design embryo {emb} has both sides in one "
                             f"library_m {multi} -- undefined")
        side_libs: dict[str, list[str]] = {}
        for lib, sides in by_lib.items():
            side_libs.setdefault(next(iter(sides)), []).append(lib)
        for side, libs in side_libs.items():
            if len(libs) < 2:
                continue
            if len(by_lib) != 2:
                raise ValueError(
                    f"design embryo {emb} repeats side {side} in {sorted(libs)} "
                    f"but spans {len(by_lib)} library_m {sorted(by_lib.index)}; "
                    f"only a two-library same-side duplicate is modelled")
            key = (frozenset(libs), side)
            rec = pairs.setdefault(key, {"libs": frozenset(libs), "side": side,
                                         "embryos": [], "stages": set()})
            rec["embryos"].append(emb)
            rec["stages"].update(grp["stage"].unique())
    out = sorted(pairs.values(), key=lambda r: sorted(r["libs"]))
    libsets = [r["libs"] for r in out]
    if len(libsets) != len(set(libsets)):
        raise ValueError("one library pair carries same-side duplicates of BOTH "
                         "sides; the duplicate side is then not design-determined")
    return out


def lane_layout_cross_check(pairs: list[dict]) -> list[dict]:
    """Assert the design duplicates match the lane-loading-record split-tube transcription.

    A duck (non-``Chick-``) lane-loading-record tube listed in two lanes with the same side
    and the same non-empty split colour is ONE tube split across two lanes.
    Returns the split-tube records [{"tube", "libs", "side", "embryos_n"}].
    """
    layout = pd.read_csv(LANE_LAYOUT_PATH, sep="\t", comment="#", dtype=str,
                     keep_default_na=False)
    duck = layout[~layout["tube"].str.startswith("Chick-")]
    splits = []
    for (tube, side), grp in duck.groupby(["tube", "side"]):
        if len(grp) < 2:
            continue
        colours = set(grp["split_color"])
        if len(grp) != 2 or len(colours) != 1 or "" in colours:
            raise ValueError(f"lane-loading-record duck tube {tube}-{side} listed in "
                             f"{grp['lane'].tolist()} with colours {colours}; "
                             f"expected one split tube (2 lanes, one colour)")
        libs = set()
        for lane in grp["lane"]:
            arc = "ARC0S-duck" if lane == "ARC0S" else lane
            if arc not in ARC_TO_D2L_LIBRARY:
                raise KeyError(f"lane-loading-record lane {lane!r} has no library mapping")
            libs.add(ARC_TO_D2L_LIBRARY[arc])
        n = {int(x) for x in grp["embryos_n"]}
        if len(n) != 1:
            raise ValueError(f"lane-loading-record tube {tube} embryos_n differs across lanes")
        splits.append({"tube": tube, "libs": frozenset(libs), "side": side,
                       "embryos_n": n.pop()})
    want = {(s["libs"], s["side"]): 0 for s in splits}
    for s in splits:
        want[(s["libs"], s["side"])] += s["embryos_n"]
    have = {(p["libs"], p["side"]): len(p["embryos"]) for p in pairs}
    if want != have:
        raise ValueError(f"design same-side duplicates {have} disagree with the "
                         f"lane-loading-record split tubes {want}")
    return splits


def _composition_by_lib(path: Path = COMPOSITION_PATH) -> pd.DataFrame:
    comp = pd.read_csv(path, sep="\t")
    for col in ("genotype_id", "library", "cell_type", "n_cells"):
        if col not in comp.columns:
            raise KeyError(f"{path.name} lacks column {col!r}")
    comp["lib_m"] = comp["library"].map(collapse_lib)
    return (comp.groupby(["genotype_id", "lib_m", "cell_type"])["n_cells"]
                .sum().unstack("cell_type", fill_value=0))


def identify_duplicate_genotypes(design: pd.DataFrame,
                                 composition_path: Path = COMPOSITION_PATH,
                                 verbose: bool = True) -> pd.DataFrame:
    """Return the per-candidate evidence table; ``is_duplicate`` marks the
    genotypes that are same-side duplicates.  Raises unless the composition
    test separates exactly k duplicates from the remaining candidates."""
    pairs = design_duplicate_pairs(design)
    splits = lane_layout_cross_check(pairs)
    cols = EVIDENCE_COLS
    if not pairs:
        return pd.DataFrame(columns=cols)
    comp = _composition_by_lib(composition_path)
    totals = comp.sum(axis=1)
    rows = []
    for pr in pairs:
        lib_a, lib_b = sorted(pr["libs"])
        k = len(pr["embryos"])
        genos = sorted(set(comp.index.get_level_values(0)))
        cand = []
        for g in genos:
            n_a = int(totals.get((g, lib_a), 0))
            n_b = int(totals.get((g, lib_b), 0))
            if n_a >= MIN_CELLS and n_b >= MIN_CELLS:
                cand.append(g)
            elif (n_a > 0 and n_b > 0) and (n_a >= MIN_CELLS or n_b >= MIN_CELLS):
                raise ValueError(
                    f"genotype {g} has {n_a} nuclei in {lib_a} and {n_b} in "
                    f"{lib_b}: too few in one lane to test, too many to ignore")
        recs = []
        for g in cand:
            t = np.vstack([comp.loc[(g, lib_a)].to_numpy(),
                           comp.loc[(g, lib_b)].to_numpy()]).T
            t = t[t.sum(axis=1) > 0]
            chi2, p, dof, _ = chi2_contingency(t)
            fa, fb = t[:, 0] / t[:, 0].sum(), t[:, 1] / t[:, 1].sum()
            recs.append({"genotype": g, "lib_a": lib_a, "lib_b": lib_b,
                         "dup_side": pr["side"], "n_a": int(t[:, 0].sum()),
                         "n_b": int(t[:, 1].sum()), "chi2": float(chi2),
                         "dof": int(dof), "p_composition": float(p),
                         "tv_distance": float(0.5 * np.abs(fa - fb).sum()),
                         "dup_embryos": ",".join(pr["embryos"])})
        same = [r for r in recs if r["p_composition"] > P_SAME]
        diff = [r for r in recs if r["p_composition"] < P_DIFF]
        if len(same) != k or len(same) + len(diff) != len(recs):
            raise ValueError(
                f"tube-split test for {lib_a}+{lib_b} ({pr['side']}; design "
                f"embryos {pr['embryos']}) is not decisive: expected exactly {k} "
                f"homogeneous candidates (p>{P_SAME}) and the rest p<{P_DIFF}; "
                f"got " + ", ".join(f"{r['genotype']} p={r['p_composition']:.3g}"
                                    for r in recs))
        p_ratio = np.nan
        one_tube = [s for s in splits if s["libs"] == pr["libs"]
                    and s["side"] == pr["side"] and s["embryos_n"] == k]
        if one_tube and k >= 2:
            ratio = np.array([[r["n_a"], r["n_b"]] for r in same])
            _, p_ratio, _, _ = chi2_contingency(ratio)
            if not p_ratio > P_SAME:
                raise ValueError(
                    f"duplicate genotypes {[r['genotype'] for r in same]} share "
                    f"lane-loading-record tube {one_tube[0]['tube']} but their nuclei ratio "
                    f"differs between {lib_a} and {lib_b} (p={p_ratio:.3g})")
        same_ids = {r["genotype"] for r in same}
        for r in recs:
            r["is_duplicate"] = r["genotype"] in same_ids
            r["p_tube_ratio"] = p_ratio if r["is_duplicate"] else np.nan
            rows.append(r)
    ev = pd.DataFrame(rows)[cols]
    if verbose:
        print("  same-side duplicate test (cell-type composition, lane A vs B):")
        for r in ev.itertuples():
            tag = "DUPLICATE" if r.is_duplicate else "distinct halves"
            print(f"    {r.genotype:22s} {r.lib_a}/{r.lib_b} n={r.n_a}/{r.n_b} "
                  f"TV={r.tv_distance:.3f} p={r.p_composition:.3g}  -> {tag}")
    return ev


def load_collapsed_design(path: Path = DESIGN_PATH) -> pd.DataFrame:
    """design/embryo_lr_design.tsv with ``lib_m`` (tech reps collapsed), deduped on
    (embryo, lib_m, side) -- the frame d2l / d2n build before calling into here."""
    design = pd.read_csv(path, sep="\t")
    design["lib_m"] = design["library"].map(collapse_lib)
    return design.drop_duplicates(["embryo", "lib_m", "side"]).reset_index(drop=True)


def load_duplicate_reference(design: pd.DataFrame,
                             path: Path = DUPLICATE_REFERENCE,
                             verbose: bool = True) -> pd.DataFrame:
    """The committed duplicate evidence table (no e3 read at run time).

    Returns the same frame ``identify_duplicate_genotypes`` would.  Raises unless
    the reference reconciles with the design and the lane-loading record: every design
    same-side duplicate (library pair, side, k embryos) has exactly k
    ``is_duplicate`` genotypes on that pair and side, citing those embryos; no
    reference row sits outside a design duplicate pair; one row per genotype.
    Whether the evidence still matches e3 is the validator's job
    (``validate_reference``)."""
    if not path.is_file():
        raise FileNotFoundError(f"same-side duplicate reference missing: {path}")
    ref = pd.read_csv(path, sep="\t", comment="#")
    missing = [c for c in EVIDENCE_COLS if c not in ref.columns]
    if missing:
        raise KeyError(f"{path.name} lacks column(s) {missing}")
    if ref["is_duplicate"].dtype != bool:
        raise ValueError(f"{path.name}: is_duplicate must be True/False, got "
                         f"{sorted(ref['is_duplicate'].astype(str).unique())}")
    if ref["genotype"].duplicated().any():
        raise ValueError(f"{path.name}: genotype listed twice: "
                         f"{sorted(ref.loc[ref['genotype'].duplicated(), 'genotype'])}")
    pairs = design_duplicate_pairs(design)
    lane_layout_cross_check(pairs)
    claimed = set()
    for pr in pairs:
        lib_a, lib_b = sorted(pr["libs"])
        rows = ref[(ref["lib_a"] == lib_a) & (ref["lib_b"] == lib_b)
                   & (ref["dup_side"] == pr["side"])]
        dup = rows[rows["is_duplicate"]]
        want_emb = ",".join(pr["embryos"])
        if len(dup) != len(pr["embryos"]):
            raise ValueError(
                f"{path.name}: design has {len(pr['embryos'])} same-side duplicate "
                f"embryo(s) {pr['embryos']} in {lib_a}+{lib_b} ({pr['side']}) but the "
                f"reference marks {len(dup)} duplicate genotype(s) "
                f"{sorted(dup['genotype'])}")
        bad_emb = sorted(set(rows["dup_embryos"]) - {want_emb})
        if bad_emb:
            raise ValueError(f"{path.name}: {lib_a}+{lib_b} rows cite embryos {bad_emb}, "
                             f"design says {want_emb}")
        claimed.update(rows.index)
    stray = ref.drop(index=sorted(claimed))
    if len(stray):
        raise ValueError(f"{path.name}: rows outside every design duplicate pair: "
                         f"{stray[['genotype', 'lib_a', 'lib_b', 'dup_side']].to_dict('records')}")
    ev = ref[EVIDENCE_COLS].reset_index(drop=True)
    if verbose:
        print(f"  same-side duplicates from {path.relative_to(REPO)} "
              f"(validated against e3 by validate_same_side_duplicates.py):")
        for r in ev.itertuples():
            tag = "DUPLICATE" if r.is_duplicate else "distinct halves"
            print(f"    {r.genotype:22s} {r.lib_a}/{r.lib_b} n={r.n_a}/{r.n_b} "
                  f"TV={r.tv_distance:.3f} p={r.p_composition:.3g}  -> {tag}")
    return ev


def validate_reference(design: pd.DataFrame,
                       reference_path: Path = DUPLICATE_REFERENCE,
                       composition_path: Path = COMPOSITION_PATH) -> pd.DataFrame:
    """Recompute the duplicate evidence from e3 and raise on ANY disagreement with
    the committed reference (genotype set, calls, counts exactly; chi2 / p / TV at
    rtol ``REFERENCE_RTOL``).  Returns the recomputed evidence on success."""
    ref = load_duplicate_reference(design, reference_path, verbose=False)
    fresh = identify_duplicate_genotypes(design, composition_path, verbose=False)
    key = ["lib_a", "lib_b", "dup_side", "genotype"]
    ref = ref.sort_values(key).reset_index(drop=True)
    fresh = fresh.sort_values(key).reset_index(drop=True)
    problems = []
    gr, gf = set(ref["genotype"]), set(fresh["genotype"])
    if gr != gf:
        problems.append(f"candidate genotypes differ: reference only {sorted(gr - gf)}, "
                        f"e3 only {sorted(gf - gr)}")
    else:
        for c in _EXACT_COLS:
            a, b = ref[c].astype(str), fresh[c].astype(str)
            for i in np.flatnonzero((a != b).to_numpy()):
                problems.append(f"{ref.at[i, 'genotype']} {c}: reference {a[i]!r}, "
                                f"e3 recomputes {b[i]!r}")
        for c in _FLOAT_COLS:
            a = ref[c].to_numpy(float)
            b = fresh[c].to_numpy(float)
            ok = np.isclose(a, b, rtol=REFERENCE_RTOL, atol=0.0, equal_nan=True)
            for i in np.flatnonzero(~ok):
                problems.append(f"{ref.at[i, 'genotype']} {c}: reference {a[i]:.10g}, "
                                f"e3 recomputes {b[i]:.10g}")
    if problems:
        raise ValueError(
            f"{reference_path.name} disagrees with the evidence recomputed from "
            f"{composition_path} ({len(problems)} difference(s)); if e3 changed on "
            f"purpose, review the calls and regenerate the reference with "
            f"`validate_same_side_duplicates.py --write-reference`:\n  "
            + "\n  ".join(problems))
    return fresh


def duplicate_map(evidence: pd.DataFrame) -> dict[str, tuple[frozenset, str]]:
    """genotype -> (frozenset of its two library_m, design side)."""
    d = evidence[evidence["is_duplicate"].astype(bool)]
    return {r.genotype: (frozenset((r.lib_a, r.lib_b)), r.dup_side)
            for r in d.itertuples()}


def assert_sidecar_duplicates(infdf: pd.DataFrame, dups: dict,
                              firmness_value: str) -> None:
    """Post-inference invariant: duplicate rows (and only them) carry the
    duplicate firmness and the design side."""
    is_dup_firm = infdf["firmness"] == firmness_value
    bad = infdf[is_dup_firm & ~infdf["genotype"].isin(dups)]
    if len(bad):
        raise RuntimeError(f"{firmness_value} on non-duplicate genotype(s): "
                           f"{sorted(bad['genotype'].unique())}")
    for g, (libs, side) in dups.items():
        sub = infdf[infdf["genotype"] == g]
        if sub.empty:
            continue
        if not set(sub["library_m"]) <= set(libs):
            raise RuntimeError(f"duplicate genotype {g} has rows in "
                               f"{sorted(set(sub['library_m']))}, outside {sorted(libs)}")
        if not ((sub["firmness"] == firmness_value).all()
                and (sub["inferred_side"] == side).all()):
            raise RuntimeError(f"duplicate genotype {g} rows not all "
                               f"{firmness_value}/{side}: "
                               f"{sub[['library_m', 'firmness', 'inferred_side']].to_dict('records')}")


def sidecar_duplicate_libraries(lr: pd.DataFrame) -> dict[str, str]:
    """genotype -> merged library label, from an L/R sidecar's duplicate rows
    (firmness ``same_side_dup`` or ``chromatin_same_side_dup``)."""
    firm = lr["firmness"].astype(str).str.removeprefix("chromatin_")
    d = lr[firm == SAME_SIDE_DUP]
    out = {}
    for g, grp in d.groupby("genotype"):
        libs = set(grp["library_m"])
        if len(libs) != 2 or grp["inferred_side"].nunique() != 1:
            raise RuntimeError(f"sidecar duplicate {g} malformed: "
                               f"{grp[['library_m', 'inferred_side']].to_dict('records')}")
        out[g] = merged_library_label(libs)
    return out


def pool_duplicate_libraries(obs_keep: pd.DataFrame, lr_tsv: str,
                             genotype_col: str = "unified_genotype_id",
                             log=print) -> pd.DataFrame:
    """Pseudobulk builders keyed on library_m (h1*): rewrite library_m of a
    same-side duplicate's cells to the merged label so its two libraries pool
    into ONE pseudobulk per (genotype, side, subcluster, stage) -- the same
    treatment LIBRARY_COLLAPSE gives the GEX03-1A/1B tech reps.  Raises if a
    duplicate-firmness cell sits in a library outside its duplicate pair."""
    lr = pd.read_csv(lr_tsv, sep="\t")
    merged = sidecar_duplicate_libraries(lr)
    mask = obs_keep["firmness_norm"] == SAME_SIDE_DUP
    if not mask.any():
        return obs_keep
    gids = obs_keep.loc[mask, genotype_col].astype(str)
    missing = sorted(set(gids) - set(merged))
    if missing:
        raise RuntimeError(f"cells carry {SAME_SIDE_DUP} for genotypes {missing} "
                           f"with no duplicate pair in {lr_tsv}")
    bad = [(g, l) for g, l in zip(gids, obs_keep.loc[mask, "library_m"])
           if l not in merged[g].split("+")]
    if bad:
        raise RuntimeError(f"{SAME_SIDE_DUP} cells outside their duplicate "
                           f"libraries: {sorted(set(bad))}")
    obs_keep = obs_keep.copy()
    obs_keep.loc[mask, "library_m"] = gids.map(merged).to_numpy()
    log(f"  same-side duplicates pooled across libraries: "
        f"{int(mask.sum()):,} cells, {dict(sorted(merged.items()))}")
    return obs_keep
