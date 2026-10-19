"""Shared L/R-resolution panel renderer.

Pivots a per-(genotype × library) inference table into the matrix layout
used by figS02 panel D and figS03 panel C, computes per-block design-vs-
inferred footers, and renders the PNG.  Both consumers feed in different
TSV sidecars (figS02D = `d2l_lr_inference.tsv`; figS03 panel C =
`lr_inference_multiome.tsv`) but share the rendering primitives so the
two panels read as a methodological progression rather than visually
distinct layouts.

Public API:
    pivot_by_genotype, assert_design_stages, assemble_blocks,
    compute_block_footers, load_design, render
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch

from .style import apply_pub_style, STAGE_COLORS, display_name
from .layout import W_FULL, apply_layout_fonts
from .colors import get_palette as _reg_palette
from .library_naming import TECH_REP, collapse_lib

apply_pub_style()
apply_layout_fonts()

FIGDATA = Path(__file__).resolve().parent.parent / "data"
DESIGN_PATH = (Path(__file__).resolve().parent.parent.parent
               / "design" / "embryo_lr_design.tsv")

STAGE_BUCKETS = ["HH13", "HH17/16", "HH21", "HH25"]
METRICS = ["pitx2_pos", "pitx2c_pos", "pitx2_mean"]
METRIC_LABELS = ["PITX2+", "PITX2c+", "mean"]
METRIC_UNITS = ["%", "%", "CP10K"]

# Chromatin-modality counterparts.  Used by the figS03 panel C
# `expanded_chromatin_columns=True` testing variant — that variant shows
# *both* modality column triples side-by-side per library so a reader
# can scan snRNA + snATAC simultaneously.  Cell 1 and Cell 3 are
# fragment-density metrics (`<region>_frags_per_cell`) — available at
# ALL four stages via either d2n (HH17/16) or f5 (HH13/HH21/HH25)
# augmentation in the fig03 driver — so the chromatin sub-block
# populates uniformly across stages and shows the corroboration claim
# (chromatin L>R paired Wilcoxon p≤0.035 across HH13/HH21/HH25;
# specifically: HH13=0.035, HH21=0.014, HH25=0.014 at PITX2c_promoter
# after MIN_CELLS=20 stability filter at the visualization layer).
# Cell 2 is the isoform-specific open-fraction (`pitx2c_open_frac`)
# computed by d2n at HH17/16 only; the cell renders empty at flanking
# stages, honestly indicating "open %" was not computed there.
METRICS_CHROM = ["pitx2_genebody_frags_per_cell", "pitx2c_open_frac",
                 "pitx2c_frags_per_cell"]
METRIC_CHROM_LABELS = ["genebody", "open", "PITX2c"]
METRIC_CHROM_UNITS = ["f/c", "%", "f/c"]

_LAT = _reg_palette("lateral")
SIDE_COLOR = {"L": _LAT["L"], "R": _LAT["R"], "?": "#bbb"}

EVIDENCE_LOG2_SATURATE = 2.0  # saturates at 4x fold
PSEUDOCOUNT = 1e-3
INCONCLUSIVE_GRAY = "#555555"
MIN_VISIBLE_ALPHA = 0.35

# 3-tier flag tints (chromatin-canonical-hh1716 §3 + Task 3): the d2m
# sidecar's `design_conflict` is the OR of three causes; the panel
# distinguishes them so readers can see which kind of conflict each
# row carries.
#   red  (#fdecea): block-level pc_fail (snRNA call/inference issue —
#                   the block's L/R counts cannot be reconciled with
#                   the design without revisiting the inference)
#   orange (#fff3e0): pc_fail==False but design_conflict==True
#                     (cardinality==0 or asymmetric multi-candidate
#                     with no interchangeable peer — irreducible from
#                     a single (library, side) measurement, regardless
#                     of which modality produced it)
#   none + asterisk on embryo_id: symmetric peer group (interchangeable
#                                  candidate set — assignment is
#                                  arbitrary but lossless)
TINT_PC_FAIL = "#fdecea"
TINT_DESIGN_AMBIGUITY = "#fff3e0"

# Slots with fewer cells than this threshold get a visual low-confidence
# cue (diagonal hatch on metric cells, italic values, red cell-count label).
# d2l inference is NOT retroactively edited — the cue is editorial honesty,
# not a logic change.  Catches 2 rows in the current sidecar: donor5 GEX04
# n=2 (all-zero metrics) and g21 GEX08 HH13 n=8 (legit low-n paired).
LOW_N_THRESHOLD = 10
LOW_N_HATCH_COLOR = "#a33"
LOW_N_TEXT_COLOR = "#c94438"


# ---------------------------------------------------------------------
# Data loading + pivot
# ---------------------------------------------------------------------

def assert_lr_input_schemas(d2l: pd.DataFrame, d2m: pd.DataFrame,
                            *, include_chromatin: bool = False) -> None:
    """Fail-fast schema check at TSV-load time.

    Originating bug class: a 2-month silent regression where `pc_fail`
    was missing from `d2m_genotype_to_embryo.tsv` — the renderer's
    `r.get("pc_fail", False)` defaulted everything to False, so 12
    `pc_fail=True` rows rendered orange instead of red.  The
    accompanying renderer change drops the `.get(default)` shape so
    column access raises KeyError on missing schema; this helper turns
    that into a friendly error pointing to the regen script before any
    rendering work happens.
    """
    d2m_required = {"pc_fail", "design_conflict"}
    missing_d2m = d2m_required - set(d2m.columns)
    if missing_d2m:
        raise RuntimeError(
            f"d2m TSV missing column(s) {sorted(missing_d2m)} — re-run "
            f"`pipeline/d2m_genotype_to_embryo.py` (or "
            f"`pipeline/d2m_multiome_genotype_to_embryo.py`) before "
            f"rendering.  These columns were silently absent in a "
            f"2-month regression."
        )

    d2l_required = set(METRICS)
    if include_chromatin:
        # cohort_indirect_side carries the d2o `apply_overrides` sidecar
        # that pivot_by_genotype reads to paint the design-ambiguity
        # orange "=L?" / "=R?" treatment.  Without this column, the
        # renderer's `if "cohort_indirect_side" in r.index` guard would
        # silently no-op, regressing donor8's soft-L resolution.
        d2l_required |= set(METRICS_CHROM) | {"cohort_indirect_side"}
    missing_d2l = d2l_required - set(d2l.columns)
    if missing_d2l:
        raise RuntimeError(
            f"d2l/multiome TSV missing metric column(s) "
            f"{sorted(missing_d2l)} — re-run `pipeline/d2l_lr_inference.py` "
            f"(or `pipeline/d2o_merge_multiome_inference.py`) before rendering."
        )


def _design_stage_match(dstage: str) -> list[str]:
    if dstage == "HH17/16":
        return ["HH17", "HH16"]
    return [dstage]


def load_design() -> pd.DataFrame:
    design = pd.read_csv(DESIGN_PATH, sep="\t")
    design["lib_m"] = design["library"].map(collapse_lib)
    design = design.drop_duplicates(["embryo", "lib_m", "side"]).reset_index(drop=True)
    return design


def pivot_by_genotype(d2l: pd.DataFrame, d2m: pd.DataFrame,
                      include_chromatin: bool = False) -> pd.DataFrame:
    """One row per genotype with alphabetical library A/B slots.

    When `include_chromatin=True`, also pull the chromatin metric
    triple (`pitx2_genebody_open_frac`, `pitx2c_open_frac`,
    `pitx2c_frags_per_cell`) into matching `{slot}_<metric>` keys.
    These will be NaN at HH13/HH21/HH25 (where the multiome merger
    keeps snRNA values verbatim) and populated at HH17/16 (where the
    chromatin call replaces snRNA).  Used by the figS03 panel C
    expanded-columns testing variant.
    """
    # Fail-fast on row-set mismatch.  d2m is built by iterating d2l's
    # genotype groupby (`narrow_one_genotype` per group), so by
    # construction d2m has every genotype d2l has.  If that invariant
    # breaks (e.g. a future d2m bug filters genotypes out of its
    # narrowing loop), the row-level `d2m_row is not None` ternaries
    # below would have silently defaulted pc_fail/design_conflict to
    # False — the same regression class as the originating .get(default)
    # bug, just shifted from "missing column" to "missing row".
    # Surface it.
    d2l_genos = set(d2l["genotype"].unique())
    d2m_genos = set(d2m["genotype"].unique())
    missing_in_d2m = d2l_genos - d2m_genos
    if missing_in_d2m:
        raise RuntimeError(
            f"d2m TSV is missing {len(missing_in_d2m)} genotype(s) "
            f"present in d2l/multiome: {sorted(missing_in_d2m)} — "
            f"re-run pipeline/d2m_genotype_to_embryo.py (or "
            f"d2m_multiome_genotype_to_embryo.py) so every inference "
            f"row has a narrowing record."
        )

    d2m_lookup = d2m.set_index("genotype")
    rows = []
    for gid, grp in d2l.groupby("genotype"):
        grp = grp.sort_values("library_m").reset_index(drop=True)
        # The A/B slot model assumes at most two libraries per
        # genotype.  A future un-collapsed tech-rep would silently
        # truncate the third library and make it invisible in the
        # panel.  Fail loudly so the caller has to extend the model.
        if len(grp) > 2:
            raise RuntimeError(
                f"genotype {gid} has {len(grp)} libraries "
                f"({grp['library_m'].tolist()}) — A/B slot model "
                f"needs extension for this case"
            )
        A = grp.iloc[0]
        B = grp.iloc[1] if len(grp) > 1 else None

        # Row-set invariant enforced above; direct lookup is safe.
        d2m_row = d2m_lookup.loc[gid]

        row: dict = {
            "genotype":           gid,
            "stage":              A["stage"],
            "design_stage":       A["design_stage"],
            "firmness":           A["firmness"],
            "pair_constraint_ok": A["pair_constraint_ok"],
            "embryo_id_display":  d2m_row["embryo_id_display"],
            "candidate_embryos":  d2m_row["candidate_embryos"],
            "cardinality":        int(d2m_row["cardinality"]),
            "interchangeable_with": (d2m_row["interchangeable_with"]
                                     if pd.notna(d2m_row["interchangeable_with"])
                                     else ""),
            "design_conflict":    bool(d2m_row["design_conflict"]),
            "pc_fail":            bool(d2m_row["pc_fail"]),
        }

        def unpack(r, slot):
            if r is None:
                base = {f"{slot}_{k}": None for k in
                        ["lib", "side", "n"] + METRICS}
                base[f"{slot}_cohort_indirect"] = False
                if include_chromatin:
                    base.update({f"{slot}_{m}": None for m in METRICS_CHROM})
                return base
            side = r["inferred_side"]
            # Cohort-indirect soft-call override: applied at d2o time
            # from references/genotype_overrides.tsv, surfaced here as
            # the optional `cohort_indirect_side` column on the multiome
            # sidecar.  d2l_lr_inference.tsv (figS02D) does not carry
            # the column — the `in r.index` guard makes this no-op
            # there.  When present and the row's d2l/d2n side is
            # unresolved ("?" / NaN), promote to the override side and
            # set the per-slot flag so the renderer paints the
            # design-ambiguity-orange "=L?" / "=R?" treatment.  d2m_
            # multiome's narrowing keys on `inferred_side` (untouched),
            # so cardinality and design_conflict honest-disclosure
            # tints survive the override.
            is_ci = False
            if "cohort_indirect_side" in r.index:
                ci_side = r["cohort_indirect_side"]
                if (pd.notna(ci_side) and ci_side in ("L", "R")
                        and side not in ("L", "R")):
                    side = ci_side
                    is_ci = True
            d = {
                f"{slot}_lib":  r["library_m"],
                f"{slot}_side": side,
                f"{slot}_n":    int(r["n_inference"]),
                f"{slot}_cohort_indirect": is_ci,
                **{f"{slot}_{m}": (float(r[m]) if not pd.isna(r[m]) else None)
                   for m in METRICS},
            }
            if include_chromatin:
                d.update({f"{slot}_{m}":
                          (float(r[m]) if not pd.isna(r[m]) else None)
                          for m in METRICS_CHROM})
            return d

        row.update(unpack(A, "A"))
        row.update(unpack(B, "B"))
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Block assembly + O3 ordering
# ---------------------------------------------------------------------

def _stage_rank(s: str) -> int:
    return STAGE_BUCKETS.index(s)   # ValueError on an unknown stage (no sort-last default)


def assert_design_stages(d2l: pd.DataFrame) -> None:
    """Fail fast on unknown design_stage labels.

    Without this guard, an unknown stage label silently flows to rank=99
    and routes the row to the end of the panel — hiding a data-quality
    problem (typo, new stage bucket not added to STAGE_BUCKETS) behind
    an innocuous-looking layout quirk.
    """
    known = set(STAGE_BUCKETS)
    unknown = sorted(set(d2l["design_stage"].dropna().unique()) - known)
    if unknown:
        raise RuntimeError(
            f"d2l design_stage contains unknown labels {unknown}; "
            f"expected one of {STAGE_BUCKETS}.  Add to STAGE_BUCKETS "
            f"(and verify stage_rank ordering) if legitimate."
        )


def _side_rank(s) -> int:
    return {"L": 0, "R": 1}.get(s if isinstance(s, str) else "", 2)


def assemble_blocks(genos: pd.DataFrame) -> pd.DataFrame:
    """Add block_key + sort so rows are in O3 order."""
    genos = genos.copy()
    genos["block_key"] = genos.apply(
        lambda r: (r["design_stage"],
                   r["A_lib"] if r["A_lib"] is not None else "",
                   r["B_lib"] if r["B_lib"] is not None else ""),
        axis=1,
    )
    genos["_stage"] = genos["design_stage"].map(_stage_rank)
    genos["_sA"] = genos["A_side"].map(_side_rank)
    genos["_sB"] = genos["B_side"].map(_side_rank)
    genos = (genos.sort_values(
        ["_stage", "A_lib", "B_lib", "_sA", "_sB", "genotype"],
        na_position="last")
        .drop(columns=["_stage", "_sA", "_sB"])
        .reset_index(drop=True))
    return genos


def compute_block_footers(genos: pd.DataFrame,
                          design: pd.DataFrame) -> pd.DataFrame:
    """One row per block_key: inferred vs design (N_L, N_R) per library.

    Design counts are filtered to embryos whose FULL library set matches
    the block's library set (not just "any embryo at this library").  This
    prevents cross-library-bridge embryos (e.g. 17-Md-13 at {GEX03-1,
    GEX08}) from leaking into unrelated blocks (e.g. the {GEX07, GEX08}
    G6 block) and inflating expected counts.
    """

    def _clean_lib(v):
        """Coerce NaN/None/empty-string into None so block-libs frozensets
        only contain real library names (block_key for singleton rows has
        NaN in the B slot)."""
        if v is None or not _has_lib(v):
            return None
        return v

    records = []
    # Precompute embryo → library set per design_stage bucket so the
    # per-block filter below is a cheap set-equality check.
    all_dstages = {bk[0] for bk in genos["block_key"].unique()}
    stage_embryo_libs: dict[tuple[str, str], frozenset] = {}
    for dstage in all_dstages:
        design_stages = _design_stage_match(dstage)
        stage_sub = design[design["stage"].isin(design_stages)]
        for emb, sub in stage_sub.groupby("embryo"):
            stage_embryo_libs[(dstage, emb)] = frozenset(sub["lib_m"].unique())
    for block_key, grp in genos.groupby("block_key", sort=False):
        dstage = block_key[0]
        a_lib = _clean_lib(block_key[1])
        b_lib = _clean_lib(block_key[2])
        design_stages = _design_stage_match(dstage)
        block_libs = frozenset(x for x in [a_lib, b_lib] if x is not None)
        # Embryos whose library set exactly equals the block's library set.
        block_embryos = [
            emb for (ds, emb), libs in stage_embryo_libs.items()
            if ds == dstage and libs == block_libs
        ]
        block_design = design[
            design["stage"].isin(design_stages)
            & design["embryo"].isin(block_embryos)
        ]

        def design_counts(lib):
            if lib is None:
                return 0, 0
            sub = block_design[block_design["lib_m"] == lib]
            return (int((sub["side"] == "L").sum()),
                    int((sub["side"] == "R").sum()))

        def infer_counts(col):
            return (int((grp[col] == "L").sum()),
                    int((grp[col] == "R").sum()))

        a_L_des, a_R_des = design_counts(a_lib)
        b_L_des, b_R_des = design_counts(b_lib)
        a_L_inf, a_R_inf = infer_counts("A_side")
        b_L_inf, b_R_inf = infer_counts("B_side")

        a_ok = (a_L_inf == a_L_des) and (a_R_inf == a_R_des)
        b_ok = (b_lib is None) or ((b_L_inf == b_L_des) and (b_R_inf == b_R_des))
        ok = a_ok and b_ok

        records.append({
            "block_key": block_key,
            "design_stage": dstage,
            "A_lib": a_lib,
            "B_lib": b_lib,
            "A_n_L_inf": a_L_inf, "A_n_R_inf": a_R_inf,
            "A_n_L_des": a_L_des, "A_n_R_des": a_R_des,
            "B_n_L_inf": b_L_inf, "B_n_R_inf": b_R_inf,
            "B_n_L_des": b_L_des, "B_n_R_des": b_R_des,
            "pair_constraint_ok": ok,
            "n_rows": len(grp),
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------
# Cell coloring
# ---------------------------------------------------------------------

def paired_cell(lib_val, partner_val, side):
    if lib_val is None or pd.isna(lib_val):
        return None
    if partner_val is None or pd.isna(partner_val):
        return None
    if lib_val <= partner_val:
        return None
    lg2 = np.log2((lib_val + PSEUDOCOUNT) / (partner_val + PSEUDOCOUNT))
    alpha = min(1.0, max(0.0, lg2 / EVIDENCE_LOG2_SATURATE))
    alpha = max(alpha, MIN_VISIBLE_ALPHA)
    if side in ("L", "R"):
        return (SIDE_COLOR[side], alpha)
    return (INCONCLUSIVE_GRAY, alpha)


def singleton_cell(lib_val, col_max):
    if lib_val is None or pd.isna(lib_val) or col_max <= 0:
        return None
    frac = float(lib_val) / col_max
    if frac <= 0:
        return None
    alpha = min(1.0, max(MIN_VISIBLE_ALPHA, frac))
    return (INCONCLUSIVE_GRAY, alpha)


# ---------------------------------------------------------------------
# 3-class firmness glyph
# ---------------------------------------------------------------------

SAME_SIDE_DUP = "same_side_dup"


def firmness_glyph(firmness: str) -> tuple[str, str]:
    """Firmness glyph.  d2o writes the canonical alphabet
    {direct, placement, singleton_rank, weak, same_side_dup} after
    stripping the chromatin_ prefix from d2n's output; the modality
    column carries the call-shape disambiguation downstream.
    same_side_dup = one embryo-half tube split across two libraries
    (lane-loading record): side is design-fixed, not inferred."""
    if firmness == "direct":
        return ("●", "#333")
    if firmness in ("placement", "singleton_rank"):
        return ("◐", "#555")
    if firmness == "weak":
        return ("○", "#aaa")
    if firmness == SAME_SIDE_DUP:
        return ("=", "#333")
    raise ValueError(f"firmness_glyph: unknown firmness {firmness!r}")


# ---------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------

COL_STAGE = 0.0
COL_MODALITY = 0.30   # only drawn when show_modality_strip=True (figS03 panel C)
COL_GENO = 0.35
COL_ID = 2.5
COL_A_LIB = 4.0
COL_A_M0 = 5.15
COL_A_M1 = 5.75
COL_A_M2 = 6.35
COL_A_N = 6.95
COL_B_LIB = 7.80
COL_B_M0 = 8.95
COL_B_M1 = 9.55
COL_B_M2 = 10.15
COL_B_N = 10.75
COL_CALL = 11.60
COL_END = 12.20

# Expanded-columns layout: 6 metric cells per library (3 snRNA + 3
# chromatin) plus all the existing non-metric columns.  Cells are
# narrower (0.40 vs 0.55) to keep the panel a manageable width
# (~14.05 vs 12.20 default).  figS03 panel C `expanded_chromatin_columns=True`
# switches to this layout so a reader can compare snRNA values and
# chromatin values side-by-side per row.
COL_A_LIB_X       = 4.0
COL_A_M0_X        = 5.05
COL_A_M1_X        = 5.50
COL_A_M2_X        = 5.95
COL_A_M0_C_X      = 6.50
COL_A_M1_C_X      = 6.95
COL_A_M2_C_X      = 7.40
COL_A_N_X         = 7.95
COL_B_LIB_X       = 8.80
COL_B_M0_X        = 9.85
COL_B_M1_X        = 10.30
COL_B_M2_X        = 10.75
COL_B_M0_C_X      = 11.30
COL_B_M1_C_X      = 11.75
COL_B_M2_C_X      = 12.20
COL_B_N_X         = 12.75
COL_CALL_X        = 13.65
COL_END_X         = 14.25


def render(genos: pd.DataFrame,
           footers: pd.DataFrame,
           out_path: str | os.PathLike | None = None,
           metric_units: list[str] | None = None,
           show_modality_strip: bool = False,
           show_stage_subtitles: bool = False,
           expanded_chromatin_columns: bool = False,
           show_legend: bool = True) -> None:
    """Render the L/R-inference panel to PNG.

    NOTE on the four display toggles (`show_modality_strip`,
    `show_stage_subtitles`, `expanded_chromatin_columns`, `show_legend`):
    they are NOT independent in practice — they co-vary into exactly TWO
    supported configurations, one per consumer:
      * figS02D  (`fig02_supp/generate_lr_inference.py`): all four defaults
        — strip=F, subtitles=F, expanded=F, legend=T  (compact 3-cell
        snRNA matrix WITH the in-image decoding legend).
      * figS03 panel C  (`figS03/generate_lr_resolution.py`): strip=T,
        subtitles=T, expanded=T, legend=F  (6-cell snRNA+chromatin matrix,
        modality strip, stage subtitles; decoding key in figS03/LEGEND.md).
    Other combinations are untested.  If a third consumer or a fifth toggle
    appears, collapse these into one `variant=` enum / RenderOptions object
    rather than adding another boolean.

    Parameters
    ----------
    out_path : optional
        Output PNG path.  Defaults to the figS02D conventional location
        (`figures/fig02_supp/panel_D_lr_inference.png`).  figS03 panel C
        passes its own destination so the same render pipeline can produce
        both panels without a path constant in this module.
    metric_units : optional
        Three-element list of unit subscripts shown under each metric column
        header (mean column especially).  Defaults to `["%", "%", "CP10K"]`
        (figS02D snRNA convention).  figS03 panel C passes `["%", "%", "mean"]`
        because chromatin rows display fragments-per-cell, not CP10K, in
        the third column slot.
    show_modality_strip : optional
        When True, render a thin colored strip in the gutter between the
        stage-color tab and the genotype-name column.  HH17/16 rows get a
        green "ATAC" tag (chromatin-derived call); other stages stay
        empty.  figS03 panel C sets this True; figS02D leaves it False.

    Cohort-indirect soft-call overrides are sourced from
    `references/genotype_overrides.tsv` and stamped into
    `lr_inference_multiome.tsv` by d2o; `pivot_by_genotype` translates
    them into per-slot `cohort_indirect` flags + the visible side
    string.  The renderer reads those flags directly — no parameter
    plumbing needed.

    show_stage_subtitles : optional
        When True, render a stage-section subtitle row at each
        design_stage transition.  The subtitle declares the modality
        ("snRNA expression" or "snATAC chromatin") and the
        modality-specific column-unit triple, so a reader scanning the
        panel sees crystal-clear "this stage is ATAC" / "this stage is
        snRNA" headers without having to consult the legend or chase
        the gutter modality strip.  figS03 panel C sets this True;
        figS02D leaves the default False.
    expanded_chromatin_columns : optional
        When True, the per-library metric block doubles to 6 cells: 3
        snRNA cells (PITX2+ %, PITX2c+ %, mean CP10K) plus 3 chromatin
        cells (gene-body open %, PITX2c open %, frags/cell).  The
        snRNA cells stay populated for HH13/HH21/HH25 rows and are
        empty at HH17/16; chromatin cells are populated at HH17/16 and
        empty at the flanking stages — readers see *which* modality
        produced the call for each row by which sub-block is filled.
        Caller MUST pass the `genos` frame produced by
        `pivot_by_genotype(..., include_chromatin=True)` so the
        `{slot}_<chrom-metric>` columns exist.  Used as a testing
        layout for figS03 panel C (panel may revert to the compact 3-
        cell layout if too busy).
    """
    if out_path is None:
        raise ValueError(
            "render() requires an explicit out_path; the figS02D / fig03 \n"
            "drivers pass their own destinations.  Previously a default of \n"
            "<figS02D dir>/panel_D_lr_inference.png was used; that lived only \n"
            "in the figS02D module and broke when render() moved to _shared."
        )
    out_path = str(out_path)
    if metric_units is None:
        metric_units = list(METRIC_UNITS)

    # Layout dict — picks compact (figS02D default) or expanded (fig03
    # figS03 panel C testing) column positions in one place so the helpers
    # downstream don't each need separate branches.
    if expanded_chromatin_columns:
        L = {
            "col_a_lib": COL_A_LIB_X,
            "col_a_m_rna": (COL_A_M0_X, COL_A_M1_X, COL_A_M2_X),
            "col_a_m_chrom": (COL_A_M0_C_X, COL_A_M1_C_X, COL_A_M2_C_X),
            "col_a_n": COL_A_N_X,
            "col_b_lib": COL_B_LIB_X,
            "col_b_m_rna": (COL_B_M0_X, COL_B_M1_X, COL_B_M2_X),
            "col_b_m_chrom": (COL_B_M0_C_X, COL_B_M1_C_X, COL_B_M2_C_X),
            "col_b_n": COL_B_N_X,
            "col_call": COL_CALL_X,
            "col_end": COL_END_X,
            "cell_width": 0.40,
            "include_chromatin_cells": True,
        }
    else:
        L = {
            "col_a_lib": COL_A_LIB,
            "col_a_m_rna": (COL_A_M0, COL_A_M1, COL_A_M2),
            "col_a_m_chrom": None,
            "col_a_n": COL_A_N,
            "col_b_lib": COL_B_LIB,
            "col_b_m_rna": (COL_B_M0, COL_B_M1, COL_B_M2),
            "col_b_m_chrom": None,
            "col_b_n": COL_B_N,
            "col_call": COL_CALL,
            "col_end": COL_END,
            "cell_width": 0.55,
            "include_chromatin_cells": False,
        }
    layout_rows: list[dict] = []
    prev_stage = None
    for block_key, grp in genos.groupby("block_key", sort=False):
        cur_stage = grp.iloc[0]["design_stage"]
        if show_stage_subtitles and cur_stage != prev_stage:
            layout_rows.append({"kind": "stage_subtitle",
                                "block_key": block_key,
                                "stage": cur_stage})
        prev_stage = cur_stage
        for _, r in grp.iterrows():
            layout_rows.append({"kind": "data", "block_key": block_key,
                                "row": r})
        footer = footers[footers["block_key"] == block_key].iloc[0]
        layout_rows.append({"kind": "footer", "block_key": block_key,
                            "row": footer})

    n = len(layout_rows)
    fig_h = max(6.0, n * 0.24)
    # Scale figure width proportionally to data-unit panel width so the
    # expanded-columns variant doesn't visually squish each cell.
    fig_w = W_FULL * (L["col_end"] / COL_END)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, L["col_end"])
    ax.set_ylim(-1.5, n - 0.5)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)

    max_n_cells = 1
    metric_keys = list(METRICS) + (list(METRICS_CHROM)
                                    if L["include_chromatin_cells"] else [])
    metric_col_max = {m: 0.0 for m in metric_keys}
    for _, row in genos.iterrows():
        for slot in ("A", "B"):
            # pivot_by_genotype guarantees per-slot {n} + per-metric
            # columns exist (None for singleton-B); bracket access is
            # fail-fast on a pivot regression.
            nc = row[f"{slot}_n"]
            if nc is not None and not pd.isna(nc):
                max_n_cells = max(max_n_cells, int(nc))
            for m in metric_keys:
                v = row[f"{slot}_{m}"]
                if v is not None and not pd.isna(v):
                    metric_col_max[m] = max(metric_col_max[m], float(v))

    # Block separators + 3-tier conflict tint:
    #   - red    = pc_fail (block-level L/R reconciliation failure;
    #              snRNA call-inference issue)
    #   - orange = design_conflict without pc_fail (asymmetric multi-
    #              candidate or cardinality 0; irreducible single-(lib,
    #              side) ambiguity in the wetlab design)
    #   - footer = light-blue background (summary)
    # See chromatin_canonical_hh1716.md §3 for the rationale.
    prev_block = None
    for i, lr in enumerate(layout_rows):
        if lr["kind"] == "stage_subtitle":
            # Stage subtitle gets its own background tint (light green for
            # ATAC stages, faint cool gray for snRNA).  The text itself
            # is rendered later in the per-row pass.
            stage = lr["stage"]
            bg = "#e8f1e3" if stage == "HH17/16" else "#f1f3f6"
            ax.add_patch(Rectangle((0, i - 0.5), L["col_end"], 1,
                                   facecolor=bg, edgecolor="none", zorder=0))
            prev_block = lr["block_key"]
            continue
        row = lr["row"]
        if lr["kind"] == "footer":
            ax.add_patch(Rectangle((0, i - 0.5), L["col_end"], 1,
                                   facecolor="#eef4f9", edgecolor="none",
                                   zorder=0))
        elif row["pc_fail"]:
            ax.add_patch(Rectangle((0, i - 0.5), L["col_end"], 1,
                                   facecolor=TINT_PC_FAIL, edgecolor="none",
                                   zorder=0))
        elif row["design_conflict"]:
            ax.add_patch(Rectangle((0, i - 0.5), L["col_end"], 1,
                                   facecolor=TINT_DESIGN_AMBIGUITY,
                                   edgecolor="none", zorder=0))
        # Block boundary marker (skip before first row)
        if lr["block_key"] != prev_block and i > 0:
            ax.plot([0, L["col_end"]], [i - 0.5, i - 0.5],
                    color="#555", linewidth=0.7, zorder=0.6)
        prev_block = lr["block_key"]

    # Column headers above first row.  Header fonts shrink in the
    # expanded layout because cells are narrower (0.40 vs 0.55) and
    # the snRNA + chromatin sub-blocks crowd the available width.
    header_y = -1.15
    cw = L["cell_width"]
    half_cw = cw / 2.0
    expanded = L["include_chromatin_cells"]
    metric_h_size  = 4.2 if expanded else 5.0
    metric_u_size  = 3.6 if expanded else 4.0
    chrom_h_size   = 4.0 if expanded else 4.5
    chrom_u_size   = 3.4 if expanded else 3.8
    lib_label_size = 5.0 if expanded else 5.5
    ax.text(COL_GENO + 0.05, header_y, "Genotype",
            fontsize=5.5, ha="left", va="bottom", fontweight="bold")
    ax.text(COL_ID + 0.05, header_y, "Embryo ID",
            fontsize=5.5, ha="left", va="bottom", fontweight="bold")
    # "Library A (alph.)" trimmed to "Library A" in expanded mode so the
    # parenthetical doesn't collide with the first metric column header.
    lib_a_text = "Library A" if expanded else "Library A (alph.)"
    ax.text(L["col_a_lib"] + 0.05, header_y, lib_a_text,
            fontsize=lib_label_size, ha="left", va="bottom", fontweight="bold")
    ax.text(L["col_b_lib"] + 0.05, header_y, "Library B",
            fontsize=lib_label_size, ha="left", va="bottom", fontweight="bold")
    for side_col_m in [L["col_a_m_rna"], L["col_b_m_rna"]]:
        for x, ml, unit in zip(side_col_m, METRIC_LABELS, metric_units):
            ax.text(x + half_cw, header_y, ml, fontsize=metric_h_size,
                    ha="center", va="bottom", fontweight="bold")
            if unit:
                ax.text(x + half_cw, header_y + 0.27, unit,
                        fontsize=metric_u_size,
                        ha="center", va="bottom", color="#666", style="italic")
    if expanded:
        for side_col_m in [L["col_a_m_chrom"], L["col_b_m_chrom"]]:
            for x, ml, unit in zip(side_col_m, METRIC_CHROM_LABELS,
                                    METRIC_CHROM_UNITS):
                ax.text(x + half_cw, header_y, ml, fontsize=chrom_h_size,
                        ha="center", va="bottom", fontweight="bold",
                        color="#1f6029")
                if unit:
                    ax.text(x + half_cw, header_y + 0.27, unit,
                            fontsize=chrom_u_size, ha="center", va="bottom",
                            color="#1f6029", style="italic")
    ax.text(L["col_a_n"] + 0.4, header_y, "cells", fontsize=5,
            ha="center", va="bottom", fontweight="bold")
    ax.text(L["col_b_n"] + 0.4, header_y, "cells", fontsize=5,
            ha="center", va="bottom", fontweight="bold")
    ax.text(L["col_call"] + 0.3, header_y, "call", fontsize=5.5,
            ha="center", va="bottom", fontweight="bold")

    stage_color = {s: STAGE_COLORS[s] for s in STAGE_BUCKETS}

    # Modality-strip header label (rendered only when figS03 panel C enables it).
    # The strip itself is drawn per-row inside _render_data_row.
    if show_modality_strip:
        ax.text(COL_MODALITY + 0.025, header_y, "ATAC",
                fontsize=4, ha="left", va="bottom",
                color="#2a7d34", fontweight="bold", rotation=0)

    # Track (block_key, A_side) to draw a hairline divider where the L→R
    # transition happens within a block (O3 visibility).
    prev_key = None
    prev_side = None
    for i, lr in enumerate(layout_rows):
        if lr["kind"] == "stage_subtitle":
            _render_stage_subtitle(ax, i, lr["stage"], L)
            prev_key = lr["block_key"]
            prev_side = None
            continue
        if lr["kind"] == "footer":
            _render_footer(ax, i, lr["row"], L)
            prev_key = None
            prev_side = None
            continue
        row = lr["row"]
        a_side = row["A_side"]
        cur_side = a_side if isinstance(a_side, str) else ""
        if (prev_key == lr["block_key"] and prev_side and cur_side
                and prev_side != cur_side):
            # Subtle hairline between L-rows and R-rows within a block
            ax.plot([COL_GENO, L["col_end"]], [i - 0.5, i - 0.5],
                    color="#d8d8d8", linewidth=0.4, linestyle=":",
                    zorder=0.4)
        _render_data_row(ax, i, row, stage_color, max_n_cells,
                         metric_col_max, L,
                         show_modality_strip=show_modality_strip)
        prev_key = lr["block_key"]
        prev_side = cur_side

    # Legend
    cell_legend = [
        Patch(facecolor=SIDE_COLOR["L"], edgecolor="#888",
              label="Left-side library  (PITX2 higher)"),
        Patch(facecolor=SIDE_COLOR["R"], edgecolor="#888",
              label="Right-side library"),
        Patch(facecolor=INCONCLUSIVE_GRAY, edgecolor="#888", alpha=0.6,
              label="inconclusive / singleton magnitude"),
        Patch(facecolor="white", edgecolor="#aaa",
              label="loses on this metric (paired)"),
        Patch(facecolor="white", edgecolor=LOW_N_HATCH_COLOR,
              hatch="////",
              label=f"slot has n<{LOW_N_THRESHOLD} cells (low confidence)"),
    ]
    call_legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#333",
               markeredgecolor="none", markersize=6, linestyle="",
               label="● direct  (4-metric firm)"),
        Line2D([0], [0], marker="$◐$", color="#555", markersize=8,
               linestyle="", label="◐ indirect  (placement + singleton rank)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#fff",
               markeredgecolor="#aaa", markersize=6, linestyle="",
               markeredgewidth=1.2, label="○ weak  (no side call)"),
        Line2D([0], [0], marker="$=$", color="#333", markersize=8,
               linestyle="",
               label="= same-side duplicate  (split tube; design side)"),
        Line2D([0], [0], marker="$*$", color="#333", markersize=8,
               linestyle="",
               label="* embryo ID swappable with peer genotype"),
        Line2D([0], [0], marker="$✓$", color="#2a7d34", markersize=8,
               linestyle="", label="✓ block pair-constraint OK"),
        Line2D([0], [0], marker="$⚠$", color="#c94438", markersize=8,
               linestyle="", label="⚠ block pair-constraint mismatch"),
        Patch(facecolor=TINT_PC_FAIL, edgecolor="#c94438",
              label="red row — block pair-constraint fail; flags "
                    "ALL block-mates (single-modality rollup; "
                    # NOTE: this legend STRING is RENDERED in figS02D (show_legend=True).
                    # Its "fig03 panel E" cross-ref is stale (the matrix is now figS03
                    # panel C), but editing it changes fig02_supp's PNG — a figure this
                    # figS03 PR does not otherwise touch.  Left verbatim on purpose;
                    # fix in a fig02_supp-scoped change so fig02_supp stays byte-identical there.
                    "chromatin resolves in fig03 panel E)"),
        Patch(facecolor=TINT_DESIGN_AMBIGUITY, edgecolor="#b25e09",
              label="orange row — design ambiguity "
                    "(asym multi-candidate, no peer)"),
    ]
    # Default True: figS02D renders the full two-row decoding legend so its
    # standalone panel stays self-decoding.  figS03 passes show_legend=False
    # because its house rule keeps all descriptive text in LEGEND.md (the
    # cell-colour / glyph / tint key is documented there) — and the two-row
    # legend is the panel's busiest element.
    if show_legend:
        fig.legend(handles=cell_legend, loc="upper center",
                   bbox_to_anchor=(0.5, 1.03),
                   frameon=False, fontsize=6, ncol=4,
                   handlelength=1.2, columnspacing=1.4,
                   title="Metric cells  (intensity ∝ log₂ fold vs partner library; "
                         f"saturates at ≥{int(2**EVIDENCE_LOG2_SATURATE)}x fold).  "
                         "Cell-count bars scaled to global max across the panel.",
                   title_fontsize=6)
        fig.legend(handles=call_legend, loc="upper center",
                   bbox_to_anchor=(0.5, 1.00),
                   frameon=False, fontsize=6, ncol=5,
                   handlelength=1.2, columnspacing=1.2,
                   title="Call column", title_fontsize=6)

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Wrote: {out_path}  ({len(genos)} genotype rows, "
          f"{len(footers)} blocks)")


def _render_stage_subtitle(ax, y_idx: int, stage: str, L: dict) -> None:
    """Stage-section subtitle: declares modality + per-stage column units.

    Used by figS03 panel C to make the modality (snRNA vs snATAC) and the
    correct unit triple per stage explicit at the header level.  The
    subtitle is rendered as a single emphasized row with a light
    background tint (already drawn in the row-tint pass) — the green
    tint flags ATAC, the cool-gray tint flags snRNA.

    In the expanded-columns layout, the snRNA subtitle annotates that the
    chromatin sub-block will be empty for that stage; the ATAC subtitle
    annotates that the snRNA sub-block will be empty.  The dual-modality
    layout naturally surfaces "this row's call comes from <X>" by which
    sub-block is filled.
    """
    y = y_idx
    if stage == "HH17/16":
        modality_label = "snATAC chromatin (canonical) + snRNA"
        modality_color = "#1f6029"
        if L["include_chromatin_cells"]:
            unit_text = ("snRNA at noise floor (PITX2c median 0%); "
                         "chromatin carries the L/R call")
        else:
            unit_text = ("gene-body open %  |  PITX2c open %  |  "
                         "PITX2c frags / cell")
    else:
        # Single string across HH13/HH21/HH25 — actual paired Wilcoxon
        # p-values from f5 stats TSV after MIN_CELLS=20 filter:
        #   HH13: p=0.035  HH21: p=0.014  HH25: p=0.014  (PITX2c_promoter, paired_wilcox_p
#   in step_f5_pitx2_region_stats.tsv at MIN_CELLS=20 — recomputed 2026-08-18; the
#   0.022 this line carried was never in that table)
        # "p≤0.035" is the conservative single-string summary.  Stable
        # against future MIN_CELLS tweaks at the visualization layer
        # because the worst flanking stage is what gets reported here.
        modality_label = ("snRNA call  +  chromatin corroborates "
                          "(paired Wilcoxon p≤0.035)")
        modality_color = "#444"
        if L["include_chromatin_cells"]:
            unit_text = ("snRNA carries the L/R call; chromatin "
                         "frags/cell corroborates")
        else:
            unit_text = "PITX2+ %  |  PITX2c+ %  |  mean CP10K"

    # Modality label on the upper half of the row; unit text on the
    # lower half — separating them on different y-positions prevents
    # the label from running into the unit text horizontally when the
    # label string is long (e.g., HH13/HH21/HH25's chromatin-
    # corroborates annotation).
    ax.text(COL_GENO + 0.05, y - 0.18,
            f"── {stage} ──  {modality_label}",
            fontsize=5.8, ha="left", va="center", fontweight="bold",
            color=modality_color)
    ax.text(COL_GENO + 0.05, y + 0.20, unit_text,
            fontsize=4.6, ha="left", va="center",
            color="#555", style="italic")


def _render_data_row(ax, y_idx: int, r, stage_color, max_n_cells,
                     metric_col_max, L: dict,
                     show_modality_strip: bool = False) -> None:
    y = y_idx
    ax.add_patch(Rectangle((COL_STAGE + 0.02, y - 0.32),
                           0.25, 0.64,
                           facecolor=stage_color[r["design_stage"]],
                           edgecolor="none", clip_on=False))

    # Modality strip — figS03 panel C only.  HH17/16 rows are chromatin-
    # derived calls (d2n); the green strip flags this visually so readers
    # don't have to map "design_stage=HH17/16" to "modality=ATAC" via the
    # legend.  Other stages (snRNA-derived from d2l) get a faint gray
    # placeholder so the column has a consistent visual width across rows.
    if show_modality_strip:
        is_chrom = r["design_stage"] == "HH17/16"
        strip_color = "#2a7d34" if is_chrom else "#dddddd"
        ax.add_patch(Rectangle((COL_MODALITY - 0.005, y - 0.32),
                               0.05, 0.64,
                               facecolor=strip_color, edgecolor="none",
                               clip_on=False))

    ax.text(COL_GENO + 0.05, y, display_name(r["genotype"]),
            fontsize=6, ha="left", va="center")

    # Embryo ID (A(ii)): single embryo + asterisk for interchangeable peers.
    # Color matches the row tint:
    #   pc_fail   → red (#c94438), normal weight
    #   design ambiguity (no pc_fail) → orange-brown (#b25e09)
    #   else (incl. peer-asterisk) → dark gray, bold
    emb = r["embryo_id_display"] if isinstance(r["embryo_id_display"], str) else ""
    has_peer = bool(r["interchangeable_with"])
    asterisk = "*" if has_peer else ""
    pc_fail = r["pc_fail"]
    conflict = r["design_conflict"]
    if pc_fail:
        id_color = "#c94438"
        id_bold = False
    elif conflict:
        id_color = "#b25e09"
        id_bold = False
    else:
        id_color = "#333"
        id_bold = True
    ax.text(COL_ID + 0.05, y, f"{emb}{asterisk}",
            fontsize=6, ha="left", va="center", color=id_color,
            fontweight="bold" if id_bold else "normal")

    cw = L["cell_width"]
    half_cw = cw / 2.0
    slot_iter = [
        ("A", L["col_a_lib"], L["col_a_m_rna"], L["col_a_m_chrom"], L["col_a_n"]),
        ("B", L["col_b_lib"], L["col_b_m_rna"], L["col_b_m_chrom"], L["col_b_n"]),
    ]
    for slot, col_lib, col_m_rna, col_m_chrom, col_n in slot_iter:
        # pivot_by_genotype guarantees {slot}_lib/_side/_n/_cohort_indirect
        # are always set (to None for singleton-B), so direct bracket
        # access is fail-fast: a missing key surfaces a pivot regression
        # rather than silently defaulting (the originating .get(default)
        # bug class).
        lib = r[f"{slot}_lib"]
        side = r[f"{slot}_side"]
        n_cells = r[f"{slot}_n"]
        if lib is None or (isinstance(lib, float) and pd.isna(lib)):
            ax.text(col_lib + 0.05, y, "—",
                    fontsize=5.5, ha="left", va="center", color="#bbb")
            continue
        # Cohort-indirect soft-call override: pivot_by_genotype reads
        # the optional `cohort_indirect_side` column from the multiome
        # sidecar (stamped at d2o time from
        # references/genotype_overrides.tsv) and translates it into
        # this per-slot flag + the overridden side string.  When the
        # flag is True, paint "=L?" / "=R?" in design-ambiguity orange.
        # figS02D's d2l input has no override column, so the flag is
        # always False there and this branch never fires.
        slot_ci = bool(r[f"{slot}_cohort_indirect"])
        if slot_ci and side in ("L", "R"):
            side_str = f"={side}?"
            side_color = "#b25e09"
            side_bold = False
        else:
            side_str = f"={side}" if side in ("L", "R") else "=?"
            side_color = SIDE_COLOR[side]
            side_bold = True
        ax.text(col_lib + 0.05, y, lib,
                fontsize=5.5, ha="left", va="center")
        # Side text right-aligned at the end of the lib name slot.  In the
        # compact layout the slot is ~1.13 wide (offset 1.08); in the
        # expanded layout the slot is also ~1.05 wide so the same offset
        # works.
        ax.text(col_lib + 1.05, y, side_str,
                fontsize=5.5, ha="right", va="center",
                color=side_color,
                fontweight="bold" if side_bold else "normal",
                fontstyle="italic" if not side_bold else "normal")

        partner_slot = "B" if slot == "A" else "A"
        partner_lib = r[f"{partner_slot}_lib"]
        # A same-side duplicate's partner library is the SAME tissue: there
        # is no L/R contrast to paint, so its cells render as magnitude-only
        # (singleton style) rather than paired win/lose.
        is_singleton = (partner_lib is None or
                        (isinstance(partner_lib, float) and pd.isna(partner_lib))
                        or r["firmness"] == SAME_SIDE_DUP)
        # Low-n cue applies to THIS slot only — partner slot may be fine.
        is_low_n = (n_cells is not None and not pd.isna(n_cells)
                    and int(n_cells) < LOW_N_THRESHOLD)

        def _draw_metric_cell(mx, m, lib_val, is_singleton, side, partner_val,
                              col_max):
            ax.add_patch(Rectangle((mx, y - 0.35), cw, 0.7,
                                   facecolor="white",
                                   edgecolor="#dddddd", linewidth=0.4))
            if is_singleton:
                cell = singleton_cell(lib_val, col_max)
            else:
                cell = paired_cell(lib_val, partner_val, side)
            cell_alpha_local = 0.0
            if cell is not None:
                color, cell_alpha_local = cell
                ax.add_patch(Rectangle((mx, y - 0.35), cw, 0.7,
                                       facecolor=color, alpha=cell_alpha_local,
                                       edgecolor="#888888", linewidth=0.3))
            if is_low_n:
                ax.add_patch(Rectangle((mx, y - 0.35), cw, 0.7,
                                       facecolor="none",
                                       edgecolor=LOW_N_HATCH_COLOR,
                                       hatch="////", linewidth=0.0))
            if lib_val is None or pd.isna(lib_val):
                val_str = ""
            elif m in ("pitx2_pos", "pitx2c_pos",
                       "pitx2_genebody_open_frac", "pitx2c_open_frac"):
                val_str = (f"{100 * float(lib_val):.0f}%"
                           if float(lib_val) >= 0.005 else "0")
            else:
                val_str = f"{float(lib_val):.2f}"
            if val_str:
                text_color = "white" if cell_alpha_local > 0.55 else "#222"
                ax.text(mx + half_cw, y, val_str,
                        fontsize=4.3, ha="center", va="center",
                        color=text_color,
                        fontstyle="italic" if is_low_n else "normal")

        # snRNA metric block — always drawn.  At HH17/16 the values are
        # NaN (in expanded layout) or substituted from chromatin (compact).
        for k, (m, mx) in enumerate(zip(METRICS, col_m_rna)):
            lib_val = r[f"{slot}_{m}"]
            partner_val = (None if is_singleton
                           else r[f"{partner_slot}_{m}"])
            _draw_metric_cell(mx, m, lib_val, is_singleton, side,
                              partner_val, metric_col_max.get(m, 0.0))

        # Chromatin metric block — only in expanded layout.  Empty cells
        # at HH13/HH21/HH25 (no chromatin call); populated at HH17/16.
        if L["include_chromatin_cells"] and col_m_chrom is not None:
            for k, (m, mx) in enumerate(zip(METRICS_CHROM, col_m_chrom)):
                lib_val = r[f"{slot}_{m}"]
                partner_val = (None if is_singleton
                               else r[f"{partner_slot}_{m}"])
                _draw_metric_cell(mx, m, lib_val, is_singleton, side,
                                  partner_val, metric_col_max.get(m, 0.0))

        if n_cells is not None and not pd.isna(n_cells):
            # Fixed-width 0.60-unit bar column; bar grows right from the
            # left edge, numeric label right-aligned at the column's right
            # edge so text never drifts outside the column.
            bar_width = 0.60
            frac = float(n_cells) / max_n_cells
            bar_color = LOW_N_TEXT_COLOR if is_low_n else "#555555"
            ax.add_patch(Rectangle((col_n + 0.02, y - 0.15),
                                   bar_width * frac, 0.30,
                                   facecolor=bar_color, edgecolor="none"))
            count_label = (f"{int(n_cells):,} low" if is_low_n
                           else f"{int(n_cells):,}")
            ax.text(col_n + bar_width + 0.20, y + 0.02,
                    count_label, fontsize=4,
                    ha="right", va="center",
                    color=LOW_N_TEXT_COLOR if is_low_n else "#555",
                    fontweight="bold" if is_low_n else "normal")

    # Cohort-indirect glyph override: when ANY of this genotype's slots
    # carries the soft-call flag, replace the d2l/d2n firmness glyph
    # (typically ○ "weak" for these rows) with ◐ in the design-
    # ambiguity orange tint, matching the side-text override above.
    # The base firmness column still flows through as the d2l/d2n
    # value; this is a render-time annotation, not a data change.
    any_ci = bool(r["A_cohort_indirect"]) or bool(r["B_cohort_indirect"])
    if any_ci and r["firmness"] in ("weak", None):
        glyph, c = "◐", "#b25e09"
    else:
        glyph, c = firmness_glyph(r["firmness"])
    ax.text(L["col_call"] + 0.3, y, glyph,
            fontsize=9, ha="center", va="center", color=c)


def _has_lib(x) -> bool:
    """True iff x is a non-empty library name (not None, not NaN)."""
    if x is None:
        return False
    if isinstance(x, float) and pd.isna(x):
        return False
    return bool(x)


def _render_footer(ax, y_idx: int, footer, L: dict) -> None:
    """B(x) footer row: per-library inferred vs design L/R counts."""
    y = y_idx
    badge_ok = footer["pair_constraint_ok"]
    badge_color = "#2a7d34" if badge_ok else "#c94438"
    badge = "✓" if badge_ok else "⚠"

    # Red bar on the LEFT edge of failure footers — high-contrast flag
    if not badge_ok:
        ax.add_patch(Rectangle((COL_STAGE, y - 0.5),
                               0.15, 1,
                               facecolor="#c94438", edgecolor="none", zorder=0.8))

    a_lib = footer["A_lib"] if _has_lib(footer["A_lib"]) else ""
    b_lib = footer["B_lib"] if _has_lib(footer["B_lib"]) else ""
    if b_lib:
        label = f"{a_lib} ↔ {b_lib}  ({footer['design_stage']})"
    else:
        label = f"{a_lib} (single library, {footer['design_stage']})"
    ax.text(COL_GENO + 0.05, y, label, fontsize=5.8,
            ha="left", va="center", color="#333",
            fontweight="bold")

    a_inf = f"{footer['A_n_L_inf']}L/{footer['A_n_R_inf']}R"
    a_des = f"{footer['A_n_L_des']}L/{footer['A_n_R_des']}R"
    a_color = "#333" if a_inf == a_des else "#c94438"
    ax.text(L["col_a_lib"] + 0.05, y, f"inferred {a_inf}  design {a_des}",
            fontsize=5.2, ha="left", va="center",
            color=a_color, fontweight="bold" if a_inf != a_des else "normal")

    if b_lib:
        b_inf = f"{footer['B_n_L_inf']}L/{footer['B_n_R_inf']}R"
        b_des = f"{footer['B_n_L_des']}L/{footer['B_n_R_des']}R"
        b_color = "#333" if b_inf == b_des else "#c94438"
        ax.text(L["col_b_lib"] + 0.05, y, f"inferred {b_inf}  design {b_des}",
                fontsize=5.2, ha="left", va="center",
                color=b_color, fontweight="bold" if b_inf != b_des else "normal")

    # Badge: larger + bolder when failing
    badge_size = 12 if not badge_ok else 9
    ax.text(L["col_call"] + 0.3, y, badge,
            fontsize=badge_size, ha="center", va="center",
            color=badge_color, fontweight="bold")
