"""Shared colors, palettes, themes, and marker genes for all figures.

Color handling: this module exposes legacy dict names (``CELLTYPE_COLORS``,
``STAGE_COLORS``, ``SEX_COLORS``, etc.) that downstream figure scripts import.
The values are sourced from ``references/color_registry.tsv`` via the helper
``figures._shared.colors``. Do NOT edit hex strings here -- update the TSV
and the dicts below will pick up the new values automatically.
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import OrderedDict

try:
    # Normal case when figure scripts add figures/ to sys.path first.
    from _shared.colors import get_palette as _reg_palette
except ImportError:
    # Fallback when style.py is imported as figures._shared.style (tests, etc.)
    from figures._shared.colors import get_palette as _reg_palette  # noqa: F401
try:
    from _shared.library_naming import LANE_0S_SHORT
except ImportError:
    from figures._shared.library_naming import LANE_0S_SHORT


# ── Cluster annotations (loaded from references/cluster_annotations.tsv) ────
# Single source of truth for cluster→cell type mapping.
# The TSV is also synced to scratch for FC_08656/ export scripts.
def _load_cluster_annotations():
    """Load CLUSTER_ANNOTATIONS from the canonical TSV."""
    tsv = os.path.join(os.path.dirname(__file__), "..", "..",
                       "references", "cluster_annotations.tsv")
    if not os.path.isfile(tsv):
        raise FileNotFoundError(
            f"cluster_annotations.tsv not found at {tsv}\n"
            "This file is the single source of truth for cluster labels."
        )
    df = pd.read_csv(tsv, sep="\t")
    return OrderedDict(
        (str(row["cluster"]), row["annotation"]) for _, row in df.iterrows()
    )

# ── Library order & labels ───────────────────────────────────────────────
LIB_ORDER = [
    "GEX03-1A", "GEX03-1B", "GEX03-2", "GEX04", "GEX0S-duck",
    "GEX05", "GEX06", "GEX07", "GEX08",
    "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B",
]
LIB_LABELS = [
    "03-1A", "03-1B", "03-2", "04", LANE_0S_SHORT["duck"],
    "05", "06", "07", "08",
    "PG3", "PG4A", "PG4B",
]

# ── Color palettes ───────────────────────────────────────────────────────
# All hex values come from references/color_registry.tsv via _shared/colors.
# To change a color, edit the TSV (NOT this file).

LIB_COLORS = _reg_palette("library")

SEX_COLORS = _reg_palette("sex")

STAGE_COLORS = _reg_palette("stage")
# After inference (e1b), HH21/16 resolves to HH21 or HH16; HH17 is retained
# as design-confirmed. Inferred stages: HH25, HH21, HH17, HH17/16, HH16, HH13.

FILTER_STAGE_COLORS = _reg_palette("filter_stage")

# ── Genotype metadata (loaded from references/genotype_metadata.tsv) ─────
_GENOTYPE_TSV_REQUIRED_COLUMNS = (
    "genotype_id", "display_name", "group", "stage", "source",
)


def _read_genotype_metadata_tsv():
    """Read + validate references/genotype_metadata.tsv.

    Fail-fast on missing file, missing columns, duplicate genotype_id rows,
    or empty cells in required columns. The TSV is the single source of truth;
    silent absorption of a malformed row would defeat the SoT contract.
    """
    tsv = os.path.join(os.path.dirname(__file__), "..", "..",
                       "references", "genotype_metadata.tsv")
    if not os.path.isfile(tsv):
        raise FileNotFoundError(
            f"genotype_metadata.tsv not found at {tsv}\n"
            "This file is checked into git as the single source of truth — "
            "verify your repo checkout includes references/."
        )
    df = pd.read_csv(tsv, sep="\t", dtype=str, keep_default_na=False)
    missing_cols = [c for c in _GENOTYPE_TSV_REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"genotype_metadata.tsv missing required columns: {missing_cols}"
        )
    if not df["genotype_id"].is_unique:
        dups = df.loc[df["genotype_id"].duplicated(keep=False), "genotype_id"].tolist()
        raise ValueError(
            f"genotype_metadata.tsv has duplicate genotype_id rows: {dups}"
        )
    blanks = df.loc[:, list(_GENOTYPE_TSV_REQUIRED_COLUMNS)].eq("").any(axis=1)
    if blanks.any():
        bad = df.loc[blanks, "genotype_id"].tolist()
        raise ValueError(
            f"genotype_metadata.tsv has empty cells in required columns "
            f"for genotype_id(s): {bad}"
        )
    return df


def _load_genotype_meta():
    """Load GENOTYPE_META: genotype_id -> {group, stage, source}."""
    df = _read_genotype_metadata_tsv()
    return OrderedDict(
        (row["genotype_id"],
         {"group": row["group"], "stage": row["stage"], "source": row["source"]})
        for _, row in df.iterrows()
    )


def _load_genotype_display_names():
    """Load DISPLAY_NAMES: genotype_id -> publication-ready label.

    Includes derived aliases for known pipeline-prefix variants
    (G3_unfiltered_donorN -> G3_donorN; dp500_donor29 -> G6_dp500_donor29).
    """
    df = _read_genotype_metadata_tsv()
    names = {}
    for _, row in df.iterrows():
        gid = row["genotype_id"]
        dname = row["display_name"]
        names[gid] = dname      # identity rows too (g00 -> g00), so display_name() can be strict
        if gid.startswith("G3_unfiltered_"):
            names[gid.replace("G3_unfiltered_", "G3_")] = dname
        if gid.startswith("dp500_"):
            names["G6_" + gid] = dname
    return names


GENOTYPE_META = _load_genotype_meta()
DISPLAY_NAMES = _load_genotype_display_names()

# Display-only aliases for alternate demux outputs (not among the 52 canonical
# embryos). Kept here (not in the TSV) because they don't correspond to
# confirmed embryos — only to intermediate pipeline donor IDs.
_ALT_DISPLAY_NAMES = {
    # G9 unlinked (PILOT4A representations, no G6 counterpart)
    "G9_donor1": "G9-d1", "G9_donor2": "G9-d2", "G9_donor11": "G9-d11",
    # G2 tri_dp50 donors (stage/identity not yet resolved)
    "G2_donor0": "G2-a", "G2_donor2": "G2-b", "G2_donor3": "G2-c",
    "G2_donor7": "G2-d", "G2_donor9": "G2-e", "G2_donor10": "G2-f",
    "G2_donor11": "G2-g",
}
# Fail-fast on overlap: the TSV owns canonical display names; ALT must stay
# disjoint. If a future TSV row uses an ALT key, remove the ALT entry rather
# than letting it shadow the TSV silently.
_alt_collisions = set(DISPLAY_NAMES) & set(_ALT_DISPLAY_NAMES)
if _alt_collisions:
    raise ValueError(
        f"_ALT_DISPLAY_NAMES collides with TSV-loaded DISPLAY_NAMES: "
        f"{sorted(_alt_collisions)}. Move these entries into "
        f"references/genotype_metadata.tsv."
    )
DISPLAY_NAMES.update(_ALT_DISPLAY_NAMES)


# ── Unified donor mapping ────────────────────────────────────────────────
# The (demux_group, raw_donor_id) -> unified_id mapping is loaded from
# unified_donor_map.tsv. Resolution is DELIBERATELY asymmetric with the demux
# chain (pipeline/c2 reads the scratch chain product first): figure code has no
# live demux chain, so it prefers the committed references/ copy and falls back
# to scratch only when the repo copy is absent. No value-level fallback — a
# missing TSV raises.
def load_unified_donor_map(tsv_path=None):
    """Load the unified donor mapping from TSV (produced by b1_donor_consolidation.py).

    Returns dict[(demux_group, raw_donor_id)] -> unified_id
    for the confirmed embryos (excludes noise/duplicates).
    Raises FileNotFoundError if the TSV does not exist.
    """
    import functools
    return _load_unified_donor_map_cached(tsv_path)


@__import__("functools").lru_cache(maxsize=4)
def _load_unified_donor_map_cached(tsv_path=None):
    import pandas as pd, os
    if tsv_path is None:
        # the committed SSOT only (no scratch fallback)
        tsv_path = os.path.join(os.path.dirname(__file__), "..", "..",
                                "references", "unified_donor_map.tsv")
    if not os.path.isfile(tsv_path):
        raise FileNotFoundError(
            f"unified_donor_map.tsv not found: {tsv_path}\n"
            "Run pipeline/b1_donor_consolidation.py first to generate it."
        )
    df = pd.read_csv(tsv_path, sep="\t")
    df = df[df["unified_id"] != "EXCLUDE"]
    return {
        (row["demux_group"], row["raw_donor_id"]): row["unified_id"]
        for _, row in df.iterrows()
    }


def get_unified_id(demux_group, donor_id, donor_map=None):
    """Map (demux_group, donor_id) to unified genotype ID.

    Strict lookup — returns NaN for any pair not in the donor map.
    This guarantees the h5ad can only contain the 52 confirmed embryos.

    Args:
        donor_map: Pre-loaded map dict. If None, loads from TSV (or fallback).
    """
    import pandas as pd, numpy as np
    if pd.isna(demux_group) or pd.isna(donor_id):
        return np.nan
    if donor_map is None:
        donor_map = load_unified_donor_map()
    return donor_map.get((demux_group, donor_id), np.nan)


CLUSTER_ANNOTATIONS = _load_cluster_annotations()

# Cell type color palette — sourced from color_registry.tsv (res=0.2 cell types only;
# the res=0.5 Level-2 mesenchyme subtype color merge was retired 2026-06-13 with res=0.5).
#
# Paraxial lineage (Sclerotome, Dermomyotome) uses a cool teal/slate family
# chosen to contrast with the warm LPM browns under deuteranopia/protanopia
# simulation (see references/color_registry.tsv rationale column). All other
# hexes are unchanged from the pre-registry palette.
CELLTYPE_COLORS = {
    **_reg_palette("celltype"),
}

# Stacking order for composition plots (LPM lineage first — paper focus).
# Order matches _PREFERRED_ORDER so the legend reads LPM + Splanchnic first,
# then paraxial (Sclerotome, Dermomyotome), then neural/other lineages.
CELLTYPE_STACK_ORDER = [
    "LPM mesenchyme", "Splanchnic mesenchyme", "Sclerotome", "Dermomyotome",
    "Neural progenitor", "Neural crest", "Diff. neuron",
    "Dorsal neural prog.", "Ventral neural prog.",
    "Epithelial", "Pharyngeal endoderm", "Otic placode",
    "Endothelial", "Notochord", "Cardiomyocyte", "RBC",
]

# ── Res=0.2 annotation (primary for main figures) ────────────────────────
# cluster_annotations.tsv is the SINGLE SOURCE OF TRUTH for annotation names.
# No remap, no aliases. If a data file contains old names, it must be re-exported.

# Display order: AUTO-DERIVED from cluster_annotations.tsv.
# _PREFERRED_ORDER controls lineage grouping. Annotations not listed here
# are appended at the end and trigger a warning — this catches pipeline
# reruns that produce new cluster identities.
_PREFERRED_ORDER = [
    "LPM mesenchyme", "Splanchnic mesenchyme", "Sclerotome", "Dermomyotome",
    "Neural progenitor", "Neural crest", "Diff. neuron",
    "Dorsal neural prog.", "Ventral neural prog.",
    "Epithelial", "Pharyngeal endoderm", "Otic placode",
    "Endothelial", "Notochord", "Cardiomyocyte", "RBC",
]

def _build_celltype_02_order():
    """Build display order from cluster_annotations.tsv, sorted by preferred order.

    Raises RuntimeError if any annotation has no color entry.
    Warns if annotations appear that aren't in _PREFERRED_ORDER.
    """
    ann_types = list(dict.fromkeys(CLUSTER_ANNOTATIONS.values()))  # unique, preserve order
    ordered = [t for t in _PREFERRED_ORDER if t in ann_types]
    remaining = [t for t in ann_types if t not in ordered]

    # Hard fail: every annotation MUST have a color
    missing_colors = [t for t in ann_types if t not in CELLTYPE_COLORS]
    if missing_colors:
        raise RuntimeError(
            f"Annotations in cluster_annotations.tsv have no CELLTYPE_COLORS entry: "
            f"{missing_colors}. Add them to CELLTYPE_COLORS in style.py before proceeding."
        )

    # Hard fail: every annotation must have a display position
    if remaining:
        raise RuntimeError(
            f"Annotations not in _PREFERRED_ORDER: {remaining}. "
            f"Add them to _PREFERRED_ORDER in style.py."
        )
    return ordered

CELLTYPE_02_ORDER = _build_celltype_02_order()


def validate_celltype_coverage(df, celltype_col="celltype"):
    """Check that all cells in df have a recognized cell type annotation.

    Call this in figure scripts after mapping cluster numbers to annotations.
    Raises RuntimeError if any cells have types not in CELLTYPE_02_ORDER,
    reporting the unaccounted types and cell counts.
    """
    known = set(CELLTYPE_02_ORDER)
    types_in_data = set(df[celltype_col].dropna().unique())
    unaccounted = types_in_data - known
    if unaccounted:
        counts = df[celltype_col].value_counts()
        detail = {t: int(counts.get(t, 0)) for t in sorted(unaccounted)}
        raise RuntimeError(
            f"Unaccounted cell types in data (not in CELLTYPE_02_ORDER): {detail}. "
            f"Update cluster_annotations.tsv and/or CELLTYPE_COLORS in style.py."
        )
    n_na = df[celltype_col].isna().sum()
    if n_na > 0:
        raise RuntimeError(
            f"{n_na:,} cells have NaN cell type annotation. "
            f"Check that all cluster numbers in the data map to cluster_annotations.tsv."
        )

# ── res=0.5 mesenchyme subclustering FULLY RETIRED 2026-06-13 (git history only) ──
# Removed with the rest of res=0.5: MESENCHYME_CL05, MESENCHYME_LINEAGE, PITX2_TIERS /
# PITX2_EXCLUDED (laterality tiers), MESENCHYME_RES02, HYBRID_CELLTYPE_COLORS /
# HYBRID_ANNOTATION_ORDER, compute_hybrid_annotation, and the Level-2 mesenchyme color merge.
# The non-seed-stable global res=0.5 (cluster_05) grouping diverged from the d2g pools that
# fig02/fig03 use and manufactured a spurious HH17/16 "inconsistency". The laterality /
# subcluster grouping is the d2g stage-restricted subclustering
# (figures/_shared/subcluster_map.py `POOLS` / `get_pool_map`). Do NOT reintroduce res=0.5.


def display_name(genotype_id):
    """Map internal genotype_id to publication-ready display name.

    Strict: an id absent from references/genotype_metadata.tsv (+ its aliases and
    _ALT_DISPLAY_NAMES) raises instead of leaking the raw pipeline id onto a figure.
    """
    try:
        return DISPLAY_NAMES[genotype_id]
    except KeyError:
        raise KeyError(f"no display name for genotype_id {genotype_id!r} — add it to "
                       "references/genotype_metadata.tsv") from None


# ── W/Z sex determination ────────────────────────────────────────────────
W_THRESHOLD = 0.002  # median percent_w > 0.002 -> Female
MIN_CELLS_VIOLIN = 50
MIN_CELLS_PCA = 500


# ── Plot helpers ─────────────────────────────────────────────────────────
def apply_pub_style():
    """Set publication-style matplotlib defaults."""
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def lib_label(lib_name):
    """Short label for a library name."""
    if lib_name not in LIB_ORDER:
        raise KeyError(f"no short label for library {lib_name!r} — add it to LIB_ORDER/LIB_LABELS")
    return LIB_LABELS[LIB_ORDER.index(lib_name)]
