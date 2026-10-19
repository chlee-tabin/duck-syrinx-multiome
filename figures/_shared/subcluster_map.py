"""Shared loader for d2g stage-restricted subcluster annotations.

Single source of truth = `references/d2g_subcluster_annotations.tsv`
(owned by `pipeline/d2g_validate_subclusters.py`, analogous to d2b's
`cluster_annotations.tsv`). All figure scripts that map d2g leiden IDs to
biological population names MUST import from here rather than hardcode
their own map — otherwise the next pipeline rerun will silently shuffle
cluster IDs out of sync with the code.

Leiden cluster IDs are not seed-stable. Do not hardcode them.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd

# Locate references/d2g_subcluster_annotations.tsv from this file's position
_REPO = Path(__file__).resolve().parent.parent.parent
_PATH = _REPO / "references" / "d2g_subcluster_annotations.tsv"
_CACHE: Optional[pd.DataFrame] = None


def _load() -> pd.DataFrame:
    global _CACHE
    if _CACHE is None:
        if not _PATH.exists():
            raise RuntimeError(
                f"{_PATH} missing — run `pipeline/d2g_validate_subclusters.py` "
                f"after any d2g rerun to regenerate it."
            )
        _CACHE = pd.read_csv(_PATH, sep="\t", dtype={"cluster": int})
    return _CACHE


def get_leiden_col(stage: str) -> str:
    """Return the canonical leiden resolution column for a stage."""
    df = _load()
    cols = df.loc[df["stage"] == stage, "leiden_col"].unique()
    if len(cols) != 1:
        raise ValueError(
            f"stage={stage!r} resolved to {len(cols)} leiden cols: {list(cols)}"
        )
    return str(cols[0])


def get_domain_map(stage: str, laterality_only: bool = False) -> dict[str, str]:
    """Return {str(cluster_id): population} for a stage.

    When `laterality_only=True`, drop populations not in the inference or
    bilateral pools — i.e. keep only laterality-relevant labels.  Excluded
    set is derived from POOLS, so a future POOLS edit (e.g. moving Mesothelium
    in/out of excluded) flows here automatically.
    """
    df = _load()
    sub = df[df["stage"] == stage]
    if sub.empty:
        raise ValueError(f"stage={stage!r} not in annotations table")
    out = {str(r["cluster"]): r["population"] for _, r in sub.iterrows()}
    if laterality_only:
        keep = set(POOLS[stage]["inference"]) | set(POOLS[stage]["bilateral"])
        out = {k: v for k, v in out.items() if v in keep}
    return out


def get_population_clusters(stage: str, population: str) -> list[int]:
    """Return the list of cluster IDs assigned to a (stage, population)."""
    df = _load()
    sub = df[(df["stage"] == stage) & (df["population"] == population)]
    return sorted(sub["cluster"].astype(int).tolist())


def get_all_populations(stage: str) -> list[str]:
    """Return unique biological populations for a stage.

    Excludes the non-biological residual buckets ('ambiguous') and technical
    junk ('Contamination') so callers enumerating real populations don't pick
    them up.
    """
    df = _load()
    sub = df[(df["stage"] == stage)
             & (~df["population"].isin(["ambiguous", "Contamination"]))]
    seen, out = set(), []
    for p in sub["population"]:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


# ── Laterality pool classification (2026-04-20 redesign; 2026-05-04 update) ────
# Per-stage pool membership for L/R inference + bilateral-control panels.
#   inference      = cells pooled across these labels for 4-metric L/R inference
#   bilateral      = dorsal (expected bilateral); run same metrics, expect weak/null
#   excluded       = Skeletogenic LPM (FOXF1-low, HOX-negative, TBX1+/EYA1+
#                    pharyngeal-cranial skeletal mesenchyme — L/R-symmetric,
#                    formerly mislabeled "Somatic"), Sclerotome/Dermomyotome
#                    (paraxial), Mesothelium (coelomic surface), Venous-pole spl.
#                    (off-axis dissection-edge introgression at the cardiac
#                    venous pole), Contamination (ribosomal/erythroid technical
#                    junk), ambiguous
#
# 2026-05-04 update (per docs/lpm_mesenchyme_subcluster_biology.md):
#   - "Splanchnic" → "Venous-pole spl." rename + reclassify as off-axis
#     introgression (van Wijk 2009 PMC2861358 — TBX18+/GATA4+/ZFPM2+
#     non-myocardial splanchnic envelope of the venous pole / proepicardium,
#     anatomically caudal-medial of the heart, NOT in the syrinx-forming
#     respiratory cross-section).  Moved from inference → excluded at all
#     stages.
#   - "Epicardial mes." → "Mesothelium" rename only (still excluded — coelomic
#     surface lining is bilateral; not part of syrinx-forming LPM-lineage
#     laterality).
#
# 2026-06-18 update (annotation vet; see scripts/annotate_d2g_via_compendium.py):
#   - "Somatic" → "Skeletogenic LPM" rename.  The population is FOXF1-low,
#     HOX-negative, TBX1+/EYA1+/EYA2+/EBF2/3+ pharyngeal-cranial SKELETAL
#     mesenchyme (a pre-chondrogenic progenitor, LPM by d2g parent gating),
#     NOT body-wall somatopleure.  It is L/R-symmetric, so it stays EXCLUDED.
#     The lateralized syrinx skeletal lineage is captured separately, in the
#     HH25 inference pool.  (That population was "Chondrogenic LPM" when this
#     note was written; the 2026-07-02 vet folded it into "Syringeal" — see the
#     HH25 entry below.  "Chondrogenic LPM" is retired from the duck MULTIOME only —
#     it stays correct in the chick arm (laterality_fold_methods.md §5.5) and in the
#     old-duck cr9 arm (a real E09-E11 population).)
#   - "Contamination" added to excluded — ribosomal/low-complexity + erythroid
#     clusters reclassified out of "Skeletogenic LPM" by the compendium
#     contamination_check pre-check (marker/QC-based, no hardcoded leiden IDs).
#
# At HH13/HH17/16 the D/V axis is not yet resolved pre-TES; at HH17/16 the
# NKX6-1+ cluster is labeled "Dorsal-biased LPM" (pre-commitment) and serves
# as bilateral control. "Dorsal sub-epi." only exists at HH21/HH25.
POOLS: dict[str, dict[str, list[str]]] = {
    "HH13": {
        "inference":  ["Early LPM", "Cardiopharyngeal LPM"],
        "bilateral":  [],  # D/V unresolved; no dorsal population
        # 2026-07-02 (annotation vet): the HH13 skeletogenic label is softened to
        # "PRRX1+/FOXF1-low mes." (branchiomeric-skeletal identity not earned
        # pre-HH17) — STILL EXCLUDED (FOXF1-low, 7%). HH17/16-HH25 keep
        # "Skeletogenic LPM". Pool membership is unchanged; only the name differs.
        "excluded":   ["PRRX1+/FOXF1-low mes.", "Sclerotome", "Dermomyotome",
                       "Mesothelium", "Venous-pole spl.", "Contamination",
                       "ambiguous"],
    },
    "HH17/16": {
        # 2026-08-09 (seeded d2g rerun): the reproducible res-0.3
        # partition emits NO "Early LPM" cluster at HH17/16 — the old cl1
        # (2,848 cells) merges into Skeletogenic/Venous-pole clusters and is
        # labeled by the compendium accordingly (excluded). Inference pool is
        # therefore Cardiopharyngeal + BARX1+ only (5,416 cells vs old 7,188).
        "inference":  ["Cardiopharyngeal LPM", "BARX1+"],
        "bilateral":  ["Dorsal-biased LPM"],  # pre-commitment NKX6-1+ (n=2,564 post-seed)
        # 2026-07-02 vet: cl4 is a PRRX1+/FOXF1-low fibroblast (branchiomeric gate);
        # cl2/cl9 remain genuine Skeletogenic LPM. Both compartments excluded.
        "excluded":   ["Skeletogenic LPM", "PRRX1+/FOXF1-low mes.",
                       "Sclerotome", "Dermomyotome",
                       "Mesothelium", "Venous-pole spl.", "Contamination",
                       "ambiguous"],
    },
    "HH21": {
        "inference":  ["Ventral spl.", "Syringeal", "BARX1+"],
        "bilateral":  ["Dorsal sub-epi."],
        "excluded":   ["Skeletogenic LPM", "Sclerotome", "Dermomyotome",
                       "Mesothelium", "Venous-pole spl.", "Contamination",
                       "ambiguous"],
    },
    "HH25": {
        # 2026-07-02 vet: the HH25 syrinx precursor (IRX1+) is cl3, relabeled from
        # "Chondrogenic LPM" → "Syringeal" (it also carries RUNX2/WNT5A cartilage
        # commitment); the old "Syringeal" cl5 (no IRX1) is Ventral spl. So
        # "Chondrogenic LPM" is no longer emitted at HH25 — dropped from inference.
        "inference":  ["Ventral spl.", "Syringeal"],
        "bilateral":  ["Dorsal sub-epi."],
        "excluded":   ["Skeletogenic LPM", "Sclerotome", "Dermomyotome",
                       "Mesothelium", "Venous-pole spl.", "Contamination",
                       "ambiguous"],
    },
}
# NB "Syringeal" (HH21 + HH25 inference) names an INFERRED identity — see
# POPULATION_DISPLAY_NAMES below. Pool membership is unaffected; only how the
# label is printed for a reader.


# ── Publication display names (hedged labels) ─────────────────────────────────
# The population strings in POOLS / d2g_subcluster_annotations.tsv are CODE keys:
# pool membership, TSV values, filters, joins. A few of them assert more than the
# data proves, so manuscript-facing panels and legends render a hedged label. This
# map is DISPLAY-ONLY and purely additive — the SSOT string never changes, and no
# code may filter, join, or compare on a value from here.
#
# "Syringeal": the identity is marker- and dissection-inferred, NOT spatially validated.
# The chondrogenic program (WNT5A/RUNX2/DPP10, PTCHD4 15x) is HH25-ONLY — the HH21
# cluster carrying the same label is IRX1-defined (UEG rank 1) with no cartilage gene,
# so the hedged name is earned at HH25 and carried back to HH21 by lineage continuity
# (docs/lpm_mesenchyme_subcluster_biology.md §1.6). Not spatially validated because the
# manuscript's only molecular in-situ at this locus is PITX2 HCR.
# "Syringeal chondrocyte precursor" keeps the marker-validated part
# (chondrocyte precursor) and flags the anatomical part as inferred, while staying
# distinct from the chick's generic "Chondrogenic LPM" (which lacks this precursor
# entirely — docs/chick_control/laterality_fold_methods.md §5.5, not distributed).
# Rationale: docs/lpm_mesenchyme_subcluster_biology.md §"Naming honesty".
POPULATION_DISPLAY_NAMES: dict[str, tuple[str, str]] = {
    # SSOT label: (long form — legends & prose, short form — panel/axis text)
    "Syringeal": ("Syringeal chondrocyte precursor", "Syringeal chondro."),
}

# A typo'd key would silently never fire.  Pin it to the POOLS vocabulary at import
# (POOLS is a literal — no file I/O, so this is safe in every context).  NB this
# catches a key matching no POOLS label; it does NOT catch a key that is a valid
# POOLS label but is never emitted into the annotations TSV (e.g. Sclerotome /
# Dermomyotome, filtered upstream by pipeline/_d2g_lpm_selection.py) — such a hedge
# would be a harmless no-op.  `assert_display_keys_emitted()` below checks that
# lazily for callers who want it.
_POOL_LABELS = {p for st in POOLS.values() for kind in st.values() for p in kind}
assert set(POPULATION_DISPLAY_NAMES) <= _POOL_LABELS, (
    "POPULATION_DISPLAY_NAMES keys absent from POOLS: "
    f"{sorted(set(POPULATION_DISPLAY_NAMES) - _POOL_LABELS)}"
)
# Validate the VALUE shape too. A bare str passes a key-only check and then silently
# renders as a single CHARACTER (str indexing: "Long name"[1] -> 'o'), which is the one
# realistic authoring mistake here — there is only ever one line to copy.
_BAD_VALUES = {
    k: v for k, v in POPULATION_DISPLAY_NAMES.items()
    if not (isinstance(v, tuple) and len(v) == 2
            and all(isinstance(x, str) and x.strip() for x in v))
}
assert not _BAD_VALUES, (
    "POPULATION_DISPLAY_NAMES values must be (long, short) tuples of non-empty str; "
    f"malformed: {_BAD_VALUES}"
)


def assert_display_keys_emitted() -> None:
    """Raise if a hedge key never appears in the annotations TSV (dead entry).

    Separate from the import-time assert because it touches the SSOT file, which
    is not guaranteed present in every import context (e.g. on scratch before sync).
    """
    emitted = set(_load()["population"])
    dead = sorted(set(POPULATION_DISPLAY_NAMES) - emitted)
    if dead:
        raise RuntimeError(
            f"POPULATION_DISPLAY_NAMES entries never emitted by d2g: {dead} — the hedge "
            f"is a no-op. Remove it, or fix the key to a label the scorer produces."
        )


def display_population(population: str, short: bool = False,
                       unhedged: Optional[str] = None) -> str:
    """Reader-facing label for an SSOT population name.

    Call this ONLY where a label is rendered for a human (dotplot rows, axis text,
    facet strips, figure legends). NEVER on a value used to filter, join, group, or
    compare AGAINST SSOT DATA — those must stay on the SSOT string.

    One deliberate exception: when the thing being matched is ITSELF a rendered label,
    the key must be derived here rather than hardcoded. `figures/data/d2g_dotplot_spatial*
    .tsv`'s `group` column is a rendered row label, and
    `fig02/generate_spatial_svg.py` reverse-looks-up its HH21 callouts against it. NB the
    writer passes `unhedged=<its own axis shortening>`, so a reverse lookup matches only
    where that shortening equals the SSOT name (or the population is hedged) — see the
    `group` caveat in `figures/data/README.md`.

    `unhedged` is the caller's own label for populations that carry no hedge (a
    figure often shortens names to fit its axis, e.g. "Skeletogenic LPM" →
    "Skeletogenic"); pass it so a hedged name overrides the figure-local
    shortening while every other row keeps it. It is NOT an error fallback —
    most populations legitimately have no hedge, the map is deliberately sparse,
    and the import-time assert above catches the one real failure mode: a key
    that matches no population at all.
    """
    hedged = POPULATION_DISPLAY_NAMES.get(population)
    if hedged is None:
        return population if unhedged is None else unhedged
    return hedged[1] if short else hedged[0]


# Per-stage allow-list of pool labels that may be absent without raising.
# These are biological no-shows: e.g. Mesothelium has 0 clusters at HH13
# in our data; BARX1+ peaks at HH17-21 and can be absent at HH25; etc.
# Any other absence (e.g. typo) raises so SoT drift is caught at panel-build.
EXCLUDED_OPTIONAL_PER_STAGE: dict[str, set[str]] = {
    "HH13":    {"Mesothelium"},          # no WT1+ cluster at HH13 in current run
    "HH17/16": set(),
    "HH21":    set(),
    # 2026-07-02 vet: Mesothelium may be absent at HH25 (cl8→Ventral spl.,
    # cl10→Venous-pole spl. after the respiratory-gate fix); Venous-pole present (cl10).
    "HH25":    {"Venous-pole spl.", "Mesothelium"},
}

# Pool labels that may be absent at any stage without raising.  Two reasons
# a label lands here: (a) paraxial / structural labels that the d2g
# LPM-restricted pipeline never emits (LPM_ANNOTATIONS filter in
# pipeline/_d2g_lpm_selection.py removes them upstream); (b) the "ambiguous"
# residual category from the compendium scorer, which appears only when a
# cluster has no clear lineage winner.  Listed module-level so a future
# rename of "Sclerotome" → "Paraxial sclerotome" is a one-grep update.
_ALWAYS_OPTIONAL: frozenset[str] = frozenset({
    "Sclerotome", "Dermomyotome",  # paraxial — filtered upstream
    "ambiguous",                   # residual category
    "Contamination",               # technical junk (ribosomal/erythroid);
                                   # present only at early stages, absent at HH21/HH25
})


def get_pool_map(stage: str, kind: str = "inference") -> dict[str, str]:
    """Return {cluster_id_str: population} for one pool kind (inference/bilateral/excluded).

    Joins POOLS[stage][kind] with the current d2g annotations TSV to resolve
    cluster IDs. Raises if any labels in POOLS[stage][kind] are absent from
    the TSV (catches post-rerun label drift). For `excluded` kind, allows
    biologically-known-absent labels per `EXCLUDED_OPTIONAL_PER_STAGE`.
    """
    if stage not in POOLS:
        raise ValueError(f"stage={stage!r} not in POOLS")
    if kind not in {"inference", "bilateral", "excluded"}:
        raise ValueError(f"kind={kind!r} not in {{inference,bilateral,excluded}}")
    wanted = set(POOLS[stage][kind])
    df = _load()
    sub = df[df["stage"] == stage]
    out = {str(r["cluster"]): r["population"]
           for _, r in sub.iterrows() if r["population"] in wanted}
    # Validate: every wanted label must either be present in the TSV or
    # explicitly opted-out via EXCLUDED_OPTIONAL_PER_STAGE.  Catches typos
    # and post-rename drift in any of the three pool kinds.
    present = set(sub["population"])
    missing = wanted - present
    if kind == "excluded":
        # Filter biologically-known-absent labels for this stage, plus the
        # always-optional set (paraxial + residual ambiguous).
        optional = EXCLUDED_OPTIONAL_PER_STAGE.get(stage, set()) | _ALWAYS_OPTIONAL
        missing = missing - optional
    if missing:
        raise RuntimeError(
            f"POOLS[{stage!r}][{kind!r}] expects {sorted(missing)} but "
            f"d2g_subcluster_annotations.tsv has none — regenerate via "
            f"pipeline/d2g_validate_subclusters.py, update POOLS, or add to "
            f"EXCLUDED_OPTIONAL_PER_STAGE if biologically-absent."
        )
    return out


def get_pool_clusters(stage: str, kind: str = "inference") -> list[int]:
    """Convenience: sorted int cluster IDs for the pool."""
    return sorted(int(c) for c in get_pool_map(stage, kind))
