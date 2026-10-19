"""Chick laterality-fold POOLS — single source of truth (duck `subcluster_map.py` analog).

The chick Panel F (orientation-free PITX2 half-vs-half) and its directional control are a
faithful port of the duck `figures/fig02/generate_paired_laterality.R`, which folds the
stage-restricted subcluster populations into three SPATIAL DOMAINS via the duck
`figures/_shared/subcluster_map.py::POOLS`:

    inference  -> "Ventral splanchnic pool" (HH13 -> its own "LPM mes." facet)
    bilateral  -> "Dorsal sub-epi."         (the within-figure bilateral control)
    excluded   -> dropped                    (Skeletogenic LPM [was Somatic] / Mesothelium /
                                              Venous-pole spl. / Sclerotome / Dermomyotome / ambiguous)

This module is the CHICK analog of that POOLS dict. It exists because the original `cAN`
folding silently diverged from the duck — it kept Somatic + Mesothelium at HH13, so the chick
HH13 facet was ~56% body-wall somatic mesenchyme, exactly what the duck excludes. `CHICK_POOLS`
removes that divergence (the genuine correction).

Full rationale + duck-vs-chick comparison (not distributed): docs/chick_control/laterality_fold_methods.md;
summarised in the manuscript Methods.

The fold is a LITERAL duck port (C-strict): HH13 is its own inference facet ("LPM mes. (HH13)");
HH15-16 and HH25 split inference (Ventral) vs bilateral (Dorsal) exactly as the duck splits its
HH17/16 + HH21/HH25. We deliberately do NOT pool HH13 + HH15-16 into one facet — that ("C-broad")
would recover the sparse Group A HH15-16 embryos but is pure count-maximization: chick HH15-16 is
the same developmental situation as the duck HH17/16, which the duck splits, so there is no
chick-specific biology to justify pooling. The dropped sparse embryos are documented in
`cAN_domain_counts.tsv`; the conclusion is unchanged (see methods doc §5.1).

One annotation note: **the HH13 early cl0 is Early LPM, not Dorsal-biased.** The signed Dorsal_subepi
panel can win on this early cluster only because HAND2/GATA4 are near-zero; its UEGs (SMOC2 / EYA4 /
WNT5A) are early lateral-plate / intermediate mesoderm and there is no dorsal sub-epithelium pre-TES,
so biologically it is the chick's Early LPM (inference). The current SSoT already labels it "Early
LPM"; `CHICK_POOLS` additionally lists the prior "Dorsal-biased LPM" compendium alias under HH13
inference so the fold stays correct regardless of label timing.

HH25 inference = {Ventral spl., BARX1+, Chondrogenic LPM} (Syringeal kept for forward-compat; the
chick run currently has no Syringeal cluster). This mirrors the duck HH21/HH25 inference, per the
2026-06-10 confirmation to include BARX1+ and Chondrogenic LPM.
"""
from __future__ import annotations

import os
from typing import Optional

import pandas as pd

# ── facet names (the three Panel F columns) — match the duck labels ────────────
FACET_HH13 = "LPM mes. (HH13)"
FACET_VENTRAL = "Ventral splanchnic pool"
FACET_DORSAL = "Dorsal sub-epi."
FACET_ORDER = [FACET_HH13, FACET_VENTRAL, FACET_DORSAL]

STAGES = ("HH13", "HH15-16", "HH25")

# ── pool membership per stage (the chick analog of duck subcluster_map.POOLS) ──
# inference = laterality-showing LPM-lineage populations
# bilateral = dorsal (expected ~bilateral; same metric, expect weak/null)
# excluded  = Skeletogenic LPM (FOXF1-low pharyngeal-cranial skeletal / pre-chondrogenic
#             progenitor; was "Somatic" — duck 2026-06-18 rename, chick aligned here),
#             Mesothelium (coelomic surface lining), Venous-pole spl. (off-axis cardiac
#             venous-pole / proepicardium), Sclerotome / Dermomyotome (paraxial — different
#             lineage), ambiguous.  All EXCLUDED AT EVERY STAGE, exactly as the duck does.
#             "Somatic" is kept as a legacy alias so a not-yet-regenerated chick export
#             (pre-rename label) still folds correctly.
_EXCLUDED_ALL = ["Skeletogenic LPM", "Somatic", "Mesothelium", "Venous-pole spl.",
                 "Sclerotome", "Dermomyotome", "ambiguous"]

CHICK_POOLS: dict[str, dict[str, list[str]]] = {
    "HH13": {
        # HH13 early cl0 = Early LPM (see module docstring; "Dorsal-biased LPM" is the prior
        # compendium alias for the same population, kept for label-timing robustness — matched on
        # the population label, NOT a non-seed-stable cluster ID). Cardiopharyngeal LPM is a
        # forward-compat member that does not realize until HH15-16. No resolvable dorsal pre-TES,
        # so bilateral is empty — matches the duck HH13 (bilateral = []).
        "inference": ["Early LPM", "Cardiopharyngeal LPM", "Dorsal-biased LPM"],
        "bilateral": [],
        "excluded": list(_EXCLUDED_ALL),
    },
    "HH15-16": {
        # Split Ventral/Dorsal exactly as the duck splits HH17/16 (its pre-TES middle stage).
        "inference": ["Early LPM", "BARX1+", "Cardiopharyngeal LPM"],
        "bilateral": ["Dorsal-biased LPM"],   # pre-commitment NKX6-1+/dorsal-biased control
        "excluded": list(_EXCLUDED_ALL),
    },
    "HH25": {
        "inference": ["Ventral spl.", "BARX1+", "Chondrogenic LPM", "Syringeal"],
        "bilateral": ["Dorsal sub-epi."],
        "excluded": list(_EXCLUDED_ALL),
    },
}


def facet_of(stage: str, population: str) -> Optional[str]:
    """Map a (firm stage, compendium population) to its Panel F facet, or None if dropped.

    Literal duck `fold_for_laterality` port (C-strict): HH13 inference -> "LPM mes. (HH13)";
    HH15-16/HH25 inference -> "Ventral splanchnic pool"; any-stage bilateral -> "Dorsal sub-epi.".
    Anything in the excluded pool — or any population not named in CHICK_POOLS[stage] — is
    dropped (returns None), so a stray/new label fails closed rather than leaking in.
    """
    pools = CHICK_POOLS.get(stage)
    if pools is None or population in pools["excluded"]:
        return None
    if population in pools["inference"]:
        return FACET_HH13 if stage == "HH13" else FACET_VENTRAL
    if population in pools["bilateral"]:
        return FACET_DORSAL
    return None


def is_inference(stage: str, population: str) -> bool:
    """True if the population is in the laterality (inference) pool for that stage.

    Used by the directional control (designed-LEFT concordance) — that test runs on the
    laterality-showing pool only, NOT the bilateral dorsal control.
    """
    pools = CHICK_POOLS.get(stage)
    return bool(pools) and population in pools["inference"]


def validate_against_annotations(tsv_path: str, log=print) -> dict[str, list[str]]:
    """Cross-check CHICK_POOLS vs the cAL subcluster annotation SSoT (transparency, not fatal).

    Reports (a) populations present in the data but not classified by any pool (would silently
    drop) and (b) POOLS labels absent from the data (forward-compat placeholders like Syringeal).
    Returns the per-stage list of unclassified populations so callers can assert if they want.
    """
    if not os.path.exists(tsv_path):
        log(f"[pools] annotations TSV not found ({tsv_path}); skipping validation")
        return {}
    ann = pd.read_csv(tsv_path, sep="\t")
    unclassified: dict[str, list[str]] = {}
    for stage in sorted(set(ann["stage"]) & set(CHICK_POOLS)):
        present = set(ann.loc[ann["stage"] == stage, "population"])
        classified = (set(CHICK_POOLS[stage]["inference"])
                      | set(CHICK_POOLS[stage]["bilateral"])
                      | set(CHICK_POOLS[stage]["excluded"]))
        missing_from_pools = sorted(present - classified)
        absent_placeholders = sorted(classified - present)
        if missing_from_pools:
            log(f"[pools] {stage}: populations present but NOT in CHICK_POOLS "
                f"(WOULD DROP): {missing_from_pools}")
            unclassified[stage] = missing_from_pools
        if absent_placeholders:
            log(f"[pools] {stage}: CHICK_POOLS labels absent from data "
                f"(forward-compat ok): {absent_placeholders}")
    if not unclassified:
        log("[pools] OK — every observed population is classified by CHICK_POOLS")
    return unclassified


if __name__ == "__main__":
    # quick self-report against the local SSoT
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(here))
    tsv = os.path.join(repo, "figures", "data", "chick_control",
                       "cAL_subcluster_annotations.tsv")
    print("CHICK_POOLS facet assignment per stage (C-strict, duck-faithful):")
    ann = pd.read_csv(tsv, sep="\t")
    for stage in STAGES:
        print(f"\n  {stage}:")
        for pop in sorted(set(ann.loc[ann["stage"] == stage, "population"])):
            print(f"    {pop:24s} -> {facet_of(stage, pop)}")
    print()
    validate_against_annotations(tsv)
