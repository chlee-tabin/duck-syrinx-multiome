"""Compendium-based cluster annotation for duck syrinx d2g subclusters.

Consumes the three marker-discovery TSVs produced by
pipeline/d2g_stage_restricted_subclustering.py:

    figures/data/d2g_uegs_<label>.tsv              top-50 UEGs per cluster (Stein 2026 U-method)
    figures/data/d2g_wilcoxon_<label>.tsv          top-50 Wilcoxon DE per cluster
    figures/data/d2g_lineage_scores_<label>.tsv    per-cluster score against 11 canonical panels

Replaces the threshold-rule classifier in pipeline/d2g_validate_subclusters.py
with a scoring/evidence-based marker annotation.

Two-level output hierarchy:
    cell_type_broad    LPM-lineage / Paraxial / Mesothelial / Skeletogenic / Apoptotic / Unknown
    cell_type          Fine-grained label (Ventral spl., Syringeal, Cardiopharyngeal LPM, etc.)

Contamination pre-check (analogous to apoptotic_check; runs before lineage
assignment):  reclassifies technical-junk clusters to cell_type="Contamination"
(broad "Unknown") — ribosomal/low-complexity clusters and a small erythroid/blood
cluster that the early-stage subclustering otherwise folds into "Skeletogenic LPM".
See contamination_check() for the marker/QC basis + threshold justification.

Stage-gated relabelings:
    Cardio_Syringeal    → "Cardiopharyngeal LPM" at HH13/HH17-16 (pre-TES)
                        → "Syringeal" at HH21/HH25 (post-TES)

Outputs:
    /tmp/d2g_audit/compendium_cluster_labels.tsv    per (stage, variant, resolution, cluster) label + evidence
    /tmp/d2g_audit/compendium_cell_labels.tsv       per cell (joined via d2g_lpm_subclustering TSVs)
    /tmp/d2g_audit/compendium_vs_rules_flips.tsv    per cluster: threshold-rule label vs compendium label, flip reason
"""
from __future__ import annotations
from pathlib import Path
import re
import pandas as pd
import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "figures" / "data"
# Output directory stays under /tmp during the compendium-migration phase so
# downstream consumers explicitly opt in.  Switch to REPO / "references"
# when the compendium supersedes the threshold-rule classifier.
OUT  = Path("/tmp/d2g_audit")
OUT.mkdir(parents=True, exist_ok=True)

TECH_REP_MERGE = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1",
                  "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}

# ── Stage-label coverage ─────────────────────────────────────────────────
# 4 stages × 2 variants (default + ccreg) = 8 (stage_label, d2g_file) combos.
STAGE_LABELS = [
    ("HH13",        "HH13",         "HH13"),
    ("HH13",        "HH13_ccreg",   "HH13"),
    ("HH17/16",     "HH1716_HH16",        "HH17/16"),
    ("HH17/16",     "HH1716_HH16_ccreg",  "HH17/16"),
    ("HH21",        "HH21",         "HH21"),
    ("HH21",        "HH21_ccreg",   "HH21"),
    ("HH25",        "HH25",         "HH25"),
    ("HH25",        "HH25_ccreg",   "HH25"),
]

# ── Lineage → biological label mapping ───────────────────────────────────
# Stage-gated Cardio_Syringeal handled separately.  Apoptotic mes. handled
# by a pre-check on QC metrics (not a lineage panel).

LINEAGE_TO_LABEL = {
    "Dorsal_subepi":   ("Dorsal sub-epi.",      "LPM-lineage"),
    "Ventral_spl":     ("Ventral spl.",         "LPM-lineage"),
    "BARX1_pos":       ("BARX1+",               "LPM-lineage"),
    "Splanchnic_vis":  ("Venous-pole spl.",     "LPM-offaxis"),
    # The lineage-panel KEY "Somatic" is plumbing (matches the
    # d2g_lineage_scores_*.tsv column + the d2g_stage_restricted_subclustering
    # CELLTYPE_MARKERS panel); only the DISPLAY label changed.  This population
    # is FOXF1-low, HOX-negative, TBX1+/EYA1+/EYA2+/EBF2/3+ pharyngeal-cranial
    # SKELETAL mesenchyme — a pre-chondrogenic progenitor, NOT body-wall
    # somatopleure.  Because d2g gates to LPM mesenchyme at the parent level it
    # is LPM by construction; it is L/R-symmetric (correctly EXCLUDED from the
    # PITX2 laterality pool).  Renamed "Somatic" → "Skeletogenic LPM" (2026-06-18
    # vet) to parallel the separately-captured lateralized syrinx lineage that
    # stays IN the pool — that population was "Chondrogenic LPM" at the time and
    # is now "Syringeal" at HH25 (2026-07-02 vet, see irx1_syringeal_gate).
    "Somatic":         ("Skeletogenic LPM",     "Skeletogenic"),
    "Epicardial":      ("Mesothelium",          "Mesothelial"),
    "Differentiating": ("Chondrogenic LPM",     "LPM-lineage"),
    "Early_LPM":       ("Early LPM",            "LPM-lineage"),
    "Sclerotome":      ("Sclerotome",           "Paraxial"),
    "Dermomyotome":    ("Dermomyotome",         "Paraxial"),
}

# Stage-gated Cardio_Syringeal renders as different labels at different
# stages (Cardiopharyngeal LPM at HH13/HH17-16, Syringeal at HH21/HH25);
# both are LPM-lineage.  Override paths (HH13 / WT1 / BARX1) can also flip
# cells to "ambiguous".  This flat dict is the single lookup for
# label -> broad — replaces the duplicated O(n) loops that used to scan
# LINEAGE_TO_LABEL.values() after each override.
# ── Skeletogenic panel gate (annotation vet 2026-06-24 → comprehensive vet 2026-07-02) ──
# The Somatic panel is up=[PRRX1, TWIST1, EBF2, EYA1, EYA2] — PRRX1/TWIST1 are
# PAN-MESENCHYME, so the panel over-fires on generic ECM fibroblasts that carry no
# skeletal identity.  Adversarial vetting
# found HH13 cl1 AND HH17/16 cl4 are pure TMEM132D/FLRT2/DLC1 fibroblasts with ZERO
# EYA/EBF/TBX1/COL25A1 in their top-50 DE, yet won "Skeletogenic LPM" on PRRX1/TWIST1
# alone; the genuinely skeletal clusters (HH17/16 cl9, HH21 cl6, HH25 cl6) place
# EYA1/EBF2/EBF3/TBX1/COL25A1 at top DE ranks.  branchiomeric_gate() therefore
# REQUIRES >=1 true pharyngeal-skeletal marker in the top-N DE/UEG before a cluster
# may keep "Skeletogenic LPM"; otherwise it is relabeled to the generic-fibroblast
# label below.  Stage-agnostic — SUPERSEDES the 2026-06-24 HH13-only Rule C.  The
# relabeled clusters stay EXCLUDED from the inference pool (FOXF1-low), so pool
# membership is unchanged.  ASCII "+" matches the "BARX1+" convention (no Unicode).
PRRX1_FOXF1LOW_LABEL = "PRRX1+/FOXF1-low mes."
BRANCHIOMERIC_MARKERS = frozenset({
    "EYA1", "EYA2", "EBF1", "EBF2", "EBF3", "TBX1", "COL25A1", "SIX1", "SIX4",
})
BRANCHIOMERIC_GATE_TOP_N = 50   # scan top-N Wilcoxon/UEG genes for a skeletal marker

# Syringeal identity gate (IRX1).  The syrinx precursor is IRX1+.  The
# Cardio_Syringeal panel [HAND2, ISL1, IRX1] can win the post-TES "Syringeal"
# display WITHOUT IRX1 (HH25 cl5 — a 3-way near-tie that is actually ventral-
# splanchnic), while the true IRX1+ syrinx precursor sits in a Differentiating-
# winning "Chondrogenic LPM" cluster (HH25 cl3).  irx1_syringeal_gate() moves the
# "Syringeal" label onto the IRX1-bearing cluster.  See the 2026-07-02 vet.
#
# TWO THINGS THE RESULTING LABEL DOES NOT MEAN:
#  1. "Syringeal" is an INFERRED anatomical identity — a marker program plus
#     syrinx-fated dissection provenance.  NB the cartilage program
#     (WNT5A/RUNX2/DPP10; PTCHD4 15x, the most specific; IRX1 54%) is HH25-ONLY;
#     the HH21 cluster carrying the same label is IRX1-defined (UEG rank 1) with no
#     cartilage gene.  (On the CURRENT run neither gate branch fires at HH21: HH21 cl3
#     wins Cardio_Syringeal (0.68) and carries IRX1 at UEG rank 1, so its label comes
#     from label_from_lineage's post-TES ternary.  That is a property of this run, not
#     of the gate: structurally only the Chondrogenic-LPM->Syringeal branch is HH25-only
#     (STAGE_GATED_PANELS["Differentiating"] == {"HH25"}); the Syringeal-without-IRX1
#     branch is live at HH21 too, since irx1_syringeal_gate() takes no stage argument.)
#     The "Syringeal" identity is NOT spatially validated (the only
#     in-situ at this locus is PITX2 HCR), so publication figures render it hedged
#     as "Syringeal chondrocyte precursor" in legends/prose and "Syringeal chondro." in
#     panel text, via
#     figures/_shared/subcluster_map.py::display_population.  This SSOT string stays
#     "Syringeal" — docs/lpm_mesenchyme_subcluster_biology.md §1.6.
#  2. The gate is DUCK-ONLY BY DESIGN, not an un-ported oversight.  Applied to the
#     chick it is a verified no-op, and the chick lacks the IRX1+/PTCHD4+ precursor
#     outright — docs/chick_control/laterality_fold_methods.md §5.5 (not distributed).
SYRINGEAL_MARKER = "IRX1"
SYRINGEAL_GATE_TOP_N = 50

LABEL_TO_BROAD: dict[str, str] = {
    label: broad for label, broad in LINEAGE_TO_LABEL.values()
}
LABEL_TO_BROAD["Cardiopharyngeal LPM"] = "LPM-lineage"
LABEL_TO_BROAD["Syringeal"]            = "LPM-lineage"
LABEL_TO_BROAD["ambiguous"]            = "Unknown"
LABEL_TO_BROAD["Apoptotic mes."]       = "Apoptotic"
# Contamination = technical junk reclassified out of lineage assignment by the
# contamination_check() pre-check (ribosomal/low-complexity or erythroid/blood
# clusters).  Broad "Unknown" so it groups with the residual/excluded bucket.
LABEL_TO_BROAD["Contamination"]        = "Unknown"
# Generic PRRX1+/FOXF1-low fibroblast — the branchiomeric_gate() relabel of
# Skeletogenic-panel clusters that lack a real skeletal marker (see above).
LABEL_TO_BROAD[PRRX1_FOXF1LOW_LABEL] = "Skeletogenic"
# Cardio_Syringeal is stage-gated — handled in label_from_lineage().
# Splanchnic_vis → "Venous-pole spl." per docs/lpm_mesenchyme_subcluster_biology.md
# (2026-05-04 decision): GATA4+/ZFPM2+/TBX20+/TBX4-low/WT1-low signature is
# non-myocardial splanchnic envelope of the venous pole / proepicardium
# (van Wijk 2009 PMC2861358), NOT respiratory-tract splanchnic. Broad
# category "LPM-offaxis" flags this for figure scripts as off-axis
# dissection-edge contamination (excluded from inference pool, exposed only
# in figS02).
#
# Epicardial → "Mesothelium": the WT1+TCF21+ cells in foregut/respiratory
# region are coelomic surface mesothelium, not heart-specific epicardium.

LINEAGE_KEYS = list(LINEAGE_TO_LABEL.keys()) + ["Cardio_Syringeal"]

# Stage-gate: lineage panels that should NOT score outside their target
# stage(s).  At out-of-stage, the compendium scorer overrides the lineage
# score to -inf so the panel never wins.
#
# - Ventral_spl: requires tracheo-esophageal (TES) septum closure →
#   post-TES stages HH21+ only.  The TBX5+/ZFPM2+/ISL1+ signature at
#   HH17/16 is pre-specification; fold it into Early_LPM instead.
# - Early_LPM: pan-LPM (TBX5+FOXF1+ISL1) pre-specification; valid at
#   HH13 and HH17/16.  At HH21+ the Ventral spl. signature dominates on
#   the same cells, so this panel never wins there anyway.
# - Chondrogenic LPM ("Differentiating" panel key, kept for compendium
#   column stability): RUNX2/WNT5A/DPP10 — chondro/osteoblast commitment
#   consistent with HH25 tracheal/bulla cartilage primordia.
STAGE_GATED_PANELS = {
    "Ventral_spl":     {"HH21", "HH25"},
    "Early_LPM":       {"HH13", "HH17/16"},
    "Differentiating": {"HH25"},
}


def apply_stage_gate(scores: dict, stage: str) -> dict:
    """Return a copy with out-of-stage gated panels set to -inf."""
    out = dict(scores)
    for panel, allowed in STAGE_GATED_PANELS.items():
        if panel in out and stage not in allowed:
            out[panel] = float("-inf")
    return out


def label_from_lineage(lin: str, stage: str) -> tuple[str, str]:
    if lin == "Cardio_Syringeal":
        return (("Cardiopharyngeal LPM", "LPM-lineage")
                if stage in ("HH13", "HH17/16")
                else ("Syringeal", "LPM-lineage"))
    return LINEAGE_TO_LABEL.get(lin, ("Unknown", "Unknown"))


# ── Annotation thresholds (tunable) ──────────────────────────────────────
# Minimum lineage-score separation between top and runner-up — below this,
# assignment is flagged as "weak" and a reader should read the UEG evidence.
MIN_SCORE_MARGIN = 0.05
# Minimum absolute score for the winning lineage (if all lineages score
# near zero, the cluster matches none → label as "ambiguous")
MIN_SCORE_ABSOLUTE = 0.0      # weak filter; lineage scoring is normalized

# Apoptotic mes. QC pre-check (mirrors rule 0 from the threshold classifier)
APOPTOTIC_RIBO_THRESHOLD = 10.0
APOPTOTIC_CYCLING_THRESHOLD = 85.0
APOPTOTIC_LINEAGE_CEILING = 20.0


def apoptotic_check(pct_ribo_mean: float, cycling_pct: float,
                    max_lineage_pct: float) -> bool:
    return (pct_ribo_mean > APOPTOTIC_RIBO_THRESHOLD
            and cycling_pct > APOPTOTIC_CYCLING_THRESHOLD
            and max_lineage_pct < APOPTOTIC_LINEAGE_CEILING)


# ── Contamination pre-check (technical junk) ─────────────────────────────
# Some early-stage (HH13, HH17/16) clusters are technical artifacts that the
# LPM-restricted subclustering otherwise folds into the skeletogenic-LPM
# population: (a) ribosomal/low-complexity clusters (HH13 has ~3 small ones;
# HH17/16 has ~1 of ~988 cells) and (b) a small erythroid/blood cluster
# (HH17/16, ~65 cells).  This pre-check reclassifies them to
# cell_type="Contamination" (broad "Unknown") BEFORE lineage assignment, the
# same way apoptotic_check short-circuits the lineage scorer.  Leiden IDs are
# NEVER hardcoded (they shuffle every d2g rerun); detection is purely
# marker/QC-based.
#
# RIBOSOMAL arm — the WILCOXON DE list flags a *candidate*, but it must be
# AND-gated on the ABSENCE of a coherent lineage before firing.  Ribosomal/
# housekeeping genes are ubiquitous, so they score LOW as U-method UEGs (per
# the compendium design) but DOMINATE the positive Wilcoxon DE of any
# LOW-DEPTH cluster — including REAL low-depth populations, not just junk.
# Calibration against the committed ccreg TSVs (2026-06-18) shows the
# Wilcoxon-ribosomal criterion alone is TOO LOOSE — it fires identically (10–13
# of top-15) on genuine junk AND on two real low-depth LPM clusters:
#   junk (no lineage)     → ribo-Wilcoxon 10–13/15, max specific-lineage ≤19%:
#                           HH13 cl5 (NKX6-1=8%), cl6 (FOXF1=11%), cl7 (ISL1=19%),
#                           HH17/16 cl7 (FOXF1=14%)
#   REAL low-depth LPM    → ribo-Wilcoxon 10–12/15 BUT a clear lineage signature:
#                           HH13 cl8  Early LPM   (FOXF1=47%, TBX5=42%, ISL1=26%)
#                           HH25 cl9  Ventral spl.(ZFPM2=42%, TBX5=58%, PITX2=17%)
#   real skeletogenic etc.→ 0 ribo-Wilcoxon hits (HH13 cl1, HH17/16 cl2/4/9, …)
# The discriminator is therefore the cluster's MAX fraction-positive over
# SPECIFIC lineage markers (SPECIFIC_LINEAGE_MARKERS — excludes pan-mesenchyme
# PRRX1/TWIST1/NR2F2, which are high everywhere and carry no lineage identity):
# junk sits at ≤19%, the two real low-depth clusters at ≥42%.  A gate at <30%
# cleanly separates them (margin 19→30→42).  Both ribosomal arms (Wilcoxon and
# pct_ribo) require this gate; an arm may fire ONLY when the cluster ALSO lacks
# a coherent lineage (max specific-lineage fraction < 30%).  pct_ribo is a
# corroborating signal but is NOT sufficient alone: a junk cluster runs pct_ribo
# 9–19% while a real population (HH25 Ventral-spl cl9=9.4%, HH13 cl8=7.6%) can
# also be mildly elevated — hence the same low-lineage gate on the pct_ribo arm.
# This catches a hypothetical future ribo cluster whose gene symbols don't match
# the regex (e.g. all-LOC ribosomal annotations) without false-positiving on the
# real low-depth populations.
#
# ERYTHROID arm — blood genes are cluster-defining for an erythroid cluster and
# appear in BOTH the Wilcoxon DE and the UEG list.  Calibration:
#   HH17/16 cl11 (65 cells) → TFRC,STEAP3,KEL in Wilcoxon top-15 (3) and
#                             TFRC,STEAP3,KEL,SPTB in UEG top-15 (4)
#   all real populations    → 0 in either list
# A threshold of >=3 distinct erythroid genes in the UNION of top-15 Wilcoxon +
# top-15 UEGs separates cl11 (3) from every real population (0).
RIBO_HK_REGEX = re.compile(
    r"^RP[SL]\d|^MRP[SL]|EEF1A1|EEF2|NPM1|NCL|PTMA|GAPDH|ACTB|HSP90AA1|TMSB4X"
)
ERYTHROID_GENES = frozenset({
    "HBA1", "HBB", "HBAD", "ALAS2", "SLC4A1", "SPTB", "KEL",
    "TFRC", "STEAP3", "GYPA", "HEMGN", "RHAG",
})
# Ribosomal arm 1: Wilcoxon-DE dominance.
CONTAM_RIBO_WILCOXON_MIN = 8       # >= this many of top-15 Wilcoxon → ribosomal
CONTAM_RIBO_TOP_N = 15
# Ribosomal arm 2: pct_ribo fingerprint (independent fallback for un-symbol'd
# ribosomal annotations).  Tuned ABOVE the highest real-population pct_ribo
# (HH13 Early-LPM cl0=5.2%, cl8=7.6%) with margin; only the junk cl5 (19%)
# would trip it on current data, and cl5 is already caught by the Wilcoxon arm.
CONTAM_RIBO_PCT_MIN = 15.0
# SPECIFIC-lineage gate shared by BOTH ribosomal arms: a ribosomal candidate
# fires ONLY if it also lacks a coherent lineage (max fraction-positive over
# SPECIFIC_LINEAGE_MARKERS below this).  Separates genuine junk (≤19%) from
# real low-depth LPM clusters (HH13 cl8 FOXF1=47%, HH25 cl9 ZFPM2=42% → both
# ≫30 → un-flagged).  The erythroid arm is intentionally NOT gated (blood is
# contamination regardless of any mesenchymal lineage signal).
CONTAM_RIBO_MAX_LINEAGE_GATE = 30.0
# SPECIFIC lineage markers (NOT pan-mesenchyme PRRX1/TWIST1/NR2F2, which are
# high everywhere and don't indicate lineage identity).  The per-cluster max
# fraction-positive over the subset of these present as <gene>_cp10k columns in
# the d2g per-cell TSV gates the ribosomal arms.  Symbols not exported (EYA1,
# SOX9, COL2A1) are absent from the current cp10k export and simply skipped;
# NKX6-1 (LOC113843492) is exported as NKX6_1_cp10k.
SPECIFIC_LINEAGE_MARKERS = (
    "NKX6_1", "FOXF1", "GATA4", "ZFPM2", "TBX5", "TBX4", "BARX1",
    "ISL1", "HAND2", "WT1", "TBX1", "EYA1", "SOX9", "COL2A1",
)
# Erythroid arm: distinct blood genes in union(top-15 Wilcoxon, top-15 UEG).
CONTAM_ERY_MIN = 3
CONTAM_ERY_TOP_N = 15


def contamination_check(top_uegs: pd.DataFrame, top_wilcoxon: pd.DataFrame,
                        pct_ribo_mean: float,
                        max_specific_lineage_pct: float) -> tuple[bool, str]:
    """Return (is_contamination, reason) for a single cluster.

    Detects technical-junk clusters (ribosomal/low-complexity or
    erythroid/blood) from their marker profile + QC, with NO hardcoded leiden
    IDs.  See module-level comment for threshold justification.

    max_specific_lineage_pct is the cluster's MAX fraction-positive (cp10k>0,
    %) over SPECIFIC_LINEAGE_MARKERS (pan-mesenchyme genes excluded).  Both
    ribosomal arms are AND-gated on this being < CONTAM_RIBO_MAX_LINEAGE_GATE:
    a ribosomal candidate fires only when the cluster ALSO lacks a coherent
    lineage, so a real low-depth LPM cluster whose top Wilcoxon DE is dominated
    by ubiquitous ribosomal genes (HH13 cl8 Early LPM, HH25 cl9 Ventral spl.)
    is NOT mis-flagged.  The erythroid arm is NOT gated (blood is contamination
    regardless of any mesenchymal signal).

    Fails LOUD (KeyError / column access) if the expected 'gene' column is
    absent — a malformed marker TSV must surface, not silently disable the
    check.
    """
    wil_genes = top_wilcoxon["gene"].head(CONTAM_RIBO_TOP_N).astype(str).tolist()
    ueg_genes = top_uegs["gene"].head(CONTAM_ERY_TOP_N).astype(str).tolist()

    # Coherent-lineage gate shared by both ribosomal arms: only fire when the
    # cluster lacks a clear lineage signature (real low-depth LPM clusters sit
    # well above this; genuine junk sits ≤19%).
    lacks_lineage = max_specific_lineage_pct < CONTAM_RIBO_MAX_LINEAGE_GATE

    # Ribosomal arm 1: Wilcoxon-DE dominance (primary), AND-gated on no lineage.
    ribo_hits = [g for g in wil_genes if RIBO_HK_REGEX.search(g)]
    if len(ribo_hits) >= CONTAM_RIBO_WILCOXON_MIN and lacks_lineage:
        return (True,
                f"CONTAM ribosomal: {len(ribo_hits)}/{len(wil_genes)} top "
                f"Wilcoxon genes are ribo/HK ({','.join(ribo_hits[:6])}…), "
                f"max specific-lineage={max_specific_lineage_pct:.0f}% "
                f"< {CONTAM_RIBO_MAX_LINEAGE_GATE:.0f}%")

    # Ribosomal arm 2: pct_ribo fingerprint (fallback), AND-gated on no lineage.
    if (pct_ribo_mean > CONTAM_RIBO_PCT_MIN and lacks_lineage):
        return (True,
                f"CONTAM ribosomal: pct_ribo={pct_ribo_mean:.1f}% > "
                f"{CONTAM_RIBO_PCT_MIN:.0f}% with max specific-lineage="
                f"{max_specific_lineage_pct:.0f}% "
                f"< {CONTAM_RIBO_MAX_LINEAGE_GATE:.0f}%")

    # Erythroid arm: blood genes in union of top Wilcoxon + UEG (NOT gated).
    ery_hits = sorted({g for g in (wil_genes + ueg_genes)
                       if g in ERYTHROID_GENES})
    if len(ery_hits) >= CONTAM_ERY_MIN:
        return (True,
                f"CONTAM erythroid: {len(ery_hits)} blood genes in top "
                f"markers ({','.join(ery_hits)})")

    return (False, "")


# ── Evidence composer ───────────────────────────────────────────────────

def hh13_biology_override(stage: str, cell_type: str, markers: dict) -> str:
    """HH13-specific post-hoc corrections.  Panel scoring struggles with
    HH13 because cells are pre-specification — most markers are mildly
    elevated and no panel fires strongly.  Apply v1-rule-style biology
    gates only for HH13, to recover the correct label.

    markers is a dict keyed by gene symbol (upper-case) with % cells
    expressing (cp10k > 0).
    """
    if stage != "HH13":
        return cell_type
    # Bracket lookups (KeyError on missing) — main() asserts MARKERS are
    # populated before calling, so missing keys here mean a code bug.
    # Rule A: Skeletogenic LPM miscall on pan-LPM — PRRX1 below the
    # skeletogenic threshold AND TBX5+/FOXF1+ suggests Early LPM, not the
    # PRRX1-high skeletogenic-LPM progenitor.
    if (cell_type == "Skeletogenic LPM" and markers["PRRX1"] < 40
            and markers["TBX5"] >= 30 and markers["FOXF1"] >= 25):
        return "Early LPM"
    # Rule B: Mesothelium (formerly Epicardial mes.) miscall when WT1 is
    # not elevated — NR2F2-high with low WT1 is not mesothelial.  Fall
    # through to the PRRX1+/FOXF1-low fibroblast label if PRRX1 is high,
    # otherwise ambiguous.  (The stage-agnostic branchiomeric_gate() now handles
    # the direct Skeletogenic path — the 2026-06-24 HH13-only Rule C is retired.)
    if cell_type == "Mesothelium" and markers["WT1"] < 25:
        if (markers["PRRX1"] >= 30 and markers["FOXF1"] < 25
                and markers["GATA4"] < 30):
            return PRRX1_FOXF1LOW_LABEL
        return "ambiguous"
    return cell_type


def _top_gene_set(top_wilcoxon: pd.DataFrame, top_uegs: pd.DataFrame,
                  top_n: int) -> set:
    """Union of the top-`top_n` gene symbols from the Wilcoxon-DE and UEG lists.
    Fails LOUD (KeyError) if the 'gene' column is absent — a malformed marker TSV
    must surface, not silently empty the marker gates that read it."""
    wil = set(top_wilcoxon["gene"].head(top_n).astype(str))
    ueg = set(top_uegs["gene"].head(top_n).astype(str))
    return wil | ueg


def branchiomeric_gate(cell_type: str, top_wilcoxon: pd.DataFrame,
                       top_uegs: pd.DataFrame) -> tuple[str, str]:
    """Skeletogenic LPM must carry >=1 true branchiomeric marker in its top DE.

    The Somatic panel (PRRX1/TWIST1 + EBF2/EYA1/EYA2) can win on pan-mesenchyme
    PRRX1/TWIST1 alone, mislabeling generic ECM fibroblasts (HH13 cl1, HH17/16
    cl4) as skeletogenic.  Require a real pharyngeal-skeletal marker
    (BRANCHIOMERIC_MARKERS) in the top-N DE/UEG; otherwise relabel to the
    PRRX1+/FOXF1-low fibroblast label.  Stage-agnostic (supersedes the HH13-only
    Rule C).  Returns (new_label, reason_suffix); ("", suffix) is never returned.
    """
    if cell_type != "Skeletogenic LPM":
        return cell_type, ""
    if BRANCHIOMERIC_MARKERS & _top_gene_set(top_wilcoxon, top_uegs,
                                             BRANCHIOMERIC_GATE_TOP_N):
        return cell_type, ""
    return (PRRX1_FOXF1LOW_LABEL,
            f" | branchiomeric gate: Skeletogenic LPM → {PRRX1_FOXF1LOW_LABEL} "
            f"(no EYA/EBF/TBX1/COL25A1/SIX in top-{BRANCHIOMERIC_GATE_TOP_N} DE)")


def irx1_syringeal_gate(cell_type: str, runner_label: str,
                        top_wilcoxon: pd.DataFrame,
                        top_uegs: pd.DataFrame) -> tuple[str, str]:
    """Put the "Syringeal" label on the IRX1+ syrinx precursor.

    A post-TES cluster the Cardio_Syringeal panel labeled "Syringeal" but that
    LACKS IRX1 (HH25 cl5 — actually ventral-splanchnic) → its runner-up lineage
    label.  A "Chondrogenic LPM" (Differentiating-winning) cluster that HAS IRX1
    (HH25 cl3 — the real syrinx precursor) → "Syringeal".  Bounded: only the
    "Syringeal"/"Chondrogenic LPM" labels can change.  Returns
    (new_label, reason_suffix).
    """
    has_irx1 = SYRINGEAL_MARKER in _top_gene_set(top_wilcoxon, top_uegs,
                                                 SYRINGEAL_GATE_TOP_N)
    if cell_type == "Syringeal" and not has_irx1:
        return (runner_label,
                f" | IRX1 gate: Syringeal (no IRX1) → {runner_label} (runner-up)")
    if cell_type == "Chondrogenic LPM" and has_irx1:
        return ("Syringeal", " | IRX1 gate: Chondrogenic LPM (IRX1+) → Syringeal")
    return cell_type, ""


# ── Stage-agnostic post-compendium overrides (2026-05-04) ───────────────
# Encodes the colleague feedback resolution documented in
# docs/lpm_mesenchyme_subcluster_biology.md §3.3.
#
# Two rules, applied AFTER label_from_lineage() and hh13_biology_override().
# PRECEDENCE (load-bearing): WT1 rule fires FIRST; BARX1 rule only fires
# if WT1 didn't trigger (early-return).  This means a hypothetical cluster
# with both WT1>35% AND BARX1>40% is classified as Mesothelium, never as
# BARX1+.  This precedence is intentional — WT1+TCF21+ identity is a more
# specific marker for mesothelium than BARX1 is for T-E septum (BARX1 is
# also expressed in stomach mesenchyme; WT1 is mesothelium-specific in
# this anatomical region).  No current cluster trips both thresholds.
# Edge cases are flagged via a margin warning (see EDGE_WARN_MARGIN below).
#
# 1. WT1 > 35% → Mesothelium.
#    The compendium Splanchnic_vis / Ventral_spl panels can win on
#    GATA4/ZFPM2 even when the cluster is dominated by mesothelial
#    WT1+TCF21+ cells (HH25 cl8 had WT1=49% labeled "Ventral spl.").
#    A direct WT1 fraction check correctly attributes these clusters.
#    Threshold tuned against per-cluster fractions.  The
#    `cell_type != "Mesothelium"` guard means the override only acts on
#    clusters NOT already labeled by the Epicardial lineage panel:
#      fires-and-flips: HH25 cl8 (WT1=49%, was Ventral spl.) → Mesothelium
#                       HH25 cl10 (WT1=41%, was Venous-pole spl.) → Mesothelium
#      fires-but-no-op (corroborates Epicardial-panel call):
#                       HH17/16 cl6 (48), cl10 (76); HH21 cl9 (57)
#      keeps:           HH25 cl5 Syringeal (WT1=28% → below threshold)
#
# 2. BARX1 > 40% AND label == "Venous-pole spl." → BARX1+.
#    Reclaims HH21 cl2 (BARX1=58%) which the compendium scored
#    Splanchnic_vis=1.06 vs BARX1_pos=0.44 (margin 0.62 — visceral
#    side wins because cl2 is also ZFPM2=89%/GATA4=25%).  Anchored on
#    the "Venous-pole spl." label so it doesn't accidentally flip
#    Ventral spl. or Syringeal clusters that happen to have BARX1.
WT1_MESOTHELIUM_THRESHOLD = 35.0
BARX1_RECLAIM_THRESHOLD = 40.0
# Respiratory-splanchnic guard on the WT1 override (2026-07-02 vet).  Coelomic
# mesothelium is respiratory-NEGATIVE (TBX5/TBX4/ISL1 low); the WT1>35% rule
# must NOT fire on a respiratory-splanchnic cluster whose WT1 is decoupled from
# TCF21 (HH25 cl8: WT1=49% but TCF21|WT1+=10%, TBX5=65%/TBX4=54%/ISL1=47% — a
# bona-fide Ventral spl. cluster the blunt override wrongly excluded from the
# HH25 inference pool).  Genuine mesothelium sits at respiratory ≤13%, so a 30%
# ceiling cleanly separates them.
WT1_MESOTHELIUM_RESPIRATORY_MAX = 30.0
# Emit a warning (printed during scorer run) when any cluster lies within
# this margin (percentage points) of either threshold — flags fragility
# before the next rerun.  Sized to catch:
#   - HH25 cl10 at WT1=41% (6pt above 35% → caught by 8pt margin)
#   - HH25 cl5 Syringeal at WT1=28% (7pt below 35% → caught by 8pt margin)
# Both motivating cases are now in the warning band.
EDGE_WARN_MARGIN = 8.0


def wt1_barx1_override(cell_type: str, markers: dict,
                       cluster_id: int | None = None,
                       stage: str | None = None) -> tuple[str, str]:
    """Apply WT1/BARX1 overrides; return (new_label, reason_suffix).

    Precedence: WT1>35% checked first; BARX1>40% checked only if WT1 didn't
    fire.  See module-level comment for biological justification.

    Thresholds (module-level):
      WT1_MESOTHELIUM_THRESHOLD = 35.0    # → Mesothelium
      BARX1_RECLAIM_THRESHOLD   = 40.0    # → BARX1+ (Venous-pole spl. only)

    Edge-margin warning emitted if any cluster lies within EDGE_WARN_MARGIN
    of either threshold — see docs/lpm_mesenchyme_subcluster_biology.md §3.3
    "Caveat" for the fragility discussion.

    reason_suffix is "" if no override fires, else a human-readable
    string describing which override fired and on what evidence.
    """
    # Explicit lookup — KeyError raises if WT1/BARX1 are missing, rather than
    # silently disabling the override with a default-0 (which would land below
    # both thresholds and produce wrong cluster labels).
    wt1 = markers["WT1"]
    barx1 = markers["BARX1"]

    # Edge-margin diagnostic — emit BEFORE the override decision so the
    # warning shows even if no rule fires (fragile-near-threshold case).
    if abs(wt1 - WT1_MESOTHELIUM_THRESHOLD) < EDGE_WARN_MARGIN:
        side = "above" if wt1 > WT1_MESOTHELIUM_THRESHOLD else "below"
        print(f"  WARN edge-of-WT1-threshold: stage={stage} cluster={cluster_id} "
              f"WT1={wt1:.1f}% ({side} {WT1_MESOTHELIUM_THRESHOLD:.0f}% "
              f"by {abs(wt1 - WT1_MESOTHELIUM_THRESHOLD):.1f}pt) — "
              f"may flip on next d2g rerun")
    if cell_type == "Venous-pole spl." and abs(barx1 - BARX1_RECLAIM_THRESHOLD) < EDGE_WARN_MARGIN:
        side = "above" if barx1 > BARX1_RECLAIM_THRESHOLD else "below"
        print(f"  WARN edge-of-BARX1-threshold: stage={stage} cluster={cluster_id} "
              f"BARX1={barx1:.1f}% ({side} {BARX1_RECLAIM_THRESHOLD:.0f}% "
              f"by {abs(barx1 - BARX1_RECLAIM_THRESHOLD):.1f}pt) — "
              f"may flip on next d2g rerun")

    # Respiratory-splanchnic guard: block the WT1→Mesothelium flip when the
    # cluster is respiratory (TBX5/TBX4/ISL1 high) — its WT1 is decoupled from
    # TCF21 and it is really Ventral spl. (HH25 cl8).  Genuine mesothelium is
    # respiratory-negative and passes.
    respiratory = max(markers["TBX5"], markers["TBX4"], markers["ISL1"])
    if (wt1 > WT1_MESOTHELIUM_THRESHOLD and cell_type != "Mesothelium"
            and respiratory < WT1_MESOTHELIUM_RESPIRATORY_MAX):
        return ("Mesothelium",
                f" | WT1 override: {cell_type} → Mesothelium "
                f"(WT1={wt1:.0f}% > {WT1_MESOTHELIUM_THRESHOLD:.0f}%, "
                f"respiratory={respiratory:.0f}% < {WT1_MESOTHELIUM_RESPIRATORY_MAX:.0f}%)")
    if barx1 > BARX1_RECLAIM_THRESHOLD and cell_type == "Venous-pole spl.":
        return ("BARX1+",
                f" | BARX1 override: Venous-pole spl. → BARX1+ "
                f"(BARX1={barx1:.0f}% > {BARX1_RECLAIM_THRESHOLD:.0f}%)")
    return (cell_type, "")


def compose_cluster_annotation(row_scores: pd.Series, top_uegs: pd.DataFrame,
                               top_wilcoxon: pd.DataFrame,
                               stage: str, qc: dict | None = None,
                               markers: dict | None = None,
                               cluster_id: int | None = None) -> dict:
    """Assign cell_type + cell_type_broad to a single cluster.

    row_scores: pd.Series indexed by lineage-panel name with mean score values.
    top_uegs: DataFrame of top-5 UEGs for this cluster (columns: gene, u_score, p_in).
    top_wilcoxon: DataFrame of top-5 Wilcoxon DE genes (columns: gene, lfc).
    qc: optional dict with 'pct_ribo_mean', 'cycling_pct', 'max_lineage_pct'
        (apoptotic ceiling, over the broad MARKERS incl. pan-mesenchyme) and
        'max_specific_lineage_pct' (ribosomal-arm gate, over
        SPECIFIC_LINEAGE_MARKERS only) for the Contamination + Apoptotic
        pre-checks.
    cluster_id: int leiden ID — used by edge-margin warnings only (so the
        WARN line says "stage=HH25 cluster=10" not "cluster=None").
    """
    # Contamination pre-check (technical junk).  Runs FIRST — a ribosomal /
    # low-complexity or erythroid/blood cluster is a technical artifact and
    # must leave the lineage assignment regardless of any marker fraction
    # (otherwise these early-stage clusters fold into "Skeletogenic LPM").
    # Marker-based (top Wilcoxon / UEG genes) + QC (pct_ribo) + specific-lineage
    # gate, no hardcoded leiden IDs.  qc carries pct_ribo_mean +
    # max_specific_lineage_pct (always supplied by main()); KeyError would
    # surface a missing QC field rather than silently disabling the check.
    if qc is not None:
        is_contam, contam_reason = contamination_check(
            top_uegs, top_wilcoxon,
            pct_ribo_mean=qc["pct_ribo_mean"],
            max_specific_lineage_pct=qc["max_specific_lineage_pct"],
        )
        if is_contam:
            return {
                "cell_type": "Contamination",
                "cell_type_broad": "Unknown",
                "top_lineage": "contamination_rule",
                "top_score": float("nan"),
                "runner_lineage": None,
                "runner_score": float("nan"),
                "margin": float("nan"),
                "confidence": "strong",
                "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
                "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
                "reason": contam_reason,
            }

    # Apoptotic pre-check (if QC supplied).  Bypasses WT1/BARX1 + HH13 biology
    # overrides — apoptotic identity dominates lineage assignment regardless of
    # marker fractions (cycling+ribo+ pattern is a dying-cell QC fingerprint,
    # not a lineage signature; WT1+ apoptotic cell is still apoptotic).
    if qc and apoptotic_check(qc.get("pct_ribo_mean", 0),
                              qc.get("cycling_pct", 0),
                              qc.get("max_lineage_pct", 100)):
        return {
            "cell_type": "Apoptotic mes.",
            "cell_type_broad": "Apoptotic",
            "top_lineage": "apoptotic_rule",
            "top_score": float("nan"),
            "runner_lineage": None,
            "runner_score": float("nan"),
            "margin": float("nan"),
            "confidence": "strong" if qc["cycling_pct"] > 90 else "moderate",
            "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
            "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
            "reason": f"QC: ribo={qc['pct_ribo_mean']:.1f} cycling={qc['cycling_pct']:.0f} max_lin={qc['max_lineage_pct']:.0f}",
        }

    scores = {k: float(row_scores.get(k, float("nan")))
              for k in LINEAGE_KEYS if k in row_scores.index}
    scores = {k: v for k, v in scores.items() if not np.isnan(v)}
    # Stage-gate out-of-stage panels (Early_LPM, Differentiating)
    scores = apply_stage_gate(scores, stage)
    # Drop any gated-to-inf entries from ranking
    scores = {k: v for k, v in scores.items() if v != float("-inf")}
    if not scores:
        return {"cell_type": "ambiguous", "cell_type_broad": "Unknown",
                "reason": "no lineage scores available"}

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top_lin, top_score = ranked[0]
    runner_lin, runner_score = ranked[1] if len(ranked) > 1 else (None, 0.0)
    margin = top_score - runner_score

    if top_score < MIN_SCORE_ABSOLUTE:
        cell_type, broad = "ambiguous", "Unknown"
        confidence = "weak"
        reason = f"no lineage above absolute threshold (top={top_score:.3f})"
    else:
        cell_type, broad = label_from_lineage(top_lin, stage)
        if margin >= 0.10:
            confidence = "strong"
            reason = f"top={top_lin} margin={margin:.3f} over {runner_lin}"
        elif margin >= MIN_SCORE_MARGIN:
            confidence = "moderate"
            reason = f"top={top_lin} margin={margin:.3f} over {runner_lin}"
        else:
            confidence = "weak"
            reason = f"top={top_lin} near-tie with {runner_lin} (Δ={margin:.3f})"

    # HH13 biology override — panels underperform at pre-specification stage
    if markers is not None:
        override = hh13_biology_override(stage, cell_type, markers)
        if override != cell_type:
            reason += f" | HH13 override: {cell_type} → {override}"
            cell_type = override
            # Lookup broad via LABEL_TO_BROAD (covers stage-gated +
            # ambiguous + override-only labels).  KeyError raises if the
            # override target is missing from the dict — fail-fast caught
            # by.
            broad = LABEL_TO_BROAD[cell_type]

    # Branchiomeric gate: Skeletogenic LPM must carry a real skeletal marker in
    # its top DE — else the Somatic panel fired on pan-mesenchyme PRRX1/TWIST1
    # and the cluster is a generic fibroblast (2026-07-02 vet).  Reads the DE/UEG
    # lists (always supplied), so no `markers is not None` guard.
    new_label, gate_suffix = branchiomeric_gate(cell_type, top_wilcoxon, top_uegs)
    if new_label != cell_type:
        reason += gate_suffix
        cell_type = new_label
        broad = LABEL_TO_BROAD[cell_type]

    # Stage-agnostic WT1/BARX1 overrides (2026-05-04 colleague feedback)
    if markers is not None:
        new_label, reason_suffix = wt1_barx1_override(
            cell_type, markers,
            cluster_id=cluster_id,
            stage=stage,
        )
        if new_label != cell_type:
            reason += reason_suffix
            cell_type = new_label
            broad = LABEL_TO_BROAD[cell_type]

    # IRX1 syringeal gate: put "Syringeal" on the IRX1+ syrinx precursor
    # (Chondrogenic-with-IRX1 → Syringeal; Syringeal-without-IRX1 → runner-up).
    # 2026-07-02 vet.  runner_label = the runner-up lineage's display label.
    runner_label = (label_from_lineage(runner_lin, stage)[0]
                    if runner_lin else "Ventral spl.")
    new_label, gate_suffix = irx1_syringeal_gate(
        cell_type, runner_label, top_wilcoxon, top_uegs)
    if new_label != cell_type:
        reason += gate_suffix
        cell_type = new_label
        broad = LABEL_TO_BROAD[cell_type]

    return {
        "cell_type": cell_type,
        "cell_type_broad": broad,
        "top_lineage": top_lin,
        "top_score": round(top_score, 4),
        "runner_lineage": runner_lin,
        "runner_score": round(runner_score, 4) if runner_score is not None else None,
        "margin": round(margin, 4),
        "confidence": confidence,
        "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
        "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
        "reason": reason,
    }


# ── Per-cluster laterality measurement (no L/R; within-pair library contrast) ──

_LAT_CACHE = {}

def _load_laterality_context():
    if _LAT_CACHE:
        return _LAT_CACHE
    iso = pd.read_csv(DATA / "d2j_pitx2_isoform_percell.tsv.gz", sep="\t")
    mapping = pd.read_csv(DATA / "b_genotype_embryo_mapping.tsv", sep="\t")
    g2e = (mapping[["genotype_id", "embryo_id"]].drop_duplicates()
           .set_index("genotype_id")["embryo_id"].to_dict())
    _LAT_CACHE["iso"] = iso
    _LAT_CACHE["g2e"] = g2e
    return _LAT_CACHE


def _pooled_laterality(pool_df: pd.DataFrame) -> dict:
    """Same as compute_cluster_laterality but operates on an already-filtered
    DataFrame (for pooling across multiple clusters of one label)."""
    ctx = _load_laterality_context()
    iso = ctx["iso"]; g2e = ctx["g2e"]
    sub = pool_df.merge(iso[["barcode", "library", "pitx2c_umi"]],
                        on=["barcode", "library"], how="left")
    sub["pitx2c_umi"] = sub["pitx2c_umi"].fillna(0)
    sub["embryo"] = sub["genotype"].map(g2e)
    sub["lib"] = sub["library"].map(lambda l: TECH_REP_MERGE.get(l, l))
    n_det = 0; folds = []; signs = []
    for emb, esub in sub.groupby("embryo"):
        libs = sorted(esub["lib"].dropna().unique())
        if len(libs) != 2: continue
        a, b = libs
        a_umi = esub.loc[esub["lib"] == a, "pitx2c_umi"].sum()
        b_umi = esub.loc[esub["lib"] == b, "pitx2c_umi"].sum()
        if a_umi == 0 and b_umi == 0: continue
        n_det += 1
        fold = (max(a_umi, b_umi) + 0.5) / (min(a_umi, b_umi) + 0.5)
        folds.append(fold)
        signs.append(+1 if a_umi >= b_umi else -1)
    if not folds:
        return {"n_embryos_detected": 0, "median_fold": None,
                "frac_strong": None, "n_concordant_sign": None}
    folds = np.asarray(folds)
    strong = int((np.log2(folds) >= 1).sum())
    pos = sum(1 for s in signs if s == +1)
    neg = sum(1 for s in signs if s == -1)
    return {
        "n_embryos_detected": n_det,
        "median_fold": round(float(np.median(folds)), 2),
        "max_fold": round(float(np.max(folds)), 2),
        "frac_strong": round(float(strong / len(folds)), 2),
        "n_concordant_sign": max(pos, neg),
    }


def _recommend(pooled_lat: dict, best_cl_fold: float | None,
               member_cls: list[int]) -> str:
    """Classify a proposed pseudobulk group as KEEP, SPLIT, or BILATERAL."""
    if len(member_cls) == 1:
        return "single cluster — no merge question"
    pooled_fold = pooled_lat.get("median_fold") or 1.0
    if pooled_fold < 1.5 and (best_cl_fold or 0) < 1.5:
        return "BILATERAL (all clusters weak) — pool is fine but carries no laterality"
    if best_cl_fold and pooled_fold / best_cl_fold < 0.5:
        return f"SPLIT — pool dilutes {best_cl_fold:.1f}× → {pooled_fold:.1f}× (2×+ loss)"
    if best_cl_fold and pooled_fold / best_cl_fold < 0.7:
        return f"SPLIT (marginal) — pool dilutes {best_cl_fold:.1f}× → {pooled_fold:.1f}× (30%+ loss)"
    return "KEEP — pool preserves per-cluster signal"


def compute_cluster_laterality(d2g_df: pd.DataFrame, leiden_col: str,
                               cluster: int) -> dict:
    """Per-embryo PITX2c asymmetry (direction-agnostic; alphabetical library
    ordering per pair) for cells in one cluster.  Returns:
        n_embryos_detected:  embryos with pitx2c > 0 in both libraries
        median_fold:         median |max(A,B) / min(A,B)| across embryos
        frac_strong:         fraction of embryos with |log2 fold| >= 1
        n_concordant_sign:   max(pos, neg) across embryos (uninterpretable
                             across pairs but diagnostic of within-pair
                             signal)
    """
    ctx = _load_laterality_context()
    iso = ctx["iso"]; g2e = ctx["g2e"]
    sub = d2g_df[d2g_df[leiden_col].astype(int) == cluster].copy()
    sub = sub.merge(iso[["barcode", "library", "pitx2c_umi"]],
                    on=["barcode", "library"], how="left")
    sub["pitx2c_umi"] = sub["pitx2c_umi"].fillna(0)
    sub["embryo"] = sub["genotype"].map(g2e)
    sub["lib"] = sub["library"].map(lambda l: TECH_REP_MERGE.get(l, l))

    n_det = 0; folds = []; signs = []
    for emb, esub in sub.groupby("embryo"):
        libs = sorted(esub["lib"].dropna().unique())
        if len(libs) != 2:
            continue
        a, b = libs
        a_umi = esub.loc[esub["lib"] == a, "pitx2c_umi"].sum()
        b_umi = esub.loc[esub["lib"] == b, "pitx2c_umi"].sum()
        if a_umi == 0 and b_umi == 0:
            continue
        n_det += 1
        fold = (max(a_umi, b_umi) + 0.5) / (min(a_umi, b_umi) + 0.5)
        folds.append(fold)
        signs.append(+1 if a_umi >= b_umi else -1)
    if not folds:
        return {"n_embryos_detected": 0, "median_fold": None,
                "frac_strong": None, "n_concordant_sign": None}
    folds = np.asarray(folds)
    strong = (np.log2(folds) >= 1).sum()
    pos = sum(1 for s in signs if s == +1)
    neg = sum(1 for s in signs if s == -1)
    return {
        "n_embryos_detected": n_det,
        "median_fold": round(float(np.median(folds)), 2),
        "max_fold": round(float(np.max(folds)), 2),
        "frac_strong": round(float(strong / len(folds)), 2),
        "n_concordant_sign": max(pos, neg),
    }


# ── Main ────────────────────────────────────────────────────────────────

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cluster_rows = []

    for stage, label, d2g_stage in STAGE_LABELS:
        uegs_path  = DATA / f"d2g_uegs_{label}.tsv"
        wil_path   = DATA / f"d2g_wilcoxon_{label}.tsv"
        score_path = DATA / f"d2g_lineage_scores_{label}.tsv"
        if not (uegs_path.exists() and wil_path.exists() and score_path.exists()):
            print(f"[skip] {label}: missing one of uegs/wilcoxon/scores TSVs")
            continue

        uegs   = pd.read_csv(uegs_path, sep="\t")
        wil    = pd.read_csv(wil_path, sep="\t")
        scores = pd.read_csv(score_path, sep="\t")

        # Need per-cluster QC (pct_ribo, cycling_pct, max_lineage_pct) for
        # the Apoptotic pre-check.  Recompute from the per-cell d2g TSV.
        d2g = pd.read_csv(DATA / f"d2g_lpm_subclustering_{label}.tsv.gz", sep="\t")
        MARKERS = ["NKX6_1","FOXF1","GATA4","ZFPM2","TBX5","TBX4","BARX1",
                   "ISL1","HAND2","PRRX1","WT1","NR2F2"]

        for _, r in scores.iterrows():
            res = r["resolution"]
            cl = int(r["cluster"])
            leiden_col = f"leiden_{res}"
            sub = d2g[d2g[leiden_col].astype(int) == cl]
            # SPECIFIC lineage markers present in this export (pan-mesenchyme
            # PRRX1/NR2F2 excluded by construction — see SPECIFIC_LINEAGE_MARKERS).
            specific_avail = [g for g in SPECIFIC_LINEAGE_MARKERS
                              if f"{g}_cp10k" in sub.columns]
            # Per-cluster QC
            qc = {
                "pct_ribo_mean": float(sub["pct_ribo"].mean()) if "pct_ribo" in sub.columns else 0,
                "cycling_pct": (((sub["phase"] == "S") | (sub["phase"] == "G2M")).mean() * 100
                                if "phase" in sub.columns else 0),
                # Apoptotic ceiling: max over the broad MARKERS (incl. pan-mes).
                "max_lineage_pct": max(
                    ((sub[f"{g}_cp10k"] > 0).mean() * 100) if f"{g}_cp10k" in sub.columns else 0
                    for g in MARKERS
                ),
                # Ribosomal-arm gate: max over SPECIFIC lineage markers only.
                "max_specific_lineage_pct": max(
                    (float((sub[f"{g}_cp10k"] > 0).mean() * 100) for g in specific_avail),
                    default=0.0,
                ),
            }

            cluster_uegs = (uegs[(uegs["resolution"] == res)
                                 & (uegs["cluster"] == cl)]
                            .sort_values("rank"))
            cluster_wil = (wil[(wil["resolution"] == res)
                               & (wil["cluster"] == cl)]
                           .sort_values("rank"))
            row_scores = r.drop(["stage_label","resolution","cluster","n_cells"])

            # Per-cluster marker fractions (for HH13 biology override + WT1/BARX1
            # override).  All MARKERS are required: a missing column would
            # silently disable an override path (default-0 below threshold).
            # A future d2g exporter rewrite that drops a marker column must
            # raise here, not silently ship wrong labels.
            missing_cols = [f"{g}_cp10k" for g in MARKERS
                            if f"{g}_cp10k" not in sub.columns]
            if missing_cols:
                raise RuntimeError(
                    f"d2g_lpm_subclustering_{label}.tsv.gz is missing marker "
                    f"columns required for compendium overrides: {missing_cols}"
                )
            markers_pct = {g: float((sub[f"{g}_cp10k"] > 0).mean() * 100)
                           for g in MARKERS}
            ann = compose_cluster_annotation(row_scores, cluster_uegs, cluster_wil,
                                             stage=d2g_stage, qc=qc,
                                             markers=markers_pct,
                                             cluster_id=cl)
            # Per-cluster laterality (pitx2c, no L/R — within-embryo pair contrast)
            lat = compute_cluster_laterality(d2g, leiden_col, cl)
            cluster_rows.append({
                "stage": d2g_stage,
                "variant": "ccreg" if label.endswith("_ccreg") else "default",
                "resolution": res,
                "cluster": cl,
                "n_cells": int(r["n_cells"]),
                **ann,
                **{f"lat_{k}": v for k, v in lat.items()},
            })

    out = pd.DataFrame(cluster_rows)
    out_path = OUT / "compendium_cluster_labels.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path} ({len(out)} rows)")

    # Pseudobulk grouping recommendation: per (stage, variant, cell_type),
    # compare per-cluster laterality against pooled laterality.  If the
    # pooled signal is much weaker than the best-cluster signal, pooling
    # dilutes and the group should be split.
    pb_rows = []
    # Only consider the 'ccreg' variant at the proposed resolutions
    PROPOSED_RES = {"HH13": 0.2, "HH17/16": 0.3, "HH21": 0.3, "HH25": 0.2}
    for (stage, variant), grp in out.groupby(["stage", "variant"]):
        res = PROPOSED_RES.get(stage)
        if res is None:
            continue
        grp = grp[grp["resolution"] == res]
        if variant != "ccreg":
            continue
        # Load the per-cell subclustering TSV once per (stage, variant)
        label = {"HH13": "HH13_ccreg", "HH17/16": "HH1716_HH16_ccreg",
                 "HH21": "HH21_ccreg", "HH25": "HH25_ccreg"}[stage]
        d2g = pd.read_csv(DATA / f"d2g_lpm_subclustering_{label}.tsv.gz", sep="\t")
        leiden_col = f"leiden_{res}"
        for ct, g2 in grp.groupby("cell_type"):
            if ct in ("ambiguous", "Apoptotic mes."):
                continue
            member_cls = g2["cluster"].tolist()
            member_ns = g2["n_cells"].tolist()
            # Pooled laterality: treat all member clusters as one pool
            pooled_df = d2g[d2g[leiden_col].astype(int).isin(member_cls)].copy()
            pooled_lat = _pooled_laterality(pooled_df)
            per_cl_folds = g2["lat_median_fold"].dropna().tolist()
            best_cl_fold = max(per_cl_folds) if per_cl_folds else None
            pb_rows.append({
                "stage": stage, "variant": variant, "resolution": res,
                "cell_type": ct,
                "n_clusters": len(member_cls),
                "member_clusters": ",".join(str(c) for c in member_cls),
                "member_n_cells": ",".join(str(n) for n in member_ns),
                "total_n_cells": sum(member_ns),
                "per_cluster_median_folds": ",".join(f"{f:.1f}" for f in per_cl_folds),
                "best_cluster_fold": round(best_cl_fold, 2) if best_cl_fold else None,
                "pooled_median_fold": pooled_lat["median_fold"],
                "pooled_n_embryos": pooled_lat["n_embryos_detected"],
                "pooled_frac_strong": pooled_lat["frac_strong"],
                "dilution_ratio": (round(pooled_lat["median_fold"] / best_cl_fold, 2)
                                   if best_cl_fold and pooled_lat["median_fold"] else None),
                "recommendation": _recommend(pooled_lat, best_cl_fold, member_cls),
            })

    pb_df = pd.DataFrame(pb_rows)
    pb_path = OUT / "pseudobulk_grouping_recommendation.tsv"
    pb_df.to_csv(pb_path, sep="\t", index=False)
    print(f"wrote {pb_path} ({len(pb_df)} rows)")

    # Comparison against the threshold-rule labels for the proposed-rules
    # panel (ccreg + proposed resolution).
    # HH13, HH17/16, HH21 use leiden_0.3 (same as baseline);
    # HH25 uses leiden_0.2 (dropped resolution).
    PROPOSED_RES = {"HH13": 0.2, "HH17/16": 0.3, "HH21": 0.3, "HH25": 0.2}
    panels = pd.read_csv(OUT / "panel_cell_labels.tsv", sep="\t")
    flips = []
    for stage, res in PROPOSED_RES.items():
        cc = out[(out["stage"] == stage) & (out["variant"] == "ccreg")
                  & (out["resolution"] == res)]
        rules_labels = (panels[panels["stage"] == stage]
                        .groupby("panel_rules").size()
                        .reset_index(name="n_cells_rules"))
        # Note: this is a coarse cross-tab — proper cell-level comparison
        # requires joining on barcode, which we'd do in the next stage.
        for _, r in cc.iterrows():
            flips.append({"stage": stage, "resolution": res,
                          "cluster": r["cluster"], "n_cells": r["n_cells"],
                          "compendium_label": r["cell_type"],
                          "confidence": r["confidence"],
                          "top_uegs": r["top_uegs"]})
    flips_df = pd.DataFrame(flips)
    flips_out = OUT / "compendium_vs_rules_flips.tsv"
    flips_df.to_csv(flips_out, sep="\t", index=False)
    print(f"wrote {flips_out} ({len(flips_df)} rows)")


if __name__ == "__main__":
    main()
