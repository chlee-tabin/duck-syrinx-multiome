# LPM mesenchyme subcluster biology and pool determination

Methods note: what each d2g stage-restricted LPM subcluster label means, the anatomy behind it,
and which labels enter the left/right (L/R) laterality pools. The label strings themselves are
emitted by `scripts/annotate_d2g_via_compendium.py` into
`references/d2g_subcluster_annotations.tsv`; pool membership is `POOLS` in
`figures/_shared/subcluster_map.py` (and its R mirror). Where anything here disagrees with
those two files, those files win.

---

## TL;DR

Each d2g stage-restricted subcluster maps to a defined anatomical compartment
in the avian foregut/respiratory tract LPM. Two principles govern whether a
compartment appears in fig02 main vs. only in figS02 supplementary:

1. **On-axis radially-distinct compartments** (mesothelium = coelomic
   surface lining of the foregut at the syrinx axial level) → **expose** in
   the main figure, **exclude from the laterality inference pool** (no PITX2
   lateralization expected from coelomic surface cells).
2. **Off-axis tissue introgression** (venous-pole / proepicardium-adjacent
   splanchnic LPM, captured by dissection edges at the cardiac pole) →
   **omit** from fig02 main, **expose in figS02** for transparency, **exclude
   from the inference pool** entirely.

The retained inference pool ("Ventral splanchnic pool" in fig02 panel C) is
the LPM-lineage cells in the syrinx-forming cross-section that are spatially
positioned to carry PITX2c laterality.

---

## 1. Subcluster inventory and biological identity

Clusters come from the d2g cell-cycle-regressed (`ccreg`) variant at canonical resolutions
(HH13 = leiden_0.2; HH17/16 = HH21 = HH25 = leiden_0.3). Leiden cluster numbers are arbitrary
and are not used here; labels are assigned from marker evidence by the compendium scorer. Cell
counts are the per-stage totals in `references/d2g_subcluster_annotations.tsv`.

### 1.1 Inference pool (lateralized populations carrying PITX2c)

| Label | Stages (cells) | Markers | Biological identity |
|---|---|---|---|
| **Ventral spl.** | HH21 (7,343), HH25 (15,194) | TBX5+, TBX4+, ZFPM2+, ISL1+, NKX6-1-low, BARX1-low | Ventral sub-epithelial splanchnic mesenchyme around the trachea/bronchi. Yan et al. 2025 ventral sub-epi. mes. equivalent. **Carries Wave 2 PITX2c laterality.** |
| **Syringeal**, displayed **"Syringeal chondrocyte precursor"** (§1.6) | HH21 (2,972), HH25 (2,584) | HH21: HAND2+, ISL1+, IRX1+, TBX5-mod. HH25: WNT5A+, RUNX2+, DPP10+, **PTCHD4+ (highly specific)**, IRX1+ | **Inferred** tracheo-bronchial junction / bulla cartilage precursor. Stage-gated label (= "Cardiopharyngeal LPM" pre-TES). The chondrocyte-precursor part is marker-validated at HH25; the *syringeal* (anatomical) part is inferred from dissection provenance and is **not spatially validated** (§1.6). |
| **BARX1+** | HH17/16 (3,156), HH21 (5,272) | BARX1+, LRIG1+, SMOC1+, NR2F2+ | Tracheo-esophageal septum (Woo et al. 2011; absent in the Yan et al. atlas). Transient HH17–HH21 population. |
| **Early LPM** | HH13 (1,196) | TBX5+, FOXF1+, ISL1+ | Pre-D/V-specification pan-LPM. Stage-gated. |
| **Cardiopharyngeal LPM** | HH13 (226), HH17/16 (2,260) | HAND2+, ISL1+, IRX1+ | Pre-syringeal cardiopharyngeal field. |

`Chondrogenic LPM` (RUNX2+, WNT5A+, DPP10+; tracheal/bulla cartilage primordia) is **not emitted
by the duck multiome**: the IRX1 gate (§3.3) assigns the IRX1+ cartilage-committed HH25 cluster to
`Syringeal`. The label is correct in two other arms: the chick control, which has generically
chondrogenic mesenchyme but no IRX1+/PTCHD4+ precursor, and the E6/E9/E11 scRNA-seq dataset,
where it is a genuine E09–E11 population (figS07/figS08/fig06). A few step scripts and sidecar
tables in this repository still carry the string `Chondrogenic LPM` for HH25; read it as the
earlier name of `Syringeal`.

### 1.2 Bilateral control

| Label | Stages (cells) | Markers | Biological identity |
|---|---|---|---|
| **Dorsal sub-epi.** | HH21 (1,996), HH25 (2,645) | NKX6-1+, PITX2+, WNT2B+, PTCH1+, GATA4-low, HAND2-low | Dorsal sub-epithelial mesenchyme between esophagus and trachea. Bilateral pan-PITX2 (shared exons / 3′UTR, not PITX2c-asymmetric). Yan et al. 2025 dorsal sub-epi. mes. equivalent. |
| **Dorsal-biased LPM** | HH17/16 (2,564) | NKX6-1+, PITX2+ (modest), WNT2B+ | Pre-TES dorsal-leaning LPM. Supports a dorsal lineage, but commitment is incomplete pre-TES. |

### 1.3 Excluded: on-axis, radially distinct

| Label | Stages (cells) | Markers | Biological identity |
|---|---|---|---|
| **Mesothelium** | HH17/16 (277), HH21 (310) | **WT1+, TCF21+**, FOXF1-low, GATA4-low, TBX5-low | Coelomic surface mesothelium lining the foregut/respiratory tract. Yan et al. 2025 mesothelium equivalent. Anatomically the **outermost radial layer** wrapping the splanchnic mesenchyme. |

### 1.4 Excluded: off-axis tissue introgression

| Label | Stages (cells) | Markers | Biological identity |
|---|---|---|---|
| **Venous-pole spl.** | HH13 (300), HH17/16 (1,735), HH21 (2,375), HH25 (2,312) | GATA4+, ZFPM2+, TBX20+, TBX4-low, BARX1-low, WT1-low | **Non-myocardial splanchnic envelope of the venous pole / vitelline veins / proepicardium**, captured by foregut-dissection-edge sampling. **Not in the syrinx-forming respiratory cross-section** (§2.4). |

### 1.5 Excluded: skeletogenic, fibroblast, paraxial, artifact

| Label | Stages (cells) | Reason |
|---|---|---|
| **Skeletogenic LPM** | HH17/16 (4,644), HH21 (3,491), HH25 (715) | FOXF1-low, HOX-negative, TBX1+/EYA1+/EYA2+/EBF2/3+ pharyngeal-cranial skeletal mesenchyme (a pre-chondrogenic progenitor, LPM by d2g parent gating, not body-wall somatopleure). L/R-symmetric, so excluded. The lateralized skeletal lineage is `Syringeal`, which is in the inference pool. |
| **PRRX1+/FOXF1-low mes.** | HH13 (354), HH17/16 (2,056) | Generic ECM fibroblast (no branchiomeric skeletal TF in the top markers; §3.3 gate 1) |
| **Contamination** | HH13 (157), HH17/16 (1,017) | Ribosomal/low-complexity or erythroid technical signal, assigned by the compendium `contamination_check` (marker/QC-based) |
| **Sclerotome / Dermomyotome** | rare | Paraxial mesoderm, not LPM |
| **ambiguous** | residual | Compendium near-tie, no clean assignment |

### 1.6 Naming honesty: what "Syringeal" does and does not claim

Every label in §1.1–§1.5 is a **transcriptional** identity: a cluster's marker profile scored
against lineage panels. Most of them name a cell *type* (chondrogenic, mesothelial, sclerotome),
which is what markers can support. **"Syringeal" is different: it names an anatomical
structure**, and the data behind it does not reach that far.

**What is established for this population.** (a) At **HH25**, a coherent, cluster-defining
chondrogenic program: WNT5A / RUNX2 / DPP10, with PTCHD4 a highly specific marker (24% of
cells vs at most 2% in any other HH25 cluster) and IRX1 at 54% vs at most 30% elsewhere. The HH21
population carrying the same label is **not** chondrogenic by markers (its `markers_hit` are
IRX1 / EYA4 / MPPED2 / LIX1 / LOC113840455, no cartilage gene); it shares IRX1, LIX1 and
LOC113840455 with the HH25 population's top-50 Wilcoxon genes, which is the lineage-continuity
bridge. The "chondrocyte precursor" reading is earned
at HH25 and carried back to HH21 by lineage continuity, not by HH21 marker evidence. (b) The
tissue was dissected as *syrinx-fated* (the whole multiome comes from the syrinx-forming
region). (c) The *precursor identity* is duck-specific: the chick's matched arm has generically
chondrogenic mesenchyme (WNT5A+) but no IRX1+/PTCHD4+ cluster.

**What is not established.** No in situ, lineage trace or trajectory places these cells in the
syrinx specifically, as opposed to general airway/tracheal cartilage. The paper's molecular
spatial evidence at this locus is **PITX2 HCR**; the multiome is used for PITX2 expression and
PITX2-promoter accessibility. The anatomical half of the name rests on dissection provenance
plus cross-species anatomy (chick Yan et al. 2025 / mouse, §2), not on this cluster's own
spatial validation.

**How that is handled.** The **code label stays `"Syringeal"`**: it is the key in
`references/d2g_subcluster_annotations.tsv`, `POOLS`, every derived table and every
filter/join. **Reader-facing labels are hedged at render time** via
`figures/_shared/subcluster_map.py::display_population` (and its R sibling):

| context | rendered as |
|---|---|
| figure legends, prose | **Syringeal chondrocyte precursor** |
| panel/axis/facet text (width-limited) | **Syringeal chondro.** |
| code, POOLS, table values, joins, filters | `Syringeal` *(unchanged)* |

The hedged name keeps the marker-validated part (chondrocyte precursor), marks the anatomical
part as inferred, and stays distinct from the chick's generic `Chondrogenic LPM`.
`display_population` is **display-only and additive**: never call it on a value used to filter,
join, group or compare.

---

## 2. Anatomical evidence

The references below anchor the biological interpretation. Figures are cited by PMC URL; none
are redistributed here.

### 2.1 LPM splitting establishes splanchnic vs somatic layers

**van den Berg & Moorman 2009** ([PMC2691808 Fig 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691808/figure/Fig2/)).

Paired ventral views + transverse sections from gastrulation through heart
looping (~HH8–HH13 / E7.5–E9.0). LPM splits into:
- **Somatopleure** (somatic mesoderm): outer, body-wall-facing
- **Splanchnopleure** (splanchnic mesoderm): inner, gut-facing

The cardiogenic crescent forms by **medial swing of splanchnic LPM** and is
contiguous with non-cardiogenic splanchnic mesoderm at the foregut. This
contiguity is why a splanchnic-LPM cluster (`Venous-pole spl.`) can share
gene programs with cardiac mesoderm without committing to cardiomyocyte fate.

### 2.2 First heart field = lateral splanchnic mesoderm

**Tzahor 2011** ([PMC3125075 Fig 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC3125075/figure/F1/)).

Chick HH8–HH16 transverse sections. The first heart field sits at
the **lateral splanchnic mesoderm**; pharyngeal mesoderm (SHF) is
medial. Both fields are continuous with non-cardiogenic splanchnic LPM; the
gene-program boundary is fuzzy, which is why GATA4+/TBX20+/ZFPM2+ cells
exist on the splanchnic-LPM side of the boundary without committing to
cardiomyocyte fate.

### 2.3 Heart elongation and venous pole

**Buckingham 2014** ([PMC4200205 Fig 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC4200205/figure/F2/)).

Mouse E7.5 → E10.5 (≈chick HH8–HH18): cardiac crescent → linear heart tube →
looped heart with arterial + venous poles + SHF added progressively. The
**venous pole** is where atria and inflow tract form, anatomically at the
caudal-medial side of the heart, immediately adjacent to the foregut at
this developmental window.

### 2.4 TBX18+ splanchnic envelope at the venous pole

**van Wijk et al. 2009**, *Circ Res* (PMID 19628790) ([PMC2861358 Fig 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC2861358/figure/F1/)).

Chick HH11, HH16, HH21: 3D reconstructions + transverse insets.
**TBX18+** splanchnic mesoderm covers the vitelline veins, then
the proepicardium, then envelopes the heart as epicardium. At HH21 this
is the splanchnic envelope of the venous pole, distinct from
myocardium and pharyngeal-mesoderm-derived OFT.

**This anatomy matches the `Venous-pole spl.` marker fingerprint:**
GATA4+/ZFPM2+/TBX20+/HAND2-variable/TBX4-low/TBX5-low/WT1-low. The
`Splanchnic_vis` lineage panel (GATA4+ZFPM2+TBX20+ up; TBX5+TBX4+BARX1+HAND2+
down) is functionally a venous-pole / proepicardium-adjacent splanchnic
detector.

**Why this is "off-axis" relative to the syrinx:** the venous pole sits
caudal-medial of the heart, not at the trachea-bronchial bifurcation. A
foregut dissection at the syrinx level that extends slightly caudal-ventral
will catch some of this tissue, but it is not part of the syrinx-forming
respiratory mesenchyme.

### 2.5 Respiratory splanchnic mesenchyme: a distinct gene program

**Whitsett et al. 2016** *Trends Dev Biol* ([PMC8547324 Fig 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8547324/figure/F2/)).

Foregut transverse section. The ventral splanchnic mesoderm wrapping the
trachea expresses **Foxf1, Islet1, Gata4/6, Tbx4, Tbx5, Wnt2, Wnt2b**.
Respiratory splanchnic GATA4 co-occurs with **TBX4 + TBX5 + WNT2/2B**; the
`Venous-pole spl.` populations are TBX4-low/TBX5-low and do not match this
signature, so they are not respiratory-tract splanchnic.

### 2.6 Trachea cross-section: the radial layers sampled

**Whitsett et al. 2016** *Trends Dev Biol* ([PMC8547324 Fig 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC8547324/figure/F3/)).
**Iber & Mederacke 2022** *Front Cell Dev Biol* ([PMC9094403 Fig 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9094403/figure/f1/)).

Mouse trachea cross-sections at E12.5–E14.5. Radial layers from lumen
outward: NKX2-1+ epithelium → SOX9+ ventral cartilage zone → ACTA2+ dorsal
smooth-muscle zone → TBX4/5+/WNT5A+ peripheral mesenchyme.

These layers map onto the retained inference + bilateral populations:
- Inner mesenchyme = `Dorsal sub-epi.` + `Ventral spl.`
- Outer surface = `Mesothelium`

No cardiac-pole / venous-pole sub-population is drawn at the trachea
cross-section level in either reference, consistent with omitting the
venous-pole splanchnic population from the main figure.

### 2.7 Heart-trachea spatial relationship over development

**Kishimoto & Morimoto 2021** *Development* ([PMC8276987 Fig 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8276987/figure/F2/)).

Mouse E7 → E18.5 timeline. The heart sits ventral and slightly caudal to
the trachea/lung primordium throughout development, and the axial separation
between heart venous pole and tracheal bifurcation increases as the embryo
elongates. Venous-pole tissue therefore enters a syrinx-level dissection
only at its caudal-ventral edge.

---

## 3. Pool determination

Pool membership lives in `figures/_shared/subcluster_map.py` `POOLS` (and
the R mirror). Three pool kinds per stage:

- **inference**: cells pooled for L/R PITX2c laterality testing (fig02
  panel C, figS02 panel D 4-metric inference)
- **bilateral**: cells pooled as bilateral control (expected weak/null
  laterality)
- **excluded**: cells pooled into nothing

### 3.1 Decision matrix

| Stage | Inference | Bilateral | Excluded |
|---|---|---|---|
| HH13 | Early LPM, Cardiopharyngeal LPM | (D/V unresolved pre-TES) | Venous-pole spl., PRRX1+/FOXF1-low mes., Mesothelium, Contamination, paraxial |
| HH17/16 | Cardiopharyngeal LPM, BARX1+ | Dorsal-biased LPM | Venous-pole spl., Skeletogenic LPM, PRRX1+/FOXF1-low mes., Mesothelium, Contamination, paraxial |
| HH21 | Ventral spl., Syringeal, BARX1+ | Dorsal sub-epi. | Venous-pole spl., Skeletogenic LPM, Mesothelium, Contamination, paraxial |
| HH25 | Ventral spl., Syringeal | Dorsal sub-epi. | Venous-pole spl., Skeletogenic LPM, Mesothelium, Contamination, paraxial |

(Mirrors `POOLS`, with `Sclerotome`/`Dermomyotome` condensed to "paraxial" and the `ambiguous`
residual omitted. A label listed as excluded at a stage where no cluster carries it is a no-op.
`Syringeal` is rendered "Syringeal chondrocyte precursor" for readers, §1.6.)

### 3.2 Pool membership rule

Changing pool membership means editing `POOLS` (not the `pipeline/d2g_validate_subclusters.py`
wrapper, which only regenerates the annotation table). Never hardcode leiden cluster IDs: they
are arbitrary labels with no cross-run meaning, and every consumer resolves populations through
`get_pool_map(stage, kind)` / `get_population_clusters(stage, population)`.

### 3.3 Override rules and marker gates (codified in `scripts/annotate_d2g_via_compendium.py`)

Two post-compendium overrides change a label based on per-cluster marker fractions:

```python
WT1_MESOTHELIUM_THRESHOLD = 35.0        # fraction of cells WT1-positive
BARX1_RECLAIM_THRESHOLD = 40.0          # fraction of cells BARX1-positive
WT1_MESOTHELIUM_RESPIRATORY_MAX = 30.0  # max(TBX5, TBX4, ISL1) ceiling for the WT1 rule

def wt1_barx1_override(cell_type, markers):   # simplified
    wt1   = markers["WT1"]    # explicit lookup — KeyError if marker absent
    barx1 = markers["BARX1"]
    if (wt1 > WT1_MESOTHELIUM_THRESHOLD and respiratory < WT1_MESOTHELIUM_RESPIRATORY_MAX
            and cell_type != "Mesothelium"):
        return "Mesothelium"
    if barx1 > BARX1_RECLAIM_THRESHOLD and cell_type == "Venous-pole spl.":
        return "BARX1+"
    return cell_type
```

**WT1>35% rationale:** above this threshold a cluster is dominated by WT1+TCF21+ mesothelial
identity regardless of which lineage panel scored highest. The respiratory guard keeps the rule
off respiratory-splanchnic clusters whose WT1 is decoupled from TCF21 (coelomic mesothelium is
TBX5/TBX4/ISL1-low; genuine mesothelium sits at ≤13% respiratory markers, so a 30% ceiling
separates the two). Where the Epicardial lineage panel already wins, the override corroborates
the panel call rather than flipping it. The 35% floor is sensitivity-adjusted from the Yan et
al. 2025 mesothelium benchmark (WT1+TCF21+ at 30–80% in scRNA-seq, lower in snRNA-seq because of
nuclear capture).

**BARX1>40% AND label == Venous-pole spl. rationale:** reclaims a cluster the compendium scores
as visceral splanchnic but whose biology is BARX1+ T-E septum. The 40% threshold avoids touching
low-BARX1 clusters, and gating on the `Venous-pole spl.` label keeps it off Ventral spl. or
Syringeal clusters with moderate BARX1. **Precedence:** the WT1 rule is checked first.

**Caveat (threshold fragility).** The scorer prints a warning for any cluster within
`EDGE_WARN_MARGIN` (8 percentage points) of either threshold, so a borderline cluster is visible
at every run. The merged HH21 BARX1+ population can mix septal and visceral-leaning BARX1+ cells.

Three **stage-agnostic, marker-based** gates (no leiden IDs) correct over-firing lineage panels:

1. **`branchiomeric_gate`**: `Skeletogenic LPM` must carry at least one true branchiomeric
   marker (EYA1/EYA2/EBF1/EBF2/EBF3/TBX1/COL25A1/SIX1) in its top-50 DE/UEG genes, else it
   becomes `PRRX1+/FOXF1-low mes.` (generic ECM fibroblast; the `Somatic` panel's
   PRRX1/TWIST1 are pan-mesenchyme and would otherwise mislabel fibroblasts as skeletogenic).
   Both labels are excluded, so pool membership does not depend on this gate.
2. **Respiratory guard on the WT1→Mesothelium override** (above).
3. **`irx1_syringeal_gate`**: the syrinx precursor is IRX1+. A post-TES `Syringeal` cluster
   lacking IRX1 takes its runner-up lineage label; a `Chondrogenic LPM` cluster with IRX1
   becomes `Syringeal`. The gate is **duck-only by design**: applied to the chick it is a
   no-op, and the chick lacks this precursor.

---

## 4. Display rules

- **fig02 main** (dotplot + HH21 schematic) omits `Venous-pole spl.` and shows `Mesothelium`
  as the outer radial layer; the schematic draws a thin outer mesothelial ring around the
  splanchnic block. Rendered row labels are shortened for axis width ("Skeletogenic LPM" →
  "Skeletogenic") and `Syringeal` is hedged by `display_population` (§1.6); the code keys are
  unchanged.
- **figS02** (extended dotplot / full-compendium view) lists every d2g compendium-assigned
  population, including `Venous-pole spl.`, for transparency; the legend states that it is
  excluded from the laterality inference pool as off-axis tissue introgression.
- **fig02 panel C** (paired laterality) shows the inference pool ("Ventral splanchnic pool"),
  the bilateral control ("Dorsal sub-epi.") and the HH13 pre-D/V pool, per §3.1. Mesothelium
  and Venous-pole spl. are not rows.
- **fig02 panel D** (PITX2 locus pileup) is per population per stage; no Venous-pole spl. or
  Mesothelium tracks are rendered.

---

## 5. Where each decision is recorded

| Decision | Source |
|---|---|
| Cluster identity | `figures/data/d2g_lpm_subclustering_*_ccreg.tsv.gz` per-cluster marker fractions |
| Lineage panel definitions | `pipeline/d2g_stage_restricted_subclustering.py` |
| Compendium scoring, overrides and gates | `scripts/annotate_d2g_via_compendium.py` |
| Labels per (stage, cluster) | `references/d2g_subcluster_annotations.tsv` |
| Pool membership and display names | `figures/_shared/subcluster_map.py` (`POOLS`, `POPULATION_DISPLAY_NAMES`) and `figures/_shared/subcluster_map.R` |

Each subcluster's identity claim traces to the marker fractions (quantitative), the compendium
lineage scores (decision logic), the Yan et al. 2025 alignment and the anatomy in §2.
