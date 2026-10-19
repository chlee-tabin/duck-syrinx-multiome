# Figure 4, panels B–E — provenance & implementation notes (repo `figures/fig02/`)

> Provenance and build notes for `figures/fig02/`. **Not part of the publication figure legend** (see `LEGEND.md` for that). This file holds data-source ownership, build steps, and the methodology detail trimmed from the legend for brevity.

## Data sources

- `references/d2g_subcluster_annotations.tsv` — stage-restricted (`d2g`) subcluster labels. Owned by `pipeline/d2g_validate_subclusters.py` (a thin wrapper over `scripts/annotate_d2g_via_compendium.py`). Never hardcode leiden IDs; figure scripts read labels through `figures/_shared/subcluster_map.py`.
- `figures/_shared/subcluster_map.py` `POOLS` — per-stage inference / bilateral-control / excluded pool membership used for the panel C domain folding.
- `figures/data/d2g_lpm_subclustering_*_ccreg.tsv.gz` — per-cell markers (panels A and C).
- `figures/data/d2j_pitx2_isoform_percell.tsv.gz` — per-cell PITX2 isoform UMIs (panel A PITX2 columns; panel C quantification).
- `figures/data/d2i_pileup_output/d2i_pitx2_rna_pileup.tsv.gz` — per-stage PITX2-locus pileup depth (panel D).

## Build / regeneration

Panels are produced by the scripts in this directory; the compositor `figures/fig02/generate.py` writes `figures/fig02.pdf`. Rerun the d2g → subcluster-map chain after any clustering rerun (`conda run -n duckfig python3 pipeline/d2g_validate_subclusters.py`) before regenerating panels.

## Methodology trimmed from the legend

### Panel A
- At coarser clustering resolutions the HH25 syrinx precursor is not separable from the surrounding ventral splanchnic mesenchyme; the split is only exposed at the d2g stage-restricted resolution. (The 2026-07-02 vet folded the former `Chondrogenic LPM` cl3 into `Syringeal`, so the duck emits no `Chondrogenic LPM` at HH25.)

### Panel B
- Stage-restricted subclustering is necessary because developmental stage dominates the global variance and masks the dorsal–ventral axis.
- A non-myocardial venous-pole splanchnic population (GATA4+ / ZFPM2+ / TBX20+ / TBX4-low) is captured by dissection-edge sampling at HH13–HH21 but is anatomically caudal-medial of the heart (van Wijk et al. 2009, *Circ Res* schematic) and is not present in this syrinx-level cross-section. It is omitted from the schematic and exposed in figure S2 panel E.
- The "Skeletogenic LPM" zone (6) was previously labeled "peripheral somatic mesenchyme"; renamed per the Somatic→Skeletogenic LPM reclassification (FOXF1-low, HOX-negative, TBX1+/EYA1+/EYA2+/EBF2/3+ pharyngeal-cranial skeletal mesenchyme; bilateral).

### Panel C
- Domain categories fold the d2g populations as: "LPM mes. (HH13)" = HH13 inference pool pre-D/V (Early LPM + Cardiopharyngeal LPM); "Ventral splanchnic pool (HH17+)" = full HH17+ inference pool, aligned with figure S2 panel D (HH17/16 = Cardiopharyngeal LPM + BARX1+ — the seeded 2026-08 res-0.3 partition emits no Early LPM at HH17/16; HH21 = Ventral spl. + Syringeal + BARX1+; HH25 = Ventral spl. + Syringeal); "Dorsal sub-epi." = bilateral control (Dorsal-biased LPM at HH17/16; Dorsal sub-epi. at HH21/HH25).
- Mesothelium (bilateral coelomic surface) and the venous-pole splanchnic population (off-axis dissection-edge introgression) are excluded from both pools.
- Library A/B labels are arbitrary (lower-numbered library = A); biological L/R comes from the 4-metric L/R inference (figure S2 panel D; `pipeline/d2l_lr_inference.py`).
- Which genotypes are drawn (`figures/_shared/paired_libraries.R`, shared with `generate_paired_laterality_isoform.R`): a genotype in two libraries is a half-vs-half point only when its library pair is an opposite-side pair in `design/embryo_lr_design.tsv` (tech reps collapsed; the 9 pairs equal the former hardcoded list). Same-side duplicates are excluded and logged: genotypes the raw d2l sidecar (`figures/data/d2l_lr_inference.tsv`) marks `same_side_dup` on both library rows — lane-loading-record tube group 17-B(2) (embryos 17-5/17-6, not repo embryo 17-B), one left half split across GEX04 and GEX0S-duck (G3_unfiltered_donor0 = g46, G3_unfiltered_donor3 = g52). Their count per library pair must equal the design's same-side embryos (17-5, 17-6), and a two-library genotype outside every design L/R pair raises. d2l's `inferred_side` / `pair_constraint_ok` do not decide pairing. Effect (2026-10-14, corrected design): HH17/16 ventral facet 13 → 11 embryos, dorsal 9 → 7; the 4D label rule now selects g36, g37, g41 (cutoff unchanged, 0.247 = g42); g52 is no longer drawn.
- Per-facet y-scales are independent (LPM mes. HH13 ≈ 0–30%; Ventral splanchnic pool HH17+ ≈ 0–40%; Dorsal sub-epi. ≈ 0–80%) so the smaller-range facets are not squished by the dorsal ceiling. Aspect ratio fixed to 1 per facet so y = x is corner-to-corner; sub-panels share facet columns so the two metrics compare column-for-column.
- Expected behavior: HH13 and HH25 off-diagonal (Wave 1 / Wave 2 PITX2c); HH17/16 near the diagonal (pre-TES noise, residual Wave 1 decline, Wave 2 not yet ramped); HH21 intermediate; Dorsal sub-epi. on-diagonal at all stages.

### Panel D
- Locus window `NC_051775.1:3,959,887–3,977,754`; per-cell pileups as a 50 bp rolling mean.
- Strand filter: PITX2-sense (+ strand) via `samtools view -F 16` (empirical cellranger GN-tag convention). Removes antisense ncRNAs LOC113843538 (inside PITX2 intron 1) and LOC113843539 (overlapping the PITX2 3′UTR).
- The pre-correction HH21/HH25 "PITX2a peak" at the dorsal 3′UTR was largely LOC113843539 (antisense) and is substantially reduced after the strand filter.
- PITX2c exon 1c = `3,971,385–3,972,064` (insulated from annotated antisense).
- ⚠ **Two of the three HH25 panel-D tracks are PRE-VET** (`dorsal_subepi_HH25` is current: 1,834 = SSOT cl7). `figures/data/d2i_barcodes/syringeal_HH25_*.txt`
  total 2,068 barcodes = the pre-2026-07-02 HH25 cl5, which the vet reassigned to
  `Ventral spl.`; the current `Syringeal` is cl3 (2,625 cells, IRX1+/PTCHD4+) — a different
  cell set. `ventral_splanchnic_HH25` is likewise the pre-vet Ventral spl. (14,273 vs 16,901
  now). The track label is therefore deliberately NOT routed through `display_population`:
  the hedged name asserts a marker identity these cells do not carry.
  HH21 tracks are current (Syringeal 2,827 = SSOT cl3, 100% barcode overlap), so the HH21
  track IS hedged while the two pre-vet HH25 tracks are not — the exemption is per-track
  (`prevet` column in the `tracks` tribble), not blanket. Regenerating barcodes + pileup
  needs the cluster BAMs (open).
- Two-tier y-axis: dorsal tracks (Dorsal sub-epithelial, Dorsal-biased LPM) at their natural per-cell depth scale; non-dorsal tracks (All LPM mesenchyme, Ventral LPM pooled, Ventral splanchnic, Syringeal) at dorsal_max / 3.5 (3.5× amplified) so their smaller peaks stay legible — a shared axis would compress them because the real (non-antisense) dorsal pan-PITX2 signal at the shared exons / 3′UTR dominates in magnitude.
- Gene model: ZJU1.0 (GCF_015476345.1). PITX2a and PITX2b share an exon at 3,966,591–3,966,646 (CDS start 3,966,601); PITX2c has an independent exon 1c with separate CDS (3,971,830). Nine intronic A-rich tracts (identified via a polyA-priming scan) are marked.

## Detail trimmed from the legend at the manuscript sync (2026-10-15)

`LEGEND.md` now carries the final manuscript Figure 4 legend (panels A–D here = Figure 4B–E; Figure 4A, the global atlas UMAP, is rendered by `figures/_manuscript_layout/render_umap_repel.py` *(not distributed)* from the same UMAP table as the uncomposited `figures/fig01/panel_C_umap.png`, with repelled labels). The legend text it replaced is kept here, unchanged except where noted.

- **Header (before the sync):** "Figure 2: Stage-restricted subclustering resolves dorsal–ventral LPM zones with domain-specific PITX2 laterality."
- **Panel A:** columns = 10 structural marker genes; color = mean CP10K (light gray at zero). The syrinx precursor's chondrogenic identity is marker-defined, while its assignment to the syrinx rather than to general airway cartilage follows from the dissected tissue and is not independently validated by in-situ hybridization. The PRRX1-marked skeletal compartment is split into the branchiomeric Skeletogenic LPM (EYA1/EBF/TBX1+, HH17/16 onward) and, at HH13/HH17-16, a distinct PRRX1+/FOXF1-low non-splanchnic mesenchyme that lacks skeletal transcription factors. (The final Figure 4 renders this row "Syringeal chondro."; the manuscript legend quotes it as "Syringeal.".)
- **Panel B, the six zones:** (1) dorsal sub-epithelial mesenchyme (NKX6-1+; bilateral total PITX2 at the shared exons / 3′UTR, panel D); (2) tracheo-esophageal septum (BARX1+; transient at HH21, absent by HH25); (3) ventral splanchnic mesenchyme (TBX5+ / TBX4+; lateralized PITX2c); (4) syringeal chondrocyte precursor / bulla (PITX2c-only, lateralized; the callouts give the adjacent ventral splanchnic zone the larger L/R fold); (5) mesothelium (WT1+ / TCF21+; the coelomic surface lining wrapping the splanchnic block, bilateral by anatomy); (6) Skeletogenic LPM (PRRX1+; bilateral pharyngeal-cranial skeletal mesenchyme).
- **Panel C:** facets = "LPM mes. (HH13)" (the HH13 inference pool, before dorsal–ventral resolution), "Ventral splanchnic pool (HH17+)" and "Dorsal sub-epi." (bilateral control); points colored by inferred stage. The pre-sync legend said HH17/16 points "cluster near the diagonal" and that HH17/16 laterality "cannot be resolved from RNA alone"; the manuscript wording is "shows no consistent separation" and "no longer reliably resolving" (2026-10-15). It did not mention the two unplotted same-side embryos (17-5, 17-6); the legend now does.
- **Panel D:** tracks span the isoform trajectory — HH13 all LPM mesenchyme (Wave 1 PITX2c); HH17/16 dorsal-biased and ventral LPM (declining PITX2c); HH21 and HH25 dorsal sub-epithelial, ventral splanchnic and syringeal populations (ventral PITX2c Wave 2). The antisense ncRNAs removed by the strand filter lie one within PITX2 intron 1 and one overlapping the 3′UTR. The dorsal signal concentrates at the shared exons (2/3/4) and the 3′UTR rather than the isoform-specific exon-1 windows (total PITX2 mRNA, not PITX2a-specific); PITX2c exon 1c is structurally insulated from annotated antisense, so the PITX2c laterality of the HH13 pan-LPM track and of the HH21/HH25 ventral splanchnic and syringeal tracks is robust to the strand filter. Track names: see the pre-vet HH25 track note above.

## Related docs

- `docs/lpm_mesenchyme_subcluster_biology.md` — subcluster biology, pool-membership rationale, the Somatic→Skeletogenic LPM rename, and the mesothelium / venous-pole reclassification.

## Print sizing (PNAS final, 2026-10, font pass)

Panels B (schematic), C (paired laterality) and D (pileup) also feed manuscript
Fig 4 (C/D/E). Fig 4 is laid out to print full width (17.8 cm, x0.667) and no
taller than 22.5 cm (native <= 955 pt; `build_fig4_rna.py`, the print-layout builder, not distributed, asserts it):
row 1 = UMAP (A, 270 pt) + schematic (C, 474 pt); row 2 = dotplot (B, ~327 pt) +
paired laterality (D, 417 x 436 pt); row 3 = pileup (E, 756 x 222 pt). Each panel is
authored at its slot size, so native text pt x 0.667 = printed pt; text is 9-10.5 pt
native (>= 6.0 pt printed, the PNAS floor), letters 16 pt (10.7 pt printed).

- Schematic (`generate_spatial_svg.py`): 950x541 landscape canvas, anatomy drawn
  through one scaled group (`ART_S = 0.68`), callouts in side columns at 19-20 SVG
  units (>= 6.3 pt printed), Left/Right rotated at the drawing's lateral midline,
  the lumens keyed "E"/"T" (expanded in the bottom-left key; the full words no
  longer fit the lumens), the swatch legend as 3 rows x 2 columns.
- **Callout values are pinned to the submitted figure.** The HH21 L/R callouts
  read `figures/data/d2g_subcluster_annotations_submitted.tsv` +
  `figures/data/d2g_dotplot_spatial_extended_submitted.tsv` (the pre-seeded-cascade
  versions), reproducing the accepted values
  (0.90x / 1.81x / 19.8x / 15.4x). `DUCK_SPATIAL_LIVE_DATA=1` renders from the live
  TSVs instead, which moves them to 0.89x / 1.82x / 20.3x / 16.2x (no manuscript
  text quotes these numbers).
- Paired laterality: theme 9.5-10 pt, ggrepel labels 3.25 mm, `max.overlaps = Inf`
  (no label dropped), legend at the bottom, the bottom row's repeated strips
  dropped, facets fill the slot (not forced square; per-facet symmetric limits keep
  y = x corner to corner). At this size ggrepel cannot fully separate every embryo
  label in the dense low-PITX2 corners (it reports unresolved overlaps after 60 s).
  Label placement (2026-10-14): each labelled embryo's label starts
  0.16 x the facet axis range from its point, in the one of 8 compass directions whose
  leader segment and label anchor keep the largest clearance from every other point and
  already-placed label of the facet (`label_nudges()`); `point.padding = 0.15` so the
  leader ends on the target ring. Before this, the dorsal PITX2+ fraction facet's g37
  label sat on the diagonal touching a neighbouring ring with no visible leader.
- Pileup: theme 9-10 pt, one labelled y tick per track (the track's top pretty
  break; baseline 0), gene-model isoform labels `label_size = 3.4`.
