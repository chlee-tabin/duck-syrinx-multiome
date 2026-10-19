# Figure S14 — provenance & implementation notes (repo `figures/figS07/`)

> Provenance and build notes for `figures/figS07/`. **Not part of the publication figure legend** (see `LEGEND.md`). Holds data-source ownership, build steps, exact statistics, methodology, and rigor caveats trimmed from the legend. Built from the **cr9** (intron-inclusive + CellBender + ccreg) re-derivation of the old E6/E9/E11 ≈ HH25/32/35 (later, growth-phase) 3′ scRNA-seq dataset, supporting the manuscript Figure 7 / Figure S14 ESR1 story. **Does not touch the multiome (FC_08656) pipeline.** The cross-platform panels (C, D, E) read committed multiome `figures/data/` for comparison only. Two-page portrait composite, 7 data-dense panels. Cell-type / population labels follow the cr9 single-source-of-truth (**Skeletogenic LPM**, **Syringeal chondrocyte**, **Syringeal fibroblast** — *no* defunct "Somatic").

## Data sources

- **Substrate (old, cr9-CB):** master `cr9_03_integrated_CB.rds` (24,816 cells) + the substage E06/E09/E11 Seurat objects + the current compendium labels. CellBender-corrected, intron-inclusive (pre-mRNA) counts, cell-cycle regressed (ccreg).
- **Substrate (multiome, read-only — cross-platform panels C/D/E):** committed `figures/data/` —
  - `d2g_signaling_{stage}.tsv.gz` (per-cell `ESR1_cp10k`),
  - `d2g_lpm_subclustering_{stage}_ccreg.tsv.gz` (subcluster membership; `ccreg` is canonical),
  - `lr_inference_multiome.tsv` (canonical L/R side = the d2l L/R-inference step at HH25),
  - `step03_sex_per_genotype.tsv` (per-genotype sex call),
  - the cross-stage MetaNeighbor (subway-map) tables for panel C.
- **Data home (shared old-duck dir, `figures/data/scrna_e6e9e11/`):**
  - `cr9_er_survey_percell.tsv.gz` — atlas A/B (UMAP_1/2, stage, broad_class, substage_population) + onset old-side rows (`log__ESR1`),
  - `substage_E09_percell.tsv.gz`, `substage_E11_percell.tsv.gz` — ESR1 level + cell-cycle `Phase` + `RNA_snn_res.0.2` + `LR.called` + `sex.called` + `donor_id`,
  - `cr9_substage_compendium_clusterlabels.tsv` — res-0.2 cluster → population map (joined on the cluster; fail-loud if a cluster ID is absent),
  - `esr1_level_lr_tidy.tsv`, `esr1_cellcycle_lr_tidy.tsv` — tidy tables exported by the panel-E / panel-F scripts,
  - `er_perembryo_ie.tsv` (gene, library, embryo, umi_exon, umi_intron), `er_intron_exon_umi.tsv` (per-gene pooled exon/intron UMIs, pct_intronic, span_kb), `esrrg_coverage_binned.tsv` (pooled cr6-BAM coverage, 200 bp bins), `esrrg_exons.bed` (16 exons).
- The xstage MetaNeighbor subway-map source figure (panel C image) lives at `figures/figS07/panel_A_pagamap.png`.

## Build / regeneration

Local panels (system Python in the `duckfig` env) + composite:
```
python3 figures/figS07/plot_atlas_umaps.py      # A (panel_atlas_stage.png), B (panel_atlas_celltype.png)
python3 figures/figS07/plot_esr1_onset.py       # D (panel_esr1_onset.png)
python3 figures/figS07/plot_esr1_level_lr.py    # E (panel_esr1_level_lr.png); exports esr1_level_lr_tidy.tsv
python3 figures/figS07/plot_esr1_cellcycle.py   # F (panel_esr1_cellcycle.png); exports esr1_cellcycle_lr_tidy.tsv
python3 figures/figS07/plot_esrrg_refined.py    # G (panel_esrrg_refined.png)
conda run -n duckfig python3 figures/figS07/generate.py   # -> figures/figS07.pdf (2-page portrait, canonical root copy)
```
Panel C is generated upstream by `figures/figS07/xstage_pagamap_mock.py` (the subway map).
Shared helpers: `_esr1_lr.py` (L/R paired engine — `paired_p`, `paired_units` — used by panels E and F); `figures/_shared/subcluster_map.py` (`get_domain_map`, `get_leiden_col`); `figures/_shared/library_naming.py` (`collapse_lib`, the tech-rep SSOT that folds PILOTGEX4A/B + GEX03-1A/B).
Compositor: `generate.py` resizes panels with PIL to harmonised widths (full rows → W; 2-up rows → equal-height pair summing to W), composes each page via `_shared/composite.py`, then merges with `pypdf` → `../figS07.pdf`. No subdirectory PDF duplicate.

## Methodology trimmed from the legend

### Panel D — ESR1 onset trajectory (exact % ESR1⁺)
- Unit = multiome genotype / old-scRNA donor, ≥10 cells/unit; positivity = `ESR1_cp10k > 0` (multiome) or `log__ESR1 > 0` (old). % detection (not CP10K) keeps the cross-platform comparison fair; lines are drawn across stages (unlike the laterality panel) specifically so the onset dynamics and junction kink are legible. Adjacent stage means are connected only between adjacent positions and only across one step — a track missing an interior stage (Syringeal has no E06) breaks rather than bridging HH25→E09 as a single step.
- **Track wiring (2026-08-21, after the seeded relabel).** The seeded relabel swapped the multiome HH25 population names' identities (the ESR1⁺/PITX2⁺ carrier that was "Syringeal" is now inside `Ventral spl.`; today's `Syringeal` is the ex-"Chondrogenic-LPM"-like population), so the tracks follow **carrier continuity, not name identity**. The pre-relabel prose ("Syringeal saturated at 57.7 %"; "Ventral 38.8 % → E06 56.5 %") described the old partition and is retired. The 2026-08-21 MetaNeighbor rerun (panel C) corroborates the carrier at the junction — HH25 Ventral spl. ↔ E06 Ventral spl. is a reciprocal-best link (AUROC 0.83), while today's Syringeal has no RBH cross-link (best affinity: the E06 Skeletogenic skeletal branch, 0.77). Note the lanes here pair populations by **ESR1⁺/PITX2⁺ carrier continuity**, not by global-transcriptome nearest neighbour (within the old data the global fan runs E06 Ventral → E09 Chondrogenic, 0.70) — panel C shows the global map; this panel tracks the marker-defined carrier.
- **all-LPM-mesenchyme aggregate track (grey dashed):** 2.2 % (HH13) → 1.8 % (HH17/16) → 3.6 % (HH21) → **31.9 % at HH25** → 44 → 56 → 58 % (E06→E09→E11).
- **Ventral spl. → Syringeal fibroblast track (carrier):** HH21 3.0 % → **HH25 44.8 %** → E06 56.5 % → E09 55.3 % → E11 59.3 % — smooth across the junction.
- **Syringeal → Syringeal chondrocyte track (sibling):** HH21 5.8 % → **HH25 9.4 %** against 59.3 / 55.4 % at its E09/E11 continuation — a later onset.
- The script fails loud if a `TRACKS` population name is absent from the *raw* (pre-MINCELL) label vocabulary at every stage (guards against a rename silently drawing an empty track), while still permitting legitimate per-stage absence (e.g. Ventral spl. only from HH21).

### Panel E — ESR1 level L/R (exact stats + estimator verification)
- **Estimator (cross-platform consistency).** CP10K = mean `expm1(log1p-cp10k)` — the estimator available on **both** platforms, since the HH25 multiome export carries only log-cp10k and forces it there. The old data additionally carries raw counts, so the count-pooled `Σcnt / ΣnCount` estimator is also computable; it gives **byte-identical L/R conclusions** (verified out-of-band, not recomputed by the committed panel script: E09 6/6 p = 0.031, E11 8/8 p = 0.0078 under both estimators). The same `expm1(log1p-cp10k)` estimator is therefore used throughout for cross-facet consistency.
- **HH25 unit construction.** Unit = multiome genotype; cells are pooled per side (a genotype's side can span >1 library — e.g. g11's two PILOTGEX4 sub-libraries are both Right — so grouping is by side, not library, giving one L and one R pseudobulk). Side from the canonical `lr_inference_multiome.tsv` (= d2l at HH25); `collapse_lib` folds tech-rep sub-libraries before the side join. Sex from `step03_sex_per_genotype.tsv`.
- **Cell floor / pairing.** Only embryos with ≥10 cells on **both** sides are kept (= the paired set used in the test, so every plotted dot has a connector); no connectors across stages (the platform junction is not a quantitative continuum); y is free per facet (within-facet L/R is the message).
- **HH25 lane mapping (2026-08-21, post- relabel):** lanes follow **ESR1⁺/PITX2⁺ carrier continuity** — current `Ventral spl.` → the *Syringeal fibroblast* lane, current `Syringeal` → the *Syringeal chondrocyte* lane (see the panel-D track-wiring note above; the earlier name-identity mapping described the pre-relabel partition). The panel-C MetaNeighbor rerun supports the junction carrier (Ventral↔Ventral RBH 0.83; Syringeal has no RBH cross-link) while noting the global-transcriptome fan within the old data pairs Ventral→Chondrogenic — the lane rule here is the marker-defined carrier, stated in the LEGEND.
- **Per-facet paired Wilcoxon (L vs R, pooled over sex):**
  - **Syringeal fibroblast — Left:** E09 ΔL−R **+2.47, p = 0.031, n = 6**; E11 **+3.15, p = 0.0078, n = 8**; HH25 (= current Ventral spl.) Left trend **+1.85, p = 0.055, n = 8**.
  - **Syringeal chondrocyte** — flat at HH25 (= current Syringeal; **+0.10, p = 0.94, n = 7**); left-leaning but not significant at E09/E11 (+0.45, p = 0.22; +0.49, p = 0.31).
  - **Dorsal sub-epi. control** — flat at HH25 (**+0.02, p = 0.94, n = 7**); the E09 control **p = 0.063** reflects the tilt-suspect old-data L/R direction.
- **Sex composition (counts; rationale in `LEGEND.md` panel E).** HH25 is **6 M / 2 F** (M/F shown descriptively only, never tested); the powered, sex-balanced contrast is the old data — **8 M / 8 F** at E09/E11, the full old dataset **12 M / 12 F** including the L/R-less E06. E11 carries no Dorsal sub-epi. cluster ("not detected").
- Fail-loud guards: missing HH25 lane names in the d2g domain map; any target-lane cell whose (genotype, library) has no/invalid (non-L/R) side; any HH25 genotype with no/invalid sex — each raises rather than silently dropping cells.

### Panel F — cycling fraction (exact stats)
- Cycling fraction = fraction of cells in S or G2M = 1 − G1, from the **Seurat `CellCycleScoring` `Phase`** call already stored on the ccreg substage objects. Unit = donor; real dissection side; E09/E11 only; 8 M / 8 F (E09 4M/4F + E11 4M/4F; full old dataset 12M/12F incl. the L/R-less E06). Shared y (cycling fraction is comparable across facets). Per-facet paired Wilcoxon (L vs R, pooled over sex) where ≥4 pairs.
- **E09 Syringeal chondrocyte:** cycles more on the Left, **ΔL−R +0.137, p = 0.031**.
- **E11 Syringeal chondrocyte (male-specific Left advantage):** male **L 0.48 vs R 0.35**; female **L 0.29 ≈ R 0.32**.
- **Descriptive / exploratory; demotable.** Cell-cycle scoring on intron-heavy single-cell data is noisier than a labelled proliferation assay (trimmed from `LEGEND.md` panel F, which now states the qualitative takeaway only); the EdU (Fig 6c) adjacency is the validated counterpart.

### Panel G — ESRRG intron-inflation control (exact counts + internal-priming method)
- **Per-embryo intron fraction (top).** One dot per embryo per gene (pseudobulk, summed across library), x = % intronic UMIs (linear), dot area ∝ total UMIs, ESRRG highlighted; pooled %-intronic drawn as a black tick.
  - **ESR1** carries the dominant exonic anchor — **81,788 pooled exon UMIs** — at a moderate intron fraction.
  - **ESRRG** is the intron-inflated outlier — **88 % intronic**, only ~12 k exon UMIs.
- **ESRRG gene model (bottom).** 16 exons, − strand, ~404 kb (locus `NC_051774.1`, GS = 20,177,932, GE = 20,581,985; 3′ terminal exon at the **low** coordinate). Coverage peaks (above the 85th-percentile depth, grouped within 2 kb, group-max as the priming site) sit in introns, not at the 3′ terminal exon — the signature of internal priming on pre-mRNA, which explains the broad intronic signal.
- **Pooled stranded cr6-BAM count** (reads, ≈2× the deduplicated UMIs above): **88.7 % of ESRRG-sense reads are intronic — 24,910 exon vs 195,044 intron**. The antisense (+ strand) intragenic LOCs `LOC110352055` / `LOC119716591` are excluded by stranded counting, so no *sense* ncRNA contaminates the count.
- **Conclusion.** The qualitative takeaway lives in `LEGEND.md` (panel G); the count evidence above — ~88 % intronic, with the pileup in introns rather than the 3′ exon — is what underpins it and the *ESRRG\** asterisk in Fig 6a.

### Panel C — MetaNeighbor subway-map mechanics
- **2026-08-22:** the old-arm E09/E11 "Chondrogenic LPM" is renamed **"Syringeal chondrocyte"** throughout (this file's live text updated; the retired pre- multiome label in the track-wiring note below is history and stays). The panel gains a **dotted magenta edge class** — direct HH25↔E09 reciprocal-best matches from the all-vs-all AUROC matrix, drawn as skip edges over the coarse E06 column (Ventral spl. ↔ Syringeal chondrocyte 0.756; Dorsal ↔ Dorsal 0.857; Skeletogenic LPM ↔ Syringeal smooth muscle 0.887; Venous-pole spl. ↔ Syringeal fibroblast 0.823; HH25 `Syringeal` has no reciprocal E09 partner). Rendered by `figures/figS07/xstage_pagamap_mock.py`.
- Rows = populations, columns = stages spanning **multiome HH13–HH25 ‖ old-scRNA E06–E11**, joined at the marked platform junction. Edge classes: black = reciprocal-best MetaNeighbor match (AUROC ≥ 0.5), thin grey = within-platform AUROC-weighted, magenta = cross-platform link (multiome↔old); node area ∝ cells within platform. AUROC ≥ 0.5 is the replicability/similarity score from MetaNeighbor (the neighbour-voting cross-validation AUROC, not a significance test). The map situates the old dataset as the continuation of the multiome series and traces Skeletogenic LPM (E06) → Syringeal chondrocyte (E09→E11) + Syringeal fibroblast.
- **Multiome-side staleness RESOLVED (2026-08-21, same day).** The MetaNeighbor chain was re-run on O2 (multiome export on the current `d2g` SSOT, job 50991311; old-cr9 export + correspondence on the reconstructed compendium labels, job 50991312) and the panel rebuilt on the refreshed `xstage_{composition,corr_ribbons}.tsv`. Junction RBH links: Dorsal sub-epi. 0.91, Ventral spl. 0.83, Skeletogenic LPM 0.82; today's `Syringeal` has no RBH cross-link (best affinity E06 Skeletogenic LPM 0.77 — the skeletal branch); within-old, E06 Ventral spl. → E09 Syringeal chondrocyte 0.70 is the strongest fan edge, and E06 Mesothelium → E09 Syringeal fibroblast 0.68 (correspondence ≠ lineage tracing; the panel states this).

### Panels A / B — atlas
- Atlas = the full cr9-CB master (**24,816 cells**); point order randomised (seed 0) to avoid overplotting bias; UMAP coordinates carry no ticks (arbitrary). Cell type (B) = the cr9 `substage_population` for LPM mesenchyme + the mapped non-mesoderm broad class (Epithelial→Endoderm, Endothelial→Endothelium, Immune→Immune), on the canonical fig02 earth-tone palette (a population keeps one colour across the whole figure). **No silent colour fallback** — a label without a palette entry raises.

## Related docs
- Companion figures built from the same cr9 dataset: `figures/figS08/` (WNT screen), and the figS09/figS10 BMP/TGF-β screens.

## References

- **MetaNeighbor** (panel C cross-stage replicability / AUROC): Crow M, Paul A, Ballouz S, Huang ZJ, Gillis J. *Characterizing the replicability of cell types defined by single cell RNA-sequencing data using MetaNeighbor.* Nat Commun. 2018;9:884. doi:10.1038/s41467-018-03282-0.
- **Cell-cycle scoring** (panel F S/G2M `Phase` calls — the scoring strategy implemented by Seurat's `CellCycleScoring`): Tirosh I, Izar B, Prakadan SM, et al. *Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.* Science. 2016;352(6282):189–196. doi:10.1126/science.aad0501.
- **Seurat v5** (panel F `CellCycleScoring` *implementation*; the cr9 pipeline runs Seurat v5): Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R. *Dictionary learning for integrative, multimodal and scalable single-cell analysis.* Nat Biotechnol. 2024;42:293–304. doi:10.1038/s41587-023-01767-y.
- **CellBender** (ambient-RNA correction, "cr9-CB"): Fleming SJ, Chaffin MD, Arduini A, Akkad A-D, Banks E, Marioni JC, Philippakis AA, Ellinor PT, Babadi M. *Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender.* Nat Methods. 2023;20(9):1323–1335. doi:10.1038/s41592-023-01943-7.
