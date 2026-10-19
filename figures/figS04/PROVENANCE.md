# Figure S12 — provenance & implementation notes (repo `figures/figS04/`)

> Provenance and build notes for `figures/figS04/`. **Not part of the publication figure legend** (see `LEGEND.md` for that). This file holds data-source ownership, build steps, and the methodology detail trimmed from the legend for brevity.

## Data sources

- `design/chick_groupA_lr_design.tsv`, `design/chick_groupB_lr_design.tsv` — the chick Group A / Group B L/R design (wetlab ground truth for the panel A schematic, the Group B designed lane→side map, and the panel F designed-direction concordance labels).
- `references/color_registry.tsv` — colours. Cell types use the shared `CELLTYPE_COLORS` registry (same hue per shared lineage as Figure 1); Left = blue / Right = red; Male = blue / Female = red; GFP⁺ = green; everything else grayscale.
- Atlas (panels B, D): the genotype-demultiplexed chick syrinx multiome, Group A + B pooled (37,362 nuclei), CellBender ambient-corrected counts (`matrix_cb`).
- Staging (panel C): per-cohort LPM-lineage pseudobulk PCA tables (Group A from cellranger-arc v2.1.0, matching the duck; Group B from standalone cellranger-7 GEX).
- Pileup (panel E): per-cell sense-strand coverage over the chick (galGal6a) PITX2 locus, Group B.
- Laterality (panel F): chick stage-restricted subcluster annotations (port of the duck stage-restricted compendium annotator), folded to spatial domains via the chick pool map; per-domain dot accounting in the panel F counts sidecar.
- Panel G: cross-species PITX2 isoform × modality L/R contrast table (duck ATAC ΔFRiP + RNA Δcp10k at the three PITX2 promoters, chick best/worst design-legal envelope).

## Build / regeneration

Panels are produced by the scripts in this directory; the compositor `figures/figS04/generate.py` writes `figures/figS04.pdf`.

- `panel_A_design.py` — panel A experimental-design schematic.
- `generate_umap.py` — panel B whole-embryo atlas UMAP.
- `generate_dotplot.py` — panel D marker dot plot (shared `_shared/dotplot.py` renderer).
- `generate_staging.py`, `generate_staging_separate.py` — panel C per-cohort LPM-lineage pseudobulk PCA (the `_separate` plotter renders the two cohort sub-panels).
- `panel_F_dv_laterality.R` — panel F orientation-free half-vs-half divergence (port of the duck paired-laterality plotter).
- `panel_G_xspecies_isoform.R` — panel G oriented cross-species isoform L/R contrast.

## Group B split doublet/unassigned exclusion (2026-10)

The resolved-12 Group B overlay (`pipeline/chick/_groupB_resolved12.py`) splits two merged vireo
donors with an N=2 re-vireo (donor11 -> donor11a/b, donor2 -> donor2a/b). Until 2026-10 the cells
the split called `doublet`/`unassigned` were `.fillna`-defaulted to `{base}a`: 163 cells (donor11
45 = 37 doublet + 8 unassigned; donor2 118 = 106 doublet + 12 unassigned). They are now EXCLUDED
(doublet removal uses vireo calls only). The list is committed as
`figures/data/chick_control/groupB_split_excluded_cells.tsv`, built by
`pipeline/chick/_groupB_split_excluded.py` from the two standby-backed split files
(`.../chick_backup/chick_demux_B/vireo_seeded/seed_qc/split_donor{11,2}/N2/donor_ids.tsv`, md5s in
the table header). The overlay now raises on any other status or on a cell missing from its split,
and asserts its excluded set equals the committed table. Committed per-cell tables that predate the
fix are filtered at read time by `drop_split_excluded()` (the hit count must be the pre-fix count
or 0, else it raises): cAM 92 cells, cAK/cAE 148, chick PITX2-locus per-cell 72 (ArchR ids).
`groupB_resolved12_percell.tsv.gz` was pruned in place (29,524 -> 29,361 rows; `--prune-resolved12`).

The cluster-side producers that fold the split files themselves (rather than through
`resolve_donor`) carried the same silent default (`sub.get(c, f"{base}a")`, plus an
`os.path.exists` guard that quietly left a donor unsplit when its split file was missing). They now
all call `fold_split()` in `_groupB_split_excluded.py`: split donor0 -> `{base}b`, donor1 ->
`{base}a`, doublet/unassigned keep that status and are dropped by the script's existing
doublet/unassigned filter; a missing split file, a base-donor cell absent from it, or an unexpected
status raises, and for Group B donor11/donor2 the non-clean set is asserted against the committed
table. Affected: `cAD_groupB_pitx2_pileup_d2i.py` (**panel E** — its committed HH13 n = 1,529
includes the 45 donor11 cells, so S12E changes only on the next cluster rerun),
`cZ_groupB_pitx2c_isoform.py`, `cP2_groupB_pitx2_refined.py`, `cB_stage_accounting.py`,
`cAB_cross_stage.py` (diagnostics, not distributed). Two side effects on rerun, both in the same direction as the
fix: cAD/cZ/cP2 also fold the legacy `split_donor8` split (not part of resolved-12), whose
doublet/unassigned cells are likewise dropped instead of defaulting to donor8a (no committed list —
counts are logged); and cAB, which group-prefixes donor labels, previously let base-run
`doublet`/`unassigned` cells through as `A_doublet`/`B_doublet` "donors" — they are now dropped.
On the cluster the scratch `pipeline/chick/` copy has no repo tree, so set
`CHICK_SPLIT_EXCLUDED_TSV` to the synced `figures/data/chick_control/groupB_split_excluded_cells.tsv`.

The table header records the paths actually read (byte-identical local copies of the standby
files, made for the build) with their md5s, and where each copy came from.

## Methodology trimmed from the legend

### Panel B
- Featured Leiden resolution = 0.5 (vs the duck's 0.2): at ~5× lower chick depth this is the resolution that separates sclerotome and dorsal-neural; argmax keeps the LPM one label regardless, so it does not over-split the displayed types.
- A "Low-quality" flag (ribo/mito-dominated clusters) is applied but trips on no cluster on the CellBender-corrected counts — the ambient hemoglobin that drove the raw atlas's single Low-quality cluster (≈1,961 cells; the RBC marker baseline drops ~0.25 → ~0.15 after correction) is removed.
- The five minor duck types that fold into parents at chick depth: Splanchnic mes. → LPM, Pharyngeal endoderm → Epithelial, Ventral neural prog. → Neural progenitor; Cardiomyocyte and Otic placode form no cluster.
- The chick res-0.2 mesenchyme splits ~99% by stage (early/broad EBF1-2/CDH11/SOX5/BNC2⁺, low TBX5/TBX4 vs HH25 TBX5/TBX4⁺ "LPM mesenchyme"), mirroring the duck cl0/cl2 stage split. Sclerotome (PAX1 logFC 6.7 + CNMD/KERA chondrogenic) is surfaced by a margin-promotion (Sclerotome score ≥ 0.30 + positive paraxial probe).

### Panel C
- PC1/PC2 variance explained on CellBender-corrected counts: PC1 (stage) = 53.5% Group A / 72.3% Group B (sharper than raw 45.9 / 60.1%); PC2 (sex) = 15.9% / 9.0%. A separate PCA per cohort is batch-free by construction (Group A and Group B are different cellranger products), so no correction is applied.
- Group A stage anchor: the design-known-LEFT library captured only the HH15-16 + HH25 embryos, so its left-half cell count = 0 for HH13 but 105–381 for HH15-16 and 28–196 for HH25; HH25 vs HH15-16 then splits on the maturation axis. Because Group A's mesenchyme PC1 tracks sex within the early cluster, the two early stages interleave on this axis — the HH13/HH15-16 assignment is anchor/composition-based, not read off the plot.
- Group B HH25 comes from the resolved-12 demultiplexing (vireo N14/seed_105 + the cS splits of two merged donors — donor11 → donor11a/b, donor2 → donor2a/b). The 8 early embryos split into HH13 vs HH15-16 on PC1 (a ~5-unit gap among the GFP⁺ embryos places donor4 alone as the lone HH15-GFP), constrained to the design counts (HH13 = 3 GFP⁺ + 1 WT, HH15-16 = 1 GFP⁺ + 3 WT). The cS splits are corroborated by the maturation axis (donor11a early; donor11b/donor2a/donor2b at HH25).

### Panel E
- Regenerated 2026-10-13 (SLURM 56118015 cb_restage + 56118016 cAD): Harmony on lane now runs directly (the earlier run's scanpy wrapper raised and silently fell back to uncorrected PCA); embryo identity is the resolved-12 set with split-donor doublet/unassigned cells dropped and donor8 unsplit; per-embryo stages come from the panel-C staging (chick_stage_alignment.tsv) instead of the script's own ranking. Track n (old -> new): HH13 LPM 1,529 -> 1,586; HH15-16 ventral 2,775 -> 2,774; HH25 dorsal sub-epi. 873 -> 1,007; HH25 ventral spl. 1,620 -> 1,321; HH25 syringeal-cardio 1,519 -> 1,831 (barcode Jaccard 0.70-0.93). The coverage pattern (5' sense signal at the a/b first exons, upstream of the PITX2c TSS) is unchanged in every track.
- snRNA 3′ bias + the unextended galGal6a 3′UTR annotation under-capture ≈ 2.3% of PITX2 UMIs (median per-cell boost 0%). By designed side the recovered UMIs split 27 LEFT / 21 RIGHT across the WT controls (0/6 laterality flips), so the laterality is unaffected.
- The antisense (−F16) control confirming the 5′ sense signal is genuine a/b transcription (antisense ≈ 0 at the sense 5′ peak) is omitted from the manuscript panel for clarity and retained in the exploratory render.

### Panel F
- Domain fold: the chick stage-restricted subclustering pools Group A + B and re-HVG/PCA/Harmony-on-library per firm stage, then ports the duck stage-restricted compendium annotator; spatial domains are folded via the duck-faithful chick pool map. "LPM mes. (HH13)" = HH13 inference pool (no dorsal split, as in the duck HH13); "Ventral splanchnic pool" = HH15-16 + HH25 inference; "Dorsal sub-epi." = HH15-16 + HH25 bilateral control. Skeletogenic LPM, Mesothelium, Venous-pole spl. and paraxial are dropped at every stage (this corrects a prior fold that retained Skeletogenic LPM + Mesothelium at HH13).
- Axes are linear with per-facet free limits so low-PITX2 domains read as low. Half assignment is arbitrary (lower lane = half 1).
- An aggregate designed-direction control was examined but dropped as lane-confounded: the designed-LEFT-concordant embryos are largely those whose higher-PITX2 lane happens to coincide with their designed-left, so it is not independent evidence of laterality.
- Lane-floor analysis: a per-population audit shows the deeper lane carries a higher PITX2 baseline in 14 of 16 stage × population cells (exceptions: HH15-16 BARX1⁺ and HH15-16 Venous-pole spl., where the shallower lane reads higher), including the excluded Skeletogenic LPM / Venous-pole domains (PITX2⁺ fraction lane→lane at HH13: inference 8.2% → 20.3%, excluded 28.9% → 43.9%; HH13 has no bilateral/dorsal pool, as in the duck). Each Group B embryo's two halves are split across the two lanes, so the half-vs-half comparison is intrinsically a between-lane comparison; lanes are balanced across sides (designed-LEFT is the deeper lane for donor1/2a/2b, the shallower lane for donor4/11b), so the apparent direction is set by which lane an embryo's designed-left landed in. A difference-in-differences lane-offset correction shows the apparent inference-pool asymmetry tracks this laterality-free floor (Pearson r = 0.79, p = 0.004, ~63% of variance), with no stage surviving correction.
- Why the floor survives CellBender: it is a sequencing-depth / detection-sensitivity difference, not droplet ambient. The deeper lane is ~39% deeper (median 3,080 vs 2,211 UMI; 1,791 vs 1,387 genes/cell); PITX2 is lowly expressed so its detection is depth-limited (PITX2⁺ fraction rises ~2.4× with depth); depth-matched the advantage vanishes (the shallower lane ≥ the deeper lane in most depth bins); empty/background-droplet PITX2 is low and not elevated in the deeper lane. CellBender removes ambient counts, not a depth-driven detection floor, so it correctly leaves it: the fraction divergence persists on corrected counts (inference |log2 CP10K| 1.65 → 1.48, not collapsed; the bilateral control rises 0.75 → 0.92, shrinking the inference-vs-control contrast 2.2× → 1.6×), while broadly-leaked ambient like hemoglobin is cleaned (RBC baseline ~0.25 → 0.15). Because PITX2⁺ fraction is itself depth-sensitive, the half-vs-half fraction divergence partly reflects the between-lane depth imbalance; the depth-normalised CP10K row is the more conservative metric.
- Dot accounting (kept per domain): LPM mes. (HH13) 8/8 (4 A + 4 B); Ventral splanchnic pool 12 (5 A + 7 B); Dorsal sub-epi. 13 (6 A + 7 B). Residual drops are 3 Group A HH15-16 embryos with a ~1-cell half (donor5/10/13) + donor8 (Group B, one physical half ~97% under-captured — a concrete data-quality failure, independent of GFP status). The most prominent off-diagonal Group A point (HH25, Dorsal CP10K row; donor16) is a 22-cell-half artifact, near-diagonal in fraction.
- At HH13 the LPM mesenchyme is not design-lateralised here because snRNA measures total PITX2, whereas the Wave-1 asymmetry is PITX2c-isoform-specific (total PITX2 is symmetric at HH13 even in the duck, and 3′-biased snRNA cannot resolve PITX2c).
- The net-readout takeaway (no clean chick PITX2 L/R laterality — the expected negative contrast with the duck) is stated in the figure summary; the panel-F prose was trimmed to the panel description.

- **Split-exclusion effect on Panel F (2026-10):** only donor11a (HH13) and donor2a (HH25) points move.
  Per-half nuclei and |log2| (mean CP10K) old -> new: donor11a LPM mes. (HH13) 62/50 -> 57/49,
  1.341 -> 1.702; donor2a ventral pool 1240/92 -> 1189/89, 0.178 -> 0.237; donor2a dorsal sub-epi.
  235/93 -> 228/92, 0.367 -> 0.357. donor11b / donor2b unchanged; still 33 points, 6 concordant /
  3 discordant design labels. donor2a directional Fisher p 0.124 -> 0.161 (still n.s.).
- **Lane-offset numbers re-checked (2026-10):** the r = 0.79 originally quoted in the legend does
  NOT reproduce from the current committed cAM. Re-running the lane-offset correction (exploratory
  script, not distributed; Pearson r of raw inference asymmetry vs lane offset over the 11 usable
  Group B embryos): June-2026 committed TSV r = 0.794 (p = 0.0035, R^2 0.63, slope 1.59); current
  inputs before the split-donor fix r = 0.451 (p = 0.16, R^2 0.20); after it r = 0.462 (p = 0.15,
  R^2 0.21). The June TSV predates later cAM regenerations (its control-pool ambient floors
  1.77/1.93 vs 3.05/2.65 now). "No stage survives correction" holds in all three (per-stage
  one-sample t on the residual: post-fix HH13 p = 0.99, HH15-16 0.52, HH25 0.33; pooled 0.21).
  The r value was therefore removed from the legend and replaced by the 14/16 lane-floor count +
  "no stage survives correction".
- **Lane-floor audit numbers (post-fix):** still 14 of 16 stage × population cells with GEX10 > GEX09
  PITX2⁺ fraction (same two HH15-16 exceptions; on linear mean CP10K it is 10 of 16, which is why the
  legend names the fraction metric); HH13 PITX2⁺ fraction lane->lane inference 8.2% -> 20.3% becomes
  8.1% -> 20.0%; excluded 28.9% -> 43.9% unchanged.

### Panel G
- Duck L/R per half is taken from the chromatin-arbitrated sides at the HH17/16 trough. Chick per-embryo orientation is design-ambiguous, so chick is shown as a best/worst envelope over the design-legal lane→side maps (best vs worst), so chick cannot be over-claimed positive. Point area ∝ min(L,R-half) cells (MINC = 20); off-scale points become edge ▲/▼ glyphs carrying the true value.

### Panel D
- Gene columns + brackets are derived from the shared duck marker panel (`figures/_shared/markers.py`, the object fig01 panel D renders) minus `KIF26B` and `NFIA` (no galGal6a annotation) — 35 of the 37 S7C genes — with the duck `LEFTY1` slot relabelled to chick `LEFTY2`. Derived, not copied.
- 2026-09-08: the chick script had carried a hand copy of fig01's dict and missed fig01's 2026-08-21 re-bracketing of `TBX5`/`TBX4` under Splanchnic, so the submitted S12D drew an "LPM" bracket over `TBX5`/`TBX4` — inverted relative to S7C. With `KIF26B`/`NFIA` absent there is no LPM bracket in the chick panel at all; `TBX5`/`TBX4` sit under Splanchnic as in S7C. The `ARX` dot is in the "Dorsal neural prog." row (the Neural progenitor marker is `PTPRZ1`); the legend was corrected to match.

## Related docs

- `docs/chick_control/ANNOTATION_DIFF.md` *(not distributed)* — the chick-vs-duck atlas annotation differences (the 11-of-16 recovery, the five fold-ins, and the "Mesenchyme (early)" stage pool).
- `docs/chick_control/laterality_fold_methods.md` *(not distributed)* — the panel F domain-fold method, pool-map rationale, and lane-floor analysis.

- **Split-donor doublets in the atlas display (S12B/S12D; declared rule, 2026-10-13, same as the demux remnants):** the 148 atlas nuclei that the donor11/donor2 N=2 splits call doublet/unassigned were base-vireo singlets of the merged donors, so they are part of the published 37,362-nucleus atlas and stay in the atlas UMAP (panel B) and marker dot plot (panel D). They are excluded from every per-embryo analysis (S12C staging, S12E pileup, S12F half-vs-half, cells-per-stage) via `groupB_split_excluded_cells.tsv`. Dropping them from the atlas display would change the published atlas count without changing any conclusion.
