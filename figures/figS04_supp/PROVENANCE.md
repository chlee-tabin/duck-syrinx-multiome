# Chick multiome QC (supports Figure S12) — provenance & implementation notes (repo `figures/figS04_supp/`)

> Provenance and build notes for `figures/figS04_supp/`. **Not part of the publication figure legend** (see `LEGEND.md` for that). This file holds data-source ownership, build steps, and the methodology detail trimmed from the legend for brevity.

## Data sources

- `pipeline/chick/cb_libqc_export.py` → `cb_libqc_perlane.tsv` — per-lane snRNA QC metrics (panel A).
- `pipeline/chick/cAK_atlas_sex_wgene.py` → `cAK_atlas_wz_percell.tsv.gz` — per-cell W/Z-gene fractions, computed from the CellBender atlas h5ad so donor IDs match Fig. S4 exactly, after the resolved-12 relabel (panel B).
- `pipeline/chick/cA_crosslane_3lane.py` → `cA_crosslane_3lane_long.tsv` — embryo × embryo cross-lane allele-frequency correlations (panel C). GEX0S-chick genotyped at the combined 283,435-variant base + projected via `arc0s/vireo_donorfile`.
- CellBender atlas cells + `chick_stage_alignment.tsv` (stage) — per-embryo L/R half cell counts (panel D).
- `design/chick_groupB_lr_design.tsv` — the designed Group-B lane→side map (panel D coloring; also referenced in panel A's depth-floor reasoning).
- `figures/data/chick_atac_metrics.tsv`, `figures/data/chick_atac_per_cell_sample.tsv`, `figures/data/chick_atac_f1_qc_summary.tsv`, `figures/data/chick_atac_f1c_amulet_summary.tsv` — per-library chick snATAC QC (panel E; producers `pipeline/chick/hQ1`–`hQ4`, staged by `hQ0_stage_atac_qc_inputs.sh`).

## Build / regeneration

- Panel scripts: `generate_libqc.py` (A), `generate_sex_violins.py` (B), `generate_demux_heatmap.py` (C), `generate_lr_cells_per_stage.py` (D). Panel D applies the designed Group-B lane→side via the orientation-class logic in `pipeline/chick/_chick_lr_orientation.py`, keyed by (group, donor) and (stage, gfp) since donor IDs are reused across groups.
- Compositor: `figures/figS04_supp/generate.py` → `figures/figS04_supp.pdf`.
- Panel E (snATAC QC) renderer: `generate_libqc_atac.py`; it and `panel_libqc_atac.png` live in `figures/figS04_supp/` and are composed as panel E by `generate.py`.

## Group B split doublet/unassigned exclusion (2026-10)

Panels B and D drop the 148 Group B cS-split doublet/unassigned cells present in the committed
`cAK_atlas_wz_percell.tsv.gz` via `drop_split_excluded()` (list:
`figures/data/chick_control/groupB_split_excluded_cells.tsv`; rationale in `figures/figS04/PROVENANCE.md`).
Only donor11a / donor2a nucleus counts change. Panels A, C, E are not affected (library QC,
pre-split genotype concordance, ATAC library QC).

## Methodology trimmed from the legend

### Panel A
- Unlike the duck, CellBender did not re-call cells for the chick — counts were corrected on the cellranger barcodes — so there is no cellranger-vs-CellBender cell-count delta to show.
- The within-pool GEX09 vs GEX10 sequencing-depth difference (not droplet ambient) is what makes the apparent PITX2 lane floor in Fig. S4F survive CellBender: PITX2 is lowly expressed, so its detection is depth-limited.
- Depth detail (moved from the legend): depth varies ≈1.6-fold across the four snRNA libraries (median UMI/nucleus 2,211–3,545). GEX09 is the shallowest; GEX10 is ≈39% deeper than its Group B partner GEX09 (3,080 vs 2,211 UMI/nucleus; 1,791 vs 1,387 genes/nucleus). The lane asymmetry is two-sided: GEX09 has ≈2× more nuclei called (25,100 vs 12,622) at lower depth and a higher doublet rate (23% vs 13%); GEX10 fewer nuclei at higher depth. (The chick snATAC QC is now panel E of this figure — formerly the standalone Figure S5-supp.)

### Panel B
- Sex violins are computed from the CellBender atlas h5ad (donor IDs therefore match Fig. S4 exactly), after the resolved-12 relabel that splits the cS donor11/donor2 merges → donor11a/b, donor2a/b.
- The within-donor M+F merge gate (>20% of cells in both the male-W and female-W modes — the duck-grade test) flags 0 of 12, confirming the cS splits resolved the donor2 and donor11 merges.
- Two sparse uncalled phantom slivers (donor5 n = 1 and donor13 n ≈ 35) are not among the 12 resolved embryos and are omitted; donor13 adjudicates female.
- Result (moved from the legend): all 12 resolved embryos are cleanly unimodal; the within-donor M+F merge gate flags 0 of 12, confirming clean sex resolution. Group B = 2 males (donor0, donor7; GFP⁺) + 10 females.

### Panel C
- GEX01 × GEX02 (12 × 12): the two complementary-half lanes were demultiplexed jointly, so this is the cross-lane concordance check — a within-lane merge of two embryos would show as a low same-embryo r. Same-embryo median r = 0.94 vs different-embryo 0.42, 0 mis-assignments.
- GEX0S-chick × GEX01 and GEX0S-chick × GEX02 (8 × 12): GEX0S-chick is the Left-only HH15-16/HH25 symmetry-breaking lane, genotyped at the identical 283,435-variant base and projected onto the already-solved Group A donors by reference-VCF vireo (2,956 cells → 1,345 assigned singlets; the remainder is the expected duck carry-over, correctly left unassigned). Because GEX0S-chick cells were assigned by matching these same genotypes, the same-embryo diagonal is expected by construction — these two grids are a projection-assignment QC, not an independent measurement. They show the assignment is clean (diagonal median r = 0.84, not eroded by contaminant cells; 0.31 off-diagonal confirms the donors stay genotypically resolvable in the GEX0S-chick data). The independent signal is the pattern of presence: exactly the 8 HH15-16 + HH25 left-half embryos light up, while the 4 HH13 embryos (d3/d4/d6/d15) carry no GEX0S-chick row at all, matching the Left-only HH15-16/HH25 design (a built-in negative control).

### Panel D
- Group B's lane→side map is designed; its design-orientable embryos (uniform or singleton design class) are colored Left/Right — including donor8 (HH25-GFP⁺, same uniform class as donor11b → GEX09 = L / GEX10 = R). The rest are grey "half" bars where the per-embryo side is genuinely not resolvable — all of Group A (ARC01/ARC02 are complementary halves of genetically identical sides) and the non-uniform Group B classes (HH13-GFP⁺, HH15-16-WT).
- donor8's Right/GEX10 half is ≈97% under-captured (32 cells), so donor8 is the one embryo excluded from the half-vs-half laterality in Fig. S4F — a concrete data-quality exclusion, independent of GFP status — even though its design side is known.
- Sampling-balance framing (moved from the legend): both physical halves are captured for every embryo, with per-half yields varying several-fold — this panel is the chick sampling-balance QC (analog of duck Fig. S2).

### Panel E (chick snATAC QC — folded in from the former Figure S5-supp)
- 5 chick ARC (multiome, galGal6a `crarc_galgal6a`) libraries: Group A `ARC01`/`ARC02` (complementary half-lanes) + `ARC0S-chick` (Left-only symmetry-breaking anchor, the single v2.0.0 library); Group B `ARC09`/`ARC10`. GFP⁺ is a pan-GFP demultiplexing label only.
- Dark = all cellranger-arc-called nuclei (`is_cell==1`); light = the `f1` snATAC-QC-passed subset (TSS enrichment ≥ 2 & ≥ 1000 fragments per nucleus; `hQ3_snapatac_f1_qc.py`). The E1 light bar counts the snATAC `f1` set (can exceed dark for ARC0S-chick); E2–E6 light violins are the cellranger∩snATAC per-nucleus distributions.
- E1 est_cells / f1-passed: ARC01 4,337/3,452; ARC02 5,732/4,410; ARC0S-chick 2,956/3,305; ARC09 21,175/16,314; ARC10 10,614/7,722. Per-nucleus median SnapATAC2 TSS enrichment 5.54/5.59/5.22(ARC09)/5.66(ARC0S-chick)/5.27(ARC10), all clear of the gate; f1 retention ~73–80% for the v2.1.0 libraries.
- E2/E4 = log10 per-nucleus raw reads / HQ fragments; E3 = duplicate rate (`dup_reads/raw_reads`, saturation proxy); E5 = per-nucleus TSS-fragment fraction; E6 = FRiP at cellranger-arc per-library peaks (stand-in until the chick ArchR peak set lands); E7 = AMULET (q < 0.01, galGal6a autosomes only) doublet rates ARC01 5.5% / ARC02 12.8% / ARC09 3.5% / ARC0S-chick 11.5% / ARC10 4.4% (the two elevated are the most sequencing-saturated). AMULET uses the duck's three patched-upstream workarounds.
- Producer → sub-panel: hQ1 (cellranger-arc metrics) → E1 dark; hQ3 (SnapATAC2 f1) → E1 light; hQ2 (per-nucleus reservoir sample, seed 42) → E2–E6; hQ4 (AMULET) → E7.

## Related docs

- Chick-control conventions: GFP⁺ embryos are pan-GFP demultiplexing labels only, treated identically to WT in every biological analysis; the only legitimate per-embryo exclusion is a concrete data-quality failure (donor8's under-captured half).
