# Figure S11 — provenance & implementation notes (repo `figures/fig04/`)

> Provenance and build notes for `figures/fig04/`. **Not part of the publication figure legend** (see `LEGEND.md` for that). This file holds data-source ownership, build steps, the methodology detail trimmed from the legend for brevity, and the honest-framing methods disclosures.

## Data sources

- `figures/data/pitx2_locus_gc_content_tad.tsv` — TAD-scope GC% / CpG-island content (panel A GC strip).
- `figures/data/chick_control/xspecies_orientfree_isoform.tsv` — per-embryo orientation-free isoform divergence + RNA divergence vs within-embryo null (panel B); built by `pipeline/chick/build_orientfree_isoform_data.py`.
- `g3f_domain_duck_motifs.tsv`, `g3f_domain_window_gc.tsv`, `g3g_domain_chick_motifs.tsv`, `g4f_domain_window_stageflags.tsv` — the `g3f`/`g3g` (PITX2-locus PWM motif scan) DOMAIN scan outputs over the curated motif windows (panel C). Scan windows built by `scripts/make_g4f_domain_windows.py`; scan run on the cluster via `pipeline/run_g3f_g3g_domain.sh`.

## Build / regeneration

Panel A is rendered by the env-parameterised locus-dynamics engine (`g4f`, the per-peak L/R accessibility engine) `figures/fig04/panel_D_lr_genetrack.R`. Panels B and C are produced by the scripts in this directory. The compositor `figures/fig04/generate.py` writes `figures/fig04.pdf`.

```
# Panel A — TAD gene-track + GC strip (env-parameterised locus-dynamics engine)
G4F_SUFFIX=_floor20 G4F_DESIGN=unpaired G4F_FILL=log2fc \
  G4F_GC_TSV=figures/data/pitx2_locus_gc_content_tad.tsv \
  G4F_OUTDIR=figures/fig04 G4F_OUTNAME=panel_A_tad_genetrack.png \
  conda run -n duckfig Rscript figures/fig04/panel_D_lr_genetrack.R

# Panel B — orientation-free 4-readout laterality magnitude (data builder → renderer)
conda run -n duckfig python3 pipeline/chick/build_orientfree_isoform_data.py
conda run -n duckfig python3 figures/fig04/panel_B_xspecies_orientfree.py

# Panel C — HH25 motif composition (consumes the g3f/g3g DOMAIN scan)
#   windows: scripts/make_g4f_domain_windows.py
#   scan (cluster, priority ~2 min): sbatch pipeline/run_g3f_g3g_domain.sh
#     -> g3f_domain_duck_motifs.tsv, g3f_domain_window_gc.tsv,
#        g3g_domain_chick_motifs.tsv, g4f_domain_window_stageflags.tsv
conda run -n duckfig python3 figures/fig04/panel_C_motif_dotplot.py

conda run -n duckfig python3 figures/fig04/generate.py   # -> figures/fig04.pdf
```

## Methodology trimmed from the legend

### Panel A
- TAD window = `NC_051775.1:3,387,000–4,476,254` (359 peaks); the upstream edge is set at C4orf32/FAM241A, the proximal end of the Pitx2 TAD. Welsh et al. 2015 report Pitx2-locus interaction enrichment consistent with a TAD — from a re-analysis of PUBLISHED mESC data (Hi-C, Dixon et al. 2012; Smc1 ChIA-PET, Dowen et al. 2014), not their own — and place C4orf32/FAM241A in the most-proximal sub-TAD with adjacent CTCF sites; their sub-TAD boundaries are given only as mm9 mESC intervals, which do not transfer to duck, so the window edge is anchored at the gene rather than at a reported boundary.
- Census = the floor-20 census (peak shown only where the pooled count clears the floor; below-floor peaks rendered grey/absent).
- Rendered from the current floor-20 `g4f` outputs, in which the split-tube halves 17-5/17-6 (G3_unfiltered_donor0 / donor3) sum their GEX04 and GEX0S-duck aliquots into one left half each: the ventral HH17/16 column reads 16 L / 11 R and carries 8 starred L peaks, and the dorsal HH17/16 column reads 12 L / 5 R (no significant peak). The HH13, HH21 and HH25 columns do not include these halves. Composite `figures/fig04.pdf` rebuilt; the supplementary booklet (`build_booklet_v5_splice.py`, untracked output) picks it up on its next build.
- The CpG island under the focal L-bias (the PITX2 promoter cluster) is ~88% GC.

### Panel B
- `gap = obs − null`, where the null is the per-embryo within-embryo cell-shuffle distribution; the L/R label is discarded so magnitude is non-circular for both species.
- ATAC rows = null-corrected `|Δlog2 FRiP|` at the PITX2c / PITX2a / PITX2b promoters; RNA row = null-corrected `|Δmean cp10k|` total PITX2.
- **Result values (moved from the legend; match the committed panel PNG):** the per-embryo divergent fraction is duck-specific at PITX2c — duck 50% vs chick 11% (p = 0.0045) — and total PITX2 RNA — duck 81% vs chick 29% (p = 2.3e-4); the within-isoform controls PITX2a (7% vs 0%) and PITX2b (3% vs 11%) show no consistent divergence in either species (n.s.).
- Circle area ∝ min-half cells, with a floor of ≥ 5 cells/half (the null-corrected gap self-calibrates for small n).
- Significance is the figure-relevant duck-vs-chick contrast: a one-sided Fisher on the per-embryo divergent fraction, where `emp_p` is the per-embryo cell-shuffle rank score (the fraction of an embryo's own within-embryo null ≥ its observed divergence) and an embryo is *divergent* at `emp_p < 0.05`. Significance is *not* obs-vs-own-null: on the ratio metric that rank test over-calls a small consistent bias at the a/b controls.
- The oriented best/worst-envelope view (chick Group-B donors placed by the designed lane→side map) is figS04 panel G.

- **Split-exclusion (2026-10):** the chick inputs drop the Group B cS-split doublet/unassigned
  cells (`groupB_split_excluded_cells.tsv`; 92 cAM cells, 72 locus ATAC nuclei) at read time in
  `build_orientfree_isoform_data.py` / `cross_species_halfvshalf.py`. Only B:donor11a and B:donor2a
  `obs` change (other chick rows shift only in the permutation null, a shared RNG stream). The
  per-readout duck-vs-chick divergent fractions and Fisher p are unchanged (PITX2c 15/30 vs 2/19,
  p = 0.0045; RNA 25/31 vs 6/21, p = 2.3e-4; PITX2a 2/30 vs 0/19; PITX2b 1/30 vs 2/19).

### Panel C
- Stage consolidation: HH17/16 *n*=8 ⊆ HH25 *n*=19 (the L-biased cis-block only grows; motif content per peak is genomic/stage-invariant), so HH25 is the superset snapshot — the memory→wave-2 establishment is the *n*-per-stage, not a second grid of identical dots.
- JASPAR2020 carries no standalone SMAD1/SMAD9 matrix; the other BMP R-SMADs share the SBE, so the BMP arm is represented by SMAD5 + the SBE row.
- The aggregated HH25 L-biased cis-block = ventral L-sig peaks (`padj_locus<0.1`, −1…+12 kb of the c-TSS; *n*=17 focal windows on the 2026-08-11 corrected floor-20 substrate — the former L+14.6 / L+17.9 kb windows are now bilateral), with PITX2c & +5.4 kb broken out as their own columns.
- The within-locus bilateral background = 36 non-lateralized peaks within ±50 kb (the fig03 / panel-A focal scope); the two peaks overlapping the PITX2a/b promoters are excluded from the null since they are plotted as their own positive-control columns.
- Cross-species lift = Smith–Waterman onto galGal6a.
- The distal R-desert is lifted into chick against a wide window covering the region orthologous to duck −460 kb (around chick FAM241A): all 4 R-desert peaks are Class-I conserved (chick %id 0.82–0.98).

**Key reads (moved from the legend):** (i) the PITX2c promoter is a motif desert — empty for every motif including the GC-immune SBE, so the emptiness is not a GC artifact — despite being the strongly L-biased element, so its laterality is chromatin-gated, not motif-encoded; (ii) the PITX2a promoter is the SMAD5 / GC-immune-SBE BMP anchor, conserved in chick; (iii) the R-desert carries generic motifs but zero SMAD5 (not a BMP target); (iv) the L-block elements are Class-I conserved in chick, so duck laterality acts on conserved elements. The distal R-desert lifts into chick with all 4 R-desert peaks Class-I conserved → the distal region is syntenic, so the R-bias is a duck artifact at conserved elements (eRNA/FAM241A-filtering interpretation), not duck-specific sequence; like duck, the chick R-desert orthologs carry generic motifs but zero SMAD5.

## Methods notes & disclosures

- **Laterality is focal, not locus-wide.** Of the 359 TAD peaks, only a small contiguous set is
  L-biased per stage (ventral, **L-biased** at `padj_locus<0.1`, floor-20 unpaired
  `g4f_lr_long_floor20.tsv`: 8 at HH17/16, 11 at HH21, 22 at HH25 — the bidirectional totals are
  8 / 11 / 45, on a ventral HH17/16 census of 16 L / 11 R; the **dorsal control is 0 L / 0 R at every stage**). The story is **spatial
  contiguity** — a coherent L-biased block at the PITX2c promoter — **not** a raw L:R peak count.
- **Count-level vs per-peak significance.** The L:R *count* asymmetry is not by itself genome-wide
  significant by permutation (HH21 label-swap L−R excess `emp_p ≈ 0.21`, `h4_summary_HH21.tsv` — a
  genome-wide *count* test, distinct from panel B's per-embryo `emp_p`; at HH25 the asymmetric
  10 L / 15 R pseudobulk null limits the count-level permutation test, so the locus-FDR / paired /
  unpaired stack — not the permutation classifier — is the HH25 evidence). The **per-peak**
  genome-wide signal, however, *is* present, strongest at the **wave-2 (HH25)** stage: the PITX2c
  promoter itself reaches genome-wide DESeq2 FDR (HH25 ventral unpaired: genome-wide
  `padj = 1.2×10⁻¹¹`, locus `padj_locus = 3.1×10⁻¹⁴`, log2FC +1.28 — read directly from
  `g4f_lr_long_floor20.tsv`, written in `%.5g` scientific notation; no reconstruction needed —
  yet a motif desert), and **13 of the 22 HH25 L-biased locus peaks** clear genome-wide
  `padj < 0.1` (2 of 11 at HH21; 0 of 8 at HH17/16). Strong claims are kept to the
  **conserved Class-I** candidates.
- **Paired design is feasible at HH21/HH25** (the embryo-half paired g4f fit returns the same
  focal block there: 10 vs 11 L-sig at HH21, 16 L vs 22 at HH25); at **HH17/16 the paired fit
  returns 0 L-sig peaks** against the unpaired 8. The paired fit has 8 complete pairs (16 halves,
  7 residual df; `g4f_lr_grid_floor20.tsv`);
  the unpaired fit adds 11 single-sided halves (8 L / 3 R, the two split-tube left halves among them). The unpaired fit is shown as the
  primary call for *n*. (This bullet used to read "2 complete pairs / 1 residual df"; that is the
  cohort of the floor-50 trough fit, `g4h_cohorts.tsv`, not of this floor-20 g4f fit.)
- **The R-biased peaks are the distal desert** (−460/−527 kb, near FAM241A/Playrr) — the
  eRNA-filtering artifact region, treated throughout as a negative landmark.
- **Panel C is descriptive, not an enrichment test.** A curated PWM scan is the right tool at this
  ~20-peak locus scale (de novo / gkm-SVM / chromVAR need ~100+ genome-wide peaks); no enrichment
  *p*-value is computed on the selected set, the background is GC-matched, and the unbiased
  genome-wide gkm-SVM logograph (homeodomain / FOX, no SMAD) agrees. All motif/TF labels are factual
  (the logo is the unbiased sequence; the annotation names the DNA-binding family). The **L-block's
  SMAD5 rate (0.35, 6/17) is not significantly enriched over the bilateral background (0.25, 9/36;
  Fisher OR 1.64, *p* = 0.52)** — the signal is the c-prom desert + a-prom anchor, not L-block SMAD5
  enrichment.
- **No combinatorial enrichment either** (`combinatorial_motif_test.py`, not distributed; re-run 2026-08-18 on the
  seeded-relabel stageflags — the previous figures here were the pre-relabel L-block of 19 peaks
  against 34 background peaks). Pre-specified per-peak metrics — motif richness, total load, SMAD
  co-occurrence (≥2 of SMAD5/SMAD2-3-4/SBE), BMP+homeodomain pairing — show the L-block (n=17) is
  **indistinguishable from the GC/locale-matched ±50 kb background** (n=36; every MWU/Fisher p ≥ 0.69
  — the minimum is 0.696, on motif richness — co-occurrence ORs 0.93 and 1.07, i.e. no direction), and
  no longer even nominally **motif-sparser**
  than the full TAD (richness 2.53 vs 3.00, MWU p = 0.16, permutation p = 0.162; the pre-relabel run
  read 2.37 vs 3.00 at perm p ≈ 0.045, which was the one nominal signal here and does not survive).
  So neither individual motifs nor their combinations are enriched — consistent with **chromatin-gated**
  (not sequence-encoded) laterality. Underpowered (n=17) for subtle grammars: report as "no enrichment
  detected (individually or combinatorially)," leading with the locale-matched comparison.
- **±50 kb background validated TAD-wide.** A robustness check scanned the full TAD bilateral set
  (`g3f_tadbil_*`; 226 duck `tadbil` + 62 chick `tadbil_ck` windows, 14 window NAMES shared between the
  two classes → 274 distinct names in the file; the duck 226 are the background): for the BMP/SMAD motifs
  (SMAD5, SMAD2/3/4, SBE)
  the full-TAD background motif rate is close to the ±50 kb proximal one (SMAD5 0.34 vs 0.25; SBE 0.93 vs
  0.89), so the ±50 kb scope is representative. Rate rule for the full-TAD side: **window-collapsed** — a
  window counts as motif-positive if any of its constituent peak rows carries the motif (SMAD5 77/226 = 0.34;
  at peak-row level it is 78/246 = 0.32, because 20 duck windows carry two peaks each under one kb-rounded
  name). The ±50 kb side has one row per window, so its 9/36 = 0.25 is rule-invariant. Only NKX6-1 / ISL1 (AT-rich) rise TAD-wide in the distal lower-GC regions the ±50 kb
  window excludes. Chick subsample 51 Class-I / 11 Class-II (82 % conserved).

## Related docs

- `figures/data/chick_control/xspecies_orientfree_isoform.tsv` provenance and the duck-vs-chick orientation-free contrast live in the chick-control documentation set (`docs/chick_control/*.md`); the oriented counterpart is figS04 panel G.

## References

- Playrr — Welsh IC, Kwak H, Chen FL, Werner M, Shopland LS, Danko CG, Lis JT, Zhang M, Martin JF, Kurpios NA. *Chromatin Architecture of the Pitx2 Locus Requires CTCF- and Pitx2-Dependent Asymmetry that Mirrors Embryonic Gut Laterality.* Cell Reports 2015;13(2):337–349. doi:10.1016/j.celrep.2015.08.075. PMID 26411685. (candidate Playrr ortholog LOC106016802, panel A)
- PANCR — Gore-Panter SR, Hsu J, Barnard J, Moravec CS, Van Wagoner DR, Chung MK, Smith JD. *PANCR, the PITX2 Adjacent Noncoding RNA, Is Expressed in Human Left Atria and Regulates PITX2c Expression.* Circ Arrhythm Electrophysiol. 2016;9(1):e003197. doi:10.1161/CIRCEP.115.003197. PMID 26783232. (candidate PANCR ortholog LOC106017192, panel A)
