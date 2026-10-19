# Figure S8: Library processing and demultiplexing validation

**Per-library and per-genotype quality metrics, clustering quality, and orthogonal validation of genotype demultiplexing for 52 duck embryos across 12 snRNA-seq libraries.**

**(A)** Per-library nuclei counts: uncorrected (dark) vs CellBender-corrected (light) paired bars.

**(B)** Per-library sequencing saturation.

**(C)** Mean reads per nucleus per library.

**(D)** Per-library mitochondrial-percentage violins (13 mitochondrially encoded genes), uncorrected (dark) vs CellBender-corrected (light), all nuclei before QC filtering; red dashed line = 5% QC threshold. Panels A–D show all 12 libraries; GEX0S-duck is a mixed duck–chick lane, so its values include chick nuclei.

**(E)** Per-library doublet rates. Circles: scDblFinder (Germain et al. 2022; transcriptome-based) on all QC-passing nuclei. Squares: vireo (Huang et al. 2019; genotype-based, homotypic-corrected) on the remaining scDblFinder singlets; arrows link the two steps. Doublets were removed using vireo (see Methods for the GEX03/GEX04 pre-screen). GEX0S-duck, a mixed duck–chick lane, has no vireo rate.

**(F-H)** Per-genotype QC on a stage-ordered axis of 52 genotypes with embryo-group brackets: **(F)** mitochondrial percentage (red dashed line = 5%), **(G)** ribosomal-gene percentage, **(H)** cell-cycle phase composition (stacked G1/S/G2M).

**(I)** PCA scree: per-PC and cumulative variance across the top 50 principal components.

**(J)** Batch integration: per-library silhouette width (lower = better mixing) for uncorrected vs Harmony-integrated (Korsunsky et al. 2019) embeddings; both sit near zero, indicating minimal residual library-batch structure, so the global atlas is clustered on uncorrected PCA (the stage-restricted LPM subclustering does use Harmony, batch = library).

**(K)** Per-cell-type mean silhouette width at the global Leiden resolution (res = 0.2).

**(L, M)** W- **(L)** and Z-chromosome **(M)** UMI fractions per genotype (blue = male ZZ, red = female ZW; dashed W-fraction threshold in L): unimodal distributions confirm each genotype is a single embryo, validating demultiplexing purity.

**(N)** Cross-group genotype distance heatmap: pairwise GT correlations across all genotype pools (top color bars = pool membership) confirm concordant calls between independent demultiplexing runs.

**(O)** Per-genotype summary on a shared stage-ordered axis. Left: sex-adjusted PC1 (PC1 of the mesenchymal pseudobulk PCA after regressing out sex) — the developmental-stage axis (cf. Figure S7d) used to resolve the age-ambiguous design pools. Right: per-genotype nucleus counts colored by inferred stage (HH13-HH25); counts range 750 to 8,400.

**(P)** Demultiplexing pipeline. Three sequencing pools were demultiplexed independently with pool-appropriate variant strategies: the main pool (GEX05-08, PG3, PG4A, PG4B; 7 libraries sharing one donor pool) — iterative sub-pool vireo at 665K universal SNP positions, resolving 33 embryos; GEX03 (3 libraries) — pool-specific variant discovery preserving Mallard/Khaki Campbell breed-discriminative positions, 10 embryos; GEX04 + GEX0S-duck (2 libraries, mixed with chick cells) — de novo variant discovery (cellSNP mode 2b, 173K pool-specific positions), 9 duck embryos. Cross-pool consolidation (pseudobulk AF correlation, barcode concordance, suffix matching) yields 52 confirmed embryos.

**(Q)** GEX03 — the study's only pool mixing two breeds — resolves into Khaki Campbell (KC) and Mallard (Md) clusters: hierarchical-clustering heatmap (average linkage on 1 − r) of pairwise Pearson r between vireo-reconstructed donor genotype profiles (12 donors; sidebars and tick labels: KC blue, Md red; greyscale r). Within-KC r ≈ 0.31-0.40 and within-Md r ≈ 0.18-0.33 (excluding one r = 0.77 pair — genotype_06 + genotype_12, two independent vireo calls of the same embryo, 17-Md-13, consolidated downstream) versus cross-breed r ≈ 0.04-0.16. Breed assignments (4 KC + 8 Md) drive the variant strategy of (P).
