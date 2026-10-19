# references/

Reference tables read by the pipeline and figure code, and the published papers
whose results specific shipped analyses rely on.

## Reference tables

| File | Contents | Read by |
|---|---|---|
| `cluster_annotations.tsv` | Global-atlas cluster → cell-type annotation (written by `pipeline/d2b_validate_clusters.py`) | `figures/_shared/style.py` and the atlas figures |
| `d2g_subcluster_annotations.tsv` | Stage-restricted LPM subcluster → population (written by `pipeline/d2g_validate_subclusters.py`) | `figures/_shared/subcluster_map.{py,R}` and every laterality figure |
| `genotype_metadata.tsv` | Genotype → display name, embryo group, stage, demultiplexing source | `figures/_shared/style.py` |
| `genotype_overrides.tsv` | Per-genotype side/stage overrides with justification | `pipeline/d2o_merge_multiome_inference.py`, `figures/_shared/lr_resolution_panel.py` |
| `unified_donor_map.tsv` | Demultiplexing donor identity map (figure-side copy) | `figures/_shared/style.py` |
| `color_registry.tsv` | Colour per cell type / population / stage | `figures/_shared/colors.{py,R}` |
| `20220910_Wgenes.tsv`, `20220910_Zgenes.tsv` | W- and Z-chromosome gene IDs (ZJU1.0) | sex calling and pseudobulk PCA |
| `duck_s_genes.txt`, `duck_g2m_genes.txt` | S-phase and G2/M gene lists (duck orthologs) | cell-cycle scoring (`pipeline/d1_normalize_pca.py`, chick atlas) |
| `pitx2_promoter_regions.bed`, `pitx2_candidate_peaks.bed`, `pitx2_locus_50kb_features.tsv` | PITX2 isoform promoters, candidate peaks and locus features (ZJU1.0) | PITX2-locus pileups, region statistics and motif scans |
| `chick_pitx2_promoter_regions.bed` | Chick PITX2 isoform promoters | chick control pileups |
| `cfc1_locus_regions.bed` | CFC1 promoter and gene-body regions | CFC1 chromatin trajectory (`f4b`/`f5`/`f6` CFC1 steps) |

## Cited papers

- **Schweickert, A., Campione, M., Steinbeisser, H. and Blum, M. (2000).** Pitx2 isoforms: involvement of Pitx2c but not Pitx2a or Pitx2b in vertebrate left–right asymmetry. *Mechanisms of Development* 90: 41–51. DOI [10.1016/s0925-4773(99)00227-0](https://doi.org/10.1016/s0925-4773(99)00227-0). Relevance: only the PITX2c isoform is expected to be asymmetric, which is why the isoform-resolved pileups test PITX2c against PITX2a/b.
- **Stein, Y., Lavon, H., Hindi Malowany, M., Arpinati, L. and Scherz-Shouval, R. (2026).** The U-method: leveraging expression probability for robust biological marker detection. *bioRxiv*. DOI [10.64898/2026.03.31.715470](https://doi.org/10.64898/2026.03.31.715470). Relevance: the uniquely-expressed-gene stream of the subcluster marker-compendium annotation (`scripts/annotate_d2g_via_compendium.py`).
- **Woo, J., Miletich, I., Kim, B.-M., Sharpe, P.T. and Shivdasani, R.A. (2011).** Barx1-mediated inhibition of Wnt signaling in the mouse thoracic foregut controls tracheo-esophageal septation and epithelial differentiation. *PLoS ONE* 6: e22493. DOI [10.1371/journal.pone.0022493](https://doi.org/10.1371/journal.pone.0022493). Relevance: the BARX1+ tracheo-esophageal septum population label.
- **Johnsgard, P.A. (1961).** The tracheal anatomy of the Anatidae and its taxonomic significance. *Wildfowl* 12: 58–69. Relevance: the presence/absence map of the asymmetric syringeal bulla across Anseriformes used by the cross-species element analysis (`scripts/stagef/`).
- **Warner, R.W. (1971).** The structural basis of the organ of voice in the genera *Anas* and *Aythya* (Aves). *Journal of Zoology* 164: 197–207. DOI [10.1111/j.1469-7998.1971.tb01306.x](https://doi.org/10.1111/j.1469-7998.1971.tb01306.x). Relevance: supporting bulla morphology for the same trait map.
