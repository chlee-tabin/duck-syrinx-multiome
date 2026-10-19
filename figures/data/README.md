# figures/data/

Tables the figure scripts read. Every figure panel renders from these files
without the sequencing data (Tier 1 in the top-level README). They were written
by the pipeline steps named in the **Produced by** column; re-running those
steps from the GEO deposit (Tier 2) regenerates them.

- `PIPELINE_STATE.json` records the pipeline run that produced the cluster tables
  (cell and genotype counts).
- `pipeline/derived_data_manifest.tsv` lists every cluster-produced table with its
  destination here and whether it depends on the master h5ad (`depends_on_h5ad`).
  When the master h5ad changes, every `depends_on_h5ad=yes` table must be
  regenerated together. `scripts/verify_derived_data.py`
  (`make check-derived-data`) cross-checks those rows against
  `references/cluster_annotations.tsv`.
- File names follow `<pipeline_step>_<description>`, without figure numbers.

**Read by** names the figure directories (`figures/<fig>/`) whose scripts read the
file, then any downstream pipeline steps. "not read by shipped code" marks
end-point results (statistics, sensitivity and robustness analyses) that no figure
script or later step reads. `HH*` and `*` in a file name stand for the stage or
library variants.

## Multiome (duck)

| File | Produced by | Read by | Description |
|---|---|---|---|
| `d2_umap_boundary_data.tsv.gz` | `pipeline/d2_cluster_markers.py` | fig01 | UMAP coordinates + Leiden cluster labels (`leiden_res0.2` etc.) for all cells |
| `e1_pca_mesenchymal.tsv` | `pipeline/e1_pseudobulk_pca.py` | fig01, `pipeline/build_genotype_embryo_table.py`, `pipeline/e1b_stage_inference.py` | Pseudobulk PCA for mesenchymal subset (res=0.2) |
| `e1_pca_neural.tsv` | `pipeline/e1_pseudobulk_pca.py` | not read by shipped code | Pseudobulk PCA for neural subset |
| `e1b_stage_assignments.tsv` | `pipeline/e1b_stage_inference.py` | fig01, fig01_supp, fig02, fig02_supp, `pipeline/_locus_chromatin_metrics.py`, `pipeline/_locus_region_counts.py`, `pipeline/d2g_signaling_export.py`, `pipeline/d2g_stage_restricted_subclustering.py`, … |  |
| `e2_pitx2_laterality.tsv` | `pipeline/e2_laterality_pitx2.py` | fig02_supp, `pipeline/build_genotype_embryo_table.py` | PITX2 expression per genotype × side × cell type |
| `e3_composition.tsv` | `pipeline/e3_lr_composition.py` | fig02_supp, `pipeline/validate_same_side_duplicates.py` | Cell type composition per genotype × library; the validator recomputes the same-side duplicate evidence from it and raises unless it matches `references/same_side_duplicate_genotypes.tsv` |
| `e3_cellcycle.tsv` | `pipeline/e3_lr_composition.py` | not read by shipped code | Cell cycle phase per genotype × library |
| `step03_sex_per_genotype.tsv` | `pipeline/c3_sex_determination.py` | _shared, fig01_supp, figS07, figS16, figS17, `FC_08656/20260329_step07c_pseudobulk_pca_subsets.py`, `FC_08656/20260330_step02f_clean_genotypes.py`, `pipeline/d1_normalize_pca.py`, `pipeline/e1_pseudobulk_pca.py`, … | Sex calls per 52 genotypes |
| `step03b_elbow.tsv` | `pipeline/d1_normalize_pca.py` | fig01_supp | PCA variance ratios (50 PCs) |
| `step03b_batch_qc.tsv` | `pipeline/d1_normalize_pca.py` | fig01_supp | Batch-mixing QC (uncorrected vs Harmony) |
| `step03_sex_per_cell.tsv.gz` | `pipeline/c3_sex_determination.py` | fig01_supp, `FC_08656/20260330_step02f_clean_genotypes.py`, `pipeline/run_pipeline.sh` | Per-cell W/Z sex evidence |
| `step02e_filter_summary.tsv` | `pipeline/c1_qc_filter.py` | `FC_08656/20260326_step02e_preCB_genotype.py` | QC-filter cascade totals (181,159 vireo singlets) |
| `step04a_cluster_counts.tsv` | `pipeline/d2_cluster_markers.py` | not read by shipped code | Per-cluster nucleus counts (res 0.2) |
| `step04b_per_genotype_metrics.tsv.gz` | `pipeline/d2b_export_genotype_metrics.py` | fig01_supp | Per-genotype QC metrics (52 genotypes) |
| `step01_cellranger_metrics.tsv` | `pipeline/d2c_parse_cellranger_metrics.py` | fig01_supp | Per-library cellranger GEX summary metrics, 12 libraries. The GEX04 and GEX0S-duck rows were transcribed from their cellranger `metrics_summary.csv`; `reads_mapped_transcriptome` is blank in all rows |
| `dotplot_res02.tsv` | generate_dotplot_export | fig01, `FC_08656/20260330_generate_dotplot.py`, `FC_08656/run_dotplot_export.sh` | 16 cell types × 37 marker fractions/means |
| `step01_library_summary.tsv`, `step01_cb_library_summary.tsv`, `step02_mito_by_library.tsv.gz`, `step02_cb_mito_by_library.tsv.gz` (4 files) | `pipeline/s8ad_regen_libqc_tables.py` (+ `.sh`) | fig01_supp | Per-library nuclei counts and QC summary (uncorrected cellranger / CellBender) and per-nucleus pct_mito, n_genes, n_counts, all 12 libraries; 13-gene duck mitochondrial list in both arms |
| `step02_doublet_comparison.tsv`, `step02b_doublet_summary.tsv` (2 files) | c1 QC exports (producer not distributed) | `FC_08656/20260320_step02b_scdblfinder.R` | Doublet summaries; kept as records. No figure reads them (fig01_supp panel E reads `step02e_doublet_comparison.tsv`) |
| `gt_distance_combined_corr.tsv` | `FC_08656/genotype_distance_heatmap.py` | fig01_supp | Combined 86×86 genotype-profile correlations, all pools |
| `genotype_distance_g2_pairwise_corr.tsv` | `FC_08656/genotype_distance_heatmap.py` | fig01_supp | G2-pool 12×12 genotype-profile correlations |
| `step03b_horn_analysis.tsv` | `pipeline/d1b_horn_pc_analysis.py` | not read by shipped code | Per-PC observed vs permuted eigenvalues (Horn's) |
| `step03b_horn_summary.tsv` | `pipeline/d1b_horn_pc_analysis.py` | not read by shipped code |  |
| `step03b_horn_umap_data.tsv.gz` | `pipeline/d1b_horn_pc_analysis.py` | not read by shipped code | UMAP coords + Leiden for baseline/Horn at 5 resolutions |
| `step03b_resolution_sweep.tsv` | `pipeline/d1b_horn_pc_analysis.py` | not read by shipped code | Cluster counts per (method, resolution) |
| `step03b_annotation_sweep.tsv` | `pipeline/d2e_annotate_resolution_sweep.py` | not read by shipped code | Auto-annotated clusters per (method, resolution) |
| `d2e_silhouette_global_res02.tsv` | `pipeline/d2e_silhouette_export.py` | fig01_supp, `pipeline/run_pipeline.sh` | Per-cluster silhouette at res=0.2 (all cells) |
| `d2e_silhouette_focused_mes.tsv` | `pipeline/d2e_silhouette_export.py` | `pipeline/run_pipeline.sh` | Per-subcluster silhouette for focused LPM |
| `d2j_pitx2_isoform_per_genotype.tsv` | `pipeline/d2j_pitx2_isoform_per_genotype.py` | `pipeline/b1_donor_consolidation.py`, `pipeline/d2j_pitx2_isoform_percell.py`, `pipeline/d2k_pitx2c_laterality_qc.py` |  |
| `d2k_laterality_qc.tsv` | `pipeline/d2k_pitx2c_laterality_qc.py` | not read by shipped code | PITX2c laterality validation vs L/R design |
| `d2g_lpm_subclustering_HH*.tsv.gz` (8 files) | `pipeline/d2g_stage_restricted_subclustering.py` | fig02, fig02_supp, figS07, figS17, `pipeline/chick/build_isoform_data.py`, `pipeline/chick/build_orientfree_isoform_data.py`, `pipeline/chick/per_stage_laterality.py`, `pipeline/d2l_lr_inference.py`, … | Per-cell markers + leiden clusters for stage-restricted LPM subclustering |
| `d2g_signaling_HH*.tsv.gz` (4 files) | `pipeline/d2g_signaling_export.py` | fig02_supp, figS07, figS17 |  |
| `d2j_pitx2_isoform_percell.tsv.gz` | `pipeline/d2j_pitx2_isoform_percell.py` | fig02, fig02_supp, figS03, `pipeline/d2l_lr_inference.py`, `scripts/annotate_d2g_via_compendium.py`, `scripts/check_biological_concordance.py` |  |
| `step_f4_pitx2_coverage.tsv.gz` | `pipeline/f4_lr_pitx2_pileup.py` | `pipeline/run_f4_hh1716_only.sh`, `scripts/merge_f4_hh1716_into_canonical.py` | Per-base PITX2 ± 50 kb ATAC coverage (aggregate L vs R per stage) |
| `step_f4_pitx2_coverage_per_genotype.tsv.gz` | f4_lr_pitx2_pileup_per_genotype | not read by shipped code | Per-base PITX2 coverage at per-genotype granularity (HH17/16 paired tracks at MIN_CELLS=30) |
| `step_f5_pitx2_region_counts.tsv` | `pipeline/f5_pitx2_region_counts.py` | fig03, figS03 | Per-(genotype × library × stage × side × region) PITX2 ATAC fragment counts at promoter / gene body / TSS regions |
| `step_f5_cfc1_region_counts.tsv` | `pipeline/f5_cfc1_region_counts.py`; HH17/16 sides then relabelled onto the corrected design (same-side duplicates pooled into one `GEX04+GEX0S-duck` sample) by a one-off script (not distributed) | figS03, `pipeline/g3_cfc1_deseq2.R` | Per-(genotype × library × stage × side × region) CFC1 ATAC fragment counts at promoter (TSS ± 300 bp) and gene body |
| `step_f5_cfc1_region_deseq2.tsv` | `pipeline/g3_cfc1_deseq2.R` | figS03 | DESeq2 L vs R per (region × stage × design) |
| `step_f4b_cfc1_coverage.tsv.gz` | `pipeline/f4b_cfc1_locus_coverage.py` | not read by shipped code | Per-base CFC1 ± 50 kb ATAC coverage (aggregate L vs R per stage) |
| `step_f4b_cfc1_cell_counts.tsv` | `pipeline/f4b_cfc1_locus_coverage.py` | not read by shipped code | Per-(group) cell-count provenance for f4b coverage normalization |
| `f6_chromatin_metrics.tsv` | `pipeline/f6_chromatin_metrics.py` | `pipeline/d2n_chromatin_lr_inference_hh1716.py`, `pipeline/d2o_merge_multiome_inference.py`, `pipeline/f6_cfc1_chromatin_metrics.py` |  |
| `f6_cfc1_chromatin_metrics.tsv` | `pipeline/f6_cfc1_chromatin_metrics.py` | not read by shipped code | Per-(genotype × library × stage) CFC1 promoter / gene-body open-fraction |
| `h3_gkmsvm_logos_strict_motifs.tsv` | `pipeline/h3_gkmsvm_logos.py` | not read by shipped code |  |
| `h3_gkmsvm_logos_strict_known.tsv` | `pipeline/h3_gkmsvm_logos.py` | not read by shipped code |  |
| `h3_gkmsvm_logos_strict_pitx2_contrib.tsv` | `pipeline/h3_gkmsvm_logos.py` | not read by shipped code | Per-PITX2-candidate-peak motif-hit table (top 10 de novo motifs × 8 candidate peaks) |
| `h3_gkmsvm_logos_strict_summary.tsv` | `pipeline/h3_gkmsvm_logos.py` | not read by shipped code | One-row run summary: tool, n_pos, n_neg, motif_lengths, n_motifs_discovered, real_top1_log10_p, null_top1_log10_p_len8, elapsed_sec |
| `step01_cellranger_arc_metrics.tsv` | `pipeline/d2c2_parse_cellranger_arc_metrics.py` | figS03 |  |
| `step01_cellranger_arc_per_cell_sample.tsv` | `pipeline/d2c2_sample_arc_per_cell_metrics.py` | figS03, `pipeline/run_f2e_array.sh` | Dual-reservoir-sampled (seed 42) 200 cells per (ARC library × qc_set) × 12 libs = up to 4800 rows |
| `step_f1_passed_barcodes.tsv` | `pipeline/d2c3_dump_f1_passed_barcodes.py` | `pipeline/d2c2_sample_arc_per_cell_metrics.py`, `pipeline/run_d2c3_dump_f1_passed_barcodes.sh` | Per-(arc_library, barcode) list of f1-passing cells (TSS ≥ 2 & fragments ≥ 1000) dumped from each library's f1 backed h5ad (snapatac2 obs_names) |
| `step_f1_qc_summary.tsv` | `pipeline/f1_coupled_atac.py` | figS03 |  |
| `step_f1c_amulet_summary.tsv` | `pipeline/f1c_amulet_doublets.py` | figS03 |  |
| `hh25_percell_counts_esr1.tsv.gz` | `scripts/export_hh25_percell_counts.py` | figS17 |  |
| `b_genotype_embryo_mapping.tsv` | `pipeline/build_genotype_embryo_table.py` | fig02, `pipeline/d2k_pitx2c_laterality_qc.py`, `pipeline/d2l_lr_inference.py`, `pipeline/d2n_chromatin_lr_inference_hh1716.py`, `scripts/annotate_d2g_via_compendium.py`, … |  |
| `d2l_lr_inference.tsv` | `pipeline/d2l_lr_inference.py` | _shared, fig02_supp, `pipeline/chick/build_isoform_data.py`, `pipeline/chick/build_orientfree_isoform_data.py`, `pipeline/chick/duck_pitx2_locus_percell_extract.py`, `pipeline/chick/per_stage_laterality.py`, … | Per-(genotype × library) 4-metric PITX2 inference + library-pair placement |
| `archr_peak_annotation.tsv.gz` | `scripts/annotate_archr_peaks.py` | figS03 |  |
| `d2l_cluster_biological_concordance.tsv` | `scripts/check_biological_concordance.py` | not read by shipped code |  |
| `d2m_genotype_to_embryo.tsv` | `pipeline/d2m_genotype_to_embryo.py` | _shared, fig02_supp | Per-genotype candidate embryo set from library + stage + d2l intersection |
| `h9_scope_arms_HH21.tsv` | `pipeline/run_h9_scope_arms.sh` | not read by shipped code |  |
| `g4e_focal_pitx2_g2c.tsv.gz`, `g4e_focal_pitx2_peaks.txt` | `pipeline/g4e_focal_pitx2_peakset.py` | `pipeline/run_rightpeak_smad_probe.sh` | Right-peak BMP-SMAD probe, Arm 2 (focal) |
| `rightpeak_smad/{dynlost,foc}_{floor,genomic,peak}__{h3b_kmer_top,h3b_summary,h9_summary}*.tsv` (18 files), `rightpeak_smad/g4d_dorsal_gained_HH17_16_to_HH21_{gained_peaks.txt,lost_peaks.txt,summary.tsv}`, `rightpeak_smad/g4d_{gained,lost}_annotated.tsv` (23 files) | `scripts/annotate_g4d_dynamic_peaks.py`, `pipeline/run_rightpeak_smad_probe.sh`, `pipeline/run_rightpeak_dynamic_lost.sh` | not read by shipped code | Right-peak BMP-SMAD probe results |
| `d2g_dotplot_spatial.tsv` | `figures/fig02/generate_dotplot_spatial.py` | not read by shipped code | Dotplot stats: 13 genes × 20 rows (intermediate) |
| `d2g_dotplot_spatial_extended.tsv` | `figures/fig02_supp/generate_dotplot_extended.py` | fig02 |  |
| `d2g_dotplot_spatial_extended_submitted.tsv`, `d2g_subcluster_annotations_submitted.tsv` (2 files) | frozen copies of `d2g_dotplot_spatial_extended.tsv` and `references/d2g_subcluster_annotations.tsv` as used for the accepted figure (not regenerated) | fig02 | Pin the fig02 schematic HH21 L/R callouts to the values in the published figure; `DUCK_SPATIAL_LIVE_DATA=1` makes `generate_spatial_svg.py` read the live tables instead |
| `grn_memory_node_promoter_scan.tsv` | grn_memnode_scan, `pipeline/run_grn_memnode_scan.sh`, `figures/figS03/grn_memory_node_scan.py` | figS03, `pipeline/run_hh21_module_v2.sh` |  |
| `hh21_hh1716_module_scores.tsv`, `hh21_hh1716_module_halfscores.tsv`, `hh21_hh1716_module_stats.tsv` | hh21_module_v2, `pipeline/run_hh21_module_v2.sh`, `figures/figS03/hh21_hh1716_module_test.py` | figS03 |  |
| `hh21_hh1716_module_labelswap_stats.tsv` | `figures/figS03/panel_hh21_hh1716_module.py` | figS03 | The formal within-embryo label-swap permutation p per (stage × module) |
| `hh21_hh1716_response_counts.tsv`, `hh21_hh1716_halfmeta.tsv` | hh21_module_v2, `pipeline/run_hh21_module_v2.sh`, `figures/figS03/hh21_hh1716_module_test.py` | figS03 |  |
| `hh21_hh1716_nb_glmm.tsv` | `figures/figS03/hh21_hh1716_nb_glmm.R` | figS03, `pipeline/run_hh21_module_v2.sh` | NB mixed-model arm test: per (stage × engine × offset × term) side/arm effect in log2 with 95% CI + one-sided p |
| `hh21_trans_promoter_map.tsv`, `hh21_trans_promoter_motifs.tsv`, `hh21_trans_motif_panel.tsv` | `figures/figS03/hh21_trans_motif.py`, `pipeline/run_hh21_trans_motif.sh`, `figures/figS03/hh21_trans_motif_scan.R` | not read by shipped code |  |
| `hh21_trans_motif_stats.tsv`, `hh21_trans_motif_groups.tsv` | `figures/figS03/hh21_trans_motif.py` | not read by shipped code |  |
| `hh21_seq_promoter_map.tsv`, `hh21_seq_promoter_motifs.tsv`, `hh21_seq1_promoter_motifs.tsv` | `figures/figS03/hh21_trans_motif_seq.py`, `pipeline/run_hh21_trans_motif_seq.sh`, `figures/figS03/hh21_trans_motif_scan.R` | figS03, `pipeline/run_hh21_smad234_dinuc_null.sh` |  |
| `hh21_seq{,1}_motif_stats.tsv`, `hh21_seq{,1}_gene_table.tsv`, `hh21_seq{,1}_da_motif_corr.tsv` (6 files) | `figures/figS03/hh21_trans_motif_seq.py` | not read by shipped code |  |
| `hh21_seq1_partition_stats.tsv` | `figures/figS03/hh21_trans_motif_seq.py` | not read by shipped code |  |
| `hh21_smad234_dinuc_null.tsv`, `hh21_smad234_trinuc_null.tsv` | `figures/figS03/hh21_smad234_dinuc_null.R`, `pipeline/run_hh21_smad234_dinuc_null.sh` | figS03 |  |
| `hh21_{dinuc,trinuc}_null_{SMAD2_3_4,PITX2}_stats.tsv` (4 files) | `figures/figS03/hh21_trans_motif_seq.py` | not read by shipped code | Null distribution of the L-up coefficient under k-let-preserving shuffling + the observed value + the empirical p (floored at 1/(K+1)). |
| `d2j_pitx2_isoform_per_genotype_nostrand.tsv` | `pipeline/d2j_pitx2_isoform_per_genotype.py`, `FC_08656/run_d2j_baseline_nostrand.sh` | not read by shipped code | PITX2 isoform UMI without the sense-strand filter |
| `d2j_pitx2_isoform_percell_nostrand.tsv.gz` | `pipeline/d2j_pitx2_isoform_percell.py`, `FC_08656/run_d2j_baseline_nostrand.sh` | not read by shipped code | Per-cell PITX2 isoform UMI without the sense-strand filter |
| `pitx2_playrr_synteny_blocks.tsv` | LASTZ duck-chick alignment (producer not distributed) | fig02_supp | Static reference data |
| `d2_umap_boundary_only.tsv.gz` | `pipeline/d2_cluster_markers.py` | not read by shipped code |  |
| `d2m_multiome_genotype_to_embryo.tsv` | `pipeline/d2m_multiome_genotype_to_embryo.py` | figS03 | Chromatin-canonical genotype→embryo narrowing |
| `d2n_chromatin_lr_inference_hh1716.tsv` | `pipeline/d2n_chromatin_lr_inference_hh1716.py` | `pipeline/chick/per_stage_laterality.py`, `pipeline/d2o_merge_multiome_inference.py` |  |
| `lr_inference_multiome.tsv` | `pipeline/d2o_merge_multiome_inference.py` | _shared, fig03, figS03, figS07, figS16, figS17, `pipeline/chick/per_stage_laterality.py`, `pipeline/d2m_multiome_genotype_to_embryo.py`, `pipeline/f4_lr_pitx2_pileup.py`, `pipeline/f5_cfc1_region_counts.py`, … |  |
| `step_f5_pitx2_region_deseq2_sweep.tsv` | `scripts/f5_threshold_sweep.sh`, `figures/figS03/beeswarm_deseq2_stats.R` | not read by shipped code |  |
| `step_f5_pitx2_region_deseq2.tsv` | `figures/figS03/beeswarm_deseq2_stats.R`, `pipeline/run_beeswarm_deseq2.sh` | figS03, `scripts/f5_threshold_sweep.sh` |  |
| `hh21_power_calibration.tsv` | `figures/figS03/hh21_power_calibration.py` | figS03 |  |
| `hh21_deg_da_concordance.tsv`, `hh21_deg_da_concordance_stats.tsv` | `figures/figS03/hh21_positive_control.py` | not read by shipped code |  |
| `step_f5_pitx2_region_stats.tsv` | `figures/figS03/region_stats_wilcoxon.R` | _shared |  |
| `step_f5_cfc1_region_stats.tsv` | `figures/figS03/region_stats_wilcoxon.R` | not read by shipped code |  |
| `validation_chromatin_summary.tsv` | `pipeline/d2m_multiome_genotype_to_embryo.py` | not read by shipped code |  |
| `validation_d2m_multiome_vs_d2m_hh1716.tsv` | `pipeline/d2m_multiome_genotype_to_embryo.py` | not read by shipped code | Per-genotype validation table comparing the chromatin-canonical multiome assignment to d2m's snRNA-only HH17/16 assignment |
| `validation_d2n_vs_d2l_hh1716.tsv` | `pipeline/d2n_chromatin_lr_inference_hh1716.py` | not read by shipped code | Per-(genotype × library) comparison of d2n chromatin call to d2l snRNA call at HH17/16 |
| `step_f4_pitx2_coverage_per_genotype_mincells100.tsv.gz` | `pipeline/f4_lr_pitx2_pileup.py` | not read by shipped code |  |
| `f4_compartment_perembryo_pitx2_coverage.tsv.gz` | `pipeline/f4_lr_pitx2_pileup.py`, `pipeline/run_f4_perembryo.sh` | fig03 | Per-base per-cell PITX2-focal-locus ATAC coverage split by (stage × compartment × genotype × side), one column per embryo-half, on the current `lr_inference_multiome.tsv` sides; a split-tube half (17-5/17-6) is one column with both aliquots summed, as in g1, f5 and h1 |
| `f4_compartment_perembryo_pitx2_cell_counts.tsv` | written with the f4 export above (`pipeline/f4_lr_pitx2_pileup.py`, `pipeline/run_f4_perembryo.sh`) | fig03 | Nuclei per f4 coverage column (`group`, `n_cells`, `coverage_sum_pitx2_locus`, `coverage_sum_per_cell`); fig03 panel C checks f4 nuclei ≥ floor-20 census nuclei for every census half |
| `f4_hh13_perembryo_pitx2_coverage.tsv.gz` | `pipeline/f4_lr_pitx2_pileup.py`, `pipeline/run_f4_perembryo.sh` | fig03 |  |
| `h0_lvr_peaks_pitx2_locus_locus_fdr.tsv` | `scripts/compute_pitx2_locus_fdr.py` | not read by shipped code |  |
| `g4f_{census,cpm_trajectory_long,lr_long,lr_grid,stage_long,stage_grid,sample_meta}.tsv`, `g4f_locus_count_slice.tsv` (8 files) | `pipeline/g4f_pitx2_locus_dynamics.py` | `pipeline/g4h_prep_trough_lr.py`, `pipeline/run_atac_relabel_campaign.sh`, `pipeline/run_hh25_relabel_chain.sh`, `pipeline/run_hh25_relabel_wrap.sh`, … | Depth-controlled per-(compartment × stage) PITX2-locus dynamics over the 359-peak TAD set |
| `g4f_dv_{census,cpm_trajectory,long,grid}{,_floor20}.tsv` (8 files) | `pipeline/g4g_pitx2_locus_dv.py`, `pipeline/run_g4g_dv.sh` | not read by shipped code | Dorsal-vs-Ventral (D/V) companion to g4f's L/R laterality |
| `g4f_depthmatch.tsv` | `scripts/depthmatch_pitx2_locus.py` | figS03 | Depth-MATCHED L/R log2(meanCPM) per NAMED locus element × ventral stage |
| `pitx2_locus_candidate_peaks.tsv` | `scripts/pitx2_locus_candidate_peaks.py` | not read by shipped code |  |
| `pitx2_locus_tad_peaks.tsv` | `scripts/derive_pitx2_locus_tad_peaks.py` | fig03, fig04, `pipeline/g4e_focal_pitx2_peakset.py`, `pipeline/g4f_pitx2_locus_dynamics.py`, `pipeline/g4g_pitx2_locus_dv.py`, `pipeline/run_atac_relabel_campaign.sh`, … |  |
| `pitx2_locus_gc_content.tsv` | `scripts/gc_content_pitx2_window.py`, `scripts/run_gc_content_pitx2.sh` | fig03 |  |
| `pitx2_locus_gc_content_tad.tsv` | `scripts/run_gc_content_pitx2.sh` | not read by shipped code |  |
| `g3f_domain_duck_motifs.tsv`, `g3f_domain_window_gc.tsv`, `g3g_domain_{chick_motifs,lift}.tsv`, `g4f_domain_window_stageflags.tsv` (5 files) | `scripts/make_g4f_domain_windows.py`, `pipeline/run_g3f_g3g_domain.sh`, `pipeline/g3f_hh21_motif_scan.R`, `pipeline/g3g_hh21_chick_lift_scan.R` | fig04 |  |
| `g3f_tadbil_duck_motifs.tsv`, `g3f_tadbil_window_gc.tsv`, `g3g_tadbil_{chick_motifs,lift}.tsv` (4 files) | `pipeline/g3f_hh21_motif_scan.R` | not read by shipped code |  |
| `g3h_bmp7_motif_windows.bed` | `scripts/make_bmp7_motif_windows.py`, `pipeline/g3g_hh21_chick_lift_scan.R`, `pipeline/run_g3h_bmp7_chick.sh` | not read by shipped code | Robustness check for the Table S3 BMP7 claim (not plotted, not a manuscript input) |
| `h4_peak_empirical_p_HH{13,17_16,21,25}.tsv.gz` (4 files) | `pipeline/h4_permutation_test.py` | figS03, `pipeline/run_hh21_trans_motif.sh`, `pipeline/run_hh21_trans_motif_seq.sh` | Per-peak empirical p from 1,000 side-shuffle permutation null at the named stage × INFERENCE_POOL pseudobulks |
| `h4_summary_HH{13,17_16,21,25}.tsv` (4 files) | `pipeline/h4_permutation_test.py` | `pipeline/run_h4_HH21.sh`, `pipeline/run_h4_HH25.sh`, `pipeline/run_h4_floor30.sh` | One-row run summary per stage: observed L/R-biased peak counts at raw p<0.05, shuffle-null medians/maxima, count-level empirical p ((b+1)/(N+1), min 1/1001) |
| `h2_deg_summary.tsv`, `h2_deg_{POOL,SUB_*}_*.tsv` (13 files) | `pipeline/h2_deseq2_snrna.R` | figS03, figS17, `pipeline/g4h_response_locus_trough_test.py`, `pipeline/g5_hh21_wave2_anticipation.py`, `pipeline/run_g4h.sh`, `pipeline/run_hh25_relabel_rna_panelB.sh`, … | Per-(stage × scope) L-vs-R RNA DEG counts at padj 0.05/0.01/0.001 + per-side sample n |
| `h2c_hh1716_trough_composition.tsv`, `h2c_hh1716_trough_composition_summary.tsv` | `pipeline/h2c_hh1716_trough_composition.py` | fig02_supp | HH17/16 inference-pool left/right composition decomposition (Cardiopharyngeal LPM vs BARX1+) of the AHR / CHRM2 right bias, per paired embryo, and its summary statistics (panel F) |
| `h1_dorsal_pseudobulk_{counts.tsv.gz,metadata.tsv}`, `h2_dorsal_{summary,POOL_*,SUB_*}.tsv` (9 files) | `pipeline/h1_dorsal_pseudobulk.py`, `pipeline/h2_deseq2_snrna.R`, `pipeline/run_pipeline.sh` | fig03, figS03, `pipeline/run_hh25_relabel_rna_panelB.sh`, `pipeline/run_hh25_relabel_rna_panelB_h2.sh` |  |
| `h3_gene_detection.tsv` | `pipeline/h3_gene_detection_export.py`, `pipeline/run_pipeline.sh` | not read by shipped code |  |
| `h8_{design_comparison,summary,subcluster_homogeneity}_HH{17_16,21,25}.tsv` (9 files) | `pipeline/h8_laterality_design_comparison.py`, `pipeline/run_h8_design_comparison.sh` | not read by shipped code | Settles the manuscript laterality convention: per stage, fits the SAME genome-wide ATAC peakset under 4 designs |
| `h6_paired_deseq2_HH{17_16,25}_balanced{,_summary,_locus}.tsv` (6 files) | `pipeline/h6_paired_deseq2.py`, `scripts/extract_h6_locus_subset.py`, `pipeline/run_h6_locus_tad_extract.sh` | figS16, `pipeline/g4h_response_locus_trough_test.py`, `pipeline/h7_unpaired_deseq2.py`, `pipeline/run_atac_relabel_campaign.sh`, `pipeline/run_h6_HH25.sh`, … |  |
| `h7_unpaired_deseq2_HH{13,17_16,21,25}{,_summary}.tsv{.gz,}` (8 files) | `pipeline/h7_unpaired_deseq2.py` | figS03, `pipeline/run_atac_relabel_campaign.sh`, `pipeline/run_h6_locus_gzip_d2n.sh`, `pipeline/run_hh25_relabel_h67.sh`, `scripts/compute_locus_fdr_h6_h7.py`, … |  |
| `h6_paired_deseq2_{HH17_16,HH25}_balanced_locus_fdr.tsv`, `h7_unpaired_deseq2_{HH13,HH17_16,HH21,HH25}_locus_fdr.tsv`, `locus_fdr_summary_h6_h7.tsv`, `locus_fdr_candidates_h6_h7.tsv` (8 files) | `scripts/compute_locus_fdr_h6_h7.py` | `scripts/xspecies_build_lbiased_set.py` |  |
| `xspecies_tad_candidate_lr.tsv` | `scripts/xspecies_tad_candidate_lr.py` | not read by shipped code |  |
| `xspecies_duck_lbiased_set_per_stage.tsv` | `scripts/xspecies_build_lbiased_set.py` | not read by shipped code |  |
| `xspecies_tad_candidate_lr_ventral.tsv` | `scripts/xspecies_ventral_subcluster_lr.py` | not read by shipped code |  |
| `xspecies_clade_conserved_divergent_elements.tsv` | `scripts/stagef/read_out.py` | not read by shipped code |  |
| `xspecies_duck_chick_2class.tsv` | `scripts/stagef/classify_2class.py` | not read by shipped code |  |
| `xspecies_tad_clade_standout.tsv` | `scripts/stagef/tad_standout.py` | not read by shipped code |  |
| `g2c_aggregated_HH{13,1716,21,25}_INFERENCE_POOL_archr.tsv`, `g2c_aggregated_summary_{archr,snapatac}.tsv` (6 files; the three HH17/16 `SUB_*` tables have their own rows below) | `pipeline/g2c_aggregated_descriptive.py` | `pipeline/g4c_motif_on_aggregated.sh`, `pipeline/h2_gkmsvm_train.py`, `pipeline/h3_jaspar_scan.R`, `scripts/compute_pitx2_locus_fdr.py`, … | Aggregated descriptive L vs R per-peak binomial test on summed pseudobulk counts per stage × pool |
| `g3_pitx2_motif_scan_p{5e-4,5e-3,1e-2}.tsv`, `g3b_candidate_peaks_motif_scan_p{5e-4,5e-3,1e-2}.tsv`, `g3b_candidate_peaks_motif_family_table.tsv`, `g3b_peaks_motif_family_{unified,wide}.tsv` (9 files) | `pipeline/g3_pitx2_motif_scan.R`, `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code | g3: motif scan over PITX2 promoter BED regions at 3 cutoffs |
| `g4_hh1716_motif_enrichment.tsv`, `g4_hh1716_motif_enrichment_aggregated_archr.tsv`, `g4_hh1716_peak_motif_matrix_archr.tsv.gz` | `pipeline/g4_hh1716_motif_enrichment.R` | figS16, `scripts/locate_lvr_peaks_pitx2_locus.py` | Motif enrichment Fisher tests on HH17/16 INFERENCE_POOL L-biased peaks (motifmatchr against JASPAR2020-vertebrate) |
| `g5_hh21_atac_anticipation_summary.tsv` | `pipeline/g5_hh21_wave2_anticipation.py` | not read by shipped code | One-row summary: HH21 anticipation analysis |
| `f2c_archr_snapatac_concordance.tsv`, `f2c_pitx2_promoter_check.tsv`, `f2c_pool_concordance.tsv` | `pipeline/f2c_archr_concordance.py` | not read by shipped code | ArchR vs SnapATAC2 peak-set concordance at the PITX2 promoter (sanity check before downstream analyses commit to ArchR peak set). |
| `archr_union_peaks.bed.gz`, `archr_per_pool_peaks.tsv.gz`, `archr_reproducibility_metadata.tsv` | `pipeline/archr_build_peaks.R`, `pipeline/chick/archr_build_peaks.sh`, `pipeline/run_archr_peaks.sh` | figS03, `pipeline/archr_pool_composition.R`, `pipeline/f2e_archr_frip_per_cell.py`, `pipeline/g4e_focal_pitx2_peakset.py`, `pipeline/g4h_response_locus_trough_test.py`, … | The ArchR 283,818-peak union (501 bp fixed width, reproducibility ≥2) across 17 callable pools |
| `archr_pool_composition.tsv` | `pipeline/archr_pool_composition.R`, `pipeline/run_archr_pool_composition.sh` | not read by shipped code |  |
| `f2f_peaks_detected_per_population.tsv` | `pipeline/f2f_peaks_detected_per_population.py`, `pipeline/run_f2f_peaks_detected.sh` | figS03 |  |
| `g2_lpm_hh1716_{INFERENCE_POOL,SUB_Early_LPM}_deseq2{,_archr}.tsv`, `g2_lpm_hh1716_coverage{,_archr}.tsv`, `g2_lpm_hh1716_summary{,_archr}.tsv` (8 files) | `pipeline/g2_lpm_hh1716_deseq2.py` | `pipeline/g4_hh1716_motif_enrichment.R`, `pipeline/g4c_motif_on_aggregated.sh`, `pipeline/g5_hh21_wave2_anticipation.py`, `pipeline/run_atac_hh1716_chain_archr.sh`, … | Pre-h6 paired DESeq2 attempts at HH17/16 LPM. Deprecated in favor of h6 / h7 |
| `step_f4_pitx2_coverage_6track_backup.tsv.gz`, `step_f4_pitx2_coverage_hh1716.tsv.gz` | `pipeline/f4_lr_pitx2_pileup.py` | `scripts/merge_f4_hh1716_into_canonical.py` |  |
| `step_f4b_pitx2_tad_coverage.tsv.gz` | `pipeline/f4b_pitx2_tad_coverage.py` | not read by shipped code | TAD-scale coverage (1.8–4.6 Mb, 5-kb bins) for the TAD-zoom multi-zoom pileup. |
| `duck_locus/locus_percell_counts.tsv.gz` | `pipeline/chick/duck_pitx2_locus_percell_extract.sh`, `pipeline/chick/duck_pitx2_locus_percell_extract.py` | `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/archr_pitx2_locus_pseudobulk.R`, `pipeline/chick/cross_species_halfvshalf.py` |  |
| `m5_esr1_lr_hh25_stats.tsv` | `figures/figS17/stats_esr1_lr_hh25.py`, `figures/figS17/s4_deseq2_refit.R`, `figures/figS17/s8_glmm.R` | figS17 |  |
| `m5_cycling_lr_hh25_stats.tsv` | `figures/figS16/stats_cycling_lr_hh25.py` | figS16 |  |
| `m5_mitogen_program_hh25.tsv` | `figures/figS16/mitogen_program_lookup.py` | figS16 | Mitogen-program + ESR1 timing table read |
| `m5_mitogen_gene_coords_zju1.tsv` | NCBI Gene coordinates for ZJU1.0 (annotation release 104) | figS16 |  |
| `m5_mitogen_loci_pitx2_motif_peaks.tsv` | `figures/figS16/mitogen_loci_pitx2_motif_check.py` | `scripts/make_bmp7_motif_windows.py` | Per-peak mitogen-locus motif/DA table |
| `m5_responding_population_e6e9.tsv` | `figures/figS16/responding_population_e6e9.py` | not read by shipped code |  |
| `bprime_locus_fdr.tsv` | `scripts/compare_bprime.py` | fig03 | The (b′) per-peak locus DESeq2 (paired + unpaired arms × floors) over the 359-peak TAD set |
| `laterality_trajectory.tsv` | `figures/figS03/export_laterality_trajectory.py`, `figures/_shared/subcluster_map.py` | figS03 | Per-(stage × population × metric) per-embryo log₂(L/R) laterality trajectory for the figS03 forest panel |
| `PIPELINE_STATE.json` | `pipeline/run_pipeline.sh` | fig01_supp |  |
| `d2e_resolution_sweep_summary.tsv` | `pipeline/d2e_resolution_umethod.py` | `FC_08656/run_d2e_silhouette.sh` |  |
| `d2e_silhouette_res_comparison.tsv` | `pipeline/d2e_resolution_umethod.py` | not read by shipped code |  |
| `d2f_canonical_markers.tsv.gz` | `pipeline/d2f_canonical_marker_export.py` | not read by shipped code |  |
| `d2f_nkx6_1.tsv.gz` | `pipeline/d2f_canonical_marker_export.py` | not read by shipped code |  |
| `d2f_somatic_markers.tsv.gz` | `pipeline/d2f_canonical_marker_export.py` | not read by shipped code |  |
| `d2g_lineage_scores_*.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2g_lineage_scores_*_ccreg.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2g_uegs_*.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2g_uegs_*_ccreg.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2g_wilcoxon_*.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2g_wilcoxon_*_ccreg.tsv` (4 files) | `pipeline/d2g_stage_restricted_subclustering.py` | `scripts/annotate_d2g_via_compendium.py` |  |
| `d2i_barcodes/` (128 files) | `figures/data/d2i_barcodes/generate_early_stage_barcodes.py` | fig02, `pipeline/chick/duck_pitx2_locus_percell_extract.py`, `pipeline/d2i_pitx2_rna_pileup.py`, `pipeline/h6_paired_deseq2.py`, `pipeline/run_h6_locus_tad_extract.sh`, … |  |
| `d2i_pileup_output/d2i_pitx2_gtf_features.tsv` | `pipeline/d2i_pitx2_rna_pileup.py` | not read by shipped code |  |
| `d2i_pileup_output/d2i_pitx2_pileup_summary.txt` | `pipeline/d2i_pitx2_rna_pileup.py` | not read by shipped code |  |
| `d2i_pileup_output/d2i_pitx2_rna_pileup.tsv.gz` | `pipeline/d2i_pitx2_rna_pileup.py` | fig02 |  |
| `d2i_pileup_output/d2i_pitx2_rna_pileup_dedup.tsv.gz` | `pipeline/d2i_pitx2_rna_pileup.py` | not read by shipped code |  |
| `duck_locus/locus_peak_coords.tsv` | `pipeline/chick/duck_pitx2_locus_percell_extract.sh` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R`, `pipeline/chick/build_isoform_data.py`, `pipeline/chick/cross_species_halfvshalf.py`, `pipeline/chick/duck_pitx2_locus_percell_extract.py` |  |
| `g2c_aggregated_HH1716_SUB_BARX1.tsv` | `pipeline/g2c_aggregated_descriptive.py` | `pipeline/_lr_binom.py`, `pipeline/g4c_motif_on_aggregated.sh`, `pipeline/h2_gkmsvm_train.py`, `pipeline/h3_jaspar_scan.R`, … |  |
| `g2c_aggregated_HH1716_SUB_Cardiopharyngeal_LPM.tsv` | `pipeline/g2c_aggregated_descriptive.py` | `pipeline/_lr_binom.py`, `pipeline/g4c_motif_on_aggregated.sh`, `pipeline/h2_gkmsvm_train.py`, `pipeline/h3_jaspar_scan.R`, … |  |
| `g2c_aggregated_HH1716_SUB_Early_LPM.tsv` | `pipeline/g2c_aggregated_descriptive.py` | `pipeline/_lr_binom.py`, `pipeline/g4c_motif_on_aggregated.sh`, `pipeline/h2_gkmsvm_train.py`, `pipeline/h3_jaspar_scan.R`, … |  |
| `g3_pitx2_motif_scan.tsv` | `pipeline/g3_pitx2_motif_scan.R` | not read by shipped code |  |
| `g3_pitx2_motif_scan_summary.tsv` | `pipeline/g3_pitx2_motif_scan.R` | not read by shipped code |  |
| `g3_pitx2_motif_scan_summary_p1e-2.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3_pitx2_motif_scan_summary_p5e-3.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3_pitx2_motif_scan_summary_p5e-4.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3b_candidate_peaks_motif_scan_summary_p1e-2.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3b_candidate_peaks_motif_scan_summary_p5e-3.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3b_candidate_peaks_motif_scan_summary_p5e-4.tsv` | `pipeline/g3b_pitx2_candidate_peaks_motif_scan.R` | not read by shipped code |  |
| `g3f_hh21_duck_motifs.tsv` | `pipeline/g3f_hh21_motif_scan.R` | not read by shipped code |  |
| `g3f_hh21_window_gc.tsv` | `pipeline/g3f_hh21_motif_scan.R` | not read by shipped code |  |
| `g3f_motif_pfms.tsv` | `pipeline/g3f_hh21_motif_scan.R` | fig04 |  |
| `g3g_hh21_chick_motifs.tsv` | `pipeline/g3g_hh21_chick_lift_scan.R` | not read by shipped code |  |
| `g3g_hh21_lift.tsv` | `pipeline/g3g_hh21_chick_lift_scan.R` | not read by shipped code |  |
| `g4f_census_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_cpm_trajectory_long_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_domain_motif_windows.bed` | `pipeline/g4f_pitx2_locus_dynamics.py` | `pipeline/g3f_hh21_motif_scan.R`, `pipeline/g3g_hh21_chick_lift_scan.R`, `pipeline/run_g3f_g3g_domain.sh`, `scripts/make_bmp7_motif_windows.py`, … |  |
| `g4f_locus_count_slice_bysubcluster.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_locus_count_slice_bysubcluster_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_locus_count_slice_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | `scripts/depthmatch_pitx2_locus.py` |  |
| `g4f_lr_grid_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_lr_long_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | fig03, `scripts/make_g4f_domain_windows.py`, `scripts/plot_focal_locus_lr_by_compartment.py` |  |
| `g4f_sample_meta_bysubcluster.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_sample_meta_bysubcluster_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_sample_meta_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | fig03, figS03, `scripts/depthmatch_pitx2_locus.py` |  |
| `g4f_stage_grid_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4f_stage_long_floor20.tsv` | `pipeline/g4f_pitx2_locus_dynamics.py` | not read by shipped code |  |
| `g4h_cohorts.tsv` | `pipeline/g4h_prep_trough_lr.py` | `pipeline/run_g4h.sh`, `pipeline/run_g4h_prep.sh`, `pipeline/run_grn_memnode_scan.sh` |  |
| `h0_lvr_peaks_pitx2_locus.tsv` | `scripts/locate_lvr_peaks_pitx2_locus.py` | `pipeline/g4e_focal_pitx2_peakset.py`, `pipeline/g4h_response_locus_trough_test.py`, `pipeline/run_g4h.sh`, `scripts/compute_locus_fdr_h6_h7.py`, … |  |
| `h1_excluded_pseudobulk_counts.tsv.gz` | `pipeline/h1_excluded_pseudobulk.py` | `pipeline/run_pipeline.sh` |  |
| `h1_excluded_pseudobulk_metadata.tsv` | `pipeline/h1_excluded_pseudobulk.py` | `pipeline/run_pipeline.sh` |  |
| `h1_pseudobulk_snrna_counts.tsv.gz` | `pipeline/h1_pseudobulk_snrna.py` | figS17, `pipeline/h2_deseq2_snrna.R`, `pipeline/run_hh25_relabel_rna_panelB.sh`, `pipeline/run_hh25_relabel_rna_panelB_h2.sh`, `pipeline/run_pipeline.sh` |  |
| `h1_pseudobulk_snrna_metadata.tsv` | `pipeline/h1_pseudobulk_snrna.py` | figS17, `pipeline/h2_deseq2_snrna.R`, `pipeline/run_hh25_relabel_rna_panelB.sh`, `pipeline/run_hh25_relabel_rna_panelB_h2.sh`, `pipeline/run_pipeline.sh` |  |
| `h2_deg_excluded_POOL_*.tsv` (4 files) | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_excluded_SUB_Contamination_HH17_16.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_excluded_SUB_PRRX1__FOXF1_low_mes__HH17_16.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_excluded_SUB_Skeletogenic_LPM_HH17_16.tsv`, `h2_deg_excluded_SUB_Skeletogenic_LPM_HH21.tsv`, `h2_deg_excluded_SUB_Skeletogenic_LPM_HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | figS03 |  |
| `h2_deg_excluded_SUB_Venous_pole_spl__HH17_16.tsv`, `h2_deg_excluded_SUB_Venous_pole_spl__HH21.tsv`, `h2_deg_excluded_SUB_Venous_pole_spl__HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_excluded_summary.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_POOL_*.tsv` (4 files) | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_SUB_BARX1__HH17_16.tsv`, `h2_deg_unpaired_SUB_BARX1__HH21.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_SUB_Cardiopharyngeal_LPM_HH17_16.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_SUB_Early_LPM_HH13.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_SUB_Syringeal_HH21.tsv`, `h2_deg_unpaired_SUB_Syringeal_HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_SUB_Ventral_spl__HH21.tsv`, `h2_deg_unpaired_SUB_Ventral_spl__HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | `pipeline/g5_hh21_wave2_anticipation.py` |  |
| `h2_deg_unpaired_excluded_POOL_*.tsv` (4 files) | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_excluded_SUB_Contamination_HH17_16.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_excluded_SUB_PRRX1__FOXF1_low_mes__HH13.tsv`, `h2_deg_unpaired_excluded_SUB_PRRX1__FOXF1_low_mes__HH17_16.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_excluded_SUB_Skeletogenic_LPM_HH17_16.tsv`, `h2_deg_unpaired_excluded_SUB_Skeletogenic_LPM_HH21.tsv`, `h2_deg_unpaired_excluded_SUB_Skeletogenic_LPM_HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_excluded_SUB_Venous_pole_spl__HH17_16.tsv`, `h2_deg_unpaired_excluded_SUB_Venous_pole_spl__HH21.tsv`, `h2_deg_unpaired_excluded_SUB_Venous_pole_spl__HH25.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_excluded_summary.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `h2_deg_unpaired_summary.tsv` | `pipeline/h2_deseq2_snrna.R` | not read by shipped code |  |
| `hh21_seq1_bound_sensitivity_PITX2.tsv` | `figures/figS03/hh21_trans_motif_seq.py` | not read by shipped code |  |
| `hh21_seq1_tail_null_SMAD2_3_4.tsv` | `figures/figS03/hh21_trans_motif_seq.py` | not read by shipped code |  |
| `step02e_doublet_comparison.tsv` | `pipeline/c1_qc_filter.py` | fig01_supp, `FC_08656/20260326_step02e_preCB_genotype.py` |  |
| `step02e_mito_by_library.tsv.gz` | `pipeline/c1_qc_filter.py` | `FC_08656/20260326_step02e_preCB_genotype.py` |  |
| `step03_sex_discordance.tsv` | `pipeline/c3_sex_determination.py` | not read by shipped code |  |
| `step04a_known_markers.tsv.gz` | `pipeline/d2_cluster_markers.py` | not read by shipped code |  |
| `step04a_markers_res0.2.tsv` | `pipeline/d2_cluster_markers.py` | `pipeline/d2b_validate_clusters.py` |  |
| `step04a_markers_res0.5.tsv` | `pipeline/d2_cluster_markers.py` | not read by shipped code |  |
| `step04a_markers_res1.0.tsv` | `pipeline/d2_cluster_markers.py` | not read by shipped code |  |
| `step04b_cluster_validation.tsv` | `pipeline/d2b_validate_clusters.py` | not read by shipped code |  |
| `step07c_genotype_celltype_composition.tsv` | generate_dotplot_export | `FC_08656/20260330_generate_dotplot.py` |  |

## Chick control

| File | Produced by | Read by | Description |
|---|---|---|---|
| `chick_atac_metrics.tsv` | `pipeline/chick/hQ1_parse_arc_metrics.py` | figS04_supp | Per-library cellranger-arc summary metrics for the 5 chick ARC libraries (galGal6a): est_cells, TSS enrichment, FRIP, duplicate %, peak counts. |
| `chick_atac_per_cell_sample.tsv` | `pipeline/chick/hQ2_sample_per_cell.py` | figS04_supp |  |
| `chick_atac_f1_qc_summary.tsv` | `pipeline/chick/hQ3_snapatac_f1_qc.py` | figS04_supp | Per-chick-library SnapATAC2 f1 QC: n_imported, n_qc_pass (TSS ≥ 2 & fragments ≥ 1000), median TSSe, median n_fragment. |
| `chick_atac_f1c_amulet_summary.tsv` | `pipeline/chick/hQ4_amulet.py` | figS04_supp |  |
| `chick_archr_coupled_metadata.tsv` | `pipeline/chick/archr_metadata_sidecar.py` | `pipeline/chick/archr_build_peaks.R`, `pipeline/chick/archr_build_peaks.sh`, `pipeline/chick/archr_pitx2_locus_pseudobulk.R` |  |
| `chick_locus/locus_percell_counts.tsv.gz` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R`, `pipeline/chick/archr_pitx2_locus_pseudobulk.sh` | `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/cross_species_halfvshalf.py`, `pipeline/chick/duck_pitx2_locus_percell_extract.py` |  |
| `chick_locus_genorepro/locus_peak_coords.tsv` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R`, `pipeline/chick/archr_pitx2_locus_pseudobulk.sh` | `pipeline/chick/build_isoform_data.py`, `pipeline/chick/cross_species_halfvshalf.py`, `pipeline/chick/duck_pitx2_locus_percell_extract.py` | Panel-G sensitivity sidecar: 69 peaks in PITX2 ±500 kb extracted from the genotype-level (vs library-level) ArchR reproducibility build |
| `chick_control/cAE_atlas_umap.tsv.gz` | `pipeline/chick/cAE_full_atlas.py` | figS04, figS04_supp, `pipeline/chick/cAF_duck_panel_dotplot.py` | Chick atlas per-cell barcode, group, library/lane, donor, stage, sex, celltype + UMAP1/2 (37,362 cells) |
| `chick_control/cAF_duck_panel_dotplot.tsv` | `pipeline/chick/cAF_duck_panel_dotplot.py` | figS04 | Chick marker dotplot export on the duck fig01 panel-D gene set (annotation, gene, frac_expressing, mean_raw) |
| `chick_control/cAJ_lpm_pseudobulk_pca.tsv` | `pipeline/chick/cAJ_lpm_pseudobulk_pca.py` | figS04 | LPM pseudobulk PCA per group×donor (PC coords + var, stage/sex/gfp/n_cells) |
| `chick_control/cAD_pitx2_pileup_depth.tsv.gz`, `cAD_pitx2_pileup_tracks.tsv`, `cAD_pitx2_gene_model.tsv`, `cAD_pitx2_polya_priming.tsv` (4 files) | `pipeline/chick/cAD_groupB_pitx2_pileup_d2i.py` | `pipeline/chick/cAD_plot_pileup.R` |  |
| `chick_control/cAN_halfvshalf.tsv`, `cAN_domain_counts.tsv`, `cAN_directional_perembryo.tsv`, `cAN_directional_summary.tsv` (4 files) | `pipeline/chick/cAN_stage_laterality.py` | figS04, figS04_supp, `pipeline/chick/_chick_subcluster_pools.py`, `pipeline/chick/build_isoform_data.py`, `pipeline/chick/per_stage_laterality.py` |  |
| `chick_control/cAK_atlas_wz_percell.tsv.gz` | `pipeline/chick/cAK_atlas_sex_wgene.py` | figS04_supp | Per-cell W/Z-gene fractions from the cAE atlas h5ad, after the resolved-12 Group B relabel (barcode, group, donor, sex, percent_w, percent_z) |
| `chick_control/groupB_resolved12_percell.tsv.gz` | `pipeline/chick/cAK_atlas_sex_wgene.py`, pruned of the 163 split-excluded cells by `pipeline/chick/_groupB_split_excluded.py --prune-resolved12` | not read by shipped code | Resolved-12 Group B per-cell subset of the cAK W/Z table (cS donor2/donor11 split record) |
| `chick_control/groupB_split_excluded_cells.tsv` | `pipeline/chick/_groupB_split_excluded.py` | figS04_supp, `pipeline/chick/_groupB_resolved12.py`, `pipeline/chick/cAN_stage_laterality.py`, `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/build_orientfree_isoform_data.py`, `pipeline/chick/cross_species_halfvshalf.py`, … | The 163 Group B nuclei that the donor11/donor2 N=2 vireo splits call doublet/unassigned (barcode, base_donor, split_status); excluded from every per-embryo analysis |
| `chick_control/chick_stage_alignment.tsv` | `pipeline/chick/cAL_stage_subclustering.py` | figS04_supp, `pipeline/chick/cAN_stage_laterality.py` | (group, donor) → firm stage for all 24 chick embryos |
| `chick_control/cA_crosslane_3lane_long.tsv`, `chick_control/cA_crosslane_3lane_block.tsv` | `pipeline/chick/cA_crosslane_3lane.py` | figS04_supp | Cross-lane allele-frequency correlation r: long form for the 3 lane-pairs, and the full (lane,donor)×(lane,donor) block matrix |
| `chick_control/cT_groupA_stage_sex.tsv` | `pipeline/chick/cT_groupA_stage_sex.py` | figS04_supp, `pipeline/chick/cAE_full_atlas.py`, `pipeline/chick/cV2b_groupA_embryo_alignment.py` | Group A per-donor stage + sex, ARC0S-triangulated 4/4/4 |
| `chick_control/chick_demux_summary_n24.tsv` | `pipeline/chick/cV2b_groupA_embryo_alignment.py` | figS04_supp, `pipeline/chick/_chick_lr_orientation.py`, `pipeline/chick/_groupB_resolved12.py`, `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/h1b_power_census.py` |  |
| `chick_control/cb_libqc_perlane.tsv` | `pipeline/chick/cb_libqc_export.py` | figS04_supp | Per-lane chick snRNA QC summary (the panel A bar data) |
| `chick_control/cAL_subcluster_annotations.tsv` | `pipeline/chick/cAM_compendium_annotate.py`, `pipeline/chick/cAL_stage_subclustering.py` | `pipeline/chick/_chick_subcluster_pools.py`, `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/cAN_stage_laterality.py` | Chick d2g subcluster SSOT: canonical (stage, resolution, cluster) → population, incl |
| `chick_control/cAM_percell_populations.tsv.gz` | `pipeline/chick/cAM_compendium_annotate.py` | `pipeline/chick/archr_metadata_sidecar.py`, `pipeline/chick/build_isoform_data.py`, `pipeline/chick/build_orientfree_isoform_data.py`, `pipeline/chick/cAN_stage_laterality.py`, … | Chick per-cell population label at the canonical resolution (group, library/lane, donor, stage, sex, leiden, PITX2_cp10k, population). Holds 92 Group B split doublet/unassigned nuclei (as donor11a/donor2a); readers drop them via `groupB_split_excluded_cells.tsv` |
| `chick_control/xspecies_pitx2_isoform_laterality.tsv` | `pipeline/chick/build_isoform_data.py` | figS04 | ORIENTED duck-vs-chick per-embryo PITX2 isoform L/R contrast (ATAC PITX2a/b/c + RNA total) |
| `chick_control/xspecies_orientfree_isoform.tsv` | `pipeline/chick/build_orientfree_isoform_data.py` | fig04, `pipeline/chick/build_isoform_data.py` | ORIENTATION-FREE per-embryo isoform divergence vs within-embryo null (modality, isoform, species, genotype, stage, n_min, obs, null, gap, emp_p) |
| `chick_control/sex_wgene/c5h_corrected_donor_table.tsv` | `pipeline/chick/c5h_validate_seed.py` | `pipeline/chick/c5h_plot_corrected.py`, `pipeline/chick/cAE_full_atlas.py`, `pipeline/chick/cAJ_lpm_pseudobulk_pca.py`, `pipeline/chick/cAK_atlas_sex_wgene.py` | Group B per-donor corrected sex + gfp calls (N14/seed105), the upstream identity table baked into cAE/cAJ/cAK on the cluster |
| `chick_control/cAJ_lpm_pseudobulk_pca_groupA.tsv` | `pipeline/chick/cAJ_lpm_pseudobulk_pca.py` | figS04 |  |
| `chick_control/cAJ_lpm_pseudobulk_pca_groupB.tsv` | `pipeline/chick/cAJ_lpm_pseudobulk_pca.py` | figS04 |  |
| `chick_locus/locus_peak_coords.tsv` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R` | `pipeline/chick/build_isoform_data.py`, `pipeline/chick/cross_species_halfvshalf.py`, `pipeline/chick/duck_pitx2_locus_percell_extract.py` |  |
| `chick_locus/locus_pseudobulk_counts.tsv` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R` | not read by shipped code |  |
| `chick_locus/locus_pseudobulk_metadata.tsv` | `pipeline/chick/archr_pitx2_locus_pseudobulk.R` | not read by shipped code |  |

## E6/E9/E11 scRNA

| File | Produced by | Read by | Description |
|---|---|---|---|
| `scrna_e6e9e11/cr9_er_survey_percell.tsv.gz`, `scrna_e6e9e11/cr9_global_cluster_markers.tsv` | `pipeline/scrna_e6e9e11/cr9_export_er_survey.R` | fig06, figS07, figS16 |  |
| `scrna_e6e9e11/substage_E09_percell.tsv.gz`, `scrna_e6e9e11/substage_E11_percell.tsv.gz`, `scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv` | `pipeline/scrna_e6e9e11/substage_subcluster_cr9.R`, `pipeline/scrna_e6e9e11/substage_compendium_annotate.R` | figS07, figS08, `pipeline/scrna_e6e9e11/cr9_screen_export.R`, `pipeline/scrna_e6e9e11/cr9_syrinx_pseudobulk_from_h5.py` |  |
| `scrna_e6e9e11/cr9_screen_E06_percell.tsv.gz`, `scrna_e6e9e11/cr9_screen_E09_percell.tsv.gz`, `scrna_e6e9e11/cr9_screen_E11_percell.tsv.gz` | `pipeline/scrna_e6e9e11/cr9_screen_export.R` | figS08, figS09, figS10 |  |
| `scrna_e6e9e11/esr1_level_lr_tidy.tsv`, `scrna_e6e9e11/esr1_cellcycle_lr_tidy.tsv` | `figures/figS07/plot_esr1_level_lr.py`, `figures/figS07/plot_esr1_cellcycle.py` | figS16 |  |
| `scrna_e6e9e11/er_perembryo_ie.tsv`, `scrna_e6e9e11/er_intron_exon_umi.tsv`, `scrna_e6e9e11/esrrg_coverage_binned.tsv`, `scrna_e6e9e11/esrrg_exons.bed` (4 files) | `pipeline/scrna_e6e9e11/cr9_er_perembryo_ie.sh`, `pipeline/scrna_e6e9e11/cr9_er_intron_exon_umi.sh`, `pipeline/scrna_e6e9e11/cr9_esrrg_pileup.sh` | figS07 |  |
| `scrna_e6e9e11/xstage_composition.tsv`, `scrna_e6e9e11/xstage_corr_ribbons.tsv`, `scrna_e6e9e11/xstage_auroc_matrix.tsv` | `pipeline/scrna_e6e9e11/xstage_correspondence.R` | figS07, figS16 |  |
| `scrna_e6e9e11/cr9_syrinx_lr_deseq2_{E09,E11}_{M,F}.tsv` (4 files) | `pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R` | figS16 |  |
| `scrna_e6e9e11/cr9_syrinx_pseudobulk_counts_E09.tsv.gz`, `scrna_e6e9e11/cr9_syrinx_pseudobulk_counts_E11.tsv.gz` | `pipeline/scrna_e6e9e11/cr9_syrinx_pseudobulk_from_h5.py` | `pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R` |  |
| `scrna_e6e9e11/cr9_syrinx_pseudobulk_meta_E09.tsv`, `scrna_e6e9e11/cr9_syrinx_pseudobulk_meta_E11.tsv` | `pipeline/scrna_e6e9e11/cr9_syrinx_pseudobulk_from_h5.py` | `pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R` |  |

## Design and reference tables (outside `figures/data/`)

| File | Contents |
|---|---|
| `design/embryo_lr_design.tsv` | Experimental left/right design: embryo halves per library |
| `design/genotype_reference.tsv` | Genotype → embryo group reference |
| `design/chick_groupA_lr_design.tsv`, `design/chick_groupB_lr_design.tsv` | Chick control left/right design (Group A and Group B lanes) |
| `references/color_registry.tsv` | Colours for cell types, populations, stages, sides, sex and libraries (read by `figures/_shared/colors.{py,R}`) |
| `references/cluster_annotations.tsv` | Global cluster → cell type (written by `pipeline/d2b_validate_clusters.py`) |
| `references/d2g_subcluster_annotations.tsv` | Stage × subcluster → population (written by `pipeline/d2g_validate_subclusters.py`; read through `figures/_shared/subcluster_map.{py,R}`) |
| `references/genotype_metadata.tsv` | Genotype → display name, group, stage and demultiplexing source for the 52 embryos (read by `figures/_shared/style.py`) |
