# Tier-2 build order — re-deriving `figures/data/` from the GEO record

The source for the README's Tier-2 section. Tier 1 (figures from committed tables) needs
none of this. Sequencing below is load-bearing: several steps hard-require another
step's output, and two figure scripts have a cross-figure order.

## Start from the deposited objects, not from a0

The deposited `multiome_step02_preCB_clean.h5ad` and
`multiome_step04_preCB_clustered.h5ad` ARE the pipeline masters (annotations,
sides, sexes, genotype→embryo in `.obs`). **Do not re-run a0–d2 to reproduce
figures**: demultiplexing and clustering are the in-kind steps (below), and a
d1/d2 re-run reshuffles Leiden IDs against `references/cluster_annotations.tsv`.
The from-raw path exists (deposited BAMs + committed `FC_08656/` + `pipeline/a*`
recipes) but is for auditing the pipeline, not for rebuilding figure inputs.

## Order of operations

1. **snRNA exports (from the deposited atlases).**
   From step04: `d2g_stage_restricted_subclustering.py` (+`d2g_signaling_export.py`,
   then `d2g_validate_subclusters.py` to regenerate the subcluster SSOT),
   `d2b_validate_clusters.py`, `d2b_export_genotype_metrics.py`,
   `d2e_silhouette_export.py`, `e1_pseudobulk_pca.py` → `e1b_stage_inference.py`,
   `e2_laterality_pitx2.py`, `e3_lr_composition.py`,
   `FC_08656/20260330_generate_dotplot.py`. From step02: `c3_sex_determination.py`,
   `h1_*` → `h2_deseq2_snrna.R` (all three arms), `h3_gene_detection_export.py`,
   `scripts/export_hh25_percell_counts.py`.
   `d2_umap_boundary_data.tsv.gz` / `step04a_cluster_counts.tsv`: the literal
   producer reads an un-deposited intermediate h5ad, but every exported column is
   in step04 — re-derive with a trivial reader, do not re-run d1/d2.
2. **BAM-derived RNA tracks.** `d2j_pitx2_isoform_{percell,per_genotype}.py` and
   `d2i_pitx2_rna_pileup.py` read the deposited `*_gex_possorted_bam.bam` + the
   deposited `ZJU1.0_PITX2ext_genes.gtf.gz`. Note: the committed
   `figures/data/d2i_barcodes/` HH25 lists are the published vintage; rebuilding
   them from the current SSOT changes the HH25 tracks.
3. **Local inference chain (laptop).** `d2l_lr_inference.py` →
   `d2m_genotype_to_embryo.py`; after step 4's f6: `d2n_…hh1716.py` →
   `d2o_merge_multiome_inference.py` (→ `lr_inference_multiome.tsv`) →
   `d2m_multiome_genotype_to_embryo.py`.
4. **ATAC.** Entry point is **`pipeline/f1_coupled_atac.py`** (NOT a
   `run_pipeline.sh` phase): deposited `*_atac_fragments.tsv.gz` + step04 h5ad +
   committed d2g SSOT → the f1 h5ads and `step_f1_coupled_atac_metadata.tsv`
   (hard-required by g1). Then either use the deposited
   `multiome_atac_peak_matrix_archr.h5ad` + `…union_peaks_archr.bed.gz` directly,
   or rebuild them (`archr_build_bsgenome.R` → `archr_build_peaks.R` →
   `f2d_archr_peak_matrix.py`). Then `f4`/`f5`/`f6` (fragments + intervals in
   `references/*.bed`), `g1_pseudobulk_atac.py` (→ the `da_archr/` tree), and the
   consumers: `g4f`/`g4g`, `h4`, `h6`/`h7`, `g4h_prep` → the memory-node scan,
   `run_hh21_module_v2.sh`, bprime. GC tracks: `scripts/run_gc_content_pitx2.sh`
   against the deposited genome FASTA.
5. **Chick.** Group-B GEX substrate first: rebuild the reference per
   `pipeline/chick/refbuild/README.md`, then `pipeline/rederive_cellranger_from_geo.sh
   gex GEX09 <ARC09_gex_possorted_bam.bam> <refbuild ref> <work>` (same for GEX10;
   the deposited ARC09/10 GEX BAMs are the same GEM wells). Chick ATAC needs no
   re-alignment — ARC01/02/09/10/ARC0S-chick are deposited as their galGal6a
   analysis-of-record outputs. Atlas-derived exports come from the deposited
   `chick_multiome_cAE_atlas.h5ad`; locus pseudobulks via `pipeline/chick/archr_*`.
6. **Old-duck (E6/E9/E11).** The deposited `scrna_cr9_03_integrated_CB.{rds,h5ad}`
   carries the per-cell demux metadata; the substage exports re-derive from it via
   the `pipeline/scrna_e6e9e11/` producers. Panel-G BAM tracks: regenerate BAMs with
   `cellranger count --create-bam` from the deposited FASTQs + deposited scRNA
   reference.
7. **QC-panel `outs/` metrics (only if regenerating fig01_supp A–E / figS03 A /
   figS04_supp E).** `pipeline/rederive_cellranger_from_geo.sh` per library
   (gex ×12, arc ×12 duck + 5 chick) — ~a day per arc library; nothing in the
   science figures needs this.
8. **Cross-figure order.** `figures/fig02_supp/generate_dotplot_extended.py`
   (figS02 panel D) BEFORE `figures/fig02/generate_spatial_svg.py` (fig02 panel B);
   `figures/figS03/panel_H_cis_memory_specificity.R` before the fig03 composite.

## Exact vs in-kind (state this; do not over-claim)

- **Exact given pinned versions** (deterministic): f1 QC/pairing, d2j/d2i BAM
  counting, f4/f5/f6 fragment counting, GC tracks, `mkref`, `cellranger{,-arc}
  count` from bamtofastq'd reads (10x BAMs retain all reads), the seeded local
  stats (depthmatch, labelswap).
- **In-kind only** (stochastic or label-permuting; the deposited objects and
  committed exports are the record): CellBender, scDblFinder, cellSNP/vireo
  demultiplexing (donor labels permute), the old-duck substage relabeling
  (`figures/figS08/PROVENANCE.md` — counts and Leiden IDs drift between runs), and any
  d1/d2/d2g re-clustering (IDs reshuffle; annotations flow only through the
  committed SSOT TSVs).
- Versions of record: cellranger 7.1.0, cellranger-arc 2.1.0, CellBender 0.3.2,
  the `duckfig` + `snapatac2` conda envs and the O2 R 4.4 library (captured in
  `envs/`).
