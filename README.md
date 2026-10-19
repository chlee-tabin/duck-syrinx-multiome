# Duck syrinx multiome analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23476404.svg)](https://doi.org/10.5281/zenodo.23476404)

Analysis code for the single-nucleus multiome (snRNA-seq + snATAC-seq) and
single-cell RNA-seq study of left-right asymmetry and sexual dimorphism in the
developing duck syrinx: genotype demultiplexing, snRNA/snATAC processing,
left/right inference, differential expression and accessibility, and every
computational figure panel.

## Citation

> Mishkind D, Lee C, Prabhakara C, Tir P, Kingsley EP, Clarke JA, Goller F,
> Tabin CJ. Integration of sexual dimorphism and left-right asymmetry in the
> development of the duck syrinx. *bioRxiv* (2024).
> [doi:10.1101/2024.12.10.627783](https://doi.org/10.1101/2024.12.10.627783)

The peer-reviewed version is in press at *Proc Natl Acad Sci USA*.

The code is archived at Zenodo: v1.0.0 is
[doi:10.5281/zenodo.23476405](https://doi.org/10.5281/zenodo.23476405), and
[doi:10.5281/zenodo.23476404](https://doi.org/10.5281/zenodo.23476404) always
points to the latest version.

Machine-readable citation metadata: [`CITATION.cff`](CITATION.cff). If you use
this code, please cite the paper.

## Data availability

All sequencing data are in the NCBI Gene Expression Omnibus, linked as one
SuperSeries:

| Dataset | GEO | SRA BioProject |
|---|---|---|
| Single-nucleus multiome (snRNA + snATAC): duck HH13–HH25 and the chick control cohort (36 samples = 18 libraries × GEX + ATAC) | [GSE336327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336327) | PRJNA1481253 |
| 3′ single-cell RNA-seq, duck E6/E9/E11, with MULTI-seq multiplexing | [GSE336329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336329) | PRJNA1481424 |

Besides raw reads and BAMs, the records carry processed objects sufficient to
reproduce the analyses: the annotated snRNA atlas (per-nucleus cell type, stage,
sex, embryo and inferred left/right side), the ArchR cell-by-peak matrix and its
283,818-peak union peak set, the SnapATAC2/MACS3 fragment and peak matrices, the
chick control atlas, the integrated E6/E9/E11 scRNA object, and the custom
references (duck ZJU1.0, GCF_015476345.1, with a 3′-extended PITX2 annotation;
chick Galgal6a) with their build scripts.

## Reproducing the results

There are two tiers.

### Tier 1: figures from the committed tables (laptop, no sequencing data)

Every figure panel renders from tables committed under
[`figures/data/`](figures/data/) (provenance table:
[`figures/data/README.md`](figures/data/README.md)). No cluster and no
sequencing data are needed.

```bash
git clone https://github.com/chlee-tabin/duck-syrinx-multiome.git
cd duck-syrinx-multiome
conda env create -f envs/duckfig_portable.yml   # Python 3.12 + R 4.4; Linux and macOS

# Compose a figure from its committed panels -> figures/fig03.pdf
conda run -n duckfig python3 figures/fig03/generate.py

# Re-render panels from figures/data/, then compose (example: fig03)
conda run -n duckfig Rscript figures/fig03/panel_A_pitx2c_signal.R
conda run -n duckfig Rscript figures/fig03/panel_B_lr_counts.R
conda run -n duckfig Rscript figures/fig03/panel_C_pileup.R
conda run -n duckfig python3 figures/fig03/generate.py
```

Each figure directory has a `LEGEND.md` (the figure legend) and a
`PROVENANCE.md` that lists, panel by panel, the data files read and the exact
commands that build them. The compositor `figures/<fig>/generate.py` writes
`figures/<fig>.pdf` from the panel PNG/SVG files in that directory (`fig06`, a
single panel, has no compositor). The
composites and panels as published are committed, so a rebuild can be compared
against them. Generators write their outputs in place (panels, composites, and
for the statistics scripts the corresponding `figures/data/` table), so after a
rebuild `git status` / `git diff` shows whether anything changed.

Notes:

- Two figure environment files are provided, both named `duckfig`.
  `envs/duckfig_portable.yml` covers every package the figure code imports,
  pinned to the versions the figures were built with, with no build strings and
  no platform-specific packages; conda resolves the builds for your platform. It
  includes the R packages of the figS03/figS17 mixed models (`lme4`, `glmmTMB`,
  `emmeans`) and of the figS03 motif scans (`motifmatchr`, `TFBSTools`,
  `universalmotif`, `JASPAR2020`, `Biostrings`, `BSgenome`, `BSgenomeForge`).
  Dry-run solves succeed for linux-64 and osx-64 (conda 26.1, `libmamba`, strict
  channel priority, conda-forge + bioconda only); the osx-64 environment has been
  built and reproduces the committed figS17 statistics table and panels. Use it on
  every platform.
  `envs/duckfig.yml` is the exact macOS (osx-64) export of the environment most
  panels were built with, kept as a record. It does **not** contain the
  mixed-model or motif-scan R packages above, so with it
  `figures/figS03/hh21_hh1716_nb_glmm.R`, `figures/figS17/s8_glmm.R` and the figS03
  motif scans do not run (those were run with a separate R 4.4 library); it also
  contains macOS-only packages and will not create on Linux.
- The figS03 motif scans (`figures/figS03/hh21_trans_motif_scan.R`,
  `hh21_smad234_dinuc_null.R`) also need the duck genome as a BSgenome package,
  `BSgenome.Anasplatyrhynchos.ZJU1.0`, which is not on Bioconductor. Build it once
  from the deposited ZJU1.0 FASTA with `pipeline/archr_build_bsgenome.R`
  (BSgenomeForge, available in `duckfig`): set `DUCK_REFERENCE_ROOT` so that
  `crarc_zju1_PITX2/fasta/genome.fa` (with its `.fai`) is the FASTA and
  `DUCK_SCRATCH_ROOT` to a work directory; the script builds the package and
  installs it into `~/R/4.4.2`, which both scripts already put first on
  `.libPaths()`. Every other panel renders without it.
- figS08 panel C is a wet-lab in situ hybridization image. Its source images are
  not distributed here; the built panel PNG is committed and the figure composes
  from it.

### Tier 2: the full pipeline from GEO (SLURM cluster)

The pipeline that produced `figures/data/` lives in [`pipeline/`](pipeline/)
(with per-flowcell processing scripts in [`FC_08656/`](FC_08656/)). It was run
on a SLURM cluster. It reproduces from the deposited GEO objects, but it is
documented reproduction, not a turnkey workflow: job scheduling, module loads
and storage layout are the user's to adapt.

Set the storage roots before running anything. Scripts read them from the
environment; a root left unset falls back to the original cluster location, so
set all of them:

```bash
export DUCK_SCRATCH_USER=/path/to/work                 # parent work tree; the chick and
                                                       # E6/E9/E11 work trees hang off it
export DUCK_SCRATCH_ROOT=$DUCK_SCRATCH_USER/FC_08656   # multiome pipeline outputs
export DUCK_REFERENCE_ROOT=/path/to/references         # genome references (FASTA, GTF)
export DUCK_HOME=$HOME                                 # conda installs (a few R scripts use ~/R/4.4.2 directly)
export DUCK_STANDBY_DATA_ROOT=/path/to/rawdata         # raw-data tree (a few staging scripts)
export DUCK_STANDBY_REPO=$DUCK_STANDBY_DATA_ROOT/duck-syrinx-multiome
```

The roots are defined in [`pipeline/_paths.py`](pipeline/_paths.py) (Python),
[`scripts/_paths.sh`](scripts/_paths.sh) (shell) and
[`scripts/_paths.R`](scripts/_paths.R) (R); see
[`docs/audit_hardcoded_paths.md`](docs/audit_hardcoded_paths.md) for the
rules, and run `make check-portable-paths` to confirm that no
script carries an unoverridable path. `#SBATCH -o/-e` log paths are read by
SLURM before any shell runs, so a direct `sbatch` needs `--output=` and
`--error=` on the command line.

The build order is in
[`docs/geo_tier2_build_order.md`](docs/geo_tier2_build_order.md). In short:

0. **Stage the GEO files into the layout the scripts read.** The scripts read
   cellranger-style paths under the roots above. Download the files from
   GSE336327 / GSE336329 and place (or symlink) them as follows; `<ARC>` is the
   ARC library name in the GEO file name (for example `ARC03-1A`, `ARC05`,
   `PILOTARC3`, `ARC0S-duck`):

   | GEO file | Path the scripts read |
   |---|---|
   | `multiome_step02_preCB_clean.h5ad` | `$DUCK_SCRATCH_ROOT/phase2/step02_preCB_clean.h5ad` |
   | `multiome_step04_preCB_clustered.h5ad` | `$DUCK_SCRATCH_ROOT/phase2/step04_preCB_clustered.h5ad` |
   | `<ARC>_atac_fragments.tsv.gz` (+ `.tbi`) | `$DUCK_SCRATCH_ROOT/<ARC>/outs/atac_fragments.tsv.gz` (+ `.tbi`) |
   | `<ARC>_{raw,filtered}_feature_bc_matrix.h5` | `$DUCK_SCRATCH_ROOT/<ARC>/outs/{raw,filtered}_feature_bc_matrix.h5` |
   | `<ARC>_gex_possorted_bam.bam` | `$DUCK_SCRATCH_ROOT/staged_v210_gex_bams/<ARC>.gex_possorted_bam.bam`, then `samtools index` it (the BAM readers require the `.bai`) |
   | `<ARC>_atac_possorted_bam.bam` | input of `pipeline/rederive_cellranger_from_geo.sh` only (step 8) |
   | `multiome_atac_peak_matrix_archr.h5ad` | `$DUCK_SCRATCH_ROOT/phase2/atac/step_f2_peak_matrix_archr.h5ad` |
   | `multiome_atac_union_peaks_archr.bed.gz` | gunzip to `$DUCK_SCRATCH_ROOT/phase2/atac_archr/peaks/union_peaks.bed` |
   | `multiome_atac_peak_matrix_snapatac2.h5ad` | `$DUCK_SCRATCH_ROOT/phase2/atac/step_f2_peak_matrix.h5ad` |
   | `ZJU1.0_PITX2ext_genome.fa`, `ZJU1.0_PITX2ext_genes.gtf.gz` | `$DUCK_REFERENCE_ROOT/crarc_zju1_PITX2/fasta/genome.fa` (then `samtools faidx`), `$DUCK_REFERENCE_ROOT/crarc_zju1_PITX2/genes/genes.gtf.gz` |
   | `galGal6a_genome.fa`, `galGal6a_genes.gtf.gz` | `$DUCK_REFERENCE_ROOT/crarc_galgal6a/fasta/genome.fa`, `$DUCK_REFERENCE_ROOT/crarc_galgal6a/genes/genes.gtf.gz` |
   | `chick_multiome_cAE_atlas.h5ad` | `$DUCK_SCRATCH_USER/chick_demux/vireo_seeded/seed_qc/cAE_atlas.h5ad` |
   | `scrna_cr9_03_integrated_CB.rds` | `$DUCK_SCRATCH_USER/old_duck_cr9/master/cr9_03_integrated_CB.rds` |
   | `<CHICK_ARC>_atac_fragments.tsv.gz` (+ `.tbi`), for `ARC01`, `ARC02`, `ARC09`, `ARC10`, `ARC0S-chick` | `$DUCK_SCRATCH_USER/chick_demux/atac_qc/<CHICK_ARC>/outs/atac_fragments.tsv.gz` (+ `.tbi`), read by the chick ArchR steps (`pipeline/chick/archr_*.R`) |
   | `scrna_D{1,2,3}_cr9_filtered_feature_bc_matrix.h5` | not read by any shipped script (the E6/E9/E11 analyses start from `scrna_cr9_03_integrated_CB.rds`) |

   The cellranger-arc library CSVs (`FC_08656/nlibrar*.csv`, `pipeline/nlibraries_v210/*.csv`, `pipeline/chick/nlibraries0S-chick.csv`) write each FASTQ directory as `<FASTQ_DIR>/<original flowcell subpath>`; replace `<FASTQ_DIR>` (and the subpath, if you keep one directory per library) with where you put the FASTQs downloaded from GEO/SRA.

   Inputs GEO does not carry, and how to obtain them:

   - **Chick ATAC QC files.** `pipeline/chick/hQ0_stage_atac_qc_inputs.sh` stages,
     next to each chick fragment file, the cellranger-arc `summary.csv`,
     `per_barcode_metrics.csv` and `atac_peaks.bed` from
     `$DUCK_STANDBY_DATA_ROOT/FC_08656/<lib>_v210/outs/` (ARC01/02/09/10) and
     `$DUCK_STANDBY_DATA_ROOT/FC_08565/ARC0S-chick/outs/`. They are needed only for
     the chick ATAC QC tables (`hQ1`–`hQ4`, figS04_supp). Re-derive each library's
     `outs/` with `pipeline/rederive_cellranger_from_geo.sh arc` from its two
     deposited BAMs and place it at those paths, then run `hQ0`.
   - **Chick demultiplexing inputs.** The chick demux front half
     (`pipeline/chick/c0_stage_inputs.sh`, `c0b_stage_groupB.sh`,
     `g0s0_stage_inputs.sh`) reads cellranger `outs/` directories
     (`gex_possorted_bam.bam` / `possorted_genome_bam.bam` and
     `filtered_feature_bc_matrix/barcodes.tsv.gz`). Re-derive them with
     `rederive_cellranger_from_geo.sh` (`arc` mode for Group A; `gex` mode against
     the `pipeline/chick/refbuild/` reference for Group B). This step is only
     needed to redo the chick demultiplexing; its result is stored in the
     deposited `chick_multiome_cAE_atlas.h5ad`, from which every chick figure
     table is rebuilt.
   - **E6/E9/E11 BAMs.** The ESR1/ESRRG pileup and intron/exon steps
     (`pipeline/scrna_e6e9e11/cr9_esrrg_pileup.sh`, `cr9_er_intron_exon_umi.sh`,
     `cr9_er_perembryo_ie.sh`) read possorted BAMs under
     `$DUCK_SCRATCH_USER/old_duck_review/bams/` (`D1_zju1_FC07361_possorted.bam`,
     `D2_zju1_FC07344_possorted.bam`, `D3_zju1_FC07344_possorted.bam`). GSE336329
     carries the FASTQs; recreate the BAMs with `cellranger count --create-bam`
     against the scRNA reference built by `build_scrna_duck_ref.sh` (deposited
     with the GEO reference files, not in this repository).
   - **E6/E9/E11 CellBender matrices.** `pipeline/scrna_e6e9e11/cr9_syrinx_pseudobulk_from_h5.py`
     (the growth-phase DESeq2 input for figS16 panel C) sums raw counts from
     per-library CellBender outputs,
     `$DUCK_SCRATCH_USER/old_duck_cr9/cellbender/<lib>_cellbender_filtered.h5`,
     which GEO does not carry. The integrated object
     `scrna_cr9_03_integrated_CB.rds` was built from these counts, so its RNA
     counts layer holds the same values; otherwise re-run CellBender on the
     cellranger output recreated from the FASTQs.

   The per-step input paths are also listed at the top of each script and in the
   figure `PROVENANCE.md` files.
1. **Start from the deposited atlases**, not from raw reads. The deposited
   `multiome_step02_preCB_clean.h5ad` and `multiome_step04_preCB_clustered.h5ad`
   are the pipeline's master objects; demultiplexing and clustering results are
   stored in their `.obs`. Re-running demultiplexing or clustering changes
   donor and cluster labels (see below), so the figure inputs are rebuilt from
   these objects.
2. **snRNA exports** from the step02/step04 atlases (stage-restricted
   subclustering, cluster validation, pseudobulk staging, sex calls, pseudobulk
   DESeq2).
3. **BAM-derived RNA tracks** (PITX2 isoform counts and pileups) from the
   deposited GEX BAMs and the deposited PITX2-extended GTF.
4. **Left/right inference from snRNA and genotype-to-embryo mapping**
   (`d2l_lr_inference.py` → `d2m_genotype_to_embryo.py`), run locally on the
   exported tables.
5. **ATAC.** The entry point is `pipeline/f1_coupled_atac.py`, run over the
   deposited fragment files, the step04 atlas and the committed subcluster
   annotation table. It writes `step_f1_coupled_atac_metadata.tsv`, which the
   ATAC pseudobulk step requires. That table is the one pipeline intermediate
   that is neither committed nor deposited, so this step must run first. Then
   use the deposited ArchR peak matrix and union peak set (or rebuild them with
   the `archr_*` scripts), followed by the fragment-counting (`f4`/`f5`/`f6`),
   pseudobulk and differential-accessibility steps.
   **After `f6`**, finish the left/right chain with the chromatin calls:
   `d2n_chromatin_lr_inference_hh1716.py` → `d2o_merge_multiome_inference.py`
   (writes `lr_inference_multiome.tsv`, the canonical per-embryo sides) →
   `d2m_multiome_genotype_to_embryo.py`.
6. **Chick control**, from the deposited chick atlas (`chick_multiome_cAE_atlas.h5ad`)
   and the chick fragment files staged in step 0; the chick ATAC QC files and the
   demultiplexing inputs are re-derived as described there.
7. **E6/E9/E11 scRNA**, from the deposited integrated object via
   `pipeline/scrna_e6e9e11/`.
8. Library QC metrics (`cellranger` `outs/`) are needed only for the QC panels
   and are re-derived from the deposited BAMs with
   `pipeline/rederive_cellranger_from_geo.sh` (about a day per ARC library).
9. **Cross-figure order.** Render figS02 panel D
   (`figures/fig02_supp/generate_dotplot_extended.py`) before fig02 panel B
   (`figures/fig02/generate_spatial_svg.py`), and figS03 panel H
   (`figures/figS03/panel_H_cis_memory_specificity.R`) before the fig03
   composite.

Where this summary and `docs/geo_tier2_build_order.md` differ, the build-order
document governs.

**Exact versus in-kind reproduction.** Deterministic steps reproduce exactly
with the pinned versions: ATAC coupling and QC, BAM and fragment counting, GC
tracks, reference building, `cellranger`/`cellranger-arc` counting, and the
seeded statistics. Other steps reproduce only in kind, because they are
stochastic or their labels are arbitrary: CellBender, scDblFinder,
cellSNP-lite/vireo demultiplexing (donor labels permute), the E6/E9/E11
substage relabeling, and any re-clustering (Leiden cluster IDs reshuffle; cell
type names flow only through the committed annotation tables
`references/cluster_annotations.tsv` and
`references/d2g_subcluster_annotations.tsv`). For these steps the deposited
objects and the committed tables are the record.

## Manuscript figures and repository directories

Repository figure directories are numbered independently of the paper. Main
figures 1–3, 6 and 8, the panels of Figure 7 other than 7A, and supplementary
figures S1–S6, S13 and S18–S20 are wet-lab or schematic and have no code here.

| Paper | Repository source |
|---|---|
| Figure 4 (snRNA atlas, dorsal/ventral subclustering, PITX2 two waves) | [`fig01`](figures/fig01/) panel C (atlas UMAP) + [`fig02`](figures/fig02/) |
| Figure 5 (PITX2c-promoter chromatin memory) | [`fig03`](figures/fig03/) + [`figS03`](figures/figS03/) panel H |
| Figure 7A (estrogen-receptor survey) | [`fig06`](figures/fig06/) panel A |
| Figure S7 (multiome atlas) | [`fig01`](figures/fig01/) |
| Figure S8 (demultiplexing and library QC) | [`fig01_supp`](figures/fig01_supp/) |
| Figure S9 (dorsal/ventral subclusters, left/right inference) | [`fig02_supp`](figures/fig02_supp/) |
| Figure S10 (snATAC QC, chromatin left/right, depth) | [`figS03`](figures/figS03/) |
| Figure S11 (duck vs chick chromatin) | [`fig04`](figures/fig04/) |
| Figure S12 (chick multiome control) | [`figS04`](figures/figS04/) |
| Figure S14 (E6/E9/E11 scRNA overview, ESR1) | [`figS07`](figures/figS07/) |
| Figure S15 (WNT/BMP/TGF-β screen) | [`figS08`](figures/figS08/) + [`figS09`](figures/figS09/) + [`figS10`](figures/figS10/) |
| Figure S16 (proliferation) | [`figS16`](figures/figS16/) |
| Figure S17 (ESR1 left/right) | [`figS17`](figures/figS17/) |

[`figS04_supp`](figures/figS04_supp/) holds chick multiome QC that supports
Figure S12 and is not a separate figure in the paper. Figure panels in the paper
were arranged for print from these repository figures; the arrangement step is
not part of this repository, so a repository composite may differ from the
printed figure in layout and typography; the data shown are the same.

## Software

| Tool | Version | Where |
|---|---|---|
| cellranger | 7.1.0 (multiome GEX); 9.0.1 (E6/E9/E11 scRNA) | `pipeline/a1_cellranger.sh` (7.1.0); 9.0.1 per the paper's Methods (run outside this repository) |
| cellranger-arc | 2.1.0 | `pipeline/rerun_arc_v210.sh` |
| CellBender | 0.3.2 | `envs/o2_cellbender.yml`; build with `scripts/cellbender/cb_build_env.sh` (applies the checkpoint patch) |
| cellSNP-lite | 1.2.3 | `envs/o2_vireo.yml` |
| vireo (vireoSNP) | 0.5.9 | `envs/o2_vireo.yml` |
| scDblFinder | 1.20.2 | R library, `envs/o2_R442_installed_packages.tsv` |
| scanpy | 1.10.3 | snRNA analysis; floor pinned in `pipeline/environment.yml` |
| leidenalg / harmonypy | 0.11.0 / 0.2.0 | clustering |
| ArchR | 1.0.3 | `envs/o2_R442_installed_packages.tsv` |
| MACS3 | 3.0.4 | peak calling, `envs/o2_snapatac2.yml` |
| SnapATAC2 | 2.9.0 | `envs/o2_snapatac2.yml` |
| DESeq2 (R) | 1.46.0 | snRNA pseudobulk differential expression, `envs/o2_R442_installed_packages.tsv` |
| PyDESeq2 | 0.5.4 | differential accessibility, `envs/o2_snapatac2.yml` |
| Seurat | 5.4.0 | E6/E9/E11 scRNA integration |
| Python / R | 3.12 / 4.4.2 | |

[`envs/`](envs/) holds the captured environments: `duckfig_portable.yml` and `duckfig.yml` (figures), the
cluster conda environments (`o2_*.yml`), and the cluster R library
(`o2_R442_sessionInfo.txt`, `o2_R442_installed_packages.tsv`; install recipe in
`pipeline/install_o2_r_packages.R`). [`envs/README.md`](envs/README.md) maps
each environment to the pipeline steps that use it.

## Repository layout

| Path | Contents |
|---|---|
| `pipeline/` | snRNA/snATAC pipeline, Python and R; `pipeline/chick/` is the chick control, `pipeline/scrna_e6e9e11/` the E6/E9/E11 scRNA |
| `FC_08656/` | per-flowcell processing scripts that the pipeline calls (genotype map, scDblFinder, demultiplexing QC, marker dotplot export) and the cellranger-arc library sheets |
| `figures/figNN/`, `figures/figSNN/` | one directory per repository figure: panel scripts, panels, `LEGEND.md`, `PROVENANCE.md` |
| `figures/<fig>.pdf` | composite figures |
| `figures/data/` | the tables every figure renders from (see its `README.md`) |
| `figures/_shared/` | shared plotting code (style, colors, dotplot renderer, compositor) |
| `references/` | annotation tables (cluster and subcluster annotations, genotype metadata, color registry), interval BEDs, gene lists; `references/README.md` lists the cited papers by DOI |
| `design/` | experimental design: embryo-to-library and left/right sample layout |
| `scripts/` | helper scripts: path roots, gates (`verify_*.py`), table builders and assemblers named in the `PROVENANCE.md` files |
| `envs/` | captured software environments |
| `docs/` | methods background (PITX2 locus annotation, LPM subcluster labels, the HH17/16 chromatin left/right rule), the path-root rules and the Tier 2 build order |

Cell-type names in figures come from two committed tables:
`references/cluster_annotations.tsv` (global clustering) and
`references/d2g_subcluster_annotations.tsv` (stage-restricted mesenchyme
subclusters). Figure code reads labels from these tables, never from numeric
cluster IDs.

## License

MIT; see [`LICENSE`](LICENSE).

Third-party code: `scripts/cellbender/cellbender_checkpoint_v0.3.2_patched.py` is
a modified copy of `checkpoint.py` from
[CellBender](https://github.com/broadinstitute/CellBender) v0.3.2
(Copyright (c) 2019, Broad Institute, Inc.). It is distributed under CellBender's
BSD 3-Clause license, not under the MIT license; the full license text is in the
file's header.

## Contact

Questions about the code or data: tabinweb@genetics.med.harvard.edu (Tabin lab,
Harvard Medical School).
