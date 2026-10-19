# Figures

Code, panels and composites for every computational figure of the paper. Each panel
renders from the tables committed under [`data/`](data/) (provenance:
[`data/README.md`](data/README.md)); no h5ad or sequencing data is needed.

## Repository figure directories and the paper

Repository directories are numbered independently of the paper.

| Directory | Paper figure | Contents |
|---|---|---|
| `fig01/` | Figure S7 (and Figure 4A, atlas UMAP) | experimental design, pipeline validation, snRNA cell-type atlas, marker dotplot, pseudobulk staging |
| `fig01_supp/` | Figure S8 | library QC, per-genotype QC, clustering quality, sex and genotype-distance validation, demultiplexing |
| `fig02/` | Figure 4 | dorsal/ventral LPM subclustering, PITX2 isoform laterality and pileups |
| `fig02_supp/` | Figure S9 | sampling balance, per-subcluster laterality, left/right inference, extended dotplot |
| `fig03/` | Figure 5 | PITX2c-promoter chromatin memory (snRNA + snATAC) |
| `figS03/` | Figure S10 (and Figure 5, panel H) | snATAC QC, chromatin left/right resolution, depth controls, trough tests |
| `fig04/` | Figure S11 | duck vs chick chromatin and motif comparison |
| `figS04/` | Figure S12 | chick multiome control |
| `figS04_supp/` | supports Figure S12 | chick multiome QC |
| `fig06/` | Figure 7A | estrogen-receptor survey (E6/E9/E11 scRNA) |
| `figS07/` | Figure S14 | E6/E9/E11 scRNA overview, ESR1 |
| `figS08/`, `figS09/`, `figS10/` | Figure S15 (pages 1–3) | WNT/BMP/TGF-β screen |
| `figS16/` | Figure S16 | proliferation |
| `figS17/` | Figure S17 | ESR1 left/right |

## Layout and conventions

```
figures/
├── figNN.pdf, figSNN.pdf   # composites (one per repository figure)
├── figNN/, figSNN/         # panel scripts, panels, LEGEND.md, PROVENANCE.md
├── data/                   # the tables every panel renders from (see data/README.md)
└── _shared/                # shared plotting code
```

- Panels are `panel_<letter>_<description>.png` (data) or `.svg` (schematics); the
  compositor `generate.py` in each directory writes the PDF composite at `figures/`
  (`fig06/` is a single panel built by `plot_fig6a_er_survey.py`, with no composite).
- `LEGEND.md` is the figure legend; `PROVENANCE.md` lists, panel by panel, the data files
  read and the commands that build them.
- Cell-type names come from `references/cluster_annotations.tsv` (global clustering) and
  `references/d2g_subcluster_annotations.tsv` (stage-restricted LPM subclusters), read through
  `_shared/style.py` and `_shared/subcluster_map.{py,R}`; figure code never maps numeric
  cluster IDs by hand.
- Genotype metadata (display names, groups, stages) is loaded from
  `references/genotype_metadata.tsv` by `_shared/style.py` at import time.
- A missing colour, label or data key raises; figure code has no silent defaults.

## Shared code (`_shared/`)

| File | Purpose |
|---|---|
| `style.py`, `colors.py`, `colors.R` | palettes, genotype metadata, annotation loaders |
| `subcluster_map.py`, `subcluster_map.R` | stage-restricted subcluster labels and laterality pools |
| `data_io.py` | path constants and TSV loading |
| `dotplot.py` | standalone dotplot renderer (reads exported TSVs) |
| `composite.py`, `layout.py` | PDF compositor and layout helpers |
| `lr_resolution_panel.py` | per-embryo left/right resolution matrix |
| `library_naming.py`, `library_naming.R` | library display names and stage anchors |
| `markers.py`, `paired_standardised.py`, `pitx2_gene_model.R` | marker lists, paired statistics, PITX2 gene model |
