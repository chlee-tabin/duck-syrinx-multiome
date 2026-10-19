# Figure 7, panel A — provenance & implementation notes (repo `figures/fig06/`)

> Provenance and build notes for `figures/fig06/`. **Not part of the publication figure legend** (see `LEGEND.md`). Holds data-source ownership, build steps, and methodology trimmed from the legend. Built from the cr9 re-derivation of the E6/E9/E11 scRNA-seq dataset; independent of the multiome (FC_08656) pipeline.

## Data sources
- `cr9_03_integrated_CB.rds` — cr9-CB master object (24,816 cells): the CellBender-corrected, intron-inclusive (pre-mRNA-inclusive) re-derivation of the old E6/E9/E11 duck 3′ scRNA-seq. Source of the per-cell 6-ER log-normalized expression, global cluster, and substage population.
- `figures/data/scrna_e6e9e11/cr9_er_survey_percell.tsv.gz` — per-cell export from the cr9-CB master. Carries `log__{ESR1,ESR2,ESRRB,ESRRG,GPER1,RERG}`, `global_cluster`, `substage_population`, `stage`, and a per-cell `broad_class`. This is the shared old-duck data home (shared with the cr6-vs-cr9 QC figures and the cluster-pull target).
- `figures/data/scrna_e6e9e11/cr9_global_cluster_markers.tsv` — refreshed per-cluster broad-lineage marker table written by the same exporter; the plotter does not read it (it requires the per-cell `broad_class`).
- `figures/fig06/panel_A_er_survey.png` — the rendered data panel (raster).

## Build / regeneration
Per-cell export (one-off, on the cluster, R/4.4.2 against the scratch cr9-CB master):
```
# scripts in pipeline/scrna_e6e9e11/
Rscript cr9_export_er_survey.R        # -> cr9_er_survey_percell.tsv.gz (+ refreshed cr9_global_cluster_markers.tsv)
```
Local panel:
```
FIG6A_PALETTE=canonical python3 figures/fig06/plot_fig6a_er_survey.py   # -> panel_A_er_survey.png
```
- `FIG6A_PALETTE` env var selects the violin palette: `canonical` (default; the repo fig02 earth-tone subcluster palette, with muted colours for the old-duck-only populations) or `rainbow` (ggplot2/Seurat `hue_pal()`, to match the collaborators' Seurat `VlnPlot` panels in Fig 6). Colours fail loud (no grey default) on an unknown population.
- **No composite is produced** (single subpanel). The directory holds only the scRNA-seq subpanel (Fig 7a) of manuscript Figure 7; the wet-lab panels (7b, 7c) are not reconstructed here.

## Methodology trimmed from the legend
- **Substrate / per-cell `broad_class`.** The non-mesoderm categories use a per-cell `broad_class` computed as the argmax of a broad lineage marker-panel on the *current* cr9-CB master, so they never depend on a stale per-cluster markers table. Mesenchymal cells carry a `substage_population`; non-mesoderm cells (`substage_population` NA) are labelled by `broad_class` (`Epithelial` → Endoderm, `Endothelial` → Endothelium, `Immune` → Immune). If the export lacks `broad_class`, the plotter stops with an error.
- **Stage-blocking and the cell floor.** Each named LPM subpopulation is shown per stage only when it has ≥ 100 cells in that stage (`MIN_CELLS = 100`); dropped sub-threshold populations are logged. `Contamination` is never a survey category (`DROP_POPS`). Non-mesoderm lineages are pooled across stages, in fixed order (Endoderm, Endothelium, Immune).
- **Non-mesoderm purity.** Each of the three non-mesoderm reference lineages is a single, well-separated global cluster: Endoderm (the epithelial cluster, 85 % epithelial), Endothelium (87 %), Immune (81 %). (Purities are marker-based calls on the cr9-CB master, not recomputable from the committed per-cell export.)
- **The `ESRRG*` asterisk.** *ESRRG* reads broad only as an artefact of intron-inclusive counting: predominantly intronic UMIs of a 404 kb gene inflate its pre-mRNA signal. It is broad and non-lateralized — see the ESRRG-mitigation supplement (Figure S14). This is marked only on the *ESRRG\** title; the caveat itself is in the figure legend (the in-panel footnote was removed 2026-08-24).
- **Violin construction.** Per-violin median bar (black); a degenerate (e.g. all-zero) distribution is drawn as a single median tick. y-axis starts at 0; per-panel y-max is the 99.5th percentile of the plotted values (×1.08, floor 0.5).

## Related docs
The ESRRG intron-inclusive pre-mRNA inflation is documented in the ESRRG-mitigation supplement (Figure S14).

## References
- **CellBender** — Fleming SJ, Chaffin MD, Arduini A, Akkad A-D, Banks E, Marioni JC, Philippakis AA, Ellinor PT, Babadi M. *Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender.* Nat Methods. 2023;20(9):1323–1335. doi:10.1038/s41592-023-01943-7

## Print sizing (PNAS final, 2026-10, font pass)

The panel is authored at about its printed size: `figsize=(7.1, 4.9)` in at 300 dpi
(tight bbox ~506 pt wide), Arial, tick labels 7.5 pt, axis label 8.5 pt, gene titles
9.5 pt italic, stage tags 7.5 pt (placed above the axes). The Fig 7 leaf
(`manuscript/m17_fig7_esr1_v64.py`, not distributed) places it at 756 pt and the 17.8 cm print fit is
504.6/756, so printed size ~= native size (factor ~1.0; was 0.43 at 170 dpi / 16.5 in).
