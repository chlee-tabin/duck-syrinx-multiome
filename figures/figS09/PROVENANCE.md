# Figure S15, page 2 of 3 — provenance & implementation notes (repo `figures/figS09/`)

> Provenance and build notes for `figures/figS09/`. **Not part of the publication figure legend** (see `LEGEND.md`). Built from the cr9 re-derivation of the E6/E9/E11 scRNA-seq dataset.

## Data sources
- `figures/data/scrna_e6e9e11/cr9_screen_E{06,09,11}_percell.tsv.gz` — see `figures/figS08/PROVENANCE.md` for the full schema and the on-cluster `population` join.

## Build / regeneration
```
Rscript figures/figS09/generate_panel_A_smad_violins.R     # panel A (subpopulation × stage, nested arms)
Rscript figures/figS08/generate_default_pooled.R           # also writes the figS09 default-pooled fallback
conda run -n duckfig python3 figures/figS09/generate.py    # -> figures/figS09.pdf
```
Loader/helpers shared from `figures/figS08/cr9_screen_io.R`.

## Displayed panel letter = "E" (2026-08-24)
Manuscript Figure S15 spans three booklet pages (figS08 = page 1, this = page 2, figS10 = page 3) and its panel letters run **continuously**, so this page's single panel is lettered **E**, not A. The letter is a compositor overlay: it lives in the `labels` entry of `LAYOUT` in `generate.py` (`["E"]`) and is **not** baked into the panel by `generate_panel_A_smad_violins.R`. The panel PNG keeps its repo-local name `panel_A_smad_violins.png` — repo filenames track the local slot, only the rendered letter is manuscript-facing. Do not "fix" the mismatch by renaming the PNG.

## Default (faithful) fallback
`panel_A_smad_violins_default_pooled.png` (written by `figS08/generate_default_pooled.R`): the manuscript S9 layout — SMAD1/5/9 + SMAD2/3 pooled mesoderm × stage. Swap into the panel-A slot of `generate.py` if the subpopulation-resolved upgrade is not adopted.

## Methodology trimmed from the legend
- Arm grouping (BMP: SMAD1/5/9 ; TGF-β: SMAD2/3) is rendered with `ggh4x::facet_nested`; a plain `facet_grid(gene ~ stage)` fallback (gene order BMP-arm first) is in the script if ggh4x is unavailable.
- Same population set / `Skeletogenic LPM` rename / `Contamination` drop as figS08. No L/R panel: the SMADs are intracellular effectors and the laterality story is carried by the ligands (figS08/figS10 panel B); a SMAD L/R panel was not built.

## Export wrapper notes
See `figures/figS08/PROVENANCE.md` (`R_LIBS_USER`/`exit $RC` export-wrapper fixes).
