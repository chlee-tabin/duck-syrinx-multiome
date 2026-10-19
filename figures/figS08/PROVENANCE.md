# Figure S15, page 1 of 3 — provenance & implementation notes (repo `figures/figS08/`)

> Provenance and build notes for `figures/figS08/`. **Not part of the publication figure legend** (see `LEGEND.md`). Holds data-source ownership, build steps, and methodology trimmed from the legend. Built from the cr9 re-derivation of the E6/E9/E11 scRNA-seq WNT/BMP screen.

## Data sources
- `figures/data/scrna_e6e9e11/cr9_screen_E{06,09,11}_percell.tsv.gz` — per-cell export from the cr9 substage objects (CellBender-corrected, genotype-carry-forward E6/E9/E11 scRNA-seq). Carries `log__`/`cp10k__` for the full S8/S9/S10 gene panels (43 genes), `population` (joined on-cluster, AUTHORITATIVE — see below), `LR.called`, `stage`, `sex.called`, `donor_id`, `RNA_snn_res.0.{2,3}`.
- `figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv` (the byte-identical duplicate `cr9_screen_clusterlabels.tsv` was dropped 2026-08-26) — the **run-matched** compendium clusterlabels the export's `population` was joined against (provenance/verification). Distinct from the repo's `substage_compendium_clusterlabels.tsv` (a different substage run) — do NOT re-join.
- `figures/figS08/panel_C_ish_s11_male.png` — rebuilt by `make_panel_C_strip.py` from the paper's in situ source PDF (`240611_InSitu_S15.pdf`, not distributed; the authors' 2026-08-26 original of the colorimetric in situ figure — wet-lab imaging, not distributed with the public code; the built PNG ships instead; old numbering S11/S15 = current S18): male block, E9–E13, panels at native 591 ppi, gene label reads WNT5B (nomenclature). Replaces the earlier ~130 ppi crop of the compiled figure.
- `references/color_registry.tsv` — subpopulation colours (`subcluster` level; rows for Chondrogenic LPM (renamed 2026-08-22 → "Syringeal chondrocyte", its own registry row; the old row is kept for the chick arm) / Syringeal fibroblast / Syringeal smooth muscle / Mesothelium), Left = `#4393C3` (blue), Right = `#D6604D` (red), via `figures/_shared/colors.R`.

## Build / regeneration
Cluster export (one-off, on scratch substage RDS):
```
# scripts in pipeline/scrna_e6e9e11/
sbatch run_cr9_screen_export.sh            # -> cr9_screen_E{06,09,11}_percell.tsv.gz on scratch
```
Local panels + composite (all in the `envs/duckfig_portable.yml` env):
```
Rscript figures/figS08/generate_panel_A_wnt_violins.R     # panel A (subpopulation × stage)
Rscript figures/figS08/generate_panel_B_wnt_lr.R          # panel B (per-embryo paired L/R)
Rscript figures/figS08/generate_panel_D_ish_dotplot.R     # panel D (correspondence dotplot)
python3 figures/figS08/make_panel_C_strip.py              # panel C strip from the in situ source PDF (not distributed; built PNG ships)
Rscript figures/figS08/generate_default_pooled.R          # FALLBACK panels (see below)
conda run -n duckfig python3 figures/figS08/generate.py   # -> figures/figS08.pdf
```
Shared helpers: `cr9_screen_io.R` (loader + `screen_long`), `tier2_lr.R` (L/R engine; also used by figS10).

## Layout — 2×2 grid, half-page panels (2026-08-24)
`generate.py` composes **A|B on the top row, C|D on the bottom row** (was a 4-row single column). The single column made a 7.00 × 21.37 in page: fitted onto a portrait manuscript page it shrank ~2.2×, which is what made the fonts unreadable. The grid gives 7.00 × 9.03 in (aspect 0.775), so it prints at ~1:1.
- Panels A, B and D are therefore rendered **portrait at half-page width** by their producers (A `w=5.75`, B `width=4.5`, D `5.9 × 6.4` in) with `base_size` unchanged — the fonts grow relative to the panel instead of being scaled down by the compositor.
- Panel D was **transposed** for this: gene on x, subpopulation on y, stage as stacked row blocks (`facet_grid(stage ~ ., scales="free_y", space="free_y")`, `scale_y_discrete(limits=rev)` to keep `POP_ORDER` reading top-to-bottom). Its old page-wide 9 × 3.6 in shape (aspect 2.5) could not pair with the portrait ISH montage without either crushing panel C to ~1.5 in or leaving a ~4 × 2.7 in hole.
- Per-panel slot widths in `LAYOUT` (px @300 dpi): A 1150, B 900 / C 875, D 1175 — each row sums to `FULL_W` 2100 px with the 50 px (12 pt) gap, i.e. **100 % width utilization on both rows**, and every slot is ≤ its panel's native pixel width so nothing is upscaled. `generate.py` raises on either invariant being broken by a re-rendered panel.
- `legend_box="vertical"` is passed to `tier2_lr_boost` for panel B: at 4.5 in the side-by-side Side + cells/embryo guides ran off the canvas. The argument defaults to `"horizontal"`, so figS10's page-wide panel B is unaffected.
- Panel letters are drawn by the compositor and are **A–D**; manuscript Figure S15 letters continue as E (figS09) and F–I (figS10).

## Default (faithful) fallback
`generate_default_pooled.R` writes `panel_A_wnt_violins_default_pooled.png` (+ the figS09/figS10 equivalents): the manuscript layout — **all** panel genes, pooled mesoderm (all non-technical cells; Contamination + Cycling dropped), x = stage. If the subpopulation-resolved upgrade is not adopted, swap this into the panel-A slot of `generate.py`. Kept committed as the drop-in.

## Methodology trimmed from the legend
- **Population join is on-cluster; `population` in the export is AUTHORITATIVE.** The scratch substage `*_obj.rds` (E09 = 8,736 cells) is a different run than the repo's `substage_compendium_clusterlabels.tsv` (E09 = 10,160) — Leiden IDs and labels drift between runs. `cr9_screen_export.R` joins `population` against the **run-matched** clusterlabels ON THE CLUSTER (committed as `figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv`). The loader reads `population` directly and never re-joins, so the figures are reproducible from the committed export; **never re-join these exports to the repo's `substage_compendium_clusterlabels.tsv` decoy.** The loader drops `Contamination` and `Cycling` (technical, non-lineage) and renames any legacy `"Somatic"` → `Skeletogenic LPM` (a no-op compat guard — the current export already says `Skeletogenic LPM`).
- **Panel A expressed-gene filter:** a WNT gene is shown only if its max subpopulation-mean (over stage × population) > 0.1, so the figure is not a stack of flat rows; the omitted genes are logged (WNT1/3/3A/4/6/7A/7B/9A/9B/16) and retained in full in panel D / the default-pooled fallback.
- **Panel B laterality method:** per-cell CP10K averaged within each (`donor_id` × side) over the syrinx-fated pool (Ventral spl. + Syringeal chondrocyte + Syringeal fibroblast + Skeletogenic LPM — at the tested stages E09/E11 the E06-only members Ventral spl. and Skeletogenic LPM contribute no cells, so the effective pool is Chondrogenic LPM + Syringeal fibroblast); a donor-side needs ≥ 10 cells (stable mean). Pairing is by `donor_id` (the vireo individual embryo) — 8 embryos/stage, each bisected with both halves captured; paired Wilcoxon signed-rank (n = 7–8 after the cell floor). E06 has no real dissection sides (whole tissue) so it is excluded from L/R. The claim is restricted to per-cell expression *level*: cell number / fraction is not comparable across the L/R split because left and right are separately-sampled, separately-normalized libraries.
- **Panel B stats (paired Wilcoxon, L−R):** WNT5A E09 p=0.35 / E11 p=0.36 (ns); SFRP2 E09 p=0.052 / **E11 p=0.014\*** (L>R); SHISA2 E09 p=0.15 / E11 p=0.44 (ns); PITX2 **E09 p=0.022\*** (L>R) / E11 p=1 (ns, expected Wave-1 decline by the old-data stages).
- **Panel C** carries a sliver of the S11 grid as source attribution; the original caption is cropped out.

## Export wrapper notes
The per-cell TSVs were md5-verified after transfer from the cluster. Export wrapper fixes: `export R_LIBS_USER=$HOME/R/4.4.2` (so `Rscript --vanilla` finds Seurat v5) and `exit $RC` (so SLURM reflects the real Rscript exit; a trailing `echo` had masked an early failure as COMPLETED).

