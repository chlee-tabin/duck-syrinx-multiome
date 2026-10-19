# Figure S15, page 3 of 3 — provenance & implementation notes (repo `figures/figS10/`)

> Provenance and build notes for `figures/figS10/`. **Not part of the publication figure legend** (see `LEGEND.md`). Built from the cr9 re-derivation of the E6/E9/E11 scRNA-seq dataset (the screen figure cited inline: "BMP4 and BMP5 … high … at this later stage").

## Data sources
- `figures/data/scrna_e6e9e11/cr9_screen_E{06,09,11}_percell.tsv.gz` — see `figures/figS08/PROVENANCE.md` for the schema and the on-cluster `population` join.

## Build / regeneration
```
Rscript figures/figS10/generate_panel_A_bmp_violins.R        # panel A = letter F (two gene blocks; uses ../figS08/cr9_screen_io.R)
Rscript figures/figS10/generate_panel_B_bmp_lr.R             # panel B = letter G (per-embryo paired L/R; uses ../figS08/tier2_lr.R)
Rscript figures/figS10/generate_panel_C_bmp_receptors.R      # panel C = letter H (transposed receptor dotplot; also prints the full number table)
Rscript figures/figS08/generate_default_pooled.R             # also writes the figS10 default-pooled fallback
conda run -n duckfig python figures/figS10/plot_screen_onset.py   # panel_D_screen_onset.png = letter I
conda run -n duckfig python figures/figS10/generate.py   # -> figures/figS10.pdf
```
`Rscript` here is the R of the figure environment (`conda run -n duckfig Rscript …`).

## Page layout + lettering
Panel A is wide and short; B and C sit left/right on the bottom row, with C transposed.
Panel letters run continuously across manuscript S15's three pages.

- **Displayed letters are F (ligands), G (paired L/R), H (receptors), I (onset trajectories,
  added 2026-08-27)** — page 3 of manuscript S15 continues the A–D (page 1 = figS08) / E
  (page 2 = figS09) sequence. The letters are drawn by the compositor overlay from
  `generate.py`'s `LAYOUT`; the panel **PNG filenames keep their historical `panel_A/B/C/D_*`
  names** (renaming them would churn every consumer for nothing).
- **Grid**: F full page width on the top row; G (left) | H (right) top-aligned on the middle
  row; I full page width on the bottom row. Row 2's slot widths are apportioned by panel
  *aspect*, so both middle panels land at the same displayed height and the row fills the page
  exactly; a slot more than 1.25× a panel's native pixel width raises (that is the < 240 dpi
  blur threshold).
- **Panel aspects** (w/h): F 2.34 (was 0.85), G 0.81 (was 1.01), H 0.59 (was 2.50, transposed),
  I 3.52. With panel I the page is 504 × 732 pt = 7.00 × 10.17 in, aspect 0.69 — the booklet
  leaf fit (manuscript booklet assembly, not distributed; `min(564/cw, 720/ch)`) is now HEIGHT-bound at
  s ≈ 0.98 (the pre-I page was width-bound at s ≈ 1.12 and filled the leaf); panels print
  ~12% smaller than the pre-I layout.
- **Render inches ARE print size.** The compositor maps a panel's ggsave `width`/`height` onto
  its page slot, so effective print pt = nominal pt × (slot inches / render inches). Two
  consequences that drove the numbers above: (i) panel A could not simply be given a 14-in-wide
  canvas — squeezed into a 7-in page slot that halves every font; (ii) panels B and C are
  rendered at *half-page* physical width (4.8 in / 3.5 in) because they now occupy half a row.
  `dpi` is orthogonal — it only sets raster density (A at 240, C at 250 keeps both ≥ 300 dpi in
  their slots; B's engine is fixed at 200 dpi → ~243 dpi in slot, matching figS08's panel B).
- **Panel A restructure**: the 9 ligands are split into two side-by-side blocks (BMP1–BMP5 |
  BMP6–BMP15) rather than one 9-row stack. A single stack cannot reach the requested ~2.0–2.5
  aspect at legible font size: the angled subpopulation-label band is a fixed ~1.1 in cost, so
  9 rows in a 3.8 in canvas leaves ~0.3 in per gene row and the per-row y ticks collide.
  Blocks reuse the shared `figS08::screen_violin()` (single-sourced theme + canonical palette);
  since that helper renders to file, each block is captured with `last_plot()` and recombined by
  patchwork. Post-hoc on the captured objects: y breaks capped at ~3 per row (the helper's
  default asked for six on a 0–2.5 row), the right block's duplicate y-axis title dropped, and
  the inter-block gutter tightened. Expansion stays `c(0, 0.05)` (no margin below zero).
- **Panel B** is unchanged analytically — same genes, pool, pairing and statistics; only its
  render size changed (`width`/`height`/`legend_box="vertical"` passed to the shared
  `figS08/tier2_lr.R` engine, which is *not* edited from here).
- **Panel C transpose**: genes on x (BMPR1A/BMPR1B/BMPR2 | BMP4), subpopulations on y,
  `facet_grid(stage ~ ., space = "free_y")` so each stage block carries only its own
  subpopulations; the receptor/ligand separator became a `geom_vline`; guides stacked with
  `legend.location = "plot"` (at 3.5 in the panel starts ~1.5 in in, so panel-justified guides
  would overflow). Encoding, scale limits, aggregation and fail-loud guards are untouched — the
  numbers below are identical to the pre-transpose panel.
- Consumers of these panel PNGs outside this directory: none. The manuscript S15
  (page 3 of 3) uses this figure's own composite `figures/figS10.pdf`.

## Default (faithful) fallback
`panel_A_bmp_violins_default_pooled.png` (written by `figS08/generate_default_pooled.R`): manuscript S10 layout — the 9-gene BMP-ligand panel (BMP1–BMP15) pooled mesoderm × stage. Swap into the first `LAYOUT` row of `generate.py` if the subpopulation-resolved upgrade is not adopted (it is portrait, so the page will be taller again).

## Methodology trimmed from the legend
- All 9 BMP-panel genes retained (no expressed-gene filter) so the BMP4/BMP5 specificity reads against the flat remainder.
- Panel B reuses the shared `figures/figS08/tier2_lr.R` engine: per-cell CP10K averaged within (`donor_id` × side) over the syrinx-fated pool (effectively Syringeal chondrocyte + Syringeal fibroblast at E09/E11 — the pool's other members are E06-only), ≥ 10 cells/donor-side, paired by individual embryo (8/stage), paired Wilcoxon. Per-cell expression *level* only (not cell number/fraction — confounded across the separately-normalized L/R libraries).
- **Panel B stats (paired Wilcoxon, L−R):** BMP4 E09 p=0.076 / E11 p=0.53 (ns); BMP5 E09 p=0.93 / E11 p=0.14 (ns); BMP7 **E09 p=0.022\*** (L>R) / E11 p=0.94 (ns); PITX2 **E09 p=0.022\*** (control) / E11 p=1. So among the BMP ligands the per-cell L/R boost is BMP7-at-E09 specific; BMP4/5 (the abundant E06 ligands) are bilateral in the syrinx-fated pool. Note BMP7 is the in situ–validated ligand (Fig S11 / Fig 5b), while BMP4/5 are the abundance standouts of panel A — the two readouts are complementary.

## Panel C — BMP receptors
Descriptive re-plot of the source/responder analysis
(`bmp_source_responder.py`, not distributed): per
(stage, subpopulation) percent-expressing (CP10K > 0) and mean CP10K, no filtering, no
normalization beyond the export's CP10K. Reproduces that script's numbers exactly; adds E06
Mesothelium (n=339), which the script's `POPS` list omitted, so panel C's subpopulation set
matches panel A's. All subpopulations present at a stage are shown, including the two small ones
(E09 Syringeal smooth muscle n=66, E11 n=116) — kept so the subpopulation vocabulary matches
panel A (they are y-axis rows since the 2026-08-24 transpose) rather than dropped; the load-bearing claims are all on the 1.0–6.2k-cell populations.

**Key numbers** (percent expressing / mean CP10K):

| | E06 | E09 | E11 |
|---|---|---|---|
| *BMPR1A* | 52.8–74.9% / 1.45–1.89 | 48.5–60.4% / 1.13–2.10 | 40.5–52.0% / 1.09–1.86 |
| *BMPR2* | 61.6–85.7% / 2.10–3.88 | 65.9–77.3% / 3.12–3.59 | 63.4–73.3% / 2.80–2.93 |
| *BMPR1B* Syringeal chondrocyte | — | **49.9% / 4.35** | **73.1% / 8.05** |
| *BMPR1B* Syringeal fibroblast | — | 27.0% / 1.13 | 27.5% / 1.78 |
| *BMPR1B* Syringeal smooth muscle | — | 25.8% / 0.58 | 6.9% / 0.15 |
| *BMPR1B* (E06, all pops) | 36.9–65.5% / 1.18–2.42 | — | — |
| *BMP4* chondrocyte vs fibroblast | — | 14.5% / 0.46 vs **33.9% / 1.09** | 6.4% / 0.14 vs 19.7% / 0.65 |

Per-stage cell counts: E06 Skeletogenic LPM 1221, Ventral spl. 1110, Dorsal sub-epi. 1015,
Mesothelium 339; E09 chondrocyte 4346, fibroblast 2958, Dorsal sub-epi. 1366, smooth muscle 66;
E11 fibroblast 6202, chondrocyte 3701, smooth muscle 116.

**Encoding choices.** Greyscale fill on ONE shared absolute CP10K scale (0–8.05), not per-gene
standardized: the whole point is that *BMPR1B*-in-chondrocyte is the only dark cell on the
panel, which a per-row renormalization would erase by making every receptor row look
differentiated. `scale_size_area` so dot area is proportional to percent-expressing with 0 → 0.
No colour is spent (nothing here encodes a biological variable). Fail-loud guards: missing
`cp10k__` column, an incomplete (stage × population) × gene grid, an NA aggregate, and a gene
silently dropped between export and aggregate all `stop()`.

### Caveat: "bilateral" is NOT established, and the panel does not claim it
The panel pools left and right halves, and a separate check says the bilateral wording of the
original decision would have been an overstatement. Running the panel-B engine
(`figS08/tier2_lr.R`, per-embryo paired per-cell CP10K over the syrinx-fated pool, paired
Wilcoxon) on the receptors gives, as median L−R: *BMPR1A* E09 **−0.30, p=0.022\*** / E11 +0.05
p=0.36 ns; *BMPR1B* E09 −2.42 p=0.076 ns / E11 **−1.59, p=0.042\*** ; *BMPR2* ns at both
(p=0.27, p=0.23); *PITX2* control +0.36 p=0.022\* at E09. So two receptor terms carry a modest
**right**-leaning bias — the opposite direction to the left-biased PITX2 response — rather than
being symmetric. Two reasons not to read that as receptor-driven laterality, and not to put it
in the legend: (i) the pooled per-cell mean over a mixed pool is confounded by
chondrocyte:fibroblast composition across the halves, and *BMPR1B* is exactly the
composition-sensitive gene (8.05 vs 1.78 CP10K between those two populations), so a pooled
*BMPR1B* shift can be pure composition; (ii) a right-leaning receptor bias cannot explain a
left-biased response in any case. The argument the panel actually makes is the population one —
reception is broad, so response specificity is not gated by receptor availability. A claim about
laterality of reception would need a per-population (not pooled) paired
test, which is not in this figure.

## Export wrapper notes
See `figures/figS08/PROVENANCE.md` (`R_LIBS_USER`/`exit $RC` export-wrapper fixes).

### Panel I — screen-gene onset trajectories (added 2026-08-27)
- Producer: `figures/figS10/plot_screen_onset.py` → `panel_D_screen_onset.png`
  (drawn letter I; construction patterned on S14D = `figures/figS07/plot_esr1_onset.py`,
  with ONE deliberate divergence: the gray aggregate track excludes the
  Contamination label at every stage that carries it — HH13 7.0%, HH17/16 5.7%
  of cells — whereas S14D keeps those cells in its aggregate; the two panels'
  identically-labeled gray tracks therefore differ slightly at HH13/HH17-16).
- Multiome side: `figures/data/d2g_signaling_{stage}.tsv.gz` (gene list extended
  with BMP7/WNT5B/SFRP2/SHISA2 via pipeline phase `d2gs`, 2026-08-27) joined to
  `d2g_lpm_subclustering_{stage}_ccreg.tsv.gz` on (genotype, library, barcode);
  positivity = `<gene>_cp10k > 0`; unit = genotype, >=10 cells.
- Old-scRNA side: `figures/data/scrna_e6e9e11/cr9_screen_{E06,E09,E11}_percell.tsv.gz`;
  positivity = `log__<gene> > 0`; unit = donor_id.
- The S14D whole-LPM baseline runs only through E06 here: the cr9 SCREEN
  exports carry the full LPM compartment at E06 but only the syrinx-region
  populations (plus Dorsal sub-epi. at E09) at E09/E11, so a whole-LPM unit
  fraction stops being computable past E06 (the gray track ends, honestly,
  at the junction+1).
- Plumbing validated 2026-08-27 on pre-extension genes (`--genes WNT5A BMP4`).
