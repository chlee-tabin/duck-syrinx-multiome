# Figure S16 — provenance (repo `figures/figS16/`)

Statistical record: `figures/figS16/STATS.md` *(not distributed)* (P-S1; pre-specified in
`docs/m5_esr1_prolif_stats_plan.md` *(not distributed)* §2.1/§2.3/§2.4;
delta = +0.137). Claim-bearing numbers are read from
`figures/data/m5_cycling_lr_hh25_stats.tsv` and cross-checked at build time
against a recomputation through the record's own loaders — mismatch raises.
Those numbers are stated in `LEGEND.md` instead of drawn on the panels (no
HL / CI / p / n / padj text in-panel); the build-time cross-check prints to stdout.

Data sources (all committed; build is local — the cluster only ever REGENERATES
the committed inputs):
- `figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz` (per-cell phase; ccreg
  partition, collider-safe for the cycling definition phase != G1)
- `figures/data/lr_inference_multiome.tsv` (canonical d2o merged sides)
- `figures/data/step03_sex_per_genotype.tsv`
- `figures/data/m5_cycling_lr_hh25_stats.tsv` (tidy record; producer
  `figures/figS16/stats_cycling_lr_hh25.py`)
- `figures/data/m5_mitogen_program_hh25.tsv` (panel C multiome facets; h2 paired
  DESeq2 read + ESR1; producer `figures/figS16/mitogen_program_lookup.py`)
- `figures/data/scrna_e6e9e11/cr9_syrinx_lr_deseq2_{E09,E11}_{M,F}.tsv` (panel C
  growth-phase facets; producer `pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R`
  on the `cr9_syrinx_pseudobulk_*` counts)

Build (from repo root; the first two producers are needed only when their
committed TSVs are stale):
```
conda run -n duckfig python figures/figS16/mitogen_program_lookup.py   # -> m5_mitogen_program_hh25.tsv
conda run -n duckfig Rscript --vanilla pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R   # -> cr9_syrinx_lr_deseq2_*.tsv
conda run -n duckfig python figures/figS16/generate_panels.py
conda run -n duckfig python figures/figS16/generate.py   # -> figures/figS16.pdf
```

Methods notes:
- Standardisation/cohort machinery = `figures/_shared/paired_standardised.py`
  (>= 30 pool cells/side cohort rule; >= 20 cells/side stratum rule; the g00
  Syringeal stratum drop is the record's disclosed substitution).
- Panel C (expanded 2026-08-27) shows 8 genes (BMP4, BMP7, ESR1, MKI67, SFRP2,
  SHISA2, TOP2A, WNT5B) across 6 facets: HH21/HH25 (multiome) + E9/E11 × sex
  (growth-phase paired DESeq2). Row order = HH25 log2FC; one
  shared x-scale; filled = genome-wide padj < 0.05; NA padj (independent
  filtering) drawn as not-significant. WNT2B is not in the committed mitogen TSV
  and is therefore not plotted.
- Panel A left subpanel y-axis does not start at 0 (paired slope readability);
  the bar-chart zero-margin rule does not apply to slope plots.
- Panel A's delta = +13.7 pp benchmark is now an UNLABELLED dashed reference
  line (the `delta = +13.7 pp (E9 effect)` text was removed 2026-08-24); the
  `set_ylim(top=...)` in `diff_axes` is what keeps it inside the axes, so do not
  drop it. Panel C keeps significance as filled-vs-open markers with no padj
  text. Both are named in `LEGEND.md`.
