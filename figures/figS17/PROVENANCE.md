# Figure S17 — provenance (repo `figures/figS17/`)

Statistical record: `figures/figS17/STATS.md` *(not distributed)* (E-S1; pre-specified in
`docs/m5_esr1_prolif_stats_plan.md` *(not distributed)* §2.1/§2.2/§2.4). Claim-bearing numbers are
read from `figures/data/m5_esr1_lr_hh25_stats.tsv` and cross-checked at build
time against a recomputation through the record's loaders — mismatch raises.
Those numbers are stated in `LEGEND.md` instead of drawn on the panels (no
HL / CI / p / n text in-panel), and the A-vs-B pool difference is carried by the
y-axis title; the build-time cross-check prints to stdout.

Data sources (all committed, local build — no cluster):
- `figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz` (pool cells)
- `figures/data/d2g_signaling_HH25.tsv.gz` (per-cell ESR1_cp10k -> detection)
- `figures/data/lr_inference_multiome.tsv` (canonical d2o merged sides)
- `figures/data/hh25_percell_counts_esr1.tsv.gz` (panel C; per-cell raw ESR1
  UMI + depth; O2 job 50993735)
- `figures/data/h1_pseudobulk_snrna_{counts.tsv.gz,metadata.tsv}` and
  `figures/data/h2_deg_POOL_HH25.tsv` (S4 DESeq2 refit and its parity check)
- `figures/data/m5_esr1_lr_hh25_stats.tsv` (tidy record; producer
  `figures/figS17/stats_esr1_lr_hh25.py` + `s4_deseq2_refit.R` + `s8_glmm.R`)

Build (env from `envs/duckfig_portable.yml`, which carries DESeq2 for S4 and glmmTMB for S8):
```
conda run -n duckfig python figures/figS17/stats_esr1_lr_hh25.py --selfcheck   # -> m5_esr1_lr_hh25_stats.tsv
cd figures/figS17
conda run -n duckfig python generate_panels.py
conda run -n duckfig python generate.py   # -> figures/figS17.pdf
```

Methods notes:
- Panel C's LEGEND quotes the S8 `S8_depth_adjusted_embryo_slope` row (+0.34
  log-odds, p 0.0025) — the CONSERVATIVE random-slope fit, per the record; the
  two random-intercept fits (~+0.41, p ~ 1e-25) are anticonservative and must
  not be quoted as the headline. `generate_panels.py` still looks that row up
  (and prints it) purely so a record rename fails loud instead of silently
  orphaning the legend sentence.
- Panel B's two cohort-rule drops are genotypes g02/g04;
  the LEGEND reports the count only (no internal IDs in figures or legends).
- Slope-plot y-axes: panel A's left subpanel starts at 0. Panel A's difference
  column and panel C are left-minus-right plots, so both must span negative
  values — the zero-margin rule does not apply to them.
- Panel C per-embryo values behind the LEGEND's depth sentence (delta-depth in
  UMIs, delta-detection in pp): g07 -1329 / +3.45, g00 -1281 / -1.57,
  g12 -1103 / +2.87, g11 -236 / +7.69, g04 +1270 / +24.46, g02 +1314 / +13.51,
  dp500_donor29 +1533 / +19.02, g06 +1742 / +6.61. So FOUR embryos had the right
  half sequenced deeper and THREE of those four are still left-higher (the
  exception, g00, is the single right-higher embryo overall). The legend says
  "three of the four"; an earlier wording said "the three", which miscounted the
  right-deeper group.
- S4 (DESeq2 refit `~ genotype + syringeal_fraction + side`) was first run on
  2026-10-14 in the portable env (osx-64 build of `duckfig_portable.yml`; DESeq2
  1.46.0, R 4.4.3); the earlier local R could not load DESeq2. The refit pools the
  committed h1 counts exactly as h2's POOL scope and is gated on parity: the same
  pooling fitted `~ genotype + side` reproduces the committed h2 ESR1 row (log2FC,
  lfcSE, p within 1e-6 relative; 13,157 genes clearing independent filtering). The
  same run reproduced every other row of the stats table byte-for-byte, and the three
  panels and the composite are unchanged. S4 is a sensitivity arm and is not quoted
  in the legend: log2FC +0.27 (Wald 95% CI -0.01 to +0.55), p 0.062, versus +0.34,
  p 0.015 without the composition covariate.
