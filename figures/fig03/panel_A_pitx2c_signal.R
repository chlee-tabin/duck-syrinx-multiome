#!/usr/bin/env Rscript
# fig03 panel A — the PITX2-locus signal per embryo-half (RNA + isoform chromatin)
# that DEFINES the L/R sides, shown as individual pseudobulks.
#
# Four rows, four stages (HH13 / HH17/16 / HH21 / HH25); each point = one
# embryo-half; paired halves of the same bisected embryo connected by a grey
# segment; L = blue circle, R = red square; black bar = median; point size = cells
# in the half.  Rows (top -> bottom):
#   1. snRNA PITX2 level        (mean expression; lr_inference_multiome pitx2_mean)
#   2. PITX2c promoter ATAC      (fragments / cell; the side-defining signal)
#   3. PITX2a promoter ATAC      (within-locus control)
#   4. PITX2b promoter ATAC      (within-locus control)
#
# The juxtaposition is the figure's opening argument: ONLY PITX2c (RNA + its own
# promoter chromatin) bisects each embryo decisively, and it does so most clearly
# at HH17/16 — exactly where the snRNA laterality signal is weak (figS02D). PITX2a/b
# promoters stay bilateral (controls). The per-stage strip carries the side-source:
#   * HH13 / HH21 / HH25: side assigned by the snRNA 4-metric inference (d2l) — so
#     the ATAC separation is CROSS-MODAL corroboration (snRNA defines, chromatin
#     confirms).
#   * HH17/16: snRNA weak -> side assigned by the chromatin PITX2c inference (d2n);
#     panel A's PITX2c row at this stage IS the side-defining signal (no independent
#     p-value reported — quantification + per-stage stats live in figS03 panel H).
#
# Data (both already on-branch):
#   figures/data/step_f5_pitx2_region_counts.tsv  (ATAC c/a/b promoter frags/cell)
#   figures/data/lr_inference_multiome.tsv        (snRNA pitx2_mean per half)
#
# Render (duckfig):
#   conda run -n duckfig Rscript figures/fig03/panel_A_pitx2c_signal.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

HERE   <- "figures/fig03"
F      <- "figures/data"
source("figures/_shared/colors.R")

MIN_CELLS <- 20
STAGES    <- c("HH13", "HH17/16", "HH21", "HH25")
# (two-line side-source tag so the strip fits its facet at print size)
SIDE_SRC  <- c("HH13"="sides:\nsnRNA", "HH17/16"="sides:\nchromatin",
               "HH21"="sides:\nsnRNA", "HH25"="sides:\nsnRNA")
# rows: snRNA PITX2 then the three isoform promoters (c first; a/b = controls)
# (3-line horizontal row strips: rotated strips did not fit the print-size rows)
MEASURES  <- c("snRNA\nPITX2\n(mean expr.)",
               "ATAC\nPITX2c prom.\n(frags/cell)",
               "ATAC\nPITX2a prom.\n(frags/cell)",
               "ATAC\nPITX2b prom.\n(frags/cell)")
REGION_TO_MEASURE <- c(PITX2c_promoter = MEASURES[2],
                       PITX2a_promoter = MEASURES[3],
                       PITX2b_promoter = MEASURES[4])

side_cols   <- get_palette("lateral", labels = c("L", "R"))
side_shapes <- c("L" = 16, "R" = 15)

# ── ATAC rows: PITX2c/a/b promoter fragments per cell (f5) ──────────────────
atac <- read.delim(file.path(F, "step_f5_pitx2_region_counts.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(region %in% names(REGION_TO_MEASURE), n_cells >= MIN_CELLS) |>
  transmute(genotype, library = library_d2l, side, stage,
            measure = REGION_TO_MEASURE[region],
            value = frags_per_cell, npool = n_cells)

# ── snRNA row: PITX2 mean expression per the SAME embryo-halves (lr_inference) ─
lr <- read.delim(file.path(F, "lr_inference_multiome.tsv"),
                 stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(!is.na(pitx2_mean)) |>
  distinct(genotype, library_m, .keep_all = TRUE) |>
  select(genotype, library_m, pitx2_mean, n_inference)

halves <- atac |> distinct(genotype, library, side, stage)

# f5 carries its own copy of each half's side, written when f5 last ran on the
# cluster.  It must agree with the canonical sidecar this panel joins on, or the
# panel draws stale sides (e.g. an f5 from before the lane-loading-record correction
# draws the tube-17-B duplicate as an L/R pair).  Every f5 half is checked:
#   * plain library label: the sidecar row exists, is NOT same_side_dup (a
#     duplicate must arrive as one pooled "<libA>+<libB>" sample), same side;
#   * pooled "<libA>+<libB>" label: every part is a same_side_dup sidecar row on
#     the same side.
# Any disagreement stops the render; the fix is an f5 rerun, not an override.
FIRMNESS_ALPHABET <- c("direct", "placement", "singleton_rank", "weak", "same_side_dup")
side_ref <- read.delim(file.path(F, "lr_inference_multiome.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE) |>
  select(genotype, library_m, firmness, inferred_side)
unknown_firm <- setdiff(unique(side_ref$firmness), FIRMNESS_ALPHABET)
if (length(unknown_firm))
  stop("lr_inference_multiome.tsv firmness value(s) outside {",
       paste(FIRMNESS_ALPHABET, collapse = ", "), "}: ", paste(unknown_firm, collapse = ", "))
if (anyDuplicated(side_ref[, c("genotype", "library_m")]))
  stop("lr_inference_multiome.tsv has duplicate (genotype, library_m) rows")
f5_halves <- read.delim(file.path(F, "step_f5_pitx2_region_counts.tsv"),
                        stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(region %in% names(REGION_TO_MEASURE)) |>
  distinct(genotype, library = library_d2l, side)
side_problems <- character(0)
for (i in seq_len(nrow(f5_halves))) {
  g <- f5_halves$genotype[i]; lib <- f5_halves$library[i]; s <- f5_halves$side[i]
  pooled <- grepl("+", lib, fixed = TRUE)
  parts <- if (pooled) strsplit(lib, "+", fixed = TRUE)[[1]] else lib
  for (p in parts) {
    ref <- side_ref[side_ref$genotype == g & side_ref$library_m == p, ]
    what <- sprintf("%s %s (f5 side %s)", g, lib, s)
    if (nrow(ref) == 0) {
      side_problems <- c(side_problems, sprintf("%s: no sidecar row for %s", what, p))
    } else if (pooled && ref$firmness != "same_side_dup") {
      side_problems <- c(side_problems, sprintf("%s: pooled label but %s is %s, not same_side_dup",
                                                what, p, ref$firmness))
    } else if (!pooled && ref$firmness == "same_side_dup") {
      side_problems <- c(side_problems, sprintf("%s: sidecar marks it same_side_dup, so f5 must pool its libraries",
                                                what))
    } else if (ref$inferred_side != s) {
      side_problems <- c(side_problems, sprintf("%s: sidecar side for %s is %s",
                                                what, p, ref$inferred_side))
    }
  }
}
if (length(side_problems))
  stop("step_f5_pitx2_region_counts.tsv disagrees with lr_inference_multiome.tsv on ",
       length(side_problems), " embryo-half check(s); f5 predates the current sidecar, rerun f5 ",
       "on the cluster:\n  ", paste(side_problems, collapse = "\n  "))

# A same-side duplicate (lane-loading-record split tube: ONE half sequenced in two
# libraries) is ONE f5 sample labelled "<libA>+<libB>" (figures/_shared/lr_firmness.py).
# Its snRNA value is the pooled mean over both libraries' inference-pool nuclei
# (weights = snRNA nuclei per library, d2l n_inference); npool sums the sidecar counts.
merged <- unique(halves$library[grepl("+", halves$library, fixed = TRUE)])
if (length(merged)) {
  d2l_n <- read.delim(file.path(F, "d2l_lr_inference.tsv"), stringsAsFactors = FALSE,
                      check.names = FALSE) |>
    select(genotype, library_m, n_rna = n_inference)
  pooled <- halves |> filter(library %in% merged) |> distinct(genotype, library) |>
    rowwise() |>
    mutate(parts = list(strsplit(library, "+", fixed = TRUE)[[1]])) |>
    ungroup()
  pooled_rows <- lapply(seq_len(nrow(pooled)), function(i) {
    g <- pooled$genotype[i]; parts <- pooled$parts[[i]]
    comp <- lr |> filter(genotype == g, library_m %in% parts) |>
      left_join(d2l_n, by = c("genotype", "library_m"))
    if (nrow(comp) != length(parts) || anyNA(comp$n_rna))
      stop("pooled half ", g, " ", pooled$library[i],
           ": snRNA rows/nuclei missing for its libraries")
    data.frame(genotype = g, library_m = pooled$library[i],
               pitx2_mean = weighted.mean(comp$pitx2_mean, comp$n_rna),
               n_inference = sum(comp$n_inference))
  })
  lr <- bind_rows(lr, do.call(rbind, pooled_rows))
}

snrna <- halves |>
  left_join(lr, by = c("genotype", "library" = "library_m"))
if (anyNA(snrna$pitx2_mean))
  stop("embryo-halves without an snRNA PITX2 value: ",
       paste(unique(paste(snrna$genotype, snrna$library)[is.na(snrna$pitx2_mean)]),
             collapse = ", "))
snrna <- snrna |>
  transmute(genotype, library, side, stage,
            measure = MEASURES[1], value = pitx2_mean, npool = n_inference)

df <- bind_rows(snrna, atac) |>
  filter(stage %in% STAGES) |>
  mutate(measure   = factor(measure, levels = MEASURES),
         side      = factor(side, levels = c("L", "R")),
         stage_lab = factor(sprintf("%s\n%s", stage, SIDE_SRC[stage]),
                            levels = sprintf("%s\n%s", STAGES, SIDE_SRC[STAGES])))

# paired segments: a genotype with BOTH halves in a (measure, stage)
pair_seg <- df |>
  group_by(measure, stage_lab, genotype) |>
  filter(all(c("L", "R") %in% side)) |>
  summarise(L = value[side == "L"][1], R = value[side == "R"][1], .groups = "drop")

p <- ggplot(df, aes(side, value)) +
  geom_segment(data = pair_seg, aes(x = "L", xend = "R", y = L, yend = R),
               inherit.aes = FALSE, colour = "grey70",
               linewidth = 0.28, alpha = 0.7) +
  geom_point(aes(colour = side, shape = side, size = npool), alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55,
               linewidth = 0.4, colour = "black", show.legend = FALSE) +
  facet_grid(measure ~ stage_lab, scales = "free_y", switch = "y") +
  scale_colour_manual(values = side_cols, name = "Side") +
  scale_shape_manual(values = side_shapes, name = "Side") +
  scale_size_continuous(range = c(0.8, 2.8),
                        breaks = c(50, 200, 500, 1500), name = "cells in half") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = NULL) +
  # PNAS final font pass (2026-10): rendered at 300 dpi on a 2/3-size canvas
  # (the compositor assumes 300 dpi) so the Fig 5 print factor is 0.653; fonts
  # are sized for >= 7 pt printed (ticks) / ~7.5-8 pt (strips).
  theme_classic(base_size = 11) +
  theme(strip.text.x = element_text(size = 11, face = "bold", lineheight = 0.88),
        strip.text.y.left = element_text(size = 10.4, face = "bold", angle = 0,
                                         hjust = 1, lineheight = 0.88),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10.8),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(0.02, "cm"),
        legend.text = element_text(size = 10.8),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.34, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.1, "cm"),
        panel.spacing.x = unit(0.45, "lines"),
        panel.spacing.y = unit(0.75, "lines"),
        plot.margin = margin(3, 4, 2, 2))

out <- file.path(HERE, "panel_A_pitx2c_signal.png")
ggsave(out, p, width = 5.067, height = 5.25, dpi = 300)
cat("wrote", out, "\n")
