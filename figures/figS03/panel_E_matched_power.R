#!/usr/bin/env Rscript
# figS03 Panel E: matched-power CPM depth control — named PITX2-locus elements.
#
# Element-resolved (matching fig03 panel C + fig04): the a/b/c promoters, the
# +5.4 kb putative distal enhancer, and a representative upstream-desert R-drift
# peak.  Every ventral stage is subsampled to a common embryo-half count
# (n = 8 L / 6 R; bootstrap B=3000) on the CONSISTENT floor-20 pseudobulks (the
# fig03/fig04 laterality floor), and the L/R log2(mean CPM) recomputed.  The
# focal PITX2c promoter stays Left at matched power at every stage, PITX2a/b
# promoters are near-bilateral (within-locus control; PITX2b leans Left at HH13),
# the desert peak is Right at HH21/HH25 — i.e. the focal L-bias is biology, not power.
#
# Print size (2026-10-14): authored at its composite slot, the figS03 half-width slot
# (generate.py HALF_W = 1025 px at 300 dpi = 3.42 in wide), 2.6 in tall so the 2-up D/E
# row is no taller than at the freeze and figS03 page 2 keeps its ~811-pt height.
# Printed size = authored size x the booklet's fit scale, min(564/504, 720/page height)
# = 0.888 for an 811-pt page (manuscript booklet layout, not distributed), so
# the 7.3-pt text here prints at ~6.5 pt.  The legend sits below the plot (three
# columns) so no whisker runs under it, and the five elements are dodged so none touch.
# Shapes: only Left/Right-meaning elements borrow the house glyphs (PITX2c = blue
# circle; R-desert peak = red square); the PITX2a/b controls are open grey circle /
# open grey diamond, so no control reads as "Right" (PITX2b was a grey square).
#
# R/ggplot2 (matches the other figS03/fig03 data panels).  Title/subtitle +
# per-stage pseudobulk counts live in LEGEND.md (no in-image title).
#
# Reads : figures/data/g4f_depthmatch.tsv  (scripts/depthmatch_pitx2_locus.py, floor-20)
# Writes: figures/figS03/panel_E_matched_power.png
# Local : conda run -n duckfig Rscript figures/figS03/panel_E_matched_power.R

local_lib <- "~/R/4.4.2"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
    library(ggplot2); library(dplyr)
})

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
REPO <- normalizePath(file.path(HERE, "..", ".."))
DATA <- file.path(REPO, "figures", "data")
source(file.path(REPO, "figures", "_shared", "colors.R"))
lr_cols <- get_palette("lateral", labels = c("L", "R"))
L_BLUE <- unname(lr_cols["L"]); R_RED <- unname(lr_cols["R"])

STAGES <- c("HH13", "HH17/16", "HH21", "HH25")
# element order (legend top->bottom) + L/R-semantic colour + shape
ELEM <- c("PITX2c promoter", "+5.4 kb distal", "PITX2a promoter", "PITX2b promoter", "R-desert peak")
COLS <- c("PITX2c promoter" = L_BLUE, "+5.4 kb distal" = "#2166AC",
          "PITX2a promoter" = "#9AA6B2", "PITX2b promoter" = "#5A6B7B",
          "R-desert peak" = R_RED)
SHP  <- c("PITX2c promoter" = 16, "+5.4 kb distal" = 17,
          "PITX2a promoter" = 1,  "PITX2b promoter" = 5, "R-desert peak" = 15)
TXT  <- 7.3   # authored pt; prints at ~6.5 pt in the booklet (header)

d <- read.delim(file.path(DATA, "g4f_depthmatch.tsv"), check.names = FALSE)
stopifnot(all(ELEM %in% d$peakset))
d <- d %>%
    filter(peakset %in% ELEM) %>%
    mutate(stage = factor(stage, levels = STAGES),
           element = factor(peakset, levels = ELEM))

pd <- position_dodge(width = 0.6)
p <- ggplot(d, aes(stage, matched_mean, colour = element, shape = element, group = element)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = L_BLUE, alpha = 0.07) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = R_RED, alpha = 0.07) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey20") +
    geom_line(position = pd, linewidth = 0.4) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), linewidth = 0.4, position = pd) +
    geom_point(position = pd, size = 1.4, fill = "white", stroke = 0.6) +
    annotate("text", x = 0.55, y = 2.45, label = "Left", hjust = 0, vjust = 1,
             fontface = "bold", colour = L_BLUE, size = TXT / .pt) +
    annotate("text", x = 0.55, y = -1.95, label = "Right", hjust = 0, vjust = 0,
             fontface = "bold", colour = R_RED, size = TXT / .pt) +
    scale_colour_manual(values = COLS, name = NULL) +
    scale_shape_manual(values = SHP, name = NULL) +
    # y-axis carries the message: each stage subsampled to a COMMON depth/n (8L/6R)
    # so a surviving L-bias = depth-robust, not a read-depth artefact.
    scale_y_continuous("depth-matched log2(mean CPM, L / R)\n(every stage at n = 8 L / 6 R)",
                       breaks = seq(-1, 2, 1)) +
    coord_cartesian(ylim = c(-2, 2.5)) +
    scale_x_discrete(NULL) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE),
           shape  = guide_legend(ncol = 3, byrow = TRUE)) +
    theme_classic(base_size = TXT) +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"), legend.spacing.x = unit(1, "pt"),
          legend.key.height = unit(7, "pt"), legend.key.width = unit(10, "pt"),
          legend.text = element_text(size = TXT, margin = margin(l = 1, r = 3)),
          axis.text = element_text(size = TXT, colour = "grey20"),
          axis.title.y = element_text(size = TXT),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.margin = margin(2, 3, 1, 1))

ggsave(file.path(HERE, "panel_E_matched_power.png"), p, width = 3.42, height = 2.6, dpi = 300)
cat("wrote", file.path(HERE, "panel_E_matched_power.png"), "\n")
