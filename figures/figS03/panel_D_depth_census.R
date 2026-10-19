#!/usr/bin/env Rscript
# figS03 Panel D — PITX2 Pitx2-TAD locus per-EMBRYO-HALF DEPTH CENSUS (ggbeeswarm).
#
# The SOURCE is the FLOOR-20 sample meta (g4f_sample_meta_floor20.tsv),
# the fig03/panel-E laterality floor — so panel D's per-half census is the SAME
# embryo-half cohort the matched-power control (panel E) and fig04 panel A are computed
# on, and the cohort the fig03 pileup restricts its HH17/16-HH25 halves to (HH17/16 ventral
# 16 L / 11 R in both, the split-tube halves 17-5/17-6 with both aliquots summed;
# fig03/PROVENANCE.md) (floor-50 admitted far fewer halves, e.g. HH13 ventral 1L/2R; floor-20 =
# HH13 8L/6R, dorsal HH17/16 12L/5R — the panel-E pools).
# (Secondary: repo-relative path resolution + the registry L/R palette, matching C/E.)
#
# One dot = one embryo-half = SUM over its inference-pool (ventral) or dorsal-control
# subclusters per genotype x side.  The reader sees:
#   - the L>R cell/fragment imbalance that grows ventral HH17/16 -> HH25 (the depth
#     confound the DESeq2 size factors normalize out — r(libsize, SF)=0.9985 on this
#     floor-20 cohort; 0.9987 is the floor-50 per-subcluster table, a different tier),
#   - the per-observation DEARTH that makes some halves low-power (HH13 ventral 8L/6R;
#     dorsal HH25 5L/7R), motivating the matched-power control in panel E.
#
# Two metrics (cells, fragments) x two compartments (ventral inference pool, dorsal
# bilateral control) — the SAME grouping as panel C.  Excluded pops (Skeletogenic /
# Mesothelium / Venous-pole spl.) are not shown.  L = blue, R = red (registry SSOT,
# matching panels C/E).  log10 y so the low-n dearth is visible.  Numbers above each
# swarm = n EMBRYO-HALVES on that side.  No in-image title (figS03 house rule).
#
# Print size (2026-10-14): authored at its composite slot, the figS03 half-width slot
# (generate.py HALF_W = 1025 px at 300 dpi = 3.42 in wide), 2.6 in tall, so the 2-up D/E
# row is no taller than at the freeze (781 px) and figS03 page 2 keeps its ~811-pt height.
# Printed size = authored size x the booklet's fit scale, min(564/504, 720/page height)
# = 0.888 for an 811-pt page (manuscript booklet layout, not distributed), so
# the 7.3-pt text here prints at ~6.5 pt in the booklet (the freeze 8.4-in canvas
# printed at ~3 pt).  Stage labels print once, under the bottom row.
#
#   conda run -n duckfig Rscript figures/figS03/panel_D_depth_census.R
local_lib <- "~/R/4.4.2"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(ggbeeswarm); library(patchwork)
})
set.seed(1)

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
REPO <- normalizePath(file.path(HERE, "..", ".."))
F   <- file.path(REPO, "figures", "data")
source(file.path(REPO, "figures", "_shared", "colors.R"))
lr_cols <- get_palette("lateral", labels = c("L", "R"))
LBLUE <- unname(lr_cols["L"]); RRED <- unname(lr_cols["R"])   # registry: L #4393C3, R #D6604D
STAGES <- c("HH13", "HH17/16", "HH21", "HH25")
DODGE  <- 0.72

meta_raw <- read_tsv(file.path(F, "g4f_sample_meta_floor20.tsv"), show_col_types = FALSE)
# fail LOUD on an unexpected compartment/side value — else the %in% filter (or the
# factor() below) would silently drop it (the no-silent-default house rule).  "none"
# = the excluded compartment (Skeletogenic / Mesothelium / Venous-pole spl.), legitimately
# filtered out; an unmapped side ("?") or a typo'd compartment must surface, not vanish.
stopifnot(all(meta_raw$inferred_side %in% c("L", "R")),
          all(meta_raw$compartment %in% c("ventral", "dorsal", "none")))
meta <- meta_raw %>%
  filter(compartment %in% c("ventral", "dorsal")) %>%
  mutate(compartment = factor(compartment, levels = c("ventral", "dorsal"),
                              labels = c("Ventral (inference pool)", "Dorsal (bilateral control)")),
         stage = factor(d2g_stage, levels = STAGES),
         side  = factor(inferred_side, levels = c("L", "R")),
         frag_M = total_fragments / 1e6)

# per-side pseudobulk counts, with explicit 0 where a side is absent so the dearth
# shows even with no point to plot.
present <- meta %>% distinct(compartment, stage)
nlab <- present %>%
  crossing(side = factor(c("L", "R"), levels = c("L", "R"))) %>%
  left_join(count(meta, compartment, stage, side), by = c("compartment", "stage", "side")) %>%
  mutate(n = coalesce(n, 0L))

# one metric panel: dodged L/R quasirandom swarm + median crossbars + top n labels
TXT <- 7.3   # authored pt; prints at ~6.5 pt in the booklet (header)
metric_panel <- function(yvar, ylab, ybreaks, ylabels, ytop, show_strip = TRUE,
                         show_x = TRUE) {
  nlab2 <- nlab %>% mutate(y = ytop)
  meds  <- meta %>% group_by(compartment, stage, side) %>%
    summarise(m = median(.data[[yvar]]), .groups = "drop")
  ggplot(meta, aes(x = stage, y = .data[[yvar]], colour = side, shape = side)) +
    geom_crossbar(data = meds, aes(x = stage, y = m, ymin = m, ymax = m, group = side),
                  position = position_dodge(width = DODGE), width = 0.52,
                  colour = "grey15", linewidth = 0.3, inherit.aes = FALSE) +
    geom_quasirandom(dodge.width = DODGE, size = 0.7, alpha = 0.85, width = 0.16,
                     method = "quasirandom") +
    # (L label right-aligned, R left-aligned on their dodge positions so two-digit
    #  counts such as 15 / 11 do not run together at the 3.42-in slot)
    geom_text(data = nlab2, aes(x = stage, y = y, label = n, colour = side,
                                hjust = ifelse(side == "L", 0.8, 0.2)),
              position = position_dodge(width = DODGE), inherit.aes = FALSE,
              size = TXT / .pt, fontface = "bold", vjust = 1, show.legend = FALSE) +
    facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = c(L = LBLUE, R = RRED),
                        labels = c(L = "Left", R = "Right"), name = NULL) +
    # house L/R encoding: Left = blue circle, Right = red square
    scale_shape_manual(values = c(L = 16, R = 15),
                       labels = c(L = "Left", R = "Right"), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)),
           shape  = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
    scale_y_log10(breaks = ybreaks, labels = ylabels,
                  expand = expansion(mult = c(0.04, 0.13))) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = TXT) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25),
          axis.text.x = if (show_x) element_text(size = TXT, margin = margin(t = 1)) else element_blank(),
          axis.text.y = element_text(size = TXT),
          axis.title.y = element_text(size = TXT),
          strip.text = if (show_strip) element_text(size = TXT, face = "bold",
                                                    margin = margin(1, 0, 1, 0)) else element_blank(),
          panel.spacing.x = unit(4, "pt"),
          legend.position = "bottom", legend.key.height = unit(7, "pt"),
          legend.text = element_text(size = TXT),
          plot.margin = margin(1, 2, 0, 1))
}

p_cells <- metric_panel("n_cells", "nuclei / half",
                        ybreaks = c(50, 100, 200, 500, 1000, 2000),
                        ylabels = c("50", "100", "200", "500", "1000", "2000"),
                        ytop = max(meta$n_cells) * 1.55, show_x = FALSE)
p_frag  <- metric_panel("frag_M", "fragments / half (M)",
                        ybreaks = c(0.25, 0.5, 1, 2, 4, 8),
                        ylabels = c("0.25", "0.5", "1", "2", "4", "8"),
                        ytop = max(meta$frag_M) * 1.55, show_strip = FALSE)

# Title/subtitle kept OFF-figure (figS03 house rule); legend lives in LEGEND.md.
p <- (p_cells / p_frag) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "pt"))

out <- file.path(HERE, "panel_D_depth_census.png")
ggsave(out, p, width = 3.42, height = 2.6, dpi = 300)
cat("wrote", out, "\n")
