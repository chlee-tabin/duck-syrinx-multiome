#!/usr/bin/env Rscript
# figS10 panel B (page letter G) — BMP-pathway hits: per-cell CP10K expression boost L vs R,
# per-embryo paired.  Rendered at HALF-page width (it shares the bottom row with the receptor
# dotplot): the compositor maps a panel's render inches onto the page, so a page-width render
# squeezed into half a page would halve every font.  Stacked legend guides at this width.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("figures/figS08/cr9_screen_io.R")
source("figures/figS08/tier2_lr.R")
source("figures/_shared/colors.R")
d <- load_cr9_screen()
tier2_lr_boost(d, c("BMP4","BMP5","BMP7","PITX2"),
               "figures/figS10/panel_B_bmp_lr.png",
               width = 4.8, height = 5.9, legend_box = "vertical",
               annotate_p = FALSE)
