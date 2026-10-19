#!/usr/bin/env Rscript
# figS08 panel B — WNT-pathway hits: per-cell CP10K expression boost L vs R, per-embryo paired.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("figures/figS08/cr9_screen_io.R")
source("figures/figS08/tier2_lr.R")
source("figures/_shared/colors.R")
d <- load_cr9_screen()
# Portrait aspect: panel B is the RIGHT half of the figS08 top row (A|B). Dimensions are
# passed here rather than changed in tier2_lr.R — that engine is shared with figS10.
tier2_lr_boost(d, c("WNT5A","SFRP2","SHISA2","PITX2"),
               "figures/figS08/panel_B_wnt_lr.png", width = 4.5, height = 6.9,
               legend_box = "vertical", annotate_p = FALSE)
