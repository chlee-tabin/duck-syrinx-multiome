#!/usr/bin/env Rscript
# figS08 panel A — WNT ligands/modulators by D/V subpopulation (canonical colours) × stage.
# Empty subpopulation columns are dropped per stage (space="free_x"). Title/observations in LEGEND.md.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("figures/figS08/cr9_screen_io.R")
source("figures/_shared/colors.R")

WNT <- c("WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6",
         "WNT7A","WNT7B","WNT9A","WNT9B","WNT11","WNT16","SFRP2","SHISA2","WIF1")

d    <- load_cr9_screen()
long <- screen_long(d, WNT, "log__")
# show only WNT genes expressed in at least one subpopulation (max group-mean > 0.1)
expr <- long |> group_by(gene, stage, population) |> summarise(gm = mean(value), .groups = "drop") |>
        group_by(gene) |> summarise(mx = max(gm), .groups = "drop") |> filter(mx > 0.1) |> pull(gene)
dropped <- setdiff(levels(long$gene), as.character(expr))
if (length(dropped)) message("WNT flat (omitted from violins): ", paste(dropped, collapse = ", "))
long <- long |> filter(gene %in% expr) |> mutate(gene = droplevels(gene))

nr <- nlevels(long$gene)
# Portrait aspect: panel A occupies the LEFT half of the figS08 top row (A|B), so it is
# rendered at half-page width — keeping base_size fixed means the fonts survive the
# side-by-side placement instead of being shrunk by the compositor.
screen_violin(long,
  facet_grid(gene ~ stage, scales = "free", space = "free_x", switch = "y"),
  "figures/figS08/panel_A_wnt_violins.png", w = 5.75, h = 1.15 + 0.72 * nr)
cat(sprintf("(%d expressed WNT genes)\n", nr))
