#!/usr/bin/env Rscript
# figS09 panel A — BMP-arm (SMAD1/5/9) vs TGF-β-arm (SMAD2/3) by D/V subpopulation (canonical
# colours) × stage. Empty subpopulation columns dropped per stage. Title/observations in LEGEND.md.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("figures/figS08/cr9_screen_io.R")
source("figures/_shared/colors.R")

SMAD <- c("SMAD1","SMAD5","SMAD9","SMAD2","SMAD3")
ARM  <- c(SMAD1="BMP arm (SMAD1/5/9)", SMAD5="BMP arm (SMAD1/5/9)", SMAD9="BMP arm (SMAD1/5/9)",
          SMAD2="TGF-β arm (SMAD2/3)", SMAD3="TGF-β arm (SMAD2/3)")

d    <- load_cr9_screen()
long <- screen_long(d, SMAD, "log__")
long$arm <- factor(ARM[as.character(long$gene)],
                   levels = c("BMP arm (SMAD1/5/9)", "TGF-β arm (SMAD2/3)"))

# ggh4x is required (envs/duckfig_portable.yml): the nested arm/gene strips are part of the panel
if (!requireNamespace("ggh4x", quietly = TRUE)) stop("package 'ggh4x' is required for the nested facets")
facet <- ggh4x::facet_nested(arm + gene ~ stage, scales = "free", space = "free_x",
                             switch = "y", nest_line = element_line())

screen_violin(long, facet, "figures/figS09/panel_A_smad_violins.png", w = 8.5, h = 6.5)
