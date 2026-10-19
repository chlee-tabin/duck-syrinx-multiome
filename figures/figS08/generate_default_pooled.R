#!/usr/bin/env Rscript
# DEFAULT (faithful) pooled-mesoderm violins — the manuscript S8/S9/S10 layout reproduced
# on cr9: RNA expression in POOLED mesoderm (all non-Contamination cells) × stage, x = stage.
# These are the drop-in FALLBACK panels: if the subpopulation-resolved upgrade (panel A of each
# figure) is not adopted, swap these in. All panel genes kept (manuscript shows specificity).
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("figures/figS08/cr9_screen_io.R")

WNT  <- c("WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6",
          "WNT7A","WNT7B","WNT9A","WNT9B","WNT11","WNT16","SFRP2","SHISA2","WIF1")
SMAD <- c("SMAD1","SMAD5","SMAD9","SMAD2","SMAD3")
BMP  <- c("BMP1","BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP15")

d <- load_cr9_screen()   # drops Contamination; all remaining cells = pooled mesenchyme

mk <- function(genes, out, ncol, w, h) {   # title lives in LEGEND/PROVENANCE, not the panel
  long <- screen_long(d, genes, "log__")
  long$gene <- factor(long$gene, levels = genes[genes %in% unique(long$gene)])
  p <- ggplot(long, aes(stage, value)) +
    geom_violin(fill = "grey85", color = "black", linewidth = 0.25, scale = "width") +
    stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    labs(x = "Identity (stage)", y = "Expression level (log-normalized)") +  # title -> LEGEND/PROVENANCE
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(face = "italic", size = 8.5),
          strip.background = element_blank(), plot.margin = margin(16, 6, 6, 6))
  ggsave(out, p, width = w, height = h, dpi = 200)
  cat("Saved", out, "\n")
}

mk(WNT,  "figures/figS08/panel_A_wnt_violins_default_pooled.png", 3, 6.5, 9)
mk(SMAD, "figures/figS09/panel_A_smad_violins_default_pooled.png", 3, 6.5, 4)
mk(BMP,  "figures/figS10/panel_A_bmp_violins_default_pooled.png", 3, 6.5, 6)
