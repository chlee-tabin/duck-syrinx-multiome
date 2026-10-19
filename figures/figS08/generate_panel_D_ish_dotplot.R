#!/usr/bin/env Rscript
# figS08 panel D — ISH ↔ subcluster correspondence (the "make sense of the differing S11
# pattern" panel). For every gene shown by in situ in manuscript S11 (WNT5b, SFRP2, SHISA2,
# BMP7, ESR1) + the active BMP ligands (BMP4/5), a dotplot of mean expression (greyscale)
# and % cells expressing (size) across D/V subpopulations × stage. Reads the differing
# spatial domains of S11 onto specific cr9 populations (SHISA2→Chondrogenic condensations,
# SFRP2→broad Syringeal fibroblast, BMP7→restricted, …). Pairs with the S11 ISH montage.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("figures/figS08/cr9_screen_io.R")

GENES <- c("WNT5B","SFRP2","SHISA2","BMP4","BMP5","BMP7","ESR1")   # S11 ISH rows + active ligands
d    <- load_cr9_screen()
long <- screen_long(d, GENES, "log__")
summ <- long |> group_by(stage, population, gene) |>
  summarise(mean_expr = mean(value), pct = 100*mean(value > 0), n = n(), .groups="drop") |>
  filter(n >= 20)
summ$gene <- factor(summ$gene, levels = GENES)   # first gene leftmost on x

# Standard dotplot: shade = mean expression (greyscale), size = % cells expressing. Genes on
# the x-axis, subpopulation on y with stage as stacked row blocks (space="free_y" so each
# stage takes only the populations it has) — a PORTRAIT aspect, so panel D pairs with the
# portrait ISH montage (panel C) in the figS08 bottom row instead of running page-wide.
# No subcluster colour here: the colour channel carries the quantitative expression.
p <- ggplot(summ, aes(gene, population)) +
  geom_point(aes(size = pct, color = mean_expr)) +
  facet_grid(stage ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(limits = rev) +   # POP_ORDER reads top-to-bottom
  scale_color_gradient(low = "grey88", high = "black", name = "mean\nlog-expr") +
  scale_size_continuous(range = c(0.2, 6), name = "% cells\nexpressing", limits = c(0, 100)) +
  labs(x = NULL, y = NULL) +   # title/correspondence text -> LEGEND.md
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "italic", size = 8),
        axis.text.y = element_text(size = 7.5),
        strip.text.y = element_text(face = "bold", angle = -90),
        strip.background = element_blank(),
        plot.margin = margin(16, 6, 6, 6),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3))

ggsave("figures/figS08/panel_D_ish_dotplot.png", p, width = 5.9, height = 6.4, dpi = 200)
cat("Saved panel_D_ish_dotplot.png\n")
