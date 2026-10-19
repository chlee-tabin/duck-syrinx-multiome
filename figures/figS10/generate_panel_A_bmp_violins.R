#!/usr/bin/env Rscript
# figS10 panel A (page letter F) — BMP ligands by D/V subpopulation (canonical colours) × stage.
# Empty subpopulation columns dropped per stage; all 9 BMP genes kept (BMP4/5 specificity).
#
# Layout: the 9 genes are split into TWO side-by-side blocks (BMP1-5 | BMP6-BMP15) instead of
# one 9-row stack, so the panel is wide-and-short (aspect ~2.3) and sits full-width above the
# B|C bottom row.  A single 9-row stack cannot be squeezed to that aspect without crushing the
# per-gene y axes: the x-label band (angled subpopulation names) is a fixed ~1.1 in cost, so
# 9 rows in a 3.8 in canvas leaves ~0.3 in per row and the y tick labels collide.
# The blocks reuse the shared figS08 `screen_violin()` style (single-sourced theme/palette);
# it renders to file, so each block is captured with last_plot() and re-combined by patchwork.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
source("figures/figS08/cr9_screen_io.R")
source("figures/_shared/colors.R")

BMP    <- c("BMP1","BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP15")
BLOCKS <- list(BMP[1:5], BMP[6:9])
OUT    <- "figures/figS10/panel_A_bmp_violins.png"
BLOCK_W <- 4.45   # in; per block. Total render width = 2 * BLOCK_W
BLOCK_H <- 3.8    # in; shared by both blocks (aspect 8.9 / 3.8 = 2.34)

d    <- load_cr9_screen()
long <- screen_long(d, BMP, "log__")
long$gene <- factor(long$gene, levels = BMP)

block <- function(genes) {
  sub <- droplevels(long[long$gene %in% genes, ])
  if (!setequal(levels(sub$gene), genes))
    stop("panel A: gene(s) missing from the screen export: ",
         paste(setdiff(genes, levels(sub$gene)), collapse = ", "))
  screen_violin(sub, facet_grid(gene ~ stage, scales = "free", space = "free_x", switch = "y"),
                tempfile(fileext = ".png"), w = BLOCK_W, h = BLOCK_H)
  last_plot()   # ggsave() print()s the plot, so last_plot() is the fully-built block
}
blocks <- lapply(BLOCKS, block)
# Shorter gene rows fit fewer y ticks: cap the breaks (the helper's default let a 0–2.5 row
# ask for six, which collide at this row height).  Keeps the no-margin-below-zero expansion.
thin_y <- scale_y_continuous(breaks = scales::breaks_pretty(3),
                             expand = expansion(mult = c(0, 0.05)))
# One y-axis title for the pair; tighten the inter-block gutter.
blocks[[1]] <- blocks[[1]] + thin_y + theme(plot.margin = margin(16, 4, 6, 6))
blocks[[2]] <- blocks[[2]] + thin_y + theme(axis.title.y = element_blank(),
                                            plot.margin = margin(16, 6, 6, 2))

# dpi only sets raster density (the physical size, hence the print font size, is fixed by
# width/height): 240 keeps the panel >= 300 dpi in its full-page-width composite slot.
ggsave(OUT, wrap_plots(blocks, nrow = 1), width = 2 * BLOCK_W, height = BLOCK_H,
       dpi = 240, limitsize = FALSE)
cat("Saved", OUT, "\n")
