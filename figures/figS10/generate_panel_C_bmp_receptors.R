#!/usr/bin/env Rscript
# figS10 panel C (page letter H) — BMP receptor competence by D/V subpopulation × stage.
# Dot size = % of cells expressing, greyscale fill = mean CP10K.  BMP4 is kept as the
# ligand-SOURCE contrast column right of a dashed rule, so the reader sees the receptors are
# flat where the ligand is not.  Greyscale only: nothing here encodes a biological
# variable, so per the house rules no colour is spent on it.
# TRANSPOSED (2026-08-24) relative to panel A: genes on x, subpopulations on y, stage as
# stacked row blocks (space="free_y" so each stage carries only its own subpopulations).
# The portrait aspect lets it sit beside the paired-L/R panel in the bottom row.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("figures/figS08/cr9_screen_io.R")

GENES <- c("BMPR1A", "BMPR1B", "BMPR2", "BMP4")   # x axis left -> right
LIGAND <- "BMP4"                                  # drawn right of the separator rule
OUT <- "figures/figS10/panel_C_bmp_receptors.png"

d <- load_cr9_screen()
miss <- GENES[!paste0("cp10k__", GENES) %in% colnames(d)]
if (length(miss)) stop("panel C: gene(s) absent from the screen export: ", paste(miss, collapse = ", "))

agg <- screen_long(d, GENES, "cp10k__") |>
  group_by(stage, population, gene) |>
  summarise(pct = 100 * mean(value > 0), mean_cp10k = mean(value), n = n(), .groups = "drop")

# Fail loud (project convention: no silent fallbacks).  A (stage, population) that exists in the
# data must carry a dot for EVERY gene — a hole would render as a missing dot that reads
# as "not measured", and an NA would vanish with only a ggplot warning.
cells <- nrow(distinct(agg, stage, population))
if (nrow(agg) != cells * length(GENES))
  stop(sprintf("panel C: incomplete (stage x population) x gene grid — %d rows, expected %d",
               nrow(agg), cells * length(GENES)))
if (anyNA(agg$pct) || anyNA(agg$mean_cp10k))
  stop("panel C: NA percent-expressing / mean CP10K in the aggregate")
if (!all(GENES %in% levels(agg$gene)))
  stop("panel C: gene dropped between export and aggregate: ",
       paste(setdiff(GENES, levels(agg$gene)), collapse = ", "))

agg$gene <- factor(agg$gene, levels = GENES)                    # first GENES entry leftmost
agg$population <- factor(agg$population, levels = rev(levels(agg$population)))  # POP_ORDER top -> bottom
sep_x <- which(levels(agg$gene) == LIGAND) - 0.5                 # rule between receptors and ligand

p <- ggplot(agg, aes(gene, population)) +
  geom_vline(xintercept = sep_x, colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_point(aes(size = pct, fill = mean_cp10k), shape = 21, colour = "black", stroke = 0.25) +
  facet_grid(stage ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "grey93", high = "black", limits = c(0, NA),
                      name = "Mean expression (CP10K)") +
  scale_size_area(max_size = 6.5, limits = c(0, 100), breaks = c(10, 25, 50, 75),
                  name = "% cells expressing") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(face = "italic", angle = 40, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.2),
    strip.text.y = element_text(face = "bold", size = 9, angle = 0),
    strip.background = element_blank(),
    # location="plot" pins the guides to the canvas edge, not the panel edge — at this narrow
    # width the panel starts ~1.5 in in, which would otherwise push the guides off the right.
    legend.position = "top", legend.box = "vertical", legend.justification = "left",
    legend.location = "plot",
    legend.box.spacing = unit(0.3, "lines"), legend.spacing.y = unit(0.1, "lines"),
    legend.key.height = unit(0.5, "lines"), legend.key.width = unit(0.9, "lines"),
    legend.title = element_text(size = 8), legend.text = element_text(size = 7),
    panel.spacing.y = unit(0.5, "lines"),
    plot.margin = margin(16, 6, 6, 6)) +
  guides(fill = guide_colourbar(order = 1, title.position = "top", barheight = unit(0.35, "lines")),
         size = guide_legend(order = 2, title.position = "top", nrow = 1))

ggsave(OUT, p, width = 3.5, height = 5.9, dpi = 250, limitsize = FALSE)  # 250: >=300 dpi in its half-row slot
cat("Saved", OUT, "\n")
print(as.data.frame(agg |> arrange(stage, population, gene) |>
                      mutate(pct = round(pct, 1), mean_cp10k = round(mean_cp10k, 2))))
