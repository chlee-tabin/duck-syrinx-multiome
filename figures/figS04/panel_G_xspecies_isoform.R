#!/usr/bin/env Rscript
# figS04 panel G — cross-species PITX2 isoform × modality per-stage L/R beeswarm (the ORIENTED
# supplement to fig04 panel B's orientation-free magnitude: "here is how it resolves once an L/R
# orientation is assigned", with the chick best/worst design-legal envelope so chick cannot be
# over-claimed positive).  PITX2a and PITX2b shown SEPARATELY — the chick galGal6a GTF annotates
# all three isoforms and the ATAC peaks resolve a vs b, so the earlier pooled "PITX2a/b" was a
# presentation choice, not an annotation limit.
#
# Duck (d2n sides at the HH17/16 trough) vs chick best/worst design-legal envelope, ventral
# splanchnic (syrinx-precursor) pool.  FOUR rows:
#   ATAC PITX2c  — DIRECTION / laterality (the asymmetric, NODAL-primed/memory isoform)
#   ATAC PITX2a  — BMP/SMAD5 cis-anchor; expected BILATERAL (bilateral BMP input)
#   ATAC PITX2b  — bystander; expected BILATERAL (within-locus control)
#   RNA  total   — FISH/HCR-comparable LEVEL (two-wave / HH17 trough)
# Point AREA ∝ min(n_L,n_R) cells; open circles, border = group; off-scale ▲▼ with true value.
# Data: figures/data/chick_control/xspecies_pitx2_isoform_laterality.tsv (built by
#   pipeline/chick/build_isoform_data.py).
#
#   conda run -n duckfig Rscript figures/figS04/panel_G_xspecies_isoform.R
suppressPackageStartupMessages({library(tidyverse); library(ggbeeswarm)})

HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
TSV <- file.path(HERE, "..", "data", "chick_control", "xspecies_pitx2_isoform_laterality.tsv")
df <- read_tsv(TSV, show_col_types = FALSE)

df <- df %>% mutate(panel = case_when(
  modality == "ATAC" & isoform == "PITX2c" ~ "ATAC · PITX2c  (ΔFRiP ×1e-5)",
  modality == "ATAC" & isoform == "PITX2a" ~ "ATAC · PITX2a  (ΔFRiP ×1e-5)",
  modality == "ATAC" & isoform == "PITX2b" ~ "ATAC · PITX2b  (ΔFRiP ×1e-5)",
  modality == "RNA"  & isoform == "total"  ~ "RNA · total PITX2  (Δcp10k)",
  TRUE ~ NA_character_))
if (any(is.na(df$panel))) stop("unmapped (modality,isoform): ",
  paste(unique(paste(df$modality, df$isoform)[is.na(df$panel)]), collapse = ", "))

panel_lev <- c("ATAC · PITX2c  (ΔFRiP ×1e-5)", "ATAC · PITX2a  (ΔFRiP ×1e-5)",
               "ATAC · PITX2b  (ΔFRiP ×1e-5)", "RNA · total PITX2  (Δcp10k)")
panel_lev <- panel_lev[panel_lev %in% df$panel]
stage_lab <- c("HH13" = "HH13", "HH17/16" = "HH17/16\n(chick HH15-16)",
               "HH21" = "HH21\n(duck only)", "HH25" = "HH25")
grp_lev <- c("duck", "chick best", "chick worst")
df <- df %>% mutate(panel = factor(panel, levels = panel_lev),
                    stage = factor(stage, levels = names(stage_lab)),
                    group = factor(group, levels = grp_lev))

# per-PANEL clamp (each panel its own units/scale) -> a few off-scale points become edge ▲/▼.
clamp <- df %>% group_by(panel) %>%
  summarise(lo = pmin(quantile(signed, 0.02, names = FALSE), -0.05),
            hi = pmax(quantile(signed, 0.98, names = FALSE),  0.05), .groups = "drop")
df <- df %>% left_join(clamp, by = "panel") %>%
  mutate(pshape = case_when(signed > hi ~ "hi", signed < lo ~ "lo", TRUE ~ "in"),
         yplot = pmin(pmax(signed, lo), hi))
off <- df %>% filter(pshape != "in")
rng <- clamp %>% pivot_longer(c(lo, hi), values_to = "y") %>%
  crossing(stage = factor("HH13", levels = names(stage_lab)))

okabe <- c("duck" = "#000000", "chick best" = "#009E73", "chick worst" = "#E69F00")

p <- ggplot(df, aes(stage, yplot, color = group, group = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_blank(data = rng, aes(stage, y), inherit.aes = FALSE) +
  geom_quasirandom(aes(size = n_cells, shape = pshape), dodge.width = 0.8,
                   width = 0.13, alpha = 0.95, stroke = 0.8) +
  geom_text(data = off, aes(label = sprintf("%+.1f", signed)), position = position_dodge(0.8),
            vjust = ifelse(off$pshape == "hi", -0.7, 1.5), size = 2.1, show.legend = FALSE) +
  facet_grid(panel ~ ., scales = "free_y", switch = "y") +
  scale_x_discrete(labels = stage_lab) +
  scale_color_manual(values = okabe) +
  scale_shape_manual(values = c("in" = 1, "hi" = 2, "lo" = 6), guide = "none") +
  scale_size_area(max_size = 5.5, breaks = c(20, 100, 400), name = "cells (min L/R half)") +
  labs(x = NULL, y = NULL, color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3.2, shape = 1, stroke = 1.1)),
         size = guide_legend(override.aes = list(shape = 1, color = "grey30"))) +
  theme_minimal(base_size = 11) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "right", axis.text.x = element_text(size = 8))

png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
  function(filename, ...) grDevices::png(filename, ..., type = "cairo")
ggsave(file.path(HERE, "panel_G_xspecies_isoform.png"), p, width = 9.5, height = 8.6,
       units = "in", dpi = 300, device = png_dev)  # 300 = composite.py resolution (matches every sibling panel)
cat(sprintf("wrote panel_G_xspecies_isoform.png (%d panels, %d off-scale)\n",
            length(panel_lev), nrow(off)))
