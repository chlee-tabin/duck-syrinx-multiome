#!/usr/bin/env Rscript
# figS04 Panel F: orientation-free PITX2 half-vs-half by spatial domain (duck fig02C analog).
#
# R/ggplot2 port of the former generate_dv_laterality.py, styled to MATCH its duck sibling
# figures/fig02/generate_paired_laterality.R (theme_bw, grey95 facet strips, per-facet square
# free-scale facets with a y=x bilateral-reference diagonal, cowplot 2-row stack + shared legend).
# This is a THIN renderer of the pre-computed cAN_halfvshalf.tsv (the cAN step already did the
# per-embryo x domain half-vs-half folding) — it does NOT recompute fractions.
#
# Per embryo x spatial DOMAIN, PITX2 in one half vs the other half; distance from the y=x diagonal
# = PITX2 asymmetry REGARDLESS of which half is Left or Right.  THREE column-facets = spatial
# domains (folded from the stage-restricted compendium populations, exactly as the duck):
#   * "LPM mes. (HH13)"          — whole HH13 LPM-lineage origin (pre-TES; D/V unresolved)
#   * "Ventral splanchnic pool"  — HH15-16 + HH25 inference pool (where laterality lives)
#   * "Dorsal sub-epi."          — the bilateral negative control (expected ON the diagonal)
# TWO metric rows: PITX2+ FRACTION (top, featured) + mean CP10K (bottom; A vs B are different
# cellranger products so absolute CP10K is not strictly cross-comparable).
#
# EVERY embryo contributes (project convention: "Chick Control Conventions"): Group A (ARC01/ARC02) AND Group
# B (GEX09/GEX10), WT AND GFP+ (GFP = pan-GFP demux label, NO L/R relevance, identical to WT).
# CHICK-SPECIFIC ENCODINGS preserved from the matplotlib version:
#   * Group B = FILLED (fill = stage, black edge); Group A = OPEN (white fill, stage-coloured edge).
#   * size = min cells per half (the reliability cue; clipped at 120 so big halves don't dominate).
#   * donor labels ONLY on design-orientable embryos (uniform/singleton Group B classes), coloured
#     by DESIGN CONCORDANCE — is PITX2+ fraction higher in the embryo's DESIGNED-LEFT half?
#       BLACK = concordant (PITX2-high in designed-LEFT) ; RED = discordant.
#     Un-orientable embryos (Group A genetically-identical halves; non-uniform Group B classes) are
#     plotted but NOT labelled.
#
# All explanatory prose lives in LEGEND.md (no in-image title/caption per the figS04 house rule).
#
# Reads:  figures/data/chick_control/cAN_halfvshalf.tsv
# Writes: figures/figS04/panel_F_dv_laterality.png  (single combined 2-row PNG)
#   conda run -n duckfig Rscript figures/figS04/panel_F_dv_laterality.R
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel); library(cowplot)
})

HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(HERE) == 0) HERE <- "figures/figS04"
TSV <- file.path(HERE, "..", "data", "chick_control", "cAN_halfvshalf.tsv")

# ── Style (duck-consistent: same STAGE_COLORS hues as _shared/style.py / fig02) ──
DOMAIN_ORDER <- c("LPM mes. (HH13)", "Ventral splanchnic pool", "Dorsal sub-epi.")
STAGE_COLORS <- c("HH13" = "#66a61e", "HH15-16" = "#b0a8d0", "HH25" = "#1b9e77")
STAGE_ORDER  <- c("HH13", "HH15-16", "HH25")
CONC_RED <- "#c0392b"
SIZE_CAP <- 120   # clip min_cells for the size aesthetic (approximates the matplotlib; ggplot size is radius- not area-scaled)

df <- read_tsv(TSV, show_col_types = FALSE)

# ── Fail LOUD, never silent-default (project convention: figure-code rule) ──
bad_stage  <- setdiff(unique(df$stage),  names(STAGE_COLORS))
if (length(bad_stage))  stop("figS04 F: stage(s) with no STAGE_COLORS entry: ", paste(bad_stage,  collapse = ", "))
bad_domain <- setdiff(unique(df$domain), DOMAIN_ORDER)
if (length(bad_domain)) stop("figS04 F: domain(s) not in DOMAIN_ORDER: ",      paste(bad_domain, collapse = ", "))

domains <- DOMAIN_ORDER[DOMAIN_ORDER %in% df$domain]
df <- df %>%
  mutate(
    domain    = factor(domain, levels = domains),
    stage     = factor(stage,  levels = STAGE_ORDER),
    edge_col  = STAGE_COLORS[as.character(stage)],                 # Group-A edge (identity colour)
    size_clip = pmin(min_cells, SIZE_CAP),
    donor_lab = sub("donor", "", donor),
    label_col = ifelse(concordant == "yes", "black",
                ifelse(concordant == "no",  CONC_RED, NA_character_))  # NA = un-orientable → unlabelled
  )
nconc <- sum(df$concordant == "yes", na.rm = TRUE); ndisc <- sum(df$concordant == "no", na.rm = TRUE)
message(sprintf("  %d points; A=%d B=%d, GFP+=%d; design-concordant=%d discordant=%d",
                nrow(df), sum(df$group == "A"), sum(df$group == "B"),
                sum(df$gfp == "GFP+"), nconc, ndisc))

# ── One faceted metric plot (h1 vs h2), styled like the duck ──
make_plot <- function(h1col, h2col, xlab, ylab, pct = FALSE, show_legend = FALSE) {
  d <- df %>% mutate(x = .data[[h1col]], y = .data[[h2col]]) %>% filter(!is.na(x), !is.na(y))
  fmt <- if (pct) scales::percent_format(accuracy = 1) else waiver()

  # Per-facet symmetric limits via invisible geom_blank → square panels with the diagonal
  # corner-to-corner (aspect.ratio = 1 makes each facet render square under scales = "free").
  lim <- d %>% group_by(domain) %>%
    summarise(lim = max(c(x, y), na.rm = TRUE) * 1.10, .groups = "drop")
  blanks <- bind_rows(transform(lim, x = lim, y = 0), transform(lim, x = 0, y = lim))

  dA <- d %>% filter(group == "A"); dB <- d %>% filter(group == "B")
  lab <- d %>% filter(!is.na(label_col))

  ggplot() +
    facet_wrap(~ domain, nrow = 1, scales = "free") +
    geom_blank(data = blanks, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.3) +
    # Group A = open (white fill, stage-coloured edge via identity colour scale)
    geom_point(data = dA, aes(x = x, y = y, color = edge_col, size = size_clip),
               shape = 21, fill = "white", stroke = 0.7, alpha = 0.9) +
    # Group B = filled (fill = stage, black edge)
    geom_point(data = dB, aes(x = x, y = y, fill = stage, size = size_clip),
               shape = 21, color = "black", stroke = 0.3, alpha = 0.85) +
    geom_text_repel(data = lab, aes(x = x, y = y, label = donor_lab, color = label_col),
                    size = 1.7, fontface = "bold", segment.color = "grey60", segment.size = 0.2,
                    min.segment.length = 0.3, max.overlaps = 20, box.padding = 0.28,
                    point.padding = 0.15, seed = 42, show.legend = FALSE) +
    scale_fill_manual(values = STAGE_COLORS, limits = STAGE_ORDER, drop = FALSE, name = "Stage") +
    scale_color_identity() +          # Group-A edges + concordance labels use literal colours
    scale_size_continuous(range = c(1.3, 3.8), limits = c(0, SIZE_CAP),
                          breaks = c(10, 40, SIZE_CAP), labels = c("10", "40", "≥120"),
                          name = "Min. cells\n(per half)") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4), labels = fmt,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 4), labels = fmt,
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = xlab, y = ylab) +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.5, color = "black")),
           size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey60", color = "grey35"))) +
    theme_bw(base_size = 7) +
    theme(strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = "grey95"),
          axis.text = element_text(size = 6), axis.title = element_text(size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
          legend.position = if (show_legend) "right" else "none",
          legend.key.size = unit(0.28, "cm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          aspect.ratio = 1, panel.spacing.x = unit(0.6, "lines"),
          plot.margin = margin(3, 5, 3, 3))
}

# ── Manual key: Group A/B (open/filled) + concordance-label colours (mirrors the matplotlib hand-built legend) ──
key_lv <- c("Group A (ARC01/02), open", "Group B (designed L/R), filled",
            "donor label: concordant", "donor label: discordant")
key_df <- data.frame(k = factor(key_lv, levels = key_lv), x = 1, y = seq_along(key_lv))
p_key <- ggplot(key_df, aes(x, y, shape = k, color = k)) +
  geom_point(size = 2.6, stroke = 0.9, fill = "grey55") +
  scale_shape_manual(name = NULL, values = c(21, 21, 15, 15)) +
  scale_color_manual(name = NULL, values = c("grey25", "black", "black", CONC_RED)) +
  guides(shape = guide_legend(override.aes = list(
           fill = c("white", "grey55", "black", CONC_RED)))) +
  theme_void() + theme(legend.position = "right", legend.text = element_text(size = 6),
                       legend.key.size = unit(0.30, "cm"))
key_legend <- get_legend(p_key)

# ── Build both metric rows; stack with the shared legend on the right ──
p1     <- make_plot("frac_h1", "frac_h2",
                    "PITX2+ fraction\nhalf 1 (ARC01 / GEX09)",
                    "PITX2+ fraction\nhalf 2 (ARC02 / GEX10)", pct = TRUE,  show_legend = FALSE)
p1_leg <- make_plot("frac_h1", "frac_h2", "", "", pct = TRUE, show_legend = TRUE)
p2     <- make_plot("mean_h1", "mean_h2",
                    "mean PITX2 CP10K\nhalf 1 (ARC01 / GEX09)",
                    "mean PITX2 CP10K\nhalf 2 (ARC02 / GEX10)", pct = FALSE, show_legend = FALSE)

native_legend <- get_legend(p1_leg)                          # Stage (fill) + Min. cells (size)
stacked    <- plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")
legend_col <- plot_grid(native_legend, key_legend, ncol = 1, rel_heights = c(1.7, 1.0), align = "v")
combined   <- plot_grid(stacked, legend_col, ncol = 2, rel_widths = c(8.6, 1.9))

out <- file.path(HERE, "panel_F_dv_laterality.png")
ggsave(out, combined, width = 10.5, height = 5.7, dpi = 300, bg = "white")
message("Saved: ", out)
