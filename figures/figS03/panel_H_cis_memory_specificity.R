#!/usr/bin/env Rscript
# panel_H_cis_memory_specificity.R — figS03 panel H (cis-memory-node specificity).
#
# One volcano (the Ventral_spl HH21 lateral DEGs with a promoter in the union peak set — 79 DEGs /
# 78 testable genes on the 2026-08 seeded relabel, after the 2026-08-22 ±1 kb promoter rescue;
# 67 under the superseded containment-only rule): x = the HH17/16 (trough) promoter L/R effect
# (DESeq2 log2FC of the scan's PRIMARY design; gene = its most significant promoter); y = the
# PRIMARY-family BH-FDR significance (dashed line = FDR 0.1); colour = HH21 DEG side (blue L-up /
# red R-up). The scan writes `padj_primary` / `log2FC_primary` + `primary_family` /
# `primary_design` (default peak family, unpaired) — this panel reads those, so a family/design
# switch in the scan regenerates the panel consistently without edits here. Significance, not magnitude, drives the read. Points are circles except
# the +5.4 kb element (diamond = a flagged PITX2-locus distal reference, not a competitor).
# Full method + caveats: figures/figS03/LEGEND.md entry (H).
#
# No title/subtitle (those live in LEGEND.md). House style: B&W; blue/red only for the L/R side.
# Run: conda run -n duckfig Rscript figures/figS03/panel_H_cis_memory_specificity.R
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggforce); library(grid)
})

HERE <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 error = function(e) "figures/figS03")
if (is.na(HERE) || !nzchar(HERE)) HERE <- "figures/figS03"
ROOT <- dirname(dirname(HERE))
SCAN <- file.path(ROOT, "figures", "data", "grn_memory_node_promoter_scan.tsv")
OUT  <- file.path(HERE, "panel_H_cis_memory_specificity.png")
lat  <- local({ source(file.path(ROOT, "figures", "_shared", "colors.R"))
                get_palette("lateral", labels = c("L", "R")) })
PADJ <- 0.10; THR <- -log10(PADJ)
nlp <- function(p) -log10(pmax(p, 1e-300))

scan <- read_tsv(SCAN, show_col_types = FALSE)
# PITX2 promoter identities come from the data, NOT hardcoded peak-id coordinates (union-peak
# boundaries shift on a peak-set rerun): the 3 reference (is_ref) rows are PITX2 a/b/c, and the
# c promoter is flagged per-PEAK by `is_pitx2c` (identity, not "whichever row has the smallest
# q" — under the paired design PITX2a wins that race). a/b are the within-locus bilateral
# control. The +5.4 kb element = the two locus ncRNAs (pitx2_locus, non-PITX2).
stopifnot(all(c("padj_primary", "log2FC_primary", "primary_family", "primary_design",
                "is_pitx2c") %in% names(scan)))
PRIM <- sprintf("%s family, %s design", scan$primary_family[1], scan$primary_design[1])
# Under the GENE family only each gene's min-p row carries a q, so a/b can be NA there. Fall
# back to the same design's PEAK-family q for the two control markers rather than dropping
# them silently (ggplot removes NA points with only a warning).
pitx2 <- scan |> filter(is_ref) |> arrange(tss) |>
  mutate(padj_marker = ifelse(is.na(padj_primary),
                              .data[[paste0("padj_family_", scan$primary_design[1])]],
                              padj_primary),
         x = log2FC_primary, y = nlp(padj_marker), e54 = FALSE)
stopifnot(nrow(pitx2) == 3L, sum(pitx2$is_pitx2c) == 1L, !anyNA(pitx2$y))
C_row <- pitx2 |> filter(is_pitx2c)                     # PITX2c promoter, by peak identity
A_row <- pitx2 |> filter(!is_pitx2c) |> slice(1)        # PITX2a (smallest TSS of the remaining)
B_row <- pitx2 |> filter(!is_pitx2c) |> slice(2)        # PITX2b
e54_genes <- scan |> filter(pitx2_locus, !is_ref) |> distinct(gene_id) |> pull(gene_id)
stopifnot(length(e54_genes) == 2L)                    # the two +5.4 kb-element ncRNAs

reps <- scan |>
  filter(!is.na(padj_primary)) |>                       # gene family: only min-p rows carry a q
  arrange(padj_primary, desc(log2FC_primary)) |>
  group_by(gene_id) |> slice(1) |> ungroup() |>
  mutate(x = log2FC_primary, y = nlp(padj_primary),
         e54 = gene_id %in% e54_genes)
ab <- bind_rows(A_row, B_row)                         # PITX2 a/b promoters (within-locus control)
# `reps` keeps one row per gene (its min-q promoter), which under the GENE family x paired design
# is PITX2a — so the c promoter would be annotated but never drawn. Add it back when absent.
c_missing <- !any(reps$peak_id == C_row$peak_id)
base <- bind_rows(reps, ab, if (c_missing) C_row else NULL)
stopifnot(any(base$peak_id == C_row$peak_id))
# always label the nearest independent gene (data-driven; NOT gated on a y threshold — the
# LEGEND names it, so it must be on the panel whatever its q is)
near <- reps |> filter(!pitx2_locus) |> slice_min(padj_primary, n = 1, with_ties = FALSE)

# a/b fall back to the peak-family q when the gene family leaves them NA; say so on the label
# rather than presenting a peak-family value on an axis labelled "primary-family FDR".
ab_suffix <- if (any(is.na(bind_rows(A_row, B_row)$padj_primary))) " (peak-family q)" else ""
mk <- bind_rows(
  C_row |> transmute(x, y, mlab = "PITX2c promoter"),
  reps  |> filter(gene_id == e54_genes[1])   |> transmute(x, y, mlab = "PITX2 +5.4 kb element"),
  # a and b sit ~0.03 log2FC apart, so two auto-routed leaders both landed in the same blob and
  # neither could be assigned to its point. The legend treats them as ONE control, so label them
  # as one at the midpoint.
  tibble(x = mean(c(A_row$x, B_row$x)), y = mean(c(A_row$y, B_row$y)),
         mlab = paste0("PITX2a/b promoters", ab_suffix)),
  # the nearest independent gene's significance is READ, never asserted in the label
  near  |> transmute(x, y, mlab = sprintf("%s (nearest, q=%.2f)", gene_id, padj_primary)))

n_indep <- sum(!reps$pitx2_locus)
# A gene competes if ANY of its promoters is lateralized-concordant — the rule the scan
# headline and the LEGEND use. Counting only each gene's min-q representative row would
# disagree with them the moment a gene's best-q promoter is the discordant one.
n_sig <- scan |> filter(classification == "lateralized_concordant", !pitx2_locus) |>
  distinct(gene_id) |> nrow()
cnt_lab <- sprintf("%d / %d independent\ngenes above FDR 0.1", n_sig, n_indep)
# fixed label anchors (lx, ly) per marked point; hj = label justification at the anchor
lab_df <- mk |> mutate(
  lx = case_when(mlab == "PITX2c promoter"       ~ x - 0.12,
                 mlab == "PITX2 +5.4 kb element" ~ x + 0.30,
                 startsWith(mlab, "PITX2a/b")    ~ x - 0.35,
                 TRUE                            ~ x - 0.25),
  ly = case_when(mlab == "PITX2c promoter"       ~ y,
                 mlab == "PITX2 +5.4 kb element" ~ THR + 0.75,
                 startsWith(mlab, "PITX2a/b")    ~ y + 0.45,
                 TRUE                            ~ THR + 0.35),
  hj = 1)

p <- ggplot(base, aes(x, y)) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_hline(yintercept = THR, colour = "grey45", linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(colour = group, shape = e54, size = e54), stroke = 0.6, alpha = 0.9) +
  scale_colour_manual(values = c("L-up" = lat[["L"]], "R-up" = lat[["R"]]), name = "HH21 DEG side",
                      labels = c("L-up (HH21 Left)", "R-up (HH21 Right)")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 2.9), guide = "none") +
  # ggforce auto-placed labels collide with the FDR text and the axis notes on this point
  # cloud, so labels sit at fixed offsets from their
  # points in the empty regions (above the FDR line, or above the low-q cloud), so placement is
  # deterministic.  Re-check by eye if the scan changes shape.
  ggforce::geom_mark_circle(data = mk, aes(x = x, y = y, group = mlab),
    inherit.aes = FALSE, colour = "grey45", linewidth = 0.3, expand = unit(1.6, "mm")) +
  geom_segment(data = lab_df, aes(x = lx, y = ly, xend = x, yend = y), inherit.aes = FALSE,
               colour = "grey45", linewidth = 0.3) +
  geom_label(data = lab_df, aes(x = lx, y = ly, label = mlab, hjust = hj), inherit.aes = FALSE,
             size = 3.8, fontface = "bold", colour = "grey15", fill = alpha("white", 0.9),
             label.size = 0) +
  annotate("text", x = -Inf, y = THR, label = "  FDR 0.1", colour = "grey45", size = 3.7, hjust = 0, vjust = 1.4) +
  # count label anchored to the DATA range (a fixed y sat under the legend box once the y-range shrank)
  annotate("text", x = Inf, y = Inf, label = cnt_lab, hjust = 1.05, vjust = 1.3, size = 3.8, colour = "grey20", lineheight = 0.9) +
  annotate("text", x = -Inf, y = -Inf, label = paste0("scan primary: ", PRIM), hjust = -0.03, vjust = -2.2, size = 3.3, colour = "grey45") +
  annotate("text", x = Inf, y = -Inf, label = "trough Left →", hjust = 1.03, vjust = -0.6, size = 3.7, colour = "grey35") +
  annotate("text", x = -Inf, y = -Inf, label = "← trough Right", hjust = -0.03, vjust = -0.6, size = 3.7, colour = "grey35") +
  scale_y_continuous(expand = expansion(mult = c(0.22, 0.3))) +
  labs(x = "trough (HH17/16) promoter L/R effect  (log2FC)",
       y = "promoter lateralization significance\n(−log10 family-FDR)") +
  coord_cartesian(clip = "off") +
  # PNAS final font pass (2026-10): 300 dpi (compositor assumes 300) at the
  # Fig 5 left-column width, so this prints at x0.653; text >= ~7 pt printed.
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10.8),
        axis.title = element_text(size = 11.5, lineheight = 0.9),
        plot.margin = margin(4, 8, 4, 4),
        legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = "grey80", linewidth = 0.3),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10.8),
        legend.key.height = unit(0.4, "cm"), legend.margin = margin(3, 4, 3, 4))
ggsave(OUT, p, width = 5.067, height = 3.55, dpi = 300); cat("->", OUT, "\n")
