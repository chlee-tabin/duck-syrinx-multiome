# Tier 2 — laterality as a PER-CELL CP10K "expression boost" (L vs R), per-embryo paired.
# Rationale (user-directed): L and R are separately-sampled / separately-normalized
# libraries, so cell NUMBER / fraction is not comparable across the split — but the
# per-cell expression LEVEL is a clean within-cell claim. We therefore aggregate to the
# per-embryo mean per-cell CP10K (pooled over the syrinx-fated LR_POPS), pair L vs R by
# embryo (avoids pseudoreplication), and test with a paired Wilcoxon signed-rank.
# Style = project canonical: L = blue circle (#4393C3), R = red square (#D6604D),
# grey segments connect the same embryo's halves, black median crossbar.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})

MIN_CELLS_EMB <- 10   # per (embryo × side) cells needed for a stable per-cell mean

# legend_box: "horizontal" (default, guides side by side — fits a page-wide panel) or
# "vertical" (guides stacked — needed when the panel is rendered at half-page width, where
# a side-by-side Side + cells/embryo legend runs off the canvas).
# annotate_p: draw the per-facet "paired p=..." strings in-panel (TRUE, legacy) or
# leave the panel bare and let the figure legend carry the numbers (FALSE — the
# 2026-08-24 convention; the ann table still prints for legend transcription).
tier2_lr_boost <- function(d, genes, outfile,
                           stages = c("E09", "E11"), pops = LR_POPS,
                           width = 7.5, height = NULL, legend_box = "horizontal",
                           annotate_p = TRUE) {
  # every named gene (incl. the PITX2 positive control) must be present — fail loud, don't
  # silently drop a control row (screen_long only message()s on absence).
  miss <- genes[!paste0("cp10k__", genes) %in% colnames(d)]
  if (length(miss)) stop("tier2_lr_boost: requested genes absent from export: ", paste(miss, collapse = ", "))
  long <- screen_long(d, genes, "cp10k__") |>
    filter(!is.na(side), stage %in% stages, population %in% pops, !is.na(donor_id),
           donor_id != "", !donor_id %in% c("NA","Negative","Doublet","doublet","unassigned"))
  long$gene <- factor(long$gene, levels = genes[genes %in% unique(long$gene)])

  # per-embryo mean per-cell CP10K (pooled over syrinx-fated populations)
  per_emb <- long |>
    group_by(stage, gene, donor_id, side) |>
    summarise(mean_cp10k = mean(value), n_cells = n(), .groups = "drop") |>
    filter(n_cells >= MIN_CELLS_EMB)

  pk <- per_emb |> group_by(stage, gene, donor_id) |>
    summarise(paired = any(side=="Left") & any(side=="Right"), .groups="drop")
  per_emb <- per_emb |> left_join(pk, by = c("stage","gene","donor_id"))
  # paired wide table (one row per donor with both halves); ensure both side columns exist
  # even if a thin subset yields no pairs (avoids a missing-`Left`/`Right`-column crash).
  wide <- per_emb |> filter(paired) |>
    select(stage, gene, donor_id, side, mean_cp10k) |>
    pivot_wider(names_from = side, values_from = mean_cp10k)
  for (s in c("Left", "Right")) if (!s %in% names(wide)) wide[[s]] <- NA_real_
  seg <- wide

  # paired Wilcoxon per (stage, gene)
  sig <- function(p) ifelse(is.na(p), "", ifelse(p<.001, "***", ifelse(p<.01, "**", ifelse(p<.05, "*", "ns"))))  # vectorized
  ann <- wide |>
    group_by(stage, gene) |>
    summarise(n_pairs = sum(!is.na(Left) & !is.na(Right)),
              p = if (sum(!is.na(Left)&!is.na(Right)) >= 3)
                    suppressWarnings(wilcox.test(Left, Right, paired=TRUE, exact=FALSE)$p.value)
                  else NA_real_,
              med_LmR = median(Left - Right, na.rm=TRUE), .groups="drop") |>
    mutate(lab = ifelse(is.na(p), sprintf("paired n=%d", n_pairs),
                        sprintf("paired p=%.2g %s (n=%d)", p, sig(p), n_pairs)))

  # registry keys are L/R; map onto our Left/Right side labels
  side_cols <- c(Left = get_color("lateral","L"), Right = get_color("lateral","R"))
  side_shp  <- c(Left = 16, Right = 15)

  p <- ggplot(per_emb, aes(side, mean_cp10k)) +
    geom_segment(data = seg, aes(x="Left", xend="Right", y=Left, yend=Right),
                 inherit.aes=FALSE, color="gray70", linewidth=0.3, alpha=0.7) +
    geom_point(aes(color=side, shape=side, size=n_cells), alpha=0.85) +
    stat_summary(fun=median, geom="crossbar", width=0.5, linewidth=0.4,
                 color="black", show.legend=FALSE) +
    (if (annotate_p) geom_text(data=ann, aes(x=1.5, y=Inf, label=lab),
                               inherit.aes=FALSE, vjust=1.4, size=2.5, color="gray20")) +
    scale_color_manual(values=side_cols, name="Side") +
    scale_shape_manual(values=side_shp, name="Side") +
    scale_size_continuous(range=c(1.3,3.8), name="cells/embryo") +
    facet_grid(gene ~ stage, scales="free_y", switch="y") +
    scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0, if (annotate_p) 0.28 else 0.08))) +
    labs(x=NULL, y="Mean per-cell expression (CP10K)") +   # title/observations -> LEGEND.md
    theme_classic(base_size=10) +
    theme(strip.text.y.left=element_text(angle=0, face="italic"),
          strip.text.x=element_text(face="bold"), strip.background=element_blank(),
          strip.placement="outside", legend.position="top", legend.box=legend_box,
          panel.spacing=unit(0.5,"lines"), plot.margin=margin(16,6,6,6))

  ng <- length(levels(per_emb$gene))
  if (is.null(height)) height <- 1.4 + 1.5 * ng
  ggsave(outfile, p, width=width, height=height, dpi=200, limitsize=FALSE)
  cat(sprintf("Saved %s\n", outfile))
  print(ann)
  invisible(ann)
}
