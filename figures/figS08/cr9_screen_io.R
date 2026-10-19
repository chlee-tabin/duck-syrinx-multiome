# Shared loader for the cr9 WNT/BMP/TGF-β screen figures (figS08/figS09/figS10).
# Joins the per-cell screen export -> compendium population labels, drops technical
# Contamination, renames "Somatic" -> "Skeletogenic LPM" (paper-wide convention), and
# returns a tidy per-cell table carrying log + CP10K for every screen gene + side/stage/
# sex/embryo.  Fail-loud: an unmapped cluster (NA population) stops() naming the cluster.
#
# source("figures/figS08/cr9_screen_io.R") then:  d <- load_cr9_screen()
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

DATA_DIR <- "figures/data/scrna_e6e9e11"

# Display order across stages (absent ones drop out per facet via scales="free_x")
POP_ORDER <- c("Ventral spl.", "Dorsal sub-epi.", "Syringeal chondrocyte",
               "Syringeal fibroblast", "Syringeal smooth muscle",
               "Skeletogenic LPM", "Mesothelium")
# Syrinx-fated populations used for the L/R "expression boost" tier
LR_POPS  <- c("Ventral spl.", "Syringeal chondrocyte", "Syringeal fibroblast",
              "Skeletogenic LPM")

load_cr9_screen <- function() {
  # `population` is AUTHORITATIVE in the export: it was joined on-cluster against the
  # run-matched clusterlabels, committed for provenance as
  # figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv.
  # Do NOT re-join these exports to the earlier cr6-era substage_compendium_clusterlabels.tsv
  # (NO cr9_ prefix; not distributed) — it is a DIFFERENT substage run (10,160 vs 8,736
  # E09 cells) with shuffled IDs.
  frames <- lapply(c("E06", "E09", "E11"), function(s) {
    f <- file.path(DATA_DIR, sprintf("cr9_screen_%s_percell.tsv.gz", s))
    if (!file.exists(f)) stop(sprintf("missing screen export: %s (run cr9_screen_export.R)", f))
    read.delim(gzfile(f), stringsAsFactors = FALSE, check.names = FALSE, quote = "")  # producer writes quote=FALSE
  })
  d <- bind_rows(frames)
  if (!"population" %in% colnames(d)) stop("export lacks `population` column — re-run cr9_screen_export.R")
  d <- d[!d$population %in% c("Contamination", "Cycling"), , drop = FALSE]  # technical (ambient junk / cell-cycle), not lineages
  d$population[d$population == "Somatic"] <- "Skeletogenic LPM"      # paper-wide rename (compat guard; no-op on current export)
  d$population <- factor(d$population, levels = POP_ORDER)
  if (any(is.na(d$population))) stop("population factor has unexpected level not in POP_ORDER")
  d$stage <- factor(d$stage, levels = c("E06", "E09", "E11"))
  d$side  <- factor(ifelse(d$LR.called %in% c("Left","Right"), d$LR.called, NA),
                    levels = c("Left","Right"))   # real dissection sides only
  d
}

# Pivot a set of genes to long for violin plotting.  prefix = "log__" or "cp10k__".
screen_long <- function(d, genes, prefix = "log__") {
  cols <- paste0(prefix, genes)
  have <- cols[cols %in% colnames(d)]
  miss <- setdiff(cols, colnames(d))
  if (length(miss)) message("screen_long: genes absent from export (skipped): ",
                            paste(sub(prefix, "", miss), collapse=", "))
  meta <- c("cell_id","stage","population","side","sex.called","embryo_id","donor_id")
  meta <- intersect(meta, colnames(d))
  long <- d[, c(meta, have)] |>
    pivot_longer(all_of(have), names_to = "gene", values_to = "value") |>
    mutate(gene = sub(prefix, "", gene))
  long$gene <- factor(long$gene, levels = genes[paste0(prefix, genes) %in% have])
  long
}

# Shared violin panel: fill = subpopulation (canonical registry colours), facet passed in
# (must use space="free_x" so absent subpopulations take no column at a stage), no embedded
# title/subtitle (those live in LEGEND.md), top margin so the compositor's letter label has room.
# Requires ggplot2 + figures/_shared/colors.R sourced by the caller.
screen_violin <- function(long, facet, outfile, w, h) {
  pops <- levels(droplevels(long$population))
  cols <- get_palette("subcluster", labels = pops)   # fail-loud on any unregistered population
  p <- ggplot2::ggplot(long, ggplot2::aes(population, value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = population), color = "black",
                         linewidth = 0.22, scale = "width") +
    ggplot2::stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    facet +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = "Expression level (log-normalized)") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 8),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "italic", size = 8),
      strip.text.x = ggplot2::element_text(face = "bold", size = 9),
      strip.background = ggplot2::element_blank(), strip.placement = "outside",
      panel.spacing.x = ggplot2::unit(0.5, "lines"),
      panel.spacing.y = ggplot2::unit(0.25, "lines"),
      plot.margin = ggplot2::margin(16, 16, 6, 6))
  ggplot2::ggsave(outfile, p, width = w, height = h, dpi = 200, limitsize = FALSE)
  cat("Saved", outfile, "\n")
}
