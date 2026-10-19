# Shared PITX2 gene-model constants + exon-drawing layer for pileup figures.
#
# Used by figures/fig02/generate_pitx2_pileup_multistage.R and an internal
# strand-mirror diagnostic (not distributed) — keeping the constants +
# draw_pitx2_gene_model() function here ensures that a PITX2 annotation
# correction (e.g., UTR boundary revision) updates every pileup together, and
# that they stay in sync on exon boundaries and CDS starts.
#
# The tracks tibble is NOT moved here because the two call sites
# differ in strand-split behavior (mirror adds fwd/rev).  Isoform
# constants and gene-model drawing are the real duplication risk.

# ── Isoform coordinates (ZJU1.0 / crarc_zju1_PITX2, NC_051775.1) ──────────
PITX2_GENE_START       <- 3961387
PITX2_GENE_END         <- 3976254
PITX2_CDS_END          <- 3975349
PITX2_UTR3_START       <- 3975353
PITX2_UTR3_END         <- 3976254
PITX2AB_CDS_START      <- 3966601
PITX2C_CDS_START       <- 3971830

PITX2_ISOFORMS <- list(
    PITX2a = list(
        label = "PITX2a", y = 3,
        ex1 = c(3961387, 3961458),
        exons = list(c(3966591, 3966646), c(3973108, 3973313),
                     c(3974789, 3976254))
    ),
    PITX2b = list(
        label = "PITX2b", y = 2,
        ex1 = c(3964253, 3964353),
        exons = list(c(3966591, 3966646), c(3973108, 3973313),
                     c(3974789, 3976254))
    ),
    PITX2c = list(
        label = "PITX2c", y = 1,
        ex1 = c(3971385, 3972064),
        exons = list(c(3973108, 3973313), c(3974789, 3976254))
    )
)


#' Append PITX2 gene-model exon / UTR / intron geoms to a ggplot.
#'
#' Takes a base ggplot (typically with highlight rectangles already
#' layered) and returns the plot with the 3-isoform gene-model drawn
#' at y = 1 (PITX2c), 2 (PITX2b), 3 (PITX2a).  Matches the rendering
#' used by fig02 panel D and the exploratory strand-mirror plot
#' exactly, so a correction here propagates to both figures.
#'
#' @param p A ggplot object (base layer).
#' @param exon_h Half-height of CDS boxes (default 0.30).
#' @param utr_h  Half-height of UTR boxes (default 0.18).
#' @return The ggplot with gene-model annotations added.
draw_pitx2_gene_model <- function(p, exon_h = 0.30, utr_h = 0.18, label_size = 2.5) {
    # label_size: isoform-label size in mm (ggplot units). Default 2.5 = the
    # historical value, so callers that do not pass it render unchanged.
    for (iso_name in names(PITX2_ISOFORMS)) {
        iso <- PITX2_ISOFORMS[[iso_name]]
        y <- iso$y
        ex1_s <- iso$ex1[1]
        ex1_e <- iso$ex1[2]
        first_shared <- iso$exons[[1]]

        # Intron line from ex1 to first shared exon
        p <- p + annotate("segment",
                          x = ex1_e / 1e3,
                          xend = first_shared[1] / 1e3,
                          y = y, yend = y, linewidth = 0.4)

        # Intron lines between downstream exons
        for (j in seq_along(iso$exons)) {
            if (j < length(iso$exons)) {
                p <- p + annotate("segment",
                    x = iso$exons[[j]][2] / 1e3,
                    xend = iso$exons[[j + 1]][1] / 1e3,
                    y = y, yend = y, linewidth = 0.4)
            }
        }

        # First exon: PITX2c has a UTR-then-CDS split within ex1;
        # PITX2a and PITX2b are all UTR in their ex1 (CDS starts at
        # PITX2AB_CDS_START in the second exon).
        if (iso_name == "PITX2c") {
            p <- p +
                annotate("rect",
                         xmin = ex1_s / 1e3,
                         xmax = PITX2C_CDS_START / 1e3,
                         ymin = y - utr_h, ymax = y + utr_h,
                         fill = "#606060", color = "black",
                         linewidth = 0.2) +
                annotate("rect",
                         xmin = PITX2C_CDS_START / 1e3,
                         xmax = ex1_e / 1e3,
                         ymin = y - exon_h, ymax = y + exon_h,
                         fill = "#303030", color = "black",
                         linewidth = 0.2)
        } else {
            p <- p +
                annotate("rect",
                         xmin = ex1_s / 1e3, xmax = ex1_e / 1e3,
                         ymin = y - utr_h, ymax = y + utr_h,
                         fill = "#606060", color = "black",
                         linewidth = 0.2)
        }

        # Downstream exons
        for (j in seq_along(iso$exons)) {
            ex_s <- iso$exons[[j]][1]
            ex_e <- iso$exons[[j]][2]

            if (j == length(iso$exons)) {
                # Last exon: CDS then 3'UTR (shared across isoforms)
                p <- p +
                    annotate("rect",
                             xmin = ex_s / 1e3,
                             xmax = PITX2_CDS_END / 1e3,
                             ymin = y - exon_h, ymax = y + exon_h,
                             fill = "#303030", color = "black",
                             linewidth = 0.2) +
                    annotate("rect",
                             xmin = PITX2_UTR3_START / 1e3,
                             xmax = PITX2_UTR3_END / 1e3,
                             ymin = y - utr_h, ymax = y + utr_h,
                             fill = "#606060", color = "black",
                             linewidth = 0.2)
            } else if (ex_s == 3966591) {
                # PITX2a/b shared exon 2: UTR up to PITX2AB_CDS_START,
                # then CDS from PITX2AB_CDS_START to ex_e.  (PITX2c
                # skips this exon entirely.)
                p <- p +
                    annotate("rect",
                             xmin = ex_s / 1e3,
                             xmax = PITX2AB_CDS_START / 1e3,
                             ymin = y - utr_h, ymax = y + utr_h,
                             fill = "#606060", color = "black",
                             linewidth = 0.2) +
                    annotate("rect",
                             xmin = PITX2AB_CDS_START / 1e3,
                             xmax = ex_e / 1e3,
                             ymin = y - exon_h, ymax = y + exon_h,
                             fill = "#303030", color = "black",
                             linewidth = 0.2)
            } else {
                # Shared CDS-only exons
                p <- p +
                    annotate("rect",
                             xmin = ex_s / 1e3, xmax = ex_e / 1e3,
                             ymin = y - exon_h, ymax = y + exon_h,
                             fill = "#303030", color = "black",
                             linewidth = 0.2)
            }
        }

        # Isoform label + TSS arrow
        p <- p +
            annotate("text", x = ex1_s / 1e3 - 0.25, y = y,
                     label = iso$label, size = label_size,
                     fontface = "italic", hjust = 1) +
            annotate("segment",
                     x = ex1_s / 1e3, xend = ex1_s / 1e3,
                     y = y + exon_h + 0.05, yend = y + exon_h + 0.3,
                     linewidth = 0.3, color = "black") +
            annotate("segment",
                     x = ex1_s / 1e3, xend = ex1_s / 1e3 + 0.4,
                     y = y + exon_h + 0.3, yend = y + exon_h + 0.3,
                     linewidth = 0.3, color = "black",
                     arrow = arrow(length = unit(1.2, "mm"),
                                   type = "closed"))
    }

    p
}
