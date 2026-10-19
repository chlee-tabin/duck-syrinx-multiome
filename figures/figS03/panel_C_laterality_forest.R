#!/usr/bin/env Rscript
# figS03 Panel C: per-subcluster PITX2 L/R laterality TRAJECTORY (R/ggplot2).
#
# Faceted BY SUBCLUSTER (one facet per d2g stage-restricted LPM population), x =
# developmental stage (HH13 -> HH17/16 -> HH21 -> HH25), y = median per-embryo
# log2(L/R) of total PITX2 (CP10K) with a bootstrap 95% CI (resample embryos with
# replacement, B=2000).  Each facet is one lineage's laterality trajectory over the
# stage window it exists in; the connecting line is the trajectory, the point/CI
# colour is the read: a CI excluding 0 = lateralized (Left = blue circle, Right = red
# square — the house L/R glyphs), a CI spanning 0 = bilateral (grey, open circle).  hline at 0 = no laterality.
#
# The effect axis is CAPPED at +5 (coord_cartesian, NOT scale limits — hard limits
# would DROP off-scale CIs).  Early-stage log2FC inflates when R~0 (Cardiopharyngeal
# HH13 ~+8.5); those points/CIs are marked off-scale (^) at the +5 ceiling.
#
# Roster (grouped): inference pool / syrinx lineage (Early LPM, Cardiopharyngeal LPM,
# BARX1+, Ventral spl., Syringeal — the IRX1+ syrinx precursor, HH25 cl3) · Skeletogenic
# LPM (excluded) · dorsal control (Dorsal-biased LPM @HH17/16, Dorsal sub-epi. @HH21/25).
# Off-axis / non-pool pops (Contamination / Mesothelium / Venous-pole spl. /
# PRRX1+/FOXF1-low mes.) are omitted on axis/lineage grounds — they remain in the
# exported TSV (see LEGEND).
#
# The message D justifies: the syrinx lineage stays Left and becomes significant at the
# fate stages; Skeletogenic LPM (branchiomeric-skeletal, HH17/16 onward — the 2026-07-02
# vet reclassified the earlier HH13 "Skeletogenic" node as a PRRX1+/FOXF1-low fibroblast)
# is Left early then resolves to bilateral at the fate stages (excluded on fate-stage
# bilaterality, NOT blanket symmetry); the dorsal control stays bilateral.  Side basis =
# canonical d2n-merged sidecar (chromatin d2n @ HH17/16; snRNA d2l elsewhere): a
# per-embryo DIRECTION sign test.  At HH17/16 the sides are called from PITX2c-promoter
# accessibility, which tracks PITX2 expression in the same halves, so a filled HH17/16
# point is NOT independent evidence of RNA laterality (other modality, same gene).
#
# No in-image title/subtitle (figS03 house rule) — all text lives in LEGEND.md.
#
# Reads : figures/data/laterality_trajectory.tsv
#         (export_laterality_trajectory.py; metric=="PITX2_mean"; median_log2FC + log2fc_ci_lo/hi)
# Writes: figures/figS03/panel_C_laterality_forest.png
# Local : conda run -n duckfig Rscript figures/figS03/panel_C_laterality_forest.R

local_lib <- "~/R/4.4.2"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ggh4x)
})

HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
REPO <- normalizePath(file.path(HERE, "..", ".."))
# exporter: figures/figS03/export_laterality_trajectory.py
TSV  <- file.path(REPO, "figures", "data", "laterality_trajectory.tsv")
source(file.path(REPO, "figures", "_shared", "colors.R"))
source(file.path(REPO, "figures", "_shared", "subcluster_map.R"))  # display_population
lr_cols <- get_palette("lateral", labels = c("L", "R"))
L_BLUE <- unname(lr_cols["L"]); R_RED <- unname(lr_cols["R"])   # registry: L #4393C3, R #D6604D
NS_GREY <- "#9AA3AB"

STAGES <- c("HH13", "HH17/16", "HH21", "HH25")
CAP    <- 5      # effect-axis ceiling (off-scale marked, not dropped)
YMIN   <- -2     # floor; the one roster CI below it (Skeletogenic LPM HH25, lower bound -2.03) gets the off-scale caret

# Focused roster grouped into the three pool roles.  Order = facet fill order
# (row-major, ncol=3): inference pool rows 1-2, then Skeletogenic + dorsal on row 3.
ROSTER <- tibble::tribble(
    ~population,             ~group,
    "Early LPM",             "Inference pool",
    "Cardiopharyngeal LPM",  "Inference pool",
    "BARX1+",                "Inference pool",
    "Ventral spl.",          "Inference pool",
    "Syringeal",             "Inference pool",
    "Skeletogenic LPM",      "Skeletogenic (excluded)",
    "Dorsal-biased LPM",     "Dorsal control",
    "Dorsal sub-epi.",       "Dorsal control",
)
GROUP_LEVELS <- c("Inference pool", "Skeletogenic (excluded)", "Dorsal control")

d <- read.delim(TSV, check.names = FALSE)
d <- d[d$metric == "PITX2_mean" & d$n_emb > 0, ]
d <- merge(d, ROSTER, by = "population")          # inner join => drops off-axis pops
stopifnot(nrow(d) > 0, all(ROSTER$population %in% d$population))

# SSOT guard: the hand-assigned display `group` must agree with the `role` the
# export derived from the canonical POOLS (figures/_shared/subcluster_map.py).  A
# d2g rerun that moves a population between pools updates `role` in the TSV; this
# fails LOUD if ROSTER's display group drifts out of sync, rather than silently
# rendering a population under the wrong facet group.
ROLE2GROUP <- c(inference = "Inference pool", bilateral = "Dorsal control",
                excluded = "Skeletogenic (excluded)")
stopifnot(all(d$role %in% names(ROLE2GROUP)))
.drift <- d$group != ROLE2GROUP[d$role]
if (any(.drift)) stop("panel C group/role drift (POOLS pool membership changed?): ",
    paste(sprintf("%s role=%s but ROSTER group=%s",
                  d$population[.drift], d$role[.drift], d$group[.drift]), collapse = "; "))
# Off-scale handling now mirrors the +5 ceiling at the YMIN floor (the seeded
# 2026-08 partition's small excluded facets — n=3 embryo bootstraps — produce
# degenerate CIs below -2, which the pre-seed roster never did): medians below
# YMIN draw as a down-triangle AT the floor, clipped CI bottoms get a down caret.

# Facet labels, asserted 1:1 with the roster — factor(labels=) would otherwise merge
# two populations into a single facet without warning.
.LBL <- display_population(ROSTER$population, short = TRUE)
stopifnot(!anyDuplicated(.LBL), length(.LBL) == nrow(ROSTER))

d <- d %>%
    mutate(
        stage = factor(stage, levels = STAGES),
        # Facet strips carry the reader-facing label. SSOT levels drive ordering and the
        # join above; `labels` changes what is DRAWN — but note factor()
        # COLLAPSES duplicate labels into one level, which would silently merge two
        # populations into one facet and run geom_line through both. .LBL is asserted
        # collision-free above, so the relabel is 1:1.
        population = factor(population, levels = ROSTER$population, labels = .LBL),
        group = factor(group, levels = GROUP_LEVELS),
        excl0 = (log2fc_ci_lo > 0) | (log2fc_ci_hi < 0),
        sig   = factor(ifelse(excl0 & median_log2FC > 0, "Left",
                       ifelse(excl0 & median_log2FC < 0, "Right", "bilateral")),
                       levels = c("Left", "Right", "bilateral")),
        col   = dplyr::recode(as.character(sig), Left = L_BLUE, Right = R_RED,
                              bilateral = NS_GREY),
        # Off-scale glyph vocabulary (the +5 cap; see header): filled circle =
        # Left-lateralized in-range · filled square = Right-lateralized in-range (house
        # L = circle / R = square) · open circle = bilateral in-range · filled ▲ (shape 24)
        # = median capped at +5 · separate ▲ caret (shape 17) = median in-range but CI
        # top clipped.  The point is placed AT the ceiling (y_pt) so coord_cartesian
        # does not clip a capped point out of existence; CIs stay raw (coord clips them).
        y_pt        = pmax(pmin(median_log2FC, CAP), YMIN),
        pt_offscale = median_log2FC > CAP,                  # median itself off the top
        pt_offbot   = median_log2FC < YMIN,                 # median itself off the bottom
        ci_offtop   = log2fc_ci_hi > CAP & median_log2FC <= CAP,  # only the CI off the top
        ci_offbot   = log2fc_ci_lo < YMIN & median_log2FC >= YMIN,  # only the CI off the bottom
        # in range: square (22) when Right-lateralized, else circle (21); up-triangle (24)
        # capped at top, down-triangle (25) at the floor
        pt_shape    = ifelse(pt_offscale, 24L, ifelse(pt_offbot, 25L,
                             ifelse(sig == "Right", 22L, 21L))),
        # filled when lateralized (sig) OR off-scale; open white when bilateral
        pt_fill     = ifelse(sig == "bilateral" & !pt_offscale & !pt_offbot, "white", col)
    )
# fail LOUD if any kept row lacks a sig/colour (the recode has no silent default —
# an NA sig would otherwise draw an invisible point, violating the no-silent-drop rule)
stopifnot(!any(is.na(d$sig)), !any(is.na(d$col)))

# caret marks at the ceiling/floor for CIs clipped there (median in range)
caret <- d %>% filter(ci_offtop) %>% mutate(y_caret = CAP)
caret_dn <- d %>% filter(ci_offbot) %>% mutate(y_caret = YMIN)

p <- ggplot(d, aes(stage, y_pt, group = population)) +
    # L/R semantic background wash (faint; the colour message is on the points)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = L_BLUE, alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = R_RED, alpha = 0.05) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey25") +
    # trajectory line (neutral; significance is carried by the points, which can change
    # category along a lineage e.g. Skeletogenic Left -> bilateral)
    geom_line(linewidth = 0.5, colour = "grey55") +
    geom_linerange(aes(ymin = log2fc_ci_lo, ymax = log2fc_ci_hi, colour = col),
                   linewidth = 0.7) +
    geom_point(aes(shape = pt_shape, fill = pt_fill, colour = col),
               size = 2.7, stroke = 0.9) +
    # off-scale carets: small up-triangle at the ceiling / down-triangle at the
    # floor where a CI is clipped
    geom_point(data = caret, aes(stage, y_caret, colour = col), inherit.aes = FALSE,
               shape = 17, size = 1.5) +
    geom_point(data = caret_dn, aes(stage, y_caret, colour = col), inherit.aes = FALSE,
               shape = 25, size = 1.5) +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_shape_identity() +
    scale_y_continuous("median per-embryo  log2(L / R),  total PITX2 RNA",
                       breaks = seq(YMIN, CAP, 1),
                       labels = function(b) ifelse(b == CAP, "5⁺", as.character(b))) +
    coord_cartesian(ylim = c(YMIN, CAP), clip = "on") +
    # by_layer_x = TRUE indexes the strip elements by NESTING DEPTH (outer group
    # strip vs inner population strip), not per-strip — robust to the group strip
    # repeating across wrapped rows (Inference pool spans rows 1-2).
    facet_nested_wrap(vars(group, population), ncol = 3, nest_line = element_line(),
                      strip = strip_nested(
                          background_x = elem_list_rect(fill = c("grey88", "white")),
                          text_x = elem_list_text(face = c("bold", "plain")),
                          by_layer_x = TRUE)) +
    labs(x = "developmental stage") +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1, colour = "grey20"),
        axis.text.y = element_text(size = 8, colour = "grey20"),
        axis.title = element_text(size = 9.5),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.5, "lines"))

ggsave(file.path(HERE, "panel_C_laterality_forest.png"), p,
       width = 9.4, height = 7.2, dpi = 200)
cat("wrote", file.path(HERE, "panel_C_laterality_forest.png"), "\n")
