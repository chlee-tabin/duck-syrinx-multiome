#!/usr/bin/env Rscript
# Fig 02 Panel D: PITX2 locus pileup — developmental isoform trajectory
#
# Multi-stage pileup tracks at the PITX2 locus.  Reads are strand-
# filtered to PITX2-sense by d2i (+ strand, samtools -F 16) so the
# annotated antisense ncRNAs LOC113843538 (PITX2 intron 1) and
# LOC113843539 (PITX2 3'UTR overlap) no longer contribute to the
# pileup — see pipeline/_pitx2_strand.py and
# docs/pitx2_locus_annotation.md §1, §6.
#
# Y-axis: two-tier (dorsal tracks natural scale; non-dorsal tracks at
# dorsal_max / NON_DORSAL_AMPLIFICATION, currently 3.5x).  Post strand
# filter the dorsal tracks still carry real pan-PITX2 signal at shared
# exons / 3'UTR; a single shared axis would compress the non-dorsal
# tracks into illegibility.  See L195+ for the amplification rationale.
#
# Gene model with 3 isoform transcript structures at bottom.
# Style: B&W, minimal color.

.libPaths(c("~/R/4.4.2", .libPaths()))

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Paths ────────────────────────────────────────────────────────────────
data_dir <- "figures/data/d2i_pileup_output"
out_dir  <- "figures/fig02"
depth_file <- file.path(data_dir, "d2i_pitx2_rna_pileup.tsv.gz")

# Shared PITX2 gene-model constants + draw_pitx2_gene_model() helper.
# Moved into figures/_shared/ so the exploratory strand-mirror plot
# (not distributed) stays in lock-step with
# this published panel when PITX2 annotation coordinates are corrected.
source("figures/_shared/pitx2_gene_model.R")
source("figures/_shared/subcluster_map.R")   # display_population

cat("Loading depth data...\n")
df <- read.delim(depth_file, check.names = FALSE)
cat(sprintf("  %s positions, %d columns\n", format(nrow(df), big.mark = ","), ncol(df)))

# Alias the shared constants to the short names this script previously
# used (keeps the rest of the file readable without a find/replace).
GENE_START <- PITX2_GENE_START
GENE_END   <- PITX2_GENE_END
isoforms   <- PITX2_ISOFORMS
CDS_END    <- PITX2_CDS_END
UTR3_START <- PITX2_UTR3_START
UTR3_END   <- PITX2_UTR3_END

# ── Track selection ─────────────────────────────────────────────────────
# 9 tracks across 4 stages, organized to show developmental trajectory.
# pop_suffix names map to columns in d2i_pitx2_rna_pileup.tsv.gz; barcode
# sets are produced by figures/data/d2i_barcodes/generate_early_stage_barcodes.py
# from the current SoT (references/d2g_subcluster_annotations.tsv).
#
# 2026-05-04 SoT update:
#   - "splanchnic_HH1716" pop_suffix retained for pileup-TSV column
#     compatibility (rendered above as "Ventral LPM (pooled)"); the column
#     now points to the per-stage HH17/16 LPM inference pool (Early LPM +
#     Cardiopharyngeal LPM + BARX1+), NOT the ex-"Splanchnic" cluster
#     (which was relabeled "Venous-pole spl." and excluded from inference
#     as off-axis introgression).
#   - "splanchnic_HH13" exists as a pileup-TSV column too but is NOT
#     rendered by THIS panel D (HH13 row uses "all_lpm_HH13"); only the
#     figS02 panel D consumer references it.  Same legacy-name caveat
#     applies if anyone wires it up later.
#   - "emerging_dorsal_HH1716" track barcodes now resolve via the
#     "Dorsal-biased LPM" label (the post-override label for the NKX6-1+
#     pre-commitment cluster at HH17/16); this fixed a 0-cell silent
#     regression on 2026-05-04.
# See docs/lpm_mesenchyme_subcluster_biology.md for biology rationale.
# `prevet` = this track's barcode set predates the 2026-07-02 subcluster vet, so the
# current SSOT name (and its publication hedge) does NOT describe these cells.
# 2026-08-24: NO track is pre-vet anymore — the d2i export was regenerated from
# barcode sets rebuilt off the current SSOT (post-vet, post-2026-08-09 seeded d2g;
# ventral_splanchnic_HH25 15,194 / syringeal_HH25 2,584 / dorsal_subepi_HH25 2,645 /
# somatic_HH25 715 cells; the v210-staged GEX BAMs replaced the purged count BAMs).
# The prevet column is retained as machinery should a stale export ever return.
# History: the pre-regen HH25 tracks (ventral 14,273 / syringeal 2,068 =
# old cl5) were exempted from display_population() to avoid asserting the
# IRX1+/PTCHD4+ identity for cells that did not carry it.
tracks <- tribble(
    ~pop_suffix,                ~stage,     ~compartment,               ~track_order, ~prevet,
    "all_lpm_HH13",             "HH13",     "All LPM mesenchyme",       1,            FALSE,
    "emerging_dorsal_HH1716",   "HH17/16",  "Dorsal-biased LPM",        2,            FALSE,
    "splanchnic_HH1716",        "HH17/16",  "Ventral LPM (pooled)",     3,            FALSE,
    "dorsal_subepi_HH21",       "HH21",     "Dorsal sub-epithelial",    4,            FALSE,
    "ventral_only_HH21",        "HH21",     "Ventral splanchnic",       5,            FALSE,
    "syringeal_HH21",           "HH21",     "Syringeal",                6,            FALSE,
    "dorsal_subepi_HH25",       "HH25",     "Dorsal sub-epithelial",    7,            FALSE,
    "ventral_splanchnic_HH25",  "HH25",     "Ventral splanchnic",       8,            FALSE,
    "syringeal_HH25",           "HH25",     "Syringeal",                9,            FALSE,
)

# Rendered compartment label: the publication hedge applies to every track
# whose cells are current; pre-vet tracks keep their original name (see `prevet` above).
tracks$render_label <- ifelse(
    tracks$prevet,
    tracks$compartment,
    display_population(tracks$compartment, short = TRUE)
)
stopifnot(!anyDuplicated(paste(tracks$stage, tracks$render_label)))

# Verify all populations exist in data
pc_cols <- grep("_percell$", names(df), value = TRUE)
available_pops <- gsub("_percell$", "", pc_cols)
for (i in seq_len(nrow(tracks))) {
    if (!tracks$pop_suffix[i] %in% available_pops) {
        stop(sprintf("Population '%s' not found in data. Available: %s",
                     tracks$pop_suffix[i], paste(available_pops, collapse = ", ")))
    }
}
cat(sprintf("All %d populations found in data.\n", nrow(tracks)))

# ── Reshape depth data ──────────────────────────────────────────────────
use_cols <- paste0(tracks$pop_suffix, "_percell")

df_long <- df %>%
    select(position, all_of(use_cols)) %>%
    pivot_longer(cols = all_of(use_cols),
                 names_to = "pop_col", values_to = "depth") %>%
    mutate(pop_suffix = gsub("_percell$", "", pop_col)) %>%
    left_join(tracks, by = "pop_suffix") %>%
    mutate(pos_kb = position / 1e3)

# Smoothing window (bp) — 0 = no smoothing, 100 or 200 for smoothed
SMOOTH_WIN <- 50
cat(sprintf("Smoothing window: %d bp\n", SMOOTH_WIN))

if (SMOOTH_WIN > 0) {
    df_long <- df_long %>%
        group_by(pop_suffix) %>%
        mutate(depth_smooth = stats::filter(depth, rep(1/SMOOTH_WIN, SMOOTH_WIN), sides = 2) |>
                   as.numeric()) %>%
        ungroup()
} else {
    df_long$depth_smooth <- df_long$depth
}

# Two normalization modes:
# - "upm": within-track UMIs per million (shape comparison, hides absolute levels)
# - "percell": raw per-cell depth (shows absolute expression differences)
NORM_MODE <- "percell"   # change to "upm" for shape comparison

df_long <- df_long %>%
    group_by(pop_suffix) %>%
    mutate(
        track_total = sum(depth_smooth, na.rm = TRUE),
        depth_upm = ifelse(track_total > 0,
                           depth_smooth / track_total * 1e6, 0),
        # Per-cell depth (already normalized by cell count in the data)
        depth_percell = depth_smooth
    ) %>%
    ungroup()

# Select y-axis value based on normalization mode
df_long$depth_plot <- if (NORM_MODE == "upm") df_long$depth_upm else df_long$depth_percell
cat(sprintf("Normalization mode: %s\n", NORM_MODE))

# Track labels: "Stage\nCompartment"
df_long <- df_long %>%
    mutate(
        # Per-track hedge: `render_label` is the hedged name for current tracks and the
        # original name for pre-vet ones (see the `tracks` tribble above).
        track_label = paste0(stage, " ", render_label),   # one line (PNAS font pass)
        track_label = factor(track_label,
                             levels = paste0(tracks$stage, " ", tracks$render_label))
    )

cat("Per-population totals (raw UMIs per cell):\n")
df_long %>%
    group_by(pop_suffix, stage, compartment) %>%
    summarise(total_percell = sum(depth, na.rm = TRUE), .groups = "drop") %>%
    arrange(match(pop_suffix, tracks$pop_suffix)) %>%
    print()

# ── Zoom to gene region ─────────────────────────────────────────────────
PLOT_START <- GENE_START - 1500
PLOT_END   <- GENE_END + 1500
df_zoom <- df_long %>% filter(position >= PLOT_START, position <= PLOT_END)

# ── Highlight bands for isoform-specific first exon regions ─────────────
highlight_bands <- data.frame(
    xmin = c(3961387, 3964253, 3971385) / 1e3,
    xmax = c(3961458, 3964353, 3972064) / 1e3,
    label = c("Exon 1a", "Exon 1b", "Exon 1c"),
    stringsAsFactors = FALSE
)

# ── PolyA internal priming markers ──────────────────────────────────────
polya_markers <- data.frame(
    pos = c(3960307, 3961576, 3963604, 3965002, 3965244,
            3969889, 3970175, 3973664, 3974401) / 1e3,
    stringsAsFactors = FALSE
)

# ── Stage separator: horizontal lines between stage groups ──────────────
# We'll use strip text styling + spacing to separate stages visually

# ── Coverage tracks ─────────────────────────────────────────────────────
# Fill keyed by the `compartment` column of the tracks tibble above.
# `fill_vals` keys must match the `compartment` column of the `tracks`
# tribble exactly — `stopifnot()` below enforces.
fill_vals <- c(
    "All LPM mesenchyme"    = "gray30",
    "Dorsal-biased LPM"     = "black",
    "Ventral LPM (pooled)"  = "gray40",
    "Dorsal sub-epithelial" = "black",
    "Ventral splanchnic"    = "gray50",
    "Syringeal"             = "gray50"
)
# Guard against future compartment renames silently falling through to
# NA fill — scale_fill_manual with unmapped values drops the geom.
stopifnot(all(unique(tracks$compartment) %in% names(fill_vals)))

# ── Two-tier y-axis: dorsal at natural scale, non-dorsal amplified ──────
# Post strand filter the dorsal tracks still dominate (real pan-PITX2
# signal at shared exons / 3'UTR — NOT antisense artifact; see
# docs/pitx2_locus_annotation.md §6).  Putting all 9 tracks on a
# single shared axis compresses the ventral / syringeal tracks into
# illegibility.  Two-tier scheme: dorsal tracks at their natural max,
# every other track at dorsal_max / NON_DORSAL_AMPLIFICATION so both the
# pan-PITX2 dorsal peaks AND the lateralized PITX2c ventral / syringeal
# peaks are readable side-by-side.
#
# Amplification factor 3.5 (set 2026-05-04):
# the tallest non-dorsal peak (HH21 Ventral spl., 0.0427 — smooth50 of
# the rendered ventral_only_HH21 track per the tracks tribble above)
# lands at ~83% of the amplified ceiling, the Syringeal tracks (peak
# 0.0123) at ~24%.  At
# the prior 2x factor the tallest non-dorsal peak only filled 48% of its
# facet, leaving the panel reading as "dorsal massive, non-dorsal flat
# noise" rather than "dorsal massive, non-dorsal lateralized peaks
# present at the 3'UTR end".  4x clips the HH21 Ventral spl. peak; 3x
# leaves more whitespace than necessary.
#
# Note the rename: "Emerging dorsal" -> "Dorsal-biased LPM" per the
# 2026-04-20 compendium redesign (see docstring on the tracks tibble).
#
# IMPORTANT: keep figures/fig02/LEGEND.md "3.5x amplified" synced with
# this constant — the legend value is hand-typed there.
NON_DORSAL_AMPLIFICATION <- 3.5
dorsal_pops <- c("Dorsal sub-epithelial", "Dorsal-biased LPM")
dorsal_max <- max(df_zoom$depth_plot[df_zoom$compartment %in% dorsal_pops],
                  na.rm = TRUE)
# Fail-fast on missing dorsal data (would otherwise produce -Inf y-ceiling
# and silently broken panels).
stopifnot(is.finite(dorsal_max) && dorsal_max > 0)
amplified_max <- dorsal_max / NON_DORSAL_AMPLIFICATION
cat(sprintf("Y-axis tiers: dorsal max = %.4f, amplified max = %.4f (%.1fx)\n",
            dorsal_max, amplified_max, NON_DORSAL_AMPLIFICATION))

y_limits <- df_zoom %>%
    select(track_label, compartment) %>%
    distinct() %>%
    mutate(
        y_max = ifelse(compartment %in% dorsal_pops, dorsal_max, amplified_max),
        pos_kb = PLOT_START / 1e3
    )

p_tracks <- ggplot(df_zoom, aes(x = pos_kb, y = depth_plot)) +
    # Highlight bands for isoform-specific first exons
    geom_rect(data = highlight_bands,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
              inherit.aes = FALSE, fill = "gray92", color = NA) +
    # PolyA markers
    geom_point(data = polya_markers,
               aes(x = pos, y = Inf),
               inherit.aes = FALSE, shape = 25, size = 0.6,
               fill = "gray60", color = "gray60") +
    # Coverage area
    geom_area(aes(fill = compartment), alpha = 0.85, linewidth = 0) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    # Invisible blank anchors per-facet y-limit: dorsal tracks get
    # dorsal_max, non-dorsal tracks get amplified_max.  Requires
    # scales="free_y" so each facet respects its own y_max.
    geom_blank(data = y_limits, aes(x = pos_kb, y = y_max)) +
    facet_grid(track_label ~ ., scales = "free_y", switch = "y") +
    scale_x_continuous(
        limits = c(PLOT_START / 1e3, PLOT_END / 1e3),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),
        # PNAS final layout: short tracks -> one tick per track (the largest
        # pretty break inside the track; baseline = 0), so labels don't collide
        breaks = function(l) {
            pb <- scales::pretty_breaks(n = 2)(l)
            pb <- pb[pb <= l[2] & pb >= 0]
            if (length(pb) < 2) stop("pileup track y-range yielded < 2 breaks: ", paste(l, collapse = ","))
            max(pb)   # one labelled tick per track; the track baseline is 0
        }
    ) +
    labs(y = if (NORM_MODE == "upm") {
        bquote("UMI fraction (" %*% 10^{-6} * "," ~ .(SMOOTH_WIN) ~ "bp smooth)")
    } else if (SMOOTH_WIN > 0) {
        bquote("UMIs per cell (" * .(SMOOTH_WIN) ~ "bp smooth)")
    } else {
        "UMIs per cell (raw)"
    }) +
    # PNAS final layout (2026-10): Fig 4 places this panel at native size and
    # prints it at x0.667 -> 9-9.5 pt text prints at 6.0-6.3 pt
    theme_classic(base_size = 9.5) +
    theme(
        # Strip text bumped 6.5 → 8 since W_FULL gives much more horizontal
        # room for left-side track labels.  Y-axis title/tick text bumped
        # in tandem so the legend doesn't appear under-weighted.
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9.5,
                                         lineheight = 1.0),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        panel.spacing = unit(0.15, "lines"),
        plot.margin = margin(12, 10, 0, 16)   # room for the composite letter "E."
    )

# ── Gene model ──────────────────────────────────────────────────────────
# All exon / UTR / CDS / TSS-arrow drawing now lives in
# figures/_shared/pitx2_gene_model.R so the exploratory strand-mirror
# plot and this published panel share one source of truth.
p_gene <- ggplot() +
    geom_rect(data = highlight_bands,
              aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = 4.5),
              fill = "gray92", color = NA)
p_gene <- draw_pitx2_gene_model(p_gene, label_size = 3.4)

p_gene <- p_gene +
    annotate("text", x = (UTR3_START + UTR3_END) / 2 / 1e3,
             y = 0.3, label = "3\u2032UTR", size = 3.2, color = "gray30") +
    annotate("text", x = GENE_END / 1e3 + 0.5, y = 2,
             label = "\u2192", size = 5, color = "black") +
    geom_point(data = polya_markers,
               aes(x = pos, y = 4.2), shape = 25, size = 0.6,
               fill = "gray60", color = "gray60") +
    scale_x_continuous(
        name = "NC_051775.1 position (kb)",
        limits = c(PLOT_START / 1e3, PLOT_END / 1e3),
        expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(-0.5, 4.8)) +
    theme_void(base_size = 9.5) +
    theme(
        axis.title.x = element_text(size = 10, margin = margin(t = 2)),
        axis.text.x = element_text(size = 9.5),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.length.x = unit(1.5, "mm"),
        axis.line.x = element_line(linewidth = 0.3),
        plot.margin = margin(0, 10, 5, 16)
    )

# ── Combine ──────────────────────────────────────────────────────────────
combined <- p_tracks / p_gene +
    plot_layout(heights = c(10, 3.2))

suffix <- sprintf("_%s_smooth%d", NORM_MODE, SMOOTH_WIN)
out_path <- file.path(out_dir, paste0("panel_D_pileup_multistage", suffix, ".png"))
# Full-width panel D (W_FULL = 10.5"). Height grows only modestly (4.4 → 5.0)
# because track count is fixed at 9 and the structural change is purely
# horizontal — the genomic axis (PLOT_START..PLOT_END = ~18kb) is now
# rendered across ~9.5" of plot area instead of ~1.5", revealing the
# isoform structure and intronic polyA markers legibly.
# PNAS final layout (2026-10): authored at its Fig 4 slot (756 x 222 pt).
ggsave(out_path, combined, width = 10.5, height = 222 / 72, dpi = 300, bg = "white")
cat("Saved:", out_path, "\n")
