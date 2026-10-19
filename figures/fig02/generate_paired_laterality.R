#!/usr/bin/env Rscript
# Figure 2 panel C: Paired laterality dot plot — ALL stages
#
# For each genotype whose two libraries hold its two opposite halves (design
# L/R library pair; same-side duplicates excluded -- _shared/paired_libraries.R):
#   PITX2 metric in Library A vs Library B, faceted by spatial domain.
#   Library A/B is an arbitrary consistent assignment (lower-numbered = A).
#   Colored by developmental stage (from _shared/style.py STAGE_COLORS).
#   Stages HH16 vs HH17/16 vs HH17 split via e1b inferred_stage.
#
# Produces TWO outputs:
#   panel_C_paired_laterality.png         — PITX2+ fraction
#   panel_C_paired_laterality_meanexpr.png — mean PITX2 CP10K

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)

# ── Paths ──
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE))
script_dir <- if (length(script_path) > 0) dirname(script_path) else "figures/fig02"
data_dir <- file.path(dirname(script_dir), "data")
shared_dir <- file.path(dirname(script_dir), "_shared")

mapping_path <- file.path(data_dir, "b_genotype_embryo_mapping.tsv")
stage_path   <- file.path(data_dir, "e1b_stage_assignments.tsv")

# ── Shared subcluster map (single source of truth from d2g validator) ──
source(file.path(shared_dir, "subcluster_map.R"))

# ── Tech-rep merge map (single source of truth: _shared/library_naming.R) ──
source(file.path(shared_dir, "library_naming.R"))
tech_rep_merge <- TECH_REP

# ── Half-vs-half pairs: design opposite-side library pairs, minus same-side
# duplicates (lane-loading-record group 17-B(2) = embryos 17-5/17-6: one left half split across GEX04 and
# GEX0S-duck).  Rule + rationale: _shared/paired_libraries.R.  Library A =
# lexically lower library of the pair.
source(file.path(shared_dir, "paired_libraries.R"))
repo_root <- dirname(dirname(script_dir))
design_pairs <- design_library_pairs(file.path(repo_root, "design", "embryo_lr_design.tsv"))
same_side_dups <- d2l_same_side_duplicates(file.path(data_dir, "d2l_lr_inference.tsv"))
check_duplicates_against_design(same_side_dups, design_pairs)

# ── Stage colors (from _shared/style.py STAGE_COLORS) ──
stage_colors <- c(
  "HH13"    = "#66a61e",
  "HH16"    = "#e7298a",
  "HH17/16" = "#b0a8d0",
  "HH17"    = "#7570b3",
  "HH21"    = "#d95f02",
  "HH25"    = "#1b9e77"
)
stage_levels <- c("HH13", "HH16", "HH17/16", "HH17", "HH21", "HH25")

# ── Per-genotype stage from e1b ──
e1b <- read_tsv(stage_path, show_col_types = FALSE) %>%
  select(genotype_id, inferred_stage)
e1b_map <- setNames(e1b$inferred_stage, e1b$genotype_id)
for (uid in grep("^G3_unfiltered_donor", e1b$genotype_id, value = TRUE)) {
  e1b_map[sub("G3_unfiltered_donor", "G3_donor", uid)] <- e1b_map[uid]
}

# ── Display names ──
mapping_names <- read_tsv(mapping_path, show_col_types = FALSE) %>%
  select(genotype_id, display_name)
name_map <- setNames(mapping_names$display_name, mapping_names$genotype_id)
for (id in grep("^G3_donor", names(name_map), value = TRUE)) {
  name_map[sub("G3_donor", "G3_unfiltered_donor", id)] <- name_map[id]
}
name_map["dp500_donor29"] <- name_map["dp500_d29"]
# Genotypes absent from b_genotype_embryo_mapping.tsv (known gap: G2_genotype_01,
# internal audit record, not distributed) take their display name
# from the hand-maintained SSOT references/genotype_metadata.tsv (style.py's
# DISPLAY_NAMES source) instead of leaking the raw pipeline id onto the figure.
gm <- read_tsv(file.path(repo_root, "references", "genotype_metadata.tsv"),
               show_col_types = FALSE)
for (i in seq_len(nrow(gm))) {
  gid <- gm$genotype_id[i]
  if (!gid %in% names(name_map)) name_map[gid] <- gm$display_name[i]
}

# ── Stage configurations (clusters resolved from POOLS in subcluster_map) ──
# Pool composition is the SSOT — `figures/_shared/subcluster_map.R` `.POOLS`
# defines per-stage inference + bilateral + excluded membership; we delegate
# to `get_pool_map(stage, kind)` rather than enumerate cluster IDs here.
# Why panel C aggregates at the full inference-pool level (not just
# Ventral spl. + adjacent clusters): an n>=15 cell filter on a stricter subset
# cherry-picked the strongest paired genotypes and produced artificially clean separation.
#
# Three display facets:
#   LPM mes. (HH13)                 — HH13 inference pool (D/V unresolved pre-TES)
#   Ventral splanchnic pool (HH17+) — HH17+ inference pool, post-2026-05-04 SoT:
#                                       Venous-pole spl. excluded as off-axis;
#                                       see docs/lpm_mesenchyme_subcluster_biology.md
#   Dorsal sub-epi.                 — bilateral control (Dorsal-biased LPM at
#                                     HH17/16; Dorsal sub-epi. at HH21/HH25)
# Skeletogenic LPM / Mesothelium / Venous-pole spl. / paraxial / Contamination
# excluded via being outside the inference + bilateral pools.
fold_for_laterality <- function(stage) {
  inf_map <- get_pool_map(stage, "inference")
  bil_map <- get_pool_map(stage, "bilateral")

  out <- list()
  for (cid in names(inf_map)) {
    if (stage == "HH13") {
      out[[cid]] <- "LPM mes. (HH13)"
    } else {
      out[[cid]] <- "Ventral splanchnic pool (HH17+)"
    }
  }
  for (cid in names(bil_map)) {
    out[[cid]] <- "Dorsal sub-epi."
  }
  unlist(out)
}

stage_tsv <- c(
  "HH13"    = "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
  "HH17/16" = "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
  "HH21"    = "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
  "HH25"    = "d2g_lpm_subclustering_HH25_ccreg.tsv.gz"
)

stage_configs <- lapply(names(stage_tsv), function(st) {
  list(file = stage_tsv[[st]],
       leiden_col = get_leiden_col(st),
       domain_map = fold_for_laterality(st),
       fallback_stage = st)
})

# ── Process each stage ──
# Display order = developmental + spatial: HH13 LPM (origin) → HH17+ ventral
# splanchnic pool (where laterality lives, full figS02D inference set) →
# dorsal sub-epi. (bilateral control).
target_domains <- c("LPM mes. (HH13)",
                    "Ventral splanchnic pool (HH17+)",
                    "Dorsal sub-epi.")

process_stage <- function(config) {
  fpath <- file.path(data_dir, config$file)
  if (!file.exists(fpath)) stop("missing stage table: ", fpath)  # fail loud, never drop a stage

  df <- read_tsv(fpath, show_col_types = FALSE)
  lcol <- config$leiden_col
  message(sprintf("  %s: %s cells", config$fallback_stage,
                  format(nrow(df), big.mark = ",")))

  df <- df %>%
    select(genotype, library, cluster = all_of(lcol), PITX2_cp10k) %>%
    mutate(cluster = as.character(cluster),
           domain = config$domain_map[cluster]) %>%
    filter(!is.na(domain)) %>%
    mutate(pitx2_pos = PITX2_cp10k > 0,
           lib_merged = ifelse(library %in% names(tech_rep_merge),
                               tech_rep_merge[library], library))

  # Half-vs-half genotypes: design opposite-side library pairs; same-side
  # duplicates excluded (logged); anything unresolvable raises
  geno_libs <- df %>%
    group_by(genotype) %>%
    summarise(libs = list(unique(lib_merged)), .groups = "drop")

  paired_df <- resolve_half_pairs(geno_libs, design_pairs, same_side_dups,
                                  config$fallback_stage)
  if (is.null(paired_df) || nrow(paired_df) == 0) return(NULL)

  message(sprintf("    %d paired genotypes", nrow(paired_df)))

  dat <- df %>%
    inner_join(paired_df, by = "genotype") %>%
    mutate(lib_member = case_when(
      lib_merged == lib_A ~ "Library A",
      lib_merged == lib_B ~ "Library B",
      TRUE ~ NA_character_
    ))
  if (anyNA(dat$lib_member))  # impossible: resolve_half_pairs takes exactly two libraries
    stop("cells outside their genotype's library pair: ",
         paste(unique(dat$genotype[is.na(dat$lib_member)]), collapse = ", "))

  no_stage <- setdiff(unique(dat$genotype), names(e1b_map))
  if (length(no_stage) > 0)
    stop("No e1b stage for genotype(s): ", paste(no_stage, collapse = ", "))
  fracs <- dat %>%
    group_by(genotype, domain, lib_member) %>%
    summarise(frac_expressing = mean(pitx2_pos),
              mean_cp10k = mean(PITX2_cp10k),
              n_cells = n(), .groups = "drop") %>%
    mutate(stage = unname(e1b_map[genotype]),
           display = unname(name_map[genotype]))
  no_stage <- unique(fracs$genotype[is.na(fracs$stage)])
  if (length(no_stage) > 0)
    stop("No e1b stage for genotype(s): ", paste(no_stage, collapse = ", "),
         " (", config$fallback_stage, ") — check e1b_stage_assignments.tsv")
  missing <- unique(fracs$genotype[is.na(fracs$display)])
  if (length(missing) > 0)
    stop("No display name for genotype(s): ", paste(missing, collapse = ", "),
         " — add to references/genotype_metadata.tsv")
  return(fracs)
}

message("Processing stages:")
all_fracs <- lapply(stage_configs, process_stage) %>% bind_rows()
message(sprintf("\nStages present: %s", paste(sort(unique(all_fracs$stage)), collapse = ", ")))

# ── Pivot helper ──
pivot_and_filter <- function(fracs_long, value_col) {
  wide <- fracs_long %>%
    select(genotype, display, domain, stage, lib_member,
           value = all_of(value_col), n_cells) %>%
    pivot_wider(id_cols = c(genotype, display, domain, stage),
                names_from = lib_member,
                values_from = c(value, n_cells)) %>%
    rename(val_A = `value_Library A`, val_B = `value_Library B`,
           n_A = `n_cells_Library A`, n_B = `n_cells_Library B`) %>%
    mutate(min_cells = pmin(n_A, n_B)) %>%
    filter(!is.na(val_A), !is.na(val_B))

  low <- wide$min_cells < 15
  if (any(low)) message(sprintf("  Dropping %d pairs with < 15 cells", sum(low)))
  wide <- wide[!low, ]

  ps <- intersect(stage_levels, unique(wide$stage))
  wide %>% mutate(domain = factor(domain, levels = target_domains),
                  stage = factor(stage, levels = ps))
}

wide_frac <- pivot_and_filter(all_fracs, "frac_expressing")
wide_mean <- pivot_and_filter(all_fracs, "mean_cp10k")

# ── Summary ──
print_summary <- function(w, label) {
  message(sprintf("\n--- %s ---", label))
  w %>% mutate(d = abs(val_A - val_B)) %>%
    group_by(stage, domain) %>%
    summarise(n = n(), md = round(mean(d), 3), .groups = "drop") %>%
    {for (i in seq_len(nrow(.))) {
      r <- .[i, ]; message(sprintf("  %s × %s (n=%d): mean |diff| = %.3f",
                                    r$stage, r$domain, r$n, r$md))
    }}
}
print_summary(wide_frac, "PITX2+ fraction")
print_summary(wide_mean, "Mean CP10K")

# ── Plotting function ──
# Both sub-variants of Panel C are built with show_legend=FALSE here; the
# legend is extracted once via cowplot::get_legend and composed to the
# right of the two stacked plots below.  cowplot::plot_grid(align = "v",
# axis = "lr") forces the left/right plot-panel edges of p1 and p2 to
# match so the three facet columns (LPM mes. (HH13) / Ventral splanchnic
# pool (HH17+) / Dorsal sub-epi.) line up vertically between the two
# metric rows.
#
# Per-facet axes (2026-04-21 layout polish): switched to
# `scales = "free"` so each domain rescales independently.  The previous
# shared y-scale squished the LPM mes. (HH13) and pre-redesign Ventral +
# Syringeal facets because Dorsal sub-epi. (HH21/HH25 PITX2a/b dominant)
# fixed the top of every facet at ~80%.  Per-facet axis breaks are picked
# dynamically by `scales::breaks_extended(n = 4)` since the per-facet
# max varies an order of magnitude (HH21/HH25 dorsal ≈ 80%, HH17/16
# splanchnic ≈ 30%).
#
# `coord_fixed(ratio = 1)` is incompatible with `scales = "free"`; instead
# we pre-compute symmetric per-facet limits (max of val_A, val_B per domain)
# and emit invisible `geom_blank()` data so x and y axes share the same
# range within each facet — this keeps the y=x bilateral-reference diagonal
# corner-to-corner inside every panel even though the panel rectangle is
# now sized by the layout (we set `aspect.ratio = 1` in the theme so each
# facet's plot panel renders square anyway).
# ── 4D labels: trough embryos whose two halves agree in the ventral pool ──
# Set for the PNAS final figure (2026-10-12). Per embryo, divergence = max over the
# two metrics of |B - A| / facet axis range. The threshold is data-derived:
# the SMALLEST divergence of any wave-stage embryo (HH13 / HH21 / HH25) in its
# inference facet. A trough embryo (HH16 / HH17/16 / HH17) is labelled when its
# ventral-pool divergence is below that, i.e. it is closer to the diagonal than
# every lateralized-stage embryo -- its RNA laterality signal is lost. The same
# embryos are labelled in the dorsal control facet. Fails loud if none qualify.
TROUGH_STAGES <- c("HH16", "HH17/16", "HH17")
WAVE_STAGES <- c("HH13", "HH21", "HH25")
INFERENCE_DOMAINS <- target_domains[1:2]
VENTRAL_DOMAIN <- target_domains[2]
LABEL_DOMAINS <- target_domains[2:3]
divergence <- function(w) {
  bind_rows(lapply(split(w, w$domain, drop = TRUE), function(v)
    data.frame(display = v$display, stage = as.character(v$stage),
               domain = as.character(v$domain),
               d = abs(v$val_B - v$val_A) / max(c(v$val_A, v$val_B)))))
}
div <- bind_rows(divergence(wide_frac), divergence(wide_mean)) %>%
  group_by(display, stage, domain) %>% summarise(d = max(d), .groups = "drop")
wave <- div[div$stage %in% WAVE_STAGES & div$domain %in% INFERENCE_DOMAINS, ]
if (nrow(wave) == 0) stop("4D label rule: no wave-stage embryos in the inference facets")
diag_threshold <- min(wave$d)
trough <- div[div$stage %in% TROUGH_STAGES & div$domain == VENTRAL_DOMAIN, ]
label_set <- trough$display[trough$d < diag_threshold]
if (length(label_set) == 0) stop("4D label rule selected no trough embryos")
message(sprintf("4D labels: trough divergence < %.3f (min wave-stage, %s): %s",
                diag_threshold, wave$display[which.min(wave$d)],
                paste(sort(label_set), collapse = ", ")))

# ── 4D label placement ──
# The labelled embryos sit on the diagonal inside the crowded HH17/16 cluster,
# where ggrepel alone parked labels touching a neighbouring ring with no
# visible leader (g37, dorsal PITX2+ fraction facet). Each label therefore
# starts a fixed distance (LABEL_NUDGE x facet axis range) from its point, in
# the compass direction with the most clearance: the larger of the minimum
# distances from (a) the label anchor and (b) the leader segment to every other
# point of the facet and to the labels already placed there (so a leader never
# runs over a neighbouring point), inside the facet. ggrepel then only resolves
# residual overlaps. Raises if no in-facet direction exists.
LABEL_NUDGE <- 0.16
LABEL_ANGLES <- seq(0, 7) * pi / 4
seg_dist <- function(q, a, b) {        # distance of points q (n x 2) to segment a-b
  d <- b - a
  t <- pmin(pmax(((q[, 1] - a[1]) * d[1] + (q[, 2] - a[2]) * d[2]) / sum(d^2), 0), 1)
  sqrt((a[1] + t * d[1] - q[, 1])^2 + (a[2] + t * d[2] - q[, 2])^2)
}
label_nudges <- function(w, facet_lim) {
  w$nudge_x <- 0
  w$nudge_y <- 0
  for (dom in LABEL_DOMAINS) {
    in_dom <- which(w$domain == dom)
    idx <- in_dom[w$display[in_dom] %in% label_set]
    if (length(idx) == 0) next
    lim <- facet_lim$lim[facet_lim$domain == dom]
    obstacles <- cbind(w$val_A[in_dom], w$val_B[in_dom]) / lim
    for (i in idx) {
      p <- c(w$val_A[i], w$val_B[i]) / lim
      others <- obstacles[-match(i, in_dom), , drop = FALSE]   # not the point itself
      best <- NULL
      best_clear <- -Inf
      for (a in LABEL_ANGLES) {
        anchor <- p + LABEL_NUDGE * c(cos(a), sin(a))
        if (any(anchor < 0.06) || any(anchor > 0.94)) next
        clear <- min(seg_dist(others, p, anchor))   # covers the anchor end too
        if (clear > best_clear) {
          best_clear <- clear
          best <- anchor
        }
      }
      if (is.null(best))
        stop("4D label placement: no in-facet direction for ", w$display[i], " in ", dom)
      obstacles <- rbind(obstacles, best)
      w$nudge_x[i] <- (best[1] - p[1]) * lim
      w$nudge_y[i] <- (best[2] - p[2]) * lim
    }
  }
  w
}

make_plot <- function(wide_df, xlab, ylab, pct_fmt = FALSE,
                      show_legend = TRUE) {
  if (pct_fmt) {
    fmt <- scales::percent_format(accuracy = 1)
  } else {
    fmt <- waiver()
  }

  # Per-domain symmetric limits (invisible geom_blank to lock axis range).
  # Pad by 8% for label room — geom_text_repel pushes labels outside data.
  per_facet_lim <- wide_df %>%
    group_by(domain) %>%
    summarise(lim = max(c(val_A, val_B), na.rm = TRUE) * 1.12,
              .groups = "drop")
  blanks <- bind_rows(
    per_facet_lim %>% mutate(val_A = lim, val_B = 0),
    per_facet_lim %>% mutate(val_A = 0, val_B = lim)
  )
  wide_df <- label_nudges(wide_df, per_facet_lim)

  p <- ggplot(wide_df, aes(x = val_A, y = val_B)) +
    facet_wrap(~ domain, nrow = 1, scales = "free",
               labeller = label_wrap_gen(width = 20)) +   # 2-line strips at print size
    geom_blank(data = blanks, aes(x = val_A, y = val_B), inherit.aes = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey75", linewidth = 0.3) +
    # HH17/16 emphasized (2026-06-24): red + thicker ring marks the
    # ambiguous/unresolved "silent trough" stage (RNA laterality can't be
    # called here -> chromatin resolves it in fig03; progressive disclosure).
    # Every other stage keeps the original thin black ring (stroke 0.3).
    geom_point(aes(fill = stage, size = min_cells,
                   color = stage == "HH17/16", stroke = stage == "HH17/16"),
               shape = 21, alpha = 0.8) +
    scale_discrete_manual(aesthetics = "stroke",
                          values = c(`TRUE` = 0.9, `FALSE` = 0.3),
                          guide = "none") +
    # Labels in dark gray for legibility regardless of stage color
    # (light HH17/16 on white was previously illegible).
    # PNAS final layout (2026-10): 3.25 mm (~9.2 pt native; ~6.2 pt printed in
    # Fig 4). max.overlaps = Inf so no embryo label is silently dropped; box
    # padding + always-drawn leaders keep labels from running together; the small
    # point.padding lets each leader end on its target ring.
    # Labels: on-diagonal trough embryos (label_set above), ventral + dorsal facets.
    # All points go to the repel layer (unselected labels blank) so the few
    # labels are pushed off every point, not only their own. Each label starts
    # at its label_nudges() offset (leaders drawn back to the data point).
    geom_text_repel(aes(label = ifelse(display %in% label_set &
                                         domain %in% LABEL_DOMAINS,
                                       display, "")), color = "#222222", size = 3.2,
                    segment.color = "grey25", segment.size = 0.35,
                    min.segment.length = 0, max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.15, force = 12,
                    force_pull = 0.08, max.time = 60, max.iter = 5000000,
                    nudge_x = wide_df$nudge_x, nudge_y = wide_df$nudge_y,
                    seed = 42, show.legend = FALSE) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4),
                       labels = fmt,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 4),
                       labels = fmt,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = stage_colors, name = "Stage", drop = TRUE) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"),
                       guide = "none") +
    scale_size_continuous(range = c(1.6, 3.8), breaks = c(30, 100, 300),
                          name = "Min. cells\n(per half)") +
    labs(x = xlab, y = ylab) +
    # PNAS final layout (2026-10): Fig 4 places this panel at native size and
    # prints it at x0.667 -> ticks/legend 9.5 pt = 6.3 pt, titles 10 pt = 6.7 pt
    theme_bw(base_size = 9.5) +
    theme(strip.text = element_text(size = 10, face = "bold",
                                    margin = margin(2, 2, 2, 2)),
          strip.background = element_rect(fill = "grey95"),
          axis.text = element_text(size = 9.5),
          axis.title = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
          legend.position = if (show_legend) "bottom" else "none",
          legend.box = "vertical",
          legend.box.just = "left",
          legend.spacing.y = unit(0, "cm"),
          legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 9.5),
          legend.title = element_text(size = 9.5),
          legend.spacing.x = unit(0.05, "cm"),
          legend.margin = margin(0, 2, 0, 2),
          panel.spacing.x = unit(0.45, "lines"),
          # PNAS final layout: facets fill the Fig 4 slot (not forced square);
          # the per-facet symmetric limits still put y = x corner to corner.
          plot.margin = margin(3, 5, 3, 3))

  if (show_legend) {
    p <- p + guides(fill = guide_legend(order = 1, nrow = 1,
                                        override.aes = list(size = 2.8)),
                    size = guide_legend(order = 2, nrow = 1))
  }
  p
}

# ── Build both sub-plots without legends, stack with a shared legend ──
# cowplot::plot_grid(..., align = "v", axis = "lr") forces the left and
# right plot-panel edges of p1 and p2 to align identically, then we paste
# the shared legend into a third column on the right.  This produces a
# single combined PNG that replaces the old two-PNG + stack_vertical
# approach and guarantees the facet dividers line up between the PITX2+
# fraction and Mean CP10K sub-panels.
p1_legend <- make_plot(wide_frac, "PITX2+ fraction (Library A)",
                       "PITX2+ fraction (Library B)", pct_fmt = TRUE,
                       show_legend = TRUE)
shared_legend <- cowplot::get_legend(p1_legend)

p1 <- make_plot(wide_frac, "PITX2+ fraction (Library A)",
                "PITX2+ fraction (Library B)", pct_fmt = TRUE,
                show_legend = FALSE)
p2 <- make_plot(wide_mean, "Mean PITX2 CP10K (Library A)",
                "Mean PITX2 CP10K (Library B)",
                show_legend = FALSE) +
  # PNAS final layout: the bottom row repeats the top row's facet columns,
  # so its strips are dropped to give the scatter more room
  theme(strip.text = element_blank(), strip.background = element_blank())

stacked <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")
combined <- cowplot::plot_grid(stacked, shared_legend,
                               ncol = 1, rel_heights = c(10, 1.1))

# Full-width panel C (W_FULL = 10.5"). Height stays compact at 5.4"
# (height grew from 6.2 → 5.4 even though width nearly doubled, because
# coord_fixed lets each facet be height-constrained: with three roughly
# square facets across ~9.4" of plotting area and 2 stacked metric rows,
# each facet panel ends up ~2.0–2.2" — readable but not bloated).
out1 <- file.path(script_dir, "panel_C_paired_laterality.png")
# PNAS final layout (2026-10): authored at its Fig 4 slot (D_W x D_H pt).
D_W <- 417; D_H <- 436
ggsave(out1, combined, width = D_W / 72, height = D_H / 72, dpi = 300, bg = "white")
message("\nSaved: ", out1)
