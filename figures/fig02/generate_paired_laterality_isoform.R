#!/usr/bin/env Rscript
# Figure 2 panel C companion: Isoform-specific paired laterality
#
# Joins d2j per-cell PITX2 isoform UMI counts with subclustering data
# to produce paired laterality plots for PITX2c and PITX2a/b separately.
#
# Expected result:
#   - PITX2c fraction: LATERALIZED in ventral/syringeal (off diagonal)
#   - PITX2a/b fraction: BILATERAL across all domains (on diagonal)
#   This directly demonstrates that laterality is PITX2c-specific.

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ── Paths ──
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE))
script_dir <- if (length(script_path) > 0) dirname(script_path) else "figures/fig02"
data_dir <- file.path(dirname(script_dir), "data")
shared_dir <- file.path(dirname(script_dir), "_shared")

mapping_path <- file.path(data_dir, "b_genotype_embryo_mapping.tsv")
stage_path   <- file.path(data_dir, "e1b_stage_assignments.tsv")
d2j_path     <- file.path(data_dir, "d2j_pitx2_isoform_percell.tsv.gz")

source(file.path(shared_dir, "subcluster_map.R"))
source(file.path(shared_dir, "library_naming.R"))

stopifnot(file.exists(d2j_path))

# ── Reuse infrastructure from generate_paired_laterality.R ──

tech_rep_merge <- TECH_REP

# Half-vs-half pairs: design opposite-side library pairs, minus same-side
# duplicates (one left half split across GEX04 and GEX0S-duck).  Same shared
# rule as generate_paired_laterality.R -- _shared/paired_libraries.R.
source(file.path(shared_dir, "paired_libraries.R"))
repo_root <- dirname(dirname(script_dir))
design_pairs <- design_library_pairs(file.path(repo_root, "design", "embryo_lr_design.tsv"))
same_side_dups <- d2l_same_side_duplicates(file.path(data_dir, "d2l_lr_inference.tsv"))
check_duplicates_against_design(same_side_dups, design_pairs)

stage_colors <- c(
  "HH13" = "#66a61e", "HH16" = "#e7298a", "HH17/16" = "#b0a8d0",
  "HH17" = "#7570b3", "HH21" = "#d95f02", "HH25" = "#1b9e77"
)
stage_levels <- c("HH13", "HH16", "HH17/16", "HH17", "HH21", "HH25")

# Stage + display name lookups
e1b <- read_tsv(stage_path, show_col_types = FALSE) %>% select(genotype_id, inferred_stage)
e1b_map <- setNames(e1b$inferred_stage, e1b$genotype_id)
for (uid in grep("^G3_unfiltered_donor", e1b$genotype_id, value = TRUE))
  e1b_map[sub("G3_unfiltered_donor", "G3_donor", uid)] <- e1b_map[uid]

mapping_names <- read_tsv(mapping_path, show_col_types = FALSE) %>% select(genotype_id, display_name)
name_map <- setNames(mapping_names$display_name, mapping_names$genotype_id)
for (id in grep("^G3_donor", names(name_map), value = TRUE))
  name_map[sub("G3_donor", "G3_unfiltered_donor", id)] <- name_map[id]
name_map["dp500_donor29"] <- name_map["dp500_d29"]
# Genotypes absent from b_genotype_embryo_mapping.tsv (G2_genotype_01) take their
# display name from the SSOT references/genotype_metadata.tsv, as in
# generate_paired_laterality.R (never the raw pipeline id).
gm <- read_tsv(file.path(repo_root, "references", "genotype_metadata.tsv"),
               show_col_types = FALSE)
for (i in seq_len(nrow(gm))) {
  if (!gm$genotype_id[i] %in% names(name_map)) name_map[gm$genotype_id[i]] <- gm$display_name[i]
}

# ── Load d2j isoform data ──
d2j <- read_tsv(d2j_path, show_col_types = FALSE)
message(sprintf("d2j: %s barcode×library entries", format(nrow(d2j), big.mark = ",")))

# ── Stage configs (clusters from POOLS, 2026-04-20 redesign) ──
# Mirror of fold_for_laterality in generate_paired_laterality.R.
fold_for_laterality <- function(stage) {
  inf_map <- get_pool_map(stage, "inference")
  bil_map <- get_pool_map(stage, "bilateral")
  out <- list()
  for (cid in names(inf_map)) {
    pop <- inf_map[[cid]]
    if (stage == "HH13") {
      out[[cid]] <- "LPM mes. (HH13)"
    } else if (pop == "Syringeal") {
      out[[cid]] <- "Syringeal"
    } else {
      out[[cid]] <- "Ventral + Syringeal"
    }
  }
  for (cid in names(bil_map)) out[[cid]] <- "Dorsal sub-epi."
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

target_domains <- c("Dorsal sub-epi.", "LPM mes. (HH13)",
                    "Ventral + Syringeal", "Syringeal")

# ── Process each stage: join subclustering with d2j ──
process_stage <- function(config) {
  fpath <- file.path(data_dir, config$file)
  if (!file.exists(fpath)) stop("missing stage table: ", fpath)  # fail loud, never drop a stage

  df <- read_tsv(fpath, show_col_types = FALSE)
  lcol <- config$leiden_col
  message(sprintf("  %s: %s cells", config$fallback_stage,
                  format(nrow(df), big.mark = ",")))

  df <- df %>%
    select(genotype, library, barcode, cluster = all_of(lcol)) %>%
    mutate(cluster = as.character(cluster),
           domain = config$domain_map[cluster]) %>%
    filter(!is.na(domain))

  # Join with d2j on (barcode, library) — cells without PITX2 reads get 0
  df <- df %>%
    left_join(d2j, by = c("barcode", "library")) %>%
    mutate(pitx2ab_umi = replace_na(pitx2ab_umi, 0L),
           pitx2c_umi  = replace_na(pitx2c_umi, 0L),
           shared_umi  = replace_na(shared_umi, 0L),
           has_pitx2c  = pitx2c_umi > 0,
           has_pitx2ab = pitx2ab_umi > 0)

  n_matched <- sum(df$pitx2c_umi > 0 | df$pitx2ab_umi > 0 | df$shared_umi > 0)
  message(sprintf("    d2j matched: %s / %s cells with any isoform UMI",
                  format(n_matched, big.mark = ","),
                  format(nrow(df), big.mark = ",")))

  # Merge tech reps for library pairing
  df <- df %>%
    mutate(lib_merged = ifelse(library %in% names(tech_rep_merge),
                               tech_rep_merge[library], library))

  # Half-vs-half genotypes (same-side duplicates excluded + logged; anything
  # unresolvable raises)
  geno_libs <- df %>%
    group_by(genotype) %>%
    summarise(libs = list(unique(lib_merged)), .groups = "drop")

  paired_df <- resolve_half_pairs(geno_libs, design_pairs, same_side_dups,
                                  config$fallback_stage)
  if (is.null(paired_df) || nrow(paired_df) == 0) return(NULL)

  dat <- df %>%
    inner_join(paired_df, by = "genotype") %>%
    mutate(lib_member = case_when(
      lib_merged == lib_A ~ "Library A",
      lib_merged == lib_B ~ "Library B",
      TRUE ~ NA_character_))
  if (anyNA(dat$lib_member))  # impossible: resolve_half_pairs takes exactly two libraries
    stop("cells outside their genotype's library pair: ",
         paste(unique(dat$genotype[is.na(dat$lib_member)]), collapse = ", "))
  no_stage <- setdiff(unique(dat$genotype), names(e1b_map))
  if (length(no_stage) > 0)
    stop("No e1b stage for genotype(s): ", paste(no_stage, collapse = ", "))
  no_name <- setdiff(unique(dat$genotype), names(name_map))
  if (length(no_name) > 0)
    stop("No display name for genotype(s): ", paste(no_name, collapse = ", "))

  # Aggregate: per-cell fractions AND pseudobulk ratios
  fracs <- dat %>%
    group_by(genotype, domain, lib_member) %>%
    summarise(
      frac_pitx2c  = mean(has_pitx2c),
      frac_pitx2ab = mean(has_pitx2ab),
      # Pseudobulk: sum UMIs across all cells in group
      pb_ab    = sum(pitx2ab_umi),
      pb_c     = sum(pitx2c_umi),
      pb_sh    = sum(shared_umi),
      pb_total = sum(pitx2ab_umi + pitx2c_umi + shared_umi),
      n_cells = n(),
      .groups = "drop") %>%
    mutate(
      # Pseudobulk ratios (add pseudocount of 1 to avoid /0)
      pb_ab_frac = pb_ab / pmax(pb_total, 1),   # a/b share of total
      pb_c_frac  = pb_c / pmax(pb_total, 1),     # c share of total
      stage = unname(e1b_map[genotype]),
      display = unname(name_map[genotype]))
  return(fracs)
}

message("Processing stages:")
all_fracs <- lapply(stage_configs, process_stage) %>% bind_rows()
message(sprintf("\nStages: %s", paste(sort(unique(all_fracs$stage)), collapse = ", ")))

# ── Pivot helper ──
# cell_floor: both halves need >= cell_floor nuclei (same 15-cell floor as the main
# panel, generate_paired_laterality.R).  The argument must not share a name with the
# min_cells COLUMN: inside filter() the column masks the argument, and the former
# `filter(min_cells >= min_cells)` compared the column with itself (a no-op floor).
pivot_and_filter <- function(fracs_long, value_col, cell_floor = 15,
                             min_pb_total = 0) {
  # For pseudobulk ratios, also pivot pb_total for filtering
  has_pb <- "pb_total" %in% names(fracs_long) && min_pb_total > 0

  cols_to_pivot <- c("value", "n_cells")
  wide <- fracs_long %>%
    select(genotype, display, domain, stage, lib_member,
           value = all_of(value_col), n_cells,
           any_of("pb_total")) %>%
    pivot_wider(id_cols = c(genotype, display, domain, stage),
                names_from = lib_member,
                values_from = c(value, n_cells,
                                if (has_pb) "pb_total" else NULL)) %>%
    rename(val_A = `value_Library A`, val_B = `value_Library B`,
           n_A = `n_cells_Library A`, n_B = `n_cells_Library B`) %>%
    mutate(min_cells = pmin(n_A, n_B)) %>%
    filter(!is.na(val_A), !is.na(val_B))

  low <- wide$min_cells < cell_floor
  if (any(low)) message(sprintf("  %s: dropping %d pairs with < %d cells in a half",
                                value_col, sum(low), cell_floor))
  wide <- wide[!low, ]

  # Filter by minimum pseudobulk total UMIs (both sides must pass)
  if (has_pb && "pb_total_Library A" %in% names(wide)) {
    wide <- wide %>%
      rename(pb_A = `pb_total_Library A`, pb_B = `pb_total_Library B`) %>%
      mutate(min_pb = pmin(pb_A, pb_B))
    low_pb <- is.na(wide$min_pb) | wide$min_pb < min_pb_total
    if (any(low_pb)) message(sprintf("  %s: dropping %d pairs with < %d locus UMIs in a half",
                                     value_col, sum(low_pb), min_pb_total))
    wide <- wide[!low_pb, ]
  }

  ps <- intersect(stage_levels, unique(wide$stage))
  wide %>% mutate(domain = factor(domain, levels = target_domains),
                  stage = factor(stage, levels = ps))
}

wide_c  <- pivot_and_filter(all_fracs, "frac_pitx2c")
wide_ab <- pivot_and_filter(all_fracs, "frac_pitx2ab")

# ── Summary ──
print_summary <- function(w, label) {
  message(sprintf("\n--- %s ---", label))
  w %>% mutate(d = abs(val_A - val_B)) %>%
    group_by(stage, domain) %>%
    summarise(n = n(), md = round(mean(d), 4), .groups = "drop") %>%
    {for (i in seq_len(nrow(.)))
      message(sprintf("  %s × %s (n=%d): mean |diff| = %.4f",
                      .[i, ]$stage, .[i, ]$domain, .[i, ]$n, .[i, ]$md))}
}
print_summary(wide_c, "PITX2c+ fraction")
print_summary(wide_ab, "PITX2a/b+ fraction")

# ── Plotting ──
make_plot <- function(wide_df, xlab, ylab, axis_max_override = NULL) {
  axis_max <- max(c(wide_df$val_A, wide_df$val_B), na.rm = TRUE) * 1.15
  axis_max <- max(axis_max, 0.05)
  if (!is.null(axis_max_override)) axis_max <- axis_max_override

  ggplot(wide_df, aes(x = val_A, y = val_B)) +
    facet_wrap(~ domain, nrow = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey60", linewidth = 0.5) +
    geom_point(aes(fill = stage, size = min_cells),
               shape = 21, color = "black", stroke = 0.4, alpha = 0.8) +
    geom_text_repel(aes(label = display, color = stage), size = 2.2,
                    segment.color = "grey75", segment.size = 0.25,
                    min.segment.length = 0.2, max.overlaps = 30,
                    box.padding = 0.3, point.padding = 0.2,
                    seed = 42, show.legend = FALSE) +
    scale_x_continuous(limits = c(-0.005, axis_max),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(limits = c(-0.005, axis_max),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = stage_colors, name = "Stage", drop = TRUE) +
    scale_color_manual(values = stage_colors, guide = "none", drop = TRUE) +
    scale_size_continuous(range = c(1.5, 4.5), breaks = c(30, 100, 300),
                          name = "Min. cells\n(per half)") +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          strip.background = element_rect(fill = "grey95"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          legend.position = "right",
          plot.margin = margin(5, 10, 5, 5)) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 3)),
           size = guide_legend(order = 2))
}

# ── Per-cell fraction plots (kept for reference) ──
p_c <- make_plot(wide_c, "PITX2c+ fraction (Library A)", "PITX2c+ fraction (Library B)")
ggsave(file.path(script_dir, "panel_C_paired_laterality_pitx2c.png"),
       p_c, width = 16, height = 5.2, dpi = 300, bg = "white")

p_ab <- make_plot(wide_ab, "PITX2a/b+ fraction (Library A)", "PITX2a/b+ fraction (Library B)")
ggsave(file.path(script_dir, "panel_C_paired_laterality_pitx2ab.png"),
       p_ab, width = 16, height = 5.2, dpi = 300, bg = "white")

# ── Pseudobulk ratio plots (key results) ──
# Require ≥5 total locus UMIs per side to avoid noise from 1-2 UMI denominators
wide_pb_ab <- pivot_and_filter(all_fracs, "pb_ab_frac", min_pb_total = 5)
wide_pb_c  <- pivot_and_filter(all_fracs, "pb_c_frac", min_pb_total = 5)

print_summary(wide_pb_ab, "Pseudobulk PITX2a/b / total ratio")
print_summary(wide_pb_c, "Pseudobulk PITX2c / total ratio")

# Pseudobulk a/b ratio: high = bilateral isoform dominates
# Prediction: dorsal = high & on diagonal; ventral/syringeal = one side lower (more c)
p_pb_ab <- make_plot(wide_pb_ab,
                     "PITX2a/b share of locus UMIs (Library A)",
                     "PITX2a/b share of locus UMIs (Library B)")
out_pb_ab <- file.path(script_dir, "panel_C_paired_laterality_pb_ab_ratio.png")
ggsave(out_pb_ab, p_pb_ab, width = 16, height = 5.2, dpi = 300, bg = "white")
message("\nSaved: ", out_pb_ab)

# Pseudobulk c ratio: cap at 30% to zoom in on the relevant range
p_pb_c <- make_plot(wide_pb_c,
                    "PITX2c share of locus UMIs (Library A)",
                    "PITX2c share of locus UMIs (Library B)",
                    axis_max_override = 0.30)
out_pb_c <- file.path(script_dir, "panel_C_paired_laterality_pb_c_ratio.png")
ggsave(out_pb_c, p_pb_c, width = 16, height = 5.2, dpi = 300, bg = "white")
message("Saved: ", out_pb_c)

# ── Pseudobulk raw data table ──
message("\n--- Pseudobulk isoform counts ---")
all_fracs %>%
  arrange(domain, stage, display, lib_member) %>%
  select(display, stage, domain, lib_member, n_cells, pb_ab, pb_c, pb_sh, pb_total,
         ab_frac = pb_ab_frac, c_frac = pb_c_frac) %>%
  mutate(across(c(ab_frac, c_frac), ~ round(., 3))) %>%
  as.data.frame() %>%
  print()
