#!/usr/bin/env Rscript
# S15B sex-split readout: does the panel-B L/R comparison
# change if partitioned by sex? Mirrors tier2_lr_boost's substrate exactly (same loader,
# LR_POPS pool, MIN_CELLS >= 10 per embryo x side, per-embryo mean CP10K, paired L-R),
# then splits pairs by sex.called. Descriptive only - not a figure change.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
# repo root DERIVED from this script's own path (figures/figS08 -> up two); the
# sourced loader's DATA_DIR is repo-root-relative, so cwd must be the checkout root.
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)][1])
setwd(file.path(dirname(normalizePath(.script)), "..", ".."))
source("figures/figS08/cr9_screen_io.R")
d <- load_cr9_screen()
genes <- c("WNT5A", "SFRP2", "SHISA2", "PITX2")
long <- screen_long(d, genes, "cp10k__") |>
  filter(!is.na(side), stage %in% c("E09", "E11"), population %in% LR_POPS,
         !is.na(donor_id), donor_id != "",
         !donor_id %in% c("NA", "Negative", "Doublet", "doublet", "unassigned"))
sexmap <- d |>
  filter(!is.na(donor_id)) |>
  distinct(stage, donor_id, sex.called)
per_emb <- long |>
  group_by(stage, gene, donor_id, side) |>
  summarise(mean_cp10k = mean(value), n_cells = n(), .groups = "drop") |>
  filter(n_cells >= 10)
wide <- per_emb |>
  select(stage, gene, donor_id, side, mean_cp10k) |>
  pivot_wider(names_from = side, values_from = mean_cp10k) |>
  filter(!is.na(Left) & !is.na(Right)) |>
  left_join(sexmap, by = c("stage", "donor_id"))
stopifnot(!any(is.na(wide$sex.called)))
tab <- wide |>
  group_by(stage, gene, sex = sex.called) |>
  summarise(n_pairs = n(),
            L_higher = sum(Left > Right),
            med_LmR = round(median(Left - Right), 3),
            p = ifelse(n() >= 3,
                       suppressWarnings(round(wilcox.test(Left, Right, paired = TRUE,
                                                          exact = FALSE)$p.value, 3)),
                       NA_real_),
            .groups = "drop") |>
  arrange(gene, stage, sex)
print(as.data.frame(tab), row.names = FALSE)
cat("\nsex of each paired donor:\n")
print(as.data.frame(wide |> distinct(stage, donor_id, sex.called) |> arrange(stage, donor_id)), row.names = FALSE)
