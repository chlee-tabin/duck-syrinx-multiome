#!/usr/bin/env Rscript
# region_stats_wilcoxon.R — canonical producer of the legend/Methods-cited
# per-(region × stage) L-vs-R Wilcoxon (+ DESeq2 reference) stats tables:
#
#   GENE=pitx2  ->  figures/data/step_f5_pitx2_region_stats.tsv
#   GENE=cfc1   ->  figures/data/step_f5_cfc1_region_stats.tsv
#
# Minimal reproducible script for the quoted Wilcoxon p-values (extracted from the
# stats sections of internal beeswarm renderers, not distributed; Wilcoxon logic,
# floor and column set unchanged).
#
# DETERMINISM (review fix, same day): the originals attempted an IN-SCRIPT
# DESeq2 fit whose result depended on whether the running env carried DESeq2
# (duckfig does NOT — requireNamespace("DESeq2") is FALSE there), which is how
# the pre-refresh pitx2 table became a mixed-provenance artifact. This script
# therefore computes NO in-script DESeq2: the deseq2_* columns come from the
# cluster-side sidecar (pitx2, hard-asserted) or are NA by construction (cfc1,
# matching the committed table). Output is byte-stable in any R env.
#
#   GENE=pitx2 conda run -n duckfig Rscript figures/figS03/region_stats_wilcoxon.R
#   GENE=cfc1  conda run -n duckfig Rscript figures/figS03/region_stats_wilcoxon.R

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
})

GENE <- Sys.getenv("GENE", "pitx2")
if (GENE == "pitx2") {
    counts_path <- "figures/data/step_f5_pitx2_region_counts.tsv"
    stats_path  <- "figures/data/step_f5_pitx2_region_stats.tsv"
    region_levels <- c("PITX2c_promoter", "PITX2b_promoter",
                       "PITX2a_promoter", "PITX2_genebody")
    # PITX2: the reported DESeq2 columns come from the cluster-side model
    # (figures/figS03/beeswarm_deseq2_stats.R -> step_f5_pitx2_region_deseq2.tsv),
    # merged with hard assertions below.
    ds_merge_file <- "figures/data/step_f5_pitx2_region_deseq2.tsv"
} else if (GENE == "cfc1") {
    counts_path <- "figures/data/step_f5_cfc1_region_counts.tsv"
    stats_path  <- "figures/data/step_f5_cfc1_region_stats.tsv"
    region_levels <- c("CFC1_promoter", "CFC1_genebody")
    # cfc1: deseq2_* stay NA, matching the committed table. CFC1's DESeq2 lives
    # in its OWN committed table (step_f5_cfc1_region_deseq2.tsv, g3_cfc1_deseq2;
    # different schema) — it is deliberately NOT merged into this stats table.
    ds_merge_file <- NA_character_
} else {
    stop("GENE must be pitx2 or cfc1, got: ", GENE)
}
stage_levels <- c("HH13", "HH17/16", "HH21", "HH25")

# ── Stability floor (rationale: generate_beeswarm.R header; project standard
# at the pseudobulk-stats level is MIN=30, the stats/visualization layer uses 20)
MIN_CELLS <- 20

cat(sprintf("Loading %s\n", counts_path))
df <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
n_pre <- nrow(df)
df <- df[df$n_cells >= MIN_CELLS, , drop = FALSE]
cat(sprintf("  MIN_CELLS=%d floor: %d -> %d observations\n",
            MIN_CELLS, n_pre, nrow(df)))

df$region <- factor(df$region, levels = region_levels)
df$stage  <- factor(df$stage,  levels = stage_levels)
df$side   <- factor(df$side,   levels = c("L", "R"))

pair_key <- df %>%
    group_by(region, stage, genotype) %>%
    summarise(has_L = any(side == "L"),
              has_R = any(side == "R"),
              .groups = "drop") %>%
    mutate(is_paired = has_L & has_R)
df <- df %>%
    left_join(pair_key %>% select(region, stage, genotype, is_paired),
              by = c("region", "stage", "genotype"))

stats_rows <- list()
for (reg in region_levels) {
    for (st in stage_levels) {
        sub_all <- df %>% filter(region == reg, stage == st)
        sub_pair <- sub_all %>%
            filter(is_paired) %>%
            select(genotype, side, frags_per_cell) %>%
            pivot_wider(names_from = side, values_from = frags_per_cell)

        n_L <- sum(sub_all$side == "L")
        n_R <- sum(sub_all$side == "R")
        n_pairs <- nrow(sub_pair)

        unpair_p <- NA_real_; unpair_W <- NA_real_
        if (n_L >= 2 && n_R >= 2) {
            ut <- wilcox.test(frags_per_cell ~ side, data = sub_all,
                              exact = FALSE)
            unpair_p <- ut$p.value
            unpair_W <- unname(ut$statistic)
        }

        pair_p <- NA_real_; pair_V <- NA_real_; median_diff <- NA_real_
        if (n_pairs >= 3) {
            pt <- wilcox.test(sub_pair$L, sub_pair$R, paired = TRUE,
                              exact = FALSE)
            pair_p <- pt$p.value
            pair_V <- unname(pt$statistic)
            median_diff <- median(sub_pair$L - sub_pair$R)
        }

        m_L <- median(sub_all$frags_per_cell[sub_all$side == "L"])
        m_R <- median(sub_all$frags_per_cell[sub_all$side == "R"])

        # deseq2_* come from the sidecar merge below (pitx2) or stay NA (cfc1)
        ds_lfc <- NA_real_; ds_se <- NA_real_; ds_padj <- NA_real_

        stats_rows[[length(stats_rows) + 1]] <- tibble(
            region = reg, stage = st,
            n_L = n_L, n_R = n_R, n_pairs = n_pairs,
            unpaired_wilcox_p = unpair_p, unpaired_wilcox_W = unpair_W,
            paired_wilcox_p   = pair_p,   paired_wilcox_V   = pair_V,
            paired_median_LminusR = median_diff,
            deseq2_log2FC = ds_lfc, deseq2_lfcSE = ds_se,
            deseq2_padj = ds_padj,
            median_L = m_L, median_R = m_R
        )
    }
}
stats_df <- bind_rows(stats_rows) |>
    mutate(across(c(region, stage), as.character))

if (!is.na(ds_merge_file)) {
    # HARD requirements (review fix): the sidecar must exist, carry all three
    # DESeq2 columns, and cover every region x stage row — a partial merge is
    # exactly the mixed-provenance failure the pre-refresh table exhibited.
    if (!file.exists(ds_merge_file))
        stop("DESeq2 sidecar missing: ", ds_merge_file,
             " (regenerate via pipeline/run_beeswarm_deseq2.sh on the cluster —",
             " R 4.4.2 + DESeq2; the local duckfig env has no DESeq2)")
    ds <- read.delim(ds_merge_file, stringsAsFactors = FALSE, check.names = FALSE)
    need <- c("deseq2_log2FC", "deseq2_lfcSE", "deseq2_padj")
    if (!all(need %in% names(ds)))
        stop("DESeq2 sidecar lacks columns: ",
             paste(setdiff(need, names(ds)), collapse = ", "))
    # guard against the env-parameterized sweep sibling (beeswarm_deseq2_stats.R
    # can legitimately write other floors/duplicated keys under OUT_PATH override).
    # min_cells is REQUIRED — the producer writes it on every row; a sidecar
    # without it cannot prove its floor and must not merge.
    if (!("min_cells" %in% names(ds)))
        stop("DESeq2 sidecar lacks the min_cells column — cannot verify the floor")
    if (!isTRUE(all(ds$min_cells == MIN_CELLS)))
        stop("DESeq2 sidecar floor mismatch (or NA min_cells): min_cells != ", MIN_CELLS)
    key  <- paste(stats_df$region, stats_df$stage)
    dkey <- paste(ds$region, ds$stage)
    if (anyDuplicated(dkey))
        stop("DESeq2 sidecar has duplicated region x stage keys (sweep-variant file?)")
    idx  <- match(key, dkey)
    if (any(is.na(idx)))
        stop("DESeq2 sidecar has no row for: ",
             paste(key[is.na(idx)], collapse = "; "))
    stats_df$deseq2_log2FC <- ds$deseq2_log2FC[idx]
    stats_df$deseq2_lfcSE  <- ds$deseq2_lfcSE[idx]
    stats_df$deseq2_padj   <- ds$deseq2_padj[idx]
    cat(sprintf("Merged cluster DESeq2 stats from %s (all %d rows matched)\n",
                ds_merge_file, length(idx)))
}
write.table(stats_df, stats_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved: %s\n", stats_path))
