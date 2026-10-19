#!/usr/bin/env Rscript
# Phase G Step 3 (CFC1 sibling): DESeq2 L vs R per (region × stage) on
# the f5_cfc1 region count table, mirroring the unpaired+sex design that
# produced the mRNA padj=0.007 finding.
#
# Why this exists separately from generate_beeswarm.R: the beeswarm R
# script is a local-render artifact that runs without DESeq2 (DESeq2
# isn't in the base R distribution).  This script runs on the cluster
# (~/R/4.4.2 user lib has DESeq2 via h2's installation) and produces
# the proper sex-covariate-controlled stats that get joined back into
# the beeswarm legend.
#
# Mirrors h2_deseq2_snrna.R's design philosophy.  This script
# unconditionally fits three designs per (region × stage), gating each
# on data sufficiency (no flags — the previous --unpaired / --control-sex
# CLI flags from h2 are inlined as automatic fits):
#   1. Unpaired (~ side):           always (when n_L>=2, n_R>=2).
#   2. Unpaired+sex (~ sex + side): when >=2 sex levels remain after
#      dropping ?-sex rows AND n_L>=2, n_R>=2.  This matches the
#      unpaired+sex run that surfaced CFC1 R>L at HH13
#      (padj=0.007 vs paired padj=0.18 in the snRNA-seq sweep).
#   3. Paired (~ genotype + side):  when n_paired_genotypes>=2.  Sex is
#      collinear with genotype in the paired design (each genotype = one
#      embryo = one sex), so the paired model already controls for sex
#      via the genotype intercepts.
#
# Differences from h2:
#   - h2 has many genes; this has 1 region per call (CFC1_promoter or
#     CFC1_genebody).  DESeq2 needs >=1 row to fit; can't estimate
#     per-gene dispersion from 1 row, so we use ``fitType = "mean"``
#     (constant dispersion across rows) — matches the PITX2 region
#     DESeq2 in figures/figS03/beeswarm_deseq2_stats.R.
#   - Size factors set explicitly to ``n_cells / mean(n_cells)`` since
#     poscounts on a single row would degenerate.
#
# Inputs
# ------
#   figures/data/step_f5_cfc1_region_counts.tsv  (from f5_cfc1)
#   figures/data/e1b_stage_assignments.tsv       (sex per genotype)
#
# Output
# ------
#   figures/data/step_f5_cfc1_region_deseq2.tsv
#     One row per (region × stage × design), columns:
#       region, stage, design, n_L, n_R, n_paired_genotypes,
#       log2FoldChange, lfcSE, pvalue, padj_within_region_stage
#
# Run on O2
# ---------
#   sbatch pipeline/run_g3_cfc1_deseq2.sh

.libPaths(c("~/R/4.4.2", .libPaths()))

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tidyr)
})

# ── Args ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
    idx <- match(paste0("--", name), args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[idx + 1]
}

counts_path <- get_arg("counts",
                       "figures/data/step_f5_cfc1_region_counts.tsv")
e1b_path    <- get_arg("e1b-tsv",
                       "figures/data/e1b_stage_assignments.tsv")
out_path    <- get_arg("out-tsv",
                       "figures/data/step_f5_cfc1_region_deseq2.tsv")
min_cells   <- as.integer(get_arg("min-cells", "20"))

cat("=== g3 CFC1 DESeq2 ===\n")
cat(sprintf("  counts:    %s\n", counts_path))
cat(sprintf("  e1b:       %s\n", e1b_path))
cat(sprintf("  out:       %s\n", out_path))
cat(sprintf("  min_cells: %d\n", min_cells))

# ── Load ─────────────────────────────────────────────────────────────────
df <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Loaded %d rows from %s\n", nrow(df), counts_path))

# Apply MIN_CELLS floor (mirrors generate_beeswarm.R)
n_pre <- nrow(df)
df <- df[df$n_cells >= min_cells, , drop = FALSE]
cat(sprintf("  Applied MIN_CELLS=%d: %d -> %d rows (%d dropped)\n",
            min_cells, n_pre, nrow(df), n_pre - nrow(df)))

e1b <- read.delim(e1b_path, stringsAsFactors = FALSE, check.names = FALSE)
sex_map <- setNames(e1b$sex, as.character(e1b$genotype_id))
cat(sprintf("Loaded sex map: %d genotypes (M=%d, F=%d, ?=%d)\n",
            length(sex_map),
            sum(sex_map == "M", na.rm = TRUE),
            sum(sex_map == "F", na.rm = TRUE),
            sum(is.na(sex_map) | sex_map == "?")))

df$sex <- sex_map[as.character(df$genotype)]
n_no_sex <- sum(is.na(df$sex) | df$sex == "?")
if (n_no_sex > 0) {
    cat(sprintf(paste0("  WARN: %d rows with missing/?-sex (kept; will skip ",
                       "sex covariate if entire stage is ?-only)\n"), n_no_sex))
}

region_levels <- c("CFC1_promoter", "CFC1_genebody")
# Derive stage_levels from the data so a new pipeline stage doesn't get
# silently skipped; warn loudly if anything outside the expected set
# appears.  KNOWN_STAGES is informational, not a filter.
KNOWN_STAGES <- c("HH13", "HH17/16", "HH21", "HH25")
stage_levels <- intersect(KNOWN_STAGES, unique(df$stage))
unknown_stages <- setdiff(unique(df$stage), KNOWN_STAGES)
if (length(unknown_stages) > 0) {
    cat(sprintf("  WARN: %d unknown stage(s) in input: %s — extending stage_levels\n",
                length(unknown_stages), paste(unknown_stages, collapse=",")))
    stage_levels <- c(stage_levels, unknown_stages)
}
missing_stages <- setdiff(KNOWN_STAGES, stage_levels)
if (length(missing_stages) > 0) {
    cat(sprintf("  WARN: KNOWN_STAGES not in input: %s — these get no DESeq2 row\n",
                paste(missing_stages, collapse=",")))
}

# ── Run DESeq2 per (region × stage × design) ──────────────────────────────
results_rows <- list()

for (reg in region_levels) {
    for (st in stage_levels) {
        sub <- df[df$region == reg & df$stage == st, , drop = FALSE]
        n_L <- sum(sub$side == "L")
        n_R <- sum(sub$side == "R")
        cat(sprintf("\n--- %s × %s: n_L=%d, n_R=%d ---\n", reg, st, n_L, n_R))

        if (n_L < 2 || n_R < 2) {
            cat("  Skip: need >=2 each side\n")
            next
        }

        # n_paired_genotypes (within-embryo L+R)
        geno_sides <- table(sub$genotype, sub$side)
        n_paired_genos <- sum(rowSums(geno_sides > 0) >= 2)

        # Build counts matrix (1 region × n_samples)
        cts <- matrix(sub$n_frags, nrow = 1)
        colnames(cts) <- paste0("s", seq_len(nrow(sub)))
        rownames(cts) <- reg

        cd <- data.frame(
            sample = colnames(cts),
            genotype = factor(sub$genotype),
            sex = factor(ifelse(is.na(sub$sex) | sub$sex == "?", NA,
                                 sub$sex)),
            side = factor(sub$side, levels = c("R", "L")),
            n_cells = sub$n_cells,
            row.names = colnames(cts)
        )

        # ── Design 1: unpaired (~ side) ───────────────────────────────
        run_deseq <- function(design_formula, design_label, cd_local,
                              cts_local) {
            tryCatch({
                dds <- DESeqDataSetFromMatrix(
                    countData = cts_local, colData = cd_local,
                    design = design_formula
                )
                # Manual size factors: n_cells / mean(n_cells).  poscounts
                # on a single row degenerates because the median-of-ratios
                # has no diversity.
                sizeFactors(dds) <- cd_local$n_cells / mean(cd_local$n_cells)
                # fitType="mean" — single-row data can't fit per-gene
                # dispersion.  Empirical mean dispersion is the standard
                # workaround (also h2's fallback for low-row scopes).
                dds <- DESeq(dds, fitType = "mean", quiet = TRUE)
                res <- results(dds, name = "side_L_vs_R", alpha = 0.05)
                tibble::tibble(
                    region = reg, stage = st, design = design_label,
                    n_L = sum(cd_local$side == "L"),
                    n_R = sum(cd_local$side == "R"),
                    n_paired_genotypes = if (design_label == "paired")
                        nlevels(droplevels(cd_local$genotype)) else NA_integer_,
                    log2FoldChange = res$log2FoldChange[1],
                    lfcSE = res$lfcSE[1],
                    pvalue = res$pvalue[1],
                    padj = res$padj[1],
                    baseMean = res$baseMean[1],
                    sex_levels = paste(levels(droplevels(cd_local$sex)),
                                       collapse = ","),
                    error_message = NA_character_
                )
            }, error = function(e) {
                # Tolerated DESeq2 fit failures (most commonly rank-deficient
                # paired+sex designs on small cohorts) emit an NA-row tagged
                # with the error message so silent failures are traceable.
                # The message also goes to stdout for log-grepping.
                msg <- conditionMessage(e)
                cat(sprintf("    [%s] FAILED: %s\n", design_label, msg))
                tibble::tibble(
                    region = reg, stage = st, design = design_label,
                    n_L = sum(cd_local$side == "L"),
                    n_R = sum(cd_local$side == "R"),
                    n_paired_genotypes = NA_integer_,
                    log2FoldChange = NA_real_, lfcSE = NA_real_,
                    pvalue = NA_real_, padj = NA_real_,
                    baseMean = NA_real_, sex_levels = NA_character_,
                    error_message = msg
                )
            })
        }

        # 1. unpaired
        cat("  [unpaired ~side]\n")
        results_rows[[length(results_rows) + 1]] <-
            run_deseq(~ side, "unpaired", cd, cts)

        # 2. unpaired+sex (drop NA-sex rows so DESeq2 doesn't error;
        #    require >=2 sex levels remaining)
        cd_sex <- cd[!is.na(cd$sex), , drop = FALSE]
        cts_sex <- cts[, rownames(cd_sex), drop = FALSE]
        sex_levels_present <- levels(droplevels(cd_sex$sex))
        if (length(sex_levels_present) >= 2 &&
            sum(cd_sex$side == "L") >= 2 && sum(cd_sex$side == "R") >= 2) {
            cat("  [unpaired+sex ~sex+side]\n")
            results_rows[[length(results_rows) + 1]] <-
                run_deseq(~ sex + side, "unpaired+sex", cd_sex, cts_sex)
        } else {
            cat(sprintf("  [unpaired+sex] Skip: %d sex levels, n_L=%d, n_R=%d\n",
                        length(sex_levels_present),
                        sum(cd_sex$side == "L"), sum(cd_sex$side == "R")))
        }

        # 3. paired (~ genotype + side); requires ≥2 paired genotypes
        if (n_paired_genos >= 2) {
            paired_geno_names <- rownames(geno_sides)[
                rowSums(geno_sides > 0) >= 2]
            keep <- as.character(cd$genotype) %in% paired_geno_names
            cd_p <- cd[keep, , drop = FALSE]
            cd_p$genotype <- droplevels(cd_p$genotype)
            cts_p <- cts[, rownames(cd_p), drop = FALSE]
            cat(sprintf("  [paired ~genotype+side; n_paired_genos=%d]\n",
                        nlevels(cd_p$genotype)))
            results_rows[[length(results_rows) + 1]] <-
                run_deseq(~ genotype + side, "paired", cd_p, cts_p)
        } else {
            cat(sprintf("  [paired] Skip: only %d paired genotypes\n",
                        n_paired_genos))
        }
    }
}

stats <- bind_rows(results_rows)
# DESeq2 error messages span several lines; with quote = FALSE an embedded newline splits one
# row across physical lines and breaks every TSV reader.  Collapse whitespace runs (incl. tabs)
# to single spaces so each (region x stage x design) stays one line.
stats$error_message <- gsub("[[:space:]]+", " ", trimws(stats$error_message))
write.table(stats, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
print(stats, n = 100)
