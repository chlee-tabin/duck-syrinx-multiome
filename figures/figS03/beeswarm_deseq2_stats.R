#!/usr/bin/env Rscript
# Beeswarm — DESeq2 (size factors from n_cells) L/R stats per region x stage,
# on the SAME f5 per-(genotype, library) region counts the beeswarm plots.
# Deps-light (DESeq2 only, no ggplot) so it runs on the cluster where the
# beeswarm's viz packages are absent. Writes the parametric log2FC/padj the
# local render then annotates. Same model as the DESeq2 block in generate_beeswarm.R
# (n_cells size factors; ~genotype+side when >=2 paired genotypes, else ~side), but NOT a
# line-for-line mirror: that block has no single-region fallback and hardcodes MIN_CELLS.
# This script is the producer; the render annotates from its committed table.
.libPaths(c("~/R/4.4.2", .libPaths()))
suppressPackageStartupMessages(library(DESeq2))

counts_path <- "figures/data/step_f5_pitx2_region_counts.tsv"
# Env overrides use the repo's empty-normalizing idiom (a set-but-empty variable falls
# back to the default; two-arg Sys.getenv would return "" -> NA floor / stdout output).
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
out_path    <- .duck_env("OUT_PATH", "figures/data/step_f5_pitx2_region_deseq2.tsv")
# The reported threshold is 20; the sweep sidecar
# scripts/f5_threshold_sweep.sh re-runs this identical model at 10/20/30/50 into
# figures/data/step_f5_pitx2_region_deseq2_sweep.tsv so the PITX2c result's
# threshold-dependence is on record rather than asserted.
MIN_CELLS   <- suppressWarnings(as.integer(.duck_env("MIN_CELLS", "20")))
stopifnot(!is.na(MIN_CELLS), MIN_CELLS >= 0, nzchar(out_path))
cat(sprintf("beeswarm_deseq2_stats: MIN_CELLS=%d  OUT_PATH=%s\n", MIN_CELLS, out_path))

df <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- df[df$n_cells >= MIN_CELLS, , drop = FALSE]
region_levels <- c("PITX2c_promoter", "PITX2b_promoter", "PITX2a_promoter", "PITX2_genebody")
stage_levels  <- c("HH13", "HH17/16", "HH21", "HH25")

rows <- list()
for (reg in region_levels) {
    for (st in stage_levels) {
        sub <- df[df$region == reg & df$stage == st, , drop = FALSE]
        n_L <- sum(sub$side == "L"); n_R <- sum(sub$side == "R")
        # n_L/n_R count the halves that clear MIN_CELLS; n_L_fit/n_R_fit count the halves the
        # model ACTUALLY fitted (the paired branch drops genotypes with only one side), which
        # is what a reader checking the design needs.
        n_L_fit <- NA_integer_; n_R_fit <- NA_integer_
        lfc <- NA_real_; se <- NA_real_; padj <- NA_real_; model <- NA_character_; note <- ""
        if (n_L >= 2 && n_R >= 2) {
            cts <- matrix(as.integer(round(sub$n_frags)), nrow = 1)
            colnames(cts) <- paste0("s", seq_len(nrow(sub))); rownames(cts) <- "region"
            cd <- data.frame(sample = colnames(cts), genotype = sub$genotype,
                             side = sub$side, n_cells = sub$n_cells, row.names = colnames(cts))
            geno_sides <- table(cd$genotype, cd$side)
            both <- rowSums(geno_sides > 0) >= 2
            if (sum(both) >= 2) {
                keep <- cd$genotype %in% names(both)[both]
                cts <- cts[, keep, drop = FALSE]; cd <- cd[keep, , drop = FALSE]
                cd$genotype <- droplevels(factor(cd$genotype))
                cd$side <- factor(cd$side, levels = c("R", "L"))
                dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ genotype + side)
                model <- "paired(~genotype+side)"
            } else {
                cd$side <- factor(cd$side, levels = c("R", "L"))
                dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ side)
                model <- "unpaired(~side)"
            }
            n_L_fit <- sum(cd$side == "L"); n_R_fit <- sum(cd$side == "R")
            sizeFactors(dds) <- cd$n_cells / mean(cd$n_cells)
            # ONE region -> DESeq2's dispersion-TREND fit has nothing to fit and stops with
            # "all gene-wise dispersion estimates are within 2 orders of magnitude from the
            # minimum value ...". Until 2026-08-17 that error was swallowed (tryCatch -> NULL)
            # and the row shipped as NA with `model` still saying "paired" — narrated
            # downstream as "not estimable". It is not a cohort limit: DESeq2's
            # own message prescribes the single-feature path (gene-wise dispersion as final ->
            # Wald), which is applied here and RECORDED in `note`. Any other error is recorded
            # verbatim in `note` and the row stays NA — never a silent NULL.
            fit <- tryCatch(DESeq(dds, fitType = "mean", quiet = TRUE), error = function(e) e)
            if (inherits(fit, "error")) {
                if (grepl("dispersion estimates are within 2 orders of magnitude", conditionMessage(fit))) {
                    fit <- tryCatch({
                        d2 <- estimateDispersionsGeneEst(dds, quiet = TRUE)
                        dispersions(d2) <- mcols(d2)$dispGeneEst
                        nbinomWaldTest(d2, quiet = TRUE)
                    }, error = function(e) e)
                    note <- "single-region: gene-wise dispersion (DESeq2 trend fit not applicable)"
                }
            }
            if (inherits(fit, "error")) {
                # append (the single-region note records that the fallback was attempted) and
                # flatten: DESeq2 messages are multi-line, and this lands in a quote=FALSE TSV.
                note <- paste0(if (nzchar(note)) paste0(note, "; ") else "",
                               "fit_failed: ", conditionMessage(fit))
            } else {
                # One feature: Cook's-distance outlier filtering and independent filtering are
                # both meaningless here and would silently return padj = NA on a 3-6 replicate
                # design (reproduced on synthetic input).
                res <- results(fit, name = "side_L_vs_R", alpha = 0.05,
                               cooksCutoff = FALSE, independentFiltering = FALSE)
                lfc <- res$log2FoldChange[1]; se <- res$lfcSE[1]; padj <- res$padj[1]
                if (is.na(res$pvalue[1])) {
                    note <- paste0(if (nzchar(note)) paste0(note, "; ") else "",
                                   "DESeq2 returned NA pvalue")
                }
            }
        } else note <- sprintf("insufficient sides (n_L=%d, n_R=%d)", n_L, n_R)
        rows[[length(rows) + 1]] <- data.frame(
            region = reg, stage = st, n_L = n_L, n_R = n_R,
            n_L_fit = n_L_fit, n_R_fit = n_R_fit,
            deseq2_log2FC = lfc, deseq2_lfcSE = se, deseq2_padj = padj, model = model,
            min_cells = MIN_CELLS, note = note)
    }
}
out <- do.call(rbind, rows)
out$note <- gsub("[\t\r\n]+", " ", out$note)          # one row per row in a quote=FALSE TSV
stopifnot(!any(grepl("[\t\n]", out$note)))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %s (%d region x stage rows)\n", out_path, nrow(out)))
print(out)
cat("DONE_BEESWARM_DESEQ2\n")
