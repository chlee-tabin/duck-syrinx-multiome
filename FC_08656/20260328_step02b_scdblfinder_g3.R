#!/usr/bin/env Rscript
#
# Phase 2 Step 2b: scDblFinder doublet detection for G3 libraries
#
# IDENTICAL procedure to 20260320_step02b_scdblfinder.R but for GEX04/GEX0S-duck.
# Previously skipped — all libraries must use the same processing pipeline.
#
# Input:  filtered_feature_bc_matrix/ per library
# Output: phase2/step02b_scdblfinder_g3.tsv (per-cell doublet scores)
#
# Run via sbatch:
#   sbatch FC_08656/20260328_step02b_scdblfinder_g3.sh
#

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))

.libPaths(c("~/R/4.4.2", .libPaths()))

suppressPackageStartupMessages({
    library(Seurat)
    library(scDblFinder)
    library(SingleCellExperiment)
    library(Matrix)
    library(scuttle)
})

cat("=== scDblFinder G3 Libraries ===\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("scDblFinder version:", packageVersion("scDblFinder"), "\n\n"))

SCRATCH <- SCRATCH_ROOT
OUTDIR <- file.path(SCRATCH, "phase2")

libraries <- data.frame(
    name = c("GEX04", "GEX0S-duck"),
    cr_dir = c("GEX04_GEX", "GEX0S-duck_GEX"),
    stringsAsFactors = FALSE
)

all_results <- list()
summary_list <- list()

for (i in seq_len(nrow(libraries))) {
    lib_name <- libraries$name[i]
    cr_dir <- libraries$cr_dir[i]

    cat(paste0("\n=== ", lib_name, " ===\n"))

    filt_dir <- file.path(SCRATCH, cr_dir, "outs", "filtered_feature_bc_matrix")
    if (!dir.exists(filt_dir)) {
        cat(paste("  ERROR: filtered matrix not found:", filt_dir, "\n"))
        next
    }

    cat("  Loading filtered matrix...\n")
    counts <- Read10X(filt_dir)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    cat(paste("  Loaded:", ncol(sce), "cells x", nrow(sce), "genes\n"))

    # Run scDblFinder — IDENTICAL parameters to step02b_scdblfinder.R
    cat("  Running scDblFinder...\n")
    set.seed(42)
    sce <- scDblFinder(sce)

    barcodes <- colnames(sce)
    cb_16bp <- sub("-1$", "", barcodes)

    dbl_class <- sce$scDblFinder.class
    dbl_score <- sce$scDblFinder.score

    n_singlet <- sum(dbl_class == "singlet")
    n_doublet <- sum(dbl_class == "doublet")
    dbl_rate <- n_doublet / length(dbl_class)

    cat(paste("  Singlets:", n_singlet, "\n"))
    cat(paste("  Doublets:", n_doublet, "\n"))
    cat(paste("  Doublet rate:", round(dbl_rate * 100, 1), "%\n"))

    result_df <- data.frame(
        cb = cb_16bp,
        library = lib_name,
        scDblFinder_class = as.character(dbl_class),
        scDblFinder_score = round(dbl_score, 4),
        stringsAsFactors = FALSE
    )
    all_results[[lib_name]] <- result_df

    summary_list[[lib_name]] <- data.frame(
        library = lib_name,
        n_cells = ncol(sce),
        n_singlet = n_singlet,
        n_doublet = n_doublet,
        doublet_rate = round(dbl_rate, 4),
        median_score_singlet = round(median(dbl_score[dbl_class == "singlet"]), 4),
        median_score_doublet = round(median(dbl_score[dbl_class == "doublet"]), 4),
        stringsAsFactors = FALSE
    )

    rm(counts, sce)
    gc(verbose = FALSE)
}

# Save G3-specific results
all_df <- do.call(rbind, all_results)
out_path <- file.path(OUTDIR, "step02b_scdblfinder_g3.tsv")
write.table(all_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("\nSaved:", out_path, "\n"))
cat(paste("Total rows:", nrow(all_df), "\n"))

# Also APPEND to the main scdblfinder file so step02e can use a single file
main_path <- file.path(OUTDIR, "step02b_scdblfinder_all.tsv")
if (file.exists(main_path)) {
    existing <- read.table(main_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Remove any old G3 rows (in case of re-run)
    existing <- existing[!(existing$library %in% c("GEX04", "GEX0S-duck")), ]
    combined <- rbind(existing, all_df)
    write.table(combined, main_path, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(paste("Appended to:", main_path, " (", nrow(combined), "total rows)\n"))
} else {
    cat(paste("WARNING: main file not found:", main_path, "\n"))
    cat("G3 results saved separately; merge manually.\n")
}

summary_df <- do.call(rbind, summary_list)
cat("\n=== scDblFinder G3 Results ===\n")
print(summary_df)

cat(paste("\nDone:", Sys.time(), "\n"))
