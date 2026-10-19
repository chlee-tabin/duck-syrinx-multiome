#!/usr/bin/env Rscript
#
# Phase 2 Step 2b: scDblFinder doublet detection per library
#
# Input:  filtered_feature_bc_matrix/ per library (loads from mtx/tsv directly)
# Output: phase2/step02b_scdblfinder_all.tsv (per-cell doublet scores)
#         phase2/figures/step02b_doublet_summary.tsv (per-library summary)
#
# Run via sbatch:
#   sbatch FC_08656/20260320_step02b_scdblfinder.sh
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

cat("=== Phase 2 Step 2b: scDblFinder Doublet Detection ===\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("scDblFinder version:", packageVersion("scDblFinder"), "\n\n"))

# ============================================================
# Configuration
# ============================================================
SCRATCH <- SCRATCH_ROOT
OUTDIR <- file.path(SCRATCH, "phase2")
FIGDIR <- file.path(OUTDIR, "figures")

libraries <- data.frame(
    name = c("GEX03-1A", "GEX03-1B", "GEX03-2",
             "GEX04", "GEX0S-duck",
             "GEX05", "GEX06", "GEX07", "GEX08",
             "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B"),
    cr_dir = c("GEX03-1A_GEX", "GEX03-1B_GEX", "GEX03-2_GEX",
               "GEX04_GEX", "GEX0S-duck_GEX",
               "GEX05_GEX", "GEX06_GEX", "GEX07_GEX", "GEX08_GEX",
               "PILOTGEX3_GEX", "PILOTGEX4A_GEX", "PILOTGEX4B_GEX"),
    stringsAsFactors = FALSE
)

# ============================================================
# Process each library
# ============================================================
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

    # Load as SingleCellExperiment
    cat("  Loading filtered matrix...\n")
    counts <- Read10X(filt_dir)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    cat(paste("  Loaded:", ncol(sce), "cells x", nrow(sce), "genes\n"))

    # Run scDblFinder
    cat("  Running scDblFinder...\n")
    set.seed(42)
    sce <- scDblFinder(sce)

    # Extract results
    barcodes <- colnames(sce)
    # Strip cellranger suffix (-1) to get 16bp barcode
    cb_16bp <- sub("-1$", "", barcodes)

    dbl_class <- sce$scDblFinder.class
    dbl_score <- sce$scDblFinder.score

    n_singlet <- sum(dbl_class == "singlet")
    n_doublet <- sum(dbl_class == "doublet")
    dbl_rate <- n_doublet / length(dbl_class)

    cat(paste("  Singlets:", n_singlet, "\n"))
    cat(paste("  Doublets:", n_doublet, "\n"))
    cat(paste("  Doublet rate:", round(dbl_rate * 100, 1), "%\n"))

    # Store per-cell results
    result_df <- data.frame(
        cb = cb_16bp,
        library = lib_name,
        scDblFinder_class = as.character(dbl_class),
        scDblFinder_score = round(dbl_score, 4),
        stringsAsFactors = FALSE
    )
    all_results[[lib_name]] <- result_df

    # Store summary
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

    # Clean up
    rm(counts, sce)
    gc(verbose = FALSE)
}

# ============================================================
# Save results
# ============================================================
# Per-cell results
all_df <- do.call(rbind, all_results)
out_path <- file.path(OUTDIR, "step02b_scdblfinder_all.tsv")
write.table(all_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("\nSaved per-cell results:", out_path, "\n"))
cat(paste("Total rows:", nrow(all_df), "\n"))

# Summary
summary_df <- do.call(rbind, summary_list)
summary_path <- file.path(FIGDIR, "step02b_doublet_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Saved summary:", summary_path, "\n"))

cat("\n=== scDblFinder Results Summary ===\n")
print(summary_df)

cat(paste("\nDone:", Sys.time(), "\n"))
