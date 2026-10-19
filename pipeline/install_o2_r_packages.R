#!/usr/bin/env Rscript
# Install / verify the R packages required by O2 cluster pipeline steps.
# Target lib: ~/R/4.4.2/  (added to .libPaths() in every cluster R script).
#
# Currently captures the DESeq2 + apeglm install needed for the h2 DESeq2
# step (`pipeline/h2_deseq2_snrna.R`).  Extend this list as new packages
# become required.
#
# Idempotent: every package is checked for an existing namespace before
# install attempts; already-present packages just smoke-test load.
#
# Run via the sbatch wrapper: `sbatch pipeline/run_install_o2_r_packages.sh`.
# Compute nodes have outbound HTTPS to CRAN/Bioconductor (verified
# 2026-05-09 with the DESeq2 install completing in 3m36s on priority).

USER_LIB <- "~/R/4.4.2"
.libPaths(c(USER_LIB, .libPaths()))

cat("R session:\n")
cat("  R version:    ", R.version.string, "\n", sep = "")
cat("  Bioc version: ", as.character(BiocManager::version()), "\n", sep = "")
cat("  User lib:     ", path.expand(USER_LIB), "\n", sep = "")
cat("  libPaths:\n")
for (p in .libPaths()) cat("    -", p, "\n")
cat("\n")

# Non-interactive CRAN mirror (BiocManager falls back to CRAN for non-Bioc
# deps like locfit, gtools, etc.)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Required Bioconductor packages.  Add new entries here as cluster R
# scripts grow; keep them paired with a one-line rationale.
TARGETS <- c(
    "DESeq2",   # h2_deseq2_snrna.R — pseudobulk L vs R DESeq2 (paired
                # + unpaired+sex models); also fig03 DESeq2 in g3.
    "apeglm"    # standard companion for lfcShrink(type='apeglm');
                # not currently called but installed for parity with
                # the local duckfig env (environment.yml) so future
                # scripts can use either env interchangeably.
)

for (pkg in TARGETS) {
    cat("==== ", pkg, " ====\n", sep = "")
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("Already installed: ", pkg, " ",
            as.character(packageVersion(pkg)), "\n", sep = "")
    } else {
        cat("Installing ", pkg, " into ", path.expand(USER_LIB), " ...\n",
            sep = "")
        BiocManager::install(
            pkg,
            lib = path.expand(USER_LIB),
            update = FALSE, ask = FALSE,
            quiet = FALSE
        )
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Install of ", pkg,
                 " did not produce a loadable namespace.")
        }
        cat("Installed:        ", pkg, " ",
            as.character(packageVersion(pkg)), "\n", sep = "")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cat("Loaded successfully: ", pkg, "\n\n", sep = "")
}

# Smoke test — confirm the canonical entry points exist before exit so
# a partial install is caught here, not in a downstream pipeline run.
cat("==== Final smoke test ====\n")
cat("DESeq2::DESeq exists:           ",
    is.function(DESeq2::DESeq), "\n", sep = "")
cat("DESeqDataSetFromMatrix exists:  ",
    is.function(DESeq2::DESeqDataSetFromMatrix), "\n", sep = "")
cat("apeglm::apeglm exists:          ",
    is.function(apeglm::apeglm), "\n", sep = "")
cat("\nAll targets installed and loadable.\n")
