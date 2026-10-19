#!/usr/bin/env Rscript
# archr_export_postprocess.R — re-export per-pool TSVs from an existing
# proj_fig03_callable, without redoing Arrow build / addGroupCoverages /
# addReproduciblePeakSet.
#
# Used as a recovery tool after job 37930531 succeeded at peak calling
# (283,818 peaks in union_peaks.bed) but crashed in the per-pool TSV
# export step due to ArchR's filename sanitisation:
#     pool_counts[[pool]]   →  subscript out of bounds
# because pool was the SANITISED name (e.g. "HH17.16.Cardiopharyngeal.LPM")
# while pool_counts has the ORIGINAL keys (e.g. "HH17/16|Cardiopharyngeal LPM").
# The fix lives in archr_build_peaks.R for future runs; this script
# rescues the existing run's outputs.

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(ArchR)
})

OUT_ROOT  <- paste0(SCRATCH_ROOT, "/phase2/atac_archr")
PEAK_DIR  <- file.path(OUT_ROOT, "peaks")
PROJ_DIR  <- file.path(OUT_ROOT, "proj_fig03_callable")

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

ts("loading ArchRProject from ", PROJ_DIR)
proj <- loadArchRProject(path = PROJ_DIR, showLogo = FALSE)
cc <- as.data.frame(getCellColData(proj))

# Mirror the SKIP_POOLS / MIN_POOL_CELLS logic from archr_build_peaks.R.
SKIP_POOLS <- c("HH21|Pharyngeal endoderm", "HH25|Pharyngeal endoderm")
MIN_POOL_CELLS <- 500L
valid <- !is.na(cc$peak_pool) & !(cc$peak_pool %in% SKIP_POOLS)
pool_counts <- table(cc$peak_pool[valid])
ts("found ", length(pool_counts), " peak_pool keys in cellColData")

# ArchR sanitises pool names: any non-[a-zA-Z0-9._] becomes ".".
archr_sanitise <- function(x) gsub("[^a-zA-Z0-9._]", ".", x)
sanitised_to_orig <- setNames(names(pool_counts),
                              archr_sanitise(names(pool_counts)))

pool_rds_dir <- file.path(PROJ_DIR, "PeakCalls")
pool_files <- list.files(pool_rds_dir,
                         pattern = "-reproduciblePeaks\\.gr\\.rds$",
                         full.names = TRUE, recursive = TRUE)
ts("found ", length(pool_files), " per-pool reproducible-peak RDS files")

per_pool_rows <- list()
repro_rows <- list()
for (f in pool_files) {
    pool_san <- sub("-reproduciblePeaks\\.gr\\.rds$", "", basename(f))
    pool <- sanitised_to_orig[[pool_san]]
    if (is.null(pool) || is.na(pool)) {
        ts("WARNING: cannot map '", pool_san, "' back to original; skipping")
        next
    }
    gr <- readRDS(f)
    n_pk <- length(gr)
    n_cells <- as.integer(pool_counts[[pool]])
    n_rep <- length(unique(cc$unified_genotype_id[cc$peak_pool == pool]))
    repro_rows[[length(repro_rows) + 1]] <- data.frame(
        peak_pool = pool, n_cells = n_cells,
        n_replicates = n_rep, n_peaks = n_pk,
        stringsAsFactors = FALSE
    )
    per_pool_rows[[pool]] <- data.frame(
        peak_pool = pool,
        chrom = as.character(seqnames(gr)),
        start = start(gr) - 1L,        # 0-based BED
        end   = end(gr),
        stringsAsFactors = FALSE
    )
}

if (length(repro_rows)) {
    repro_df <- do.call(rbind, repro_rows)
    write.table(repro_df, file.path(PEAK_DIR, "reproducibility_metadata.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    ts("wrote reproducibility_metadata.tsv (", nrow(repro_df), " pools)")
    per_pool_df <- do.call(rbind, per_pool_rows)
    write.table(per_pool_df, file.path(PEAK_DIR, "per_pool_peaks.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    ts("wrote per_pool_peaks.tsv (", nrow(per_pool_df), " rows)")
} else {
    stop("no per-pool peak RDS files mapped — nothing to export")
}

ts("=== archr_export_postprocess.R done ===")
