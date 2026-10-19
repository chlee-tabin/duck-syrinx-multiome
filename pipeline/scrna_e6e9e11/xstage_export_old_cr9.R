#!/usr/bin/env Rscript
# =====================================================================
# Cross-stage lineage export — OLD side, CORRECTED cr9-CB version.
#
# Re-export of xstage_export_old.R for the CORRECTED pipeline: the cr9 (intron-
# inclusive) + CellBender master with the v2 24-embryo demux + categorical
# sex_numeric ccreg, labelled by the cr9 compendium.  REPLACES the prior non-CB
# old-side export (which was deliberately non-CB only because the OLD CB marker-
# CLASSIFIER lost the E09 Chondrogenic call — now moot: the compendium U-method on
# cr9-CB recovers E09 Chondrogenic AND surfaces the E09 Somatic progenitor pool).
# Using cr9-CB here makes BOTH sides CellBender-corrected — cleaner apples-to-apples
# vs the (CB) multiome than the old non-CB old-side.
#
# Writes to a SEPARATE link dir (xstage_link_cr9/) so the prior non-CB export is
# preserved; the multiome mm_* are symlinked in by the runner (unchanged).
# Output per stage (identical schema to xstage_export_old.R, consumed by
# xstage_correspondence.R): old_<stage>_{matrix.mtx.gz,features.tsv.gz,
# barcodes.tsv.gz,obs.tsv.gz,pb.tsv.gz,pb_meta.tsv}.
#
# Labels: the cr9 compendium per-cell population (substage_<stage>_percell_compendium
# .tsv.gz: barcode,library,cluster,population,stage) joined to the master meta by
# (library, barcode) — the master carries `barcode`(=barcode16) + `library`, so no
# bridge file is needed.  RAW (CellBender-cleaned) counts; normalization deferred
# downstream, identically to the multiome export.
# Run on o2 via SLURM (reads the 2.27 GB master rds).  R/4.4.2.  Exploratory only.
# =====================================================================
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressWarnings(suppressMessages({
  .libPaths(c(paste0(HOME_ROOT, "/R/4.4.2"), .libPaths()))
  library(Seurat); library(Matrix); library(jsonlite)
}))
set.seed(1); options(future.globals.maxSize = 32 * 1024^3)
try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
say <- function(...) cat(sprintf(...), "\n"); t0 <- Sys.time()

CR9    <- paste0(SCRATCH_USER, "/old_duck_cr9")
OBJRDS <- file.path(CR9, "master", "cr9_03_integrated_CB.rds")   # cr9-CB master
COMPD  <- file.path(CR9, "substage")                             # cr9 compendium per-cell labels
OUT    <- Sys.getenv("XSTAGE_LINK_CR9", file.path(CR9, "xstage_link_cr9"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
STAGES <- c("E06", "E09", "E11")
EXCL   <- c("Epithelial (excl.)", "Endothelial (excl.)", "Immune (excl.)", "Low-signal",
            "Contamination", "Cycling")

say("== load cr9-CB master %s", OBJRDS); stopifnot(file.exists(OBJRDS))
obj <- readRDS(OBJRDS); DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
md <- obj@meta.data
C_stage <- if ("stage.called" %in% colnames(md)) "stage.called" else "stage"
stopifnot("barcode" %in% colnames(md), "library" %in% colnames(md))
say("master: %d genes x %d cells; stage col=%s", nrow(obj), ncol(obj), C_stage)

prov <- list(step = "xstage_export_old_cr9", obj_rds = OBJRDS,
             n_genes_namespace = nrow(obj), stages = list(),
             seurat = as.character(packageVersion("Seurat")), date = as.character(Sys.time()),
             note = "RAW CellBender counts (cr9-CB master, v2 demux, sex_numeric ccreg); compendium labels by (library,barcode)")

export_stage <- function(s) {
  say("\n#### STAGE %s ####", s)
  cmpf <- file.path(COMPD, sprintf("substage_%s_percell_compendium.tsv.gz", s))
  stopifnot(file.exists(cmpf))
  cmp <- read.table(gzfile(cmpf), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
  stopifnot(all(c("barcode", "library", "population") %in% colnames(cmp)))
  # subset the master to this stage, then label each cell by (library, barcode)
  cells_stage <- rownames(md)[!is.na(md[[C_stage]]) & md[[C_stage]] == s]
  ss <- subset(obj, cells = cells_stage)
  k_ss  <- paste(ss$library, ss$barcode)
  k_cmp <- paste(cmp$library, cmp$barcode)
  pop <- cmp$population[match(k_ss, k_cmp)]
  say("[%s] compendium join: %d/%d master-stage cells labelled", s, sum(!is.na(pop)), ncol(ss))
  keep <- !is.na(pop) & !(pop %in% EXCL)
  ss <- subset(ss, cells = colnames(ss)[keep]); pop <- pop[keep]
  say("[%s] %d labelled mesenchyme cells; labels: %s", s, ncol(ss),
      paste(sort(unique(pop)), collapse = " | "))

  cnt <- GetAssayData(ss, assay = "RNA", layer = "counts")
  if (!inherits(cnt, "dgCMatrix")) cnt <- as(cnt, "CsparseMatrix")
  n_umi <- Matrix::colSums(cnt)
  Matrix::writeMM(cnt, file.path(OUT, sprintf("old_%s_matrix.mtx", s)))
  if (system(sprintf("gzip -f %s", file.path(OUT, sprintf("old_%s_matrix.mtx", s)))) != 0)
    stop(sprintf("[%s] gzip failed for old_%s_matrix.mtx", s, s))
  gz <- gzfile(file.path(OUT, sprintf("old_%s_features.tsv.gz", s)), "w"); writeLines(rownames(cnt), gz); close(gz)
  gz <- gzfile(file.path(OUT, sprintf("old_%s_barcodes.tsv.gz", s)), "w"); writeLines(colnames(cnt), gz); close(gz)
  obsdf <- data.frame(cell_id = colnames(ss), library = ss$library, label = pop,
                      n_umi = n_umi, stringsAsFactors = FALSE)
  gz <- gzfile(file.path(OUT, sprintf("old_%s_obs.tsv.gz", s)), "w")
  write.table(obsdf, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

  pops <- sort(unique(pop))
  pb <- sapply(pops, function(p) Matrix::rowSums(cnt[, pop == p, drop = FALSE]))
  pb <- as.data.frame(as.matrix(pb)); rownames(pb) <- rownames(cnt); colnames(pb) <- pops
  gz <- gzfile(file.path(OUT, sprintf("old_%s_pb.tsv.gz", s)), "w")
  write.table(data.frame(gene = rownames(pb), pb, check.names = FALSE), gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  pbmeta <- data.frame(stage = s, population = pops, n_cells = as.integer(table(pop)[pops]),
                       dataset = "old_scrna", stringsAsFactors = FALSE)
  write.table(pbmeta, file.path(OUT, sprintf("old_%s_pb_meta.tsv", s)), sep = "\t", quote = FALSE, row.names = FALSE)
  say("[%s] wrote mtx (%d x %d) + pb (%d pops): %s", s, nrow(cnt), ncol(cnt), length(pops),
      paste(sprintf("%s=%d", pops, pbmeta$n_cells), collapse = " "))
  prov$stages[[s]] <<- list(n_cells = ncol(cnt), n_genes = nrow(cnt),
                            populations = as.list(setNames(pbmeta$n_cells, pops)))
}
ok <- vapply(STAGES, function(s)
  tryCatch({ export_stage(s); TRUE },
           error = function(e) { say("[%s] ERR %s", s, conditionMessage(e)); FALSE }),
  logical(1))
writeLines(toJSON(prov, pretty = TRUE, auto_unbox = TRUE), file.path(OUT, "xstage_export_old_cr9_provenance.json"))
# fail loud: a silently dropped stage would vanish from the cross-stage map while the job
# still reports success — abort before the DONE marker so the poller catches it.
if (!all(ok)) stop(sprintf("export failed for stage(s): %s", paste(STAGES[!ok], collapse = ", ")))
say("\n== xstage_export_old_cr9 DONE in %.1f min ==", as.numeric(difftime(Sys.time(), t0, "mins")))
cat("SBATCH_DONE_OK\n")
