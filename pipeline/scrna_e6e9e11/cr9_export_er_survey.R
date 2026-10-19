#!/usr/bin/env Rscript
# cr9_export_er_survey.R — per-cell log-normalized expression of the 6 estrogen
# receptors surveyed in manuscript Fig 6a (ESR1/ESR2/ESRRB/ESRRG/GPER1/RERG),
# from the cr9-CB MASTER (all cell classes), annotated per cell with BOTH the
# global cluster AND the named substage subcluster population (so the Fig 6a
# violin can show the proper LPM subclusters + Endoderm/Endothelium/Immune,
# not raw global cluster ids).
#
# Why the master: Fig 6a's claim is "ESR1 is the only ER expressed THROUGHOUT the
# mesoderm" — a mesoderm-vs-rest contrast.  The substage objects are mesoderm-only
# and carry only a 3-ER subset, so the 6-ER global panel must come from the master.
# The substage population label is joined back by (library, barcode); cells absent
# from every substage object are non-mesoderm (Endo/Endothelium/Immune) and get a
# NA substage_population (filled downstream from the global-cluster broad class).
#
# Output (scratch): old_duck_cr9/qc_export/cr9_er_survey_percell.tsv.gz
#   + cr9_global_cluster_markers.tsv (refreshes the stale committed copy)
# Repo destination (copied back from the cluster): figures/data/scrna_e6e9e11/
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")

suppressWarnings(.libPaths(c(file.path(Sys.getenv("HOME"), "R", "4.4.2"), .libPaths())))
suppressMessages({library(Seurat); library(Matrix)})
say <- function(...) cat(sprintf(...), "\n")

CR9    <- paste0(SCRATCH_USER, "/old_duck_cr9")
MASTER <- file.path(CR9, "master", "cr9_03_integrated_CB.rds")
SUBDIR <- file.path(CR9, "substage")
LABELS <- file.path(SUBDIR, "substage_compendium_clusterlabels.tsv")
OUTDIR <- file.path(CR9, "qc_export"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUT    <- file.path(OUTDIR, "cr9_er_survey_percell.tsv.gz")
STAGES <- c("E06", "E09", "E11")
RES    <- "RNA_snn_res.0.2"                # canonical substage resolution (matches compendium)

# manuscript Fig 6a surveyed these 6 estrogen receptors (same ZJU1.0/Gnomon namespace as cr6)
ER <- c("ESR1", "ESR2", "ESRRB", "ESRRG", "GPER1", "RERG")

GCOL <- "SCT_snn_res.0.4"                   # current-master global clustering (NOT seed-stable across reruns)
# Broad lineage marker panels (identical to export_cr6_master_broad.R) — used to compute a
# per-cell broad_class (argmax over panels, per current global cluster) so the Fig 6a non-mesoderm
# categories (Endoderm/Endothelium/Immune) are derived from the CURRENT 24,816 master and never
# from the committed STALE per-cluster markers TSV (built on the older 24,521 master; IDs shuffle).
RBC   <- c("LOC101798290","LOC113839608","LOC113844783","LOC101798111","LOC101792412","LOC101800520","LOC101800713")
BROAD <- list(LPM_mesenchyme=c("PRRX1","TWIST1","FOXF1","GATA4","COL1A1","COL3A1"),
              Epithelial=c("EPCAM","CDH1","CLDN1","WNT7B"), Endothelial=c("CDH5","PECAM1","EGFL7","KDR"),
              Immune=c("PTPRC","LCP1","CD74"), Blood=RBC, Neural=c("SOX10","STMN2","TUBB3","ELAVL4"),
              Muscle=c("MYH3","MYOD1","MYF5","ACTC1"))

# ── substage population lookup: (library|barcode) -> compendium population + embryo ──
if (!file.exists(LABELS)) stop(sprintf("FATAL: compendium labels missing: %s", LABELS))
labs <- read.delim(LABELS, check.names = FALSE)
lab_col <- if ("population" %in% names(labs)) "population" else "label"
say("== substage labels: %s (%d rows, label col '%s')", LABELS, nrow(labs), lab_col)

key2pop <- character(0); key2emb <- character(0)
for (st in STAGES) {
  pf <- file.path(SUBDIR, sprintf("substage_%s_percell.tsv.gz", st))
  if (!file.exists(pf)) { say("   WARNING: %s missing — skip", pf); next }
  pc <- read.delim(pf, check.names = FALSE)
  if (!(RES %in% names(pc))) stop(sprintf("FATAL: %s lacks %s", pf, RES))
  m  <- labs[labs$stage == st, ]
  cl2pop <- setNames(as.character(m[[lab_col]]), as.character(m$cluster))
  unmapped <- setdiff(unique(as.character(pc[[RES]])), names(cl2pop))
  if (length(unmapped)) stop(sprintf("FATAL: %s clusters not in labels: %s", st, paste(unmapped, collapse=",")))
  k   <- paste(pc$library, pc$barcode, sep = "|")
  pop <- cl2pop[as.character(pc[[RES]])]
  emb <- if ("embryo_id" %in% names(pc)) as.character(pc$embryo_id) else rep(NA_character_, nrow(pc))
  key2pop[k] <- pop; key2emb[k] <- emb
  say("   %s: %d cells -> %d populations", st, nrow(pc), length(unique(pop)))
}
say("== substage lookup built: %d (library|barcode) keys", length(key2pop))

# ── master: 6-ER log-norm + meta ──
say("== load master: %s", MASTER)
m <- readRDS(MASTER); DefaultAssay(m) <- "RNA"
m <- tryCatch(JoinLayers(m), error = function(e) m)   # single counts/data matrix (for broad detection)
say("   cells=%d  RNA features=%d", ncol(m), nrow(m[["RNA"]]))

rn <- rownames(m[["RNA"]])
present <- ER[ER %in% rn]; missing <- setdiff(ER, rn)
say("   ER present (%d/6): %s", length(present), paste(present, collapse = ","))
if (length(missing)) {
  say("   ER MISSING: %s", paste(missing, collapse = ","))
  for (fam in c("ESR", "GPER", "RERG", "ESRR", "NR3")) {
    hits <- grep(fam, rn, value = TRUE, ignore.case = TRUE)
    if (length(hits)) say("     candidates ~%s: %s", fam, paste(head(hits, 15), collapse = ","))
  }
}
if (!length(present)) stop("FATAL: none of the 6 ER symbols found in master RNA assay")

m <- NormalizeData(m, assay = "RNA", normalization.method = "LogNormalize",
                   scale.factor = 1e4, verbose = FALSE)
expr <- as.matrix(GetAssayData(m, assay = "RNA", layer = "data")[present, , drop = FALSE])

# ── per-cell broad lineage: argmax of marker-panel detection per CURRENT global cluster ──
# Computed on THIS master so the Fig 6a non-mesoderm categories (Endoderm/Endothelium/Immune)
# never depend on the committed STALE markers TSV (older 24,521 master; global cluster IDs
# shuffle on every rerun → that was why only Endoderm surfaced and Immune/Endothelium were lost).
if (!(GCOL %in% colnames(m@meta.data))) stop(sprintf("FATAL: master lacks global cluster col %s", GCOL))
clv      <- as.character(m@meta.data[[GCOL]])
det      <- as(GetAssayData(m, assay = "RNA", layer = "counts") > 0, "CsparseMatrix")
clusters <- sort(unique(clv))
bm <- data.frame(global_cluster = clusters,
                 n_cells = as.integer(table(clv)[clusters]), stringsAsFactors = FALSE)
for (k in names(BROAD)) {
  g <- intersect(BROAD[[k]], rownames(det))
  if (!length(g)) { say("   WARNING: broad panel %s has 0 genes in master", k); bm[[k]] <- NA_real_; next }
  fc <- Matrix::colMeans(det[g, , drop = FALSE])               # per-cell fraction of panel-k markers detected
  bm[[k]] <- sapply(clusters, function(c) mean(fc[clv == c]))  # cluster-mean detection
}
cls_cols <- names(BROAD)
mat <- as.matrix(bm[, cls_cols]); mat[is.na(mat)] <- -1        # never let an absent panel win argmax
bm$broad_class <- cls_cols[max.col(mat, ties.method = "first")]
cl2broad <- setNames(bm$broad_class, bm$global_cluster)
broad_class_percell <- unname(cl2broad[clv])
write.table(bm, file.path(OUTDIR, "cr9_global_cluster_markers.tsv"),   # refresh stale committed copy
            sep = "\t", quote = FALSE, row.names = FALSE)
say("== broad-class per global cluster (%s); cluster->lineage:", GCOL)
for (i in order(as.integer(bm$global_cluster)))
  say("   cl%s: n=%d -> %s", bm$global_cluster[i], bm$n_cells[i], bm$broad_class[i])

md <- m@meta.data
gv <- function(x) {
  if (!x %in% colnames(md)) stop("cr9 master meta.data lacks required column '", x, "'")
  md[[x]]
}
if (!"umap" %in% names(m@reductions)) stop("cr9 master has no 'umap' reduction")
um <- Embeddings(m, "umap")
mkey <- paste(gv("library"), gv("barcode"), sep = "|")

df <- data.frame(
  barcode            = gv("barcode"),
  library            = gv("library"),
  stage              = gv("stage.called"),
  global_cluster     = as.character(gv("SCT_snn_res.0.4")),
  broad_class        = broad_class_percell,      # current-master lineage (drives non-meso labels)
  substage_population = unname(key2pop[mkey]),   # NA = non-mesoderm (not in any substage object)
  donor_id           = gv("donor_id"),
  embryo_id          = if ("embryo_id" %in% colnames(md)) gv("embryo_id") else unname(key2emb[mkey]),
  LR_called          = gv("LR.called"),
  sex_called         = gv("sex.called"),
  stage_source       = gv("stage_source"),
  scDblFinder_class  = gv("scDblFinder_class"),
  UMAP_1             = um[, 1],
  UMAP_2             = um[, 2],
  check.names = FALSE, stringsAsFactors = FALSE
)
for (g in present) df[[paste0("log__", g)]] <- expr[g, ]

gz <- gzfile(OUT, "w"); write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
say("== wrote %s (%d cells x %d cols)", OUT, nrow(df), ncol(df))

# sanity to the log
nmeso <- sum(!is.na(df$substage_population))
say("   cells with a substage population (mesoderm): %d / %d", nmeso, nrow(df))
say("   broad_class of NON-mesoderm cells (substage_population NA) — Endoderm/Endothelium/Immune should appear:")
print(sort(table(df$broad_class[is.na(df$substage_population)]), decreasing = TRUE))
say("   substage_population counts:")
print(sort(table(df$substage_population, useNA = "ifany"), decreasing = TRUE))
if ("log__ESR1" %in% colnames(df)) {
  say("   ESR1+ fraction by substage_population (NA=non-mesoderm):")
  print(round(tapply(df[["log__ESR1"]] > 0, addNA(factor(df$substage_population)), mean), 3))
}
say("ER_SURVEY_DONE_OK")
