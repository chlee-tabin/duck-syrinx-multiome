#!/usr/bin/env Rscript
# =====================================================================
# Cross-stage LINEAGE/composition correspondence backbone.
#
# Two transparent, mutually-checking correspondence measures across the full
# 7-stage duck developmental series (multiome HH13, HH17/16, HH21, HH25 ->
# old-scRNA E06, E09, E11):
#
#   (1) MetaNeighbor AUROC  — per-cell, rank-standardized, batch-robust.
#       study_id = stage (7 levels), cell_type = HARMONIZED population.
#       MetaNeighborUS (unsupervised, fast version) over the shared HVGs gives a
#       cell-type x cell-type AUROC matrix; the cross-stage off-diagonal blocks
#       (e.g. HH25-pop x E06-pop) are the correspondence the alluvial ribbons use.
#   (2) Pseudobulk Spearman — transparent cross-check.  Per (stage,population)
#       pseudobulk -> CP10K -> log1p -> Spearman over the shared HVGs.  Adjacent
#       (and all-pairs) population correlations corroborate the AUROC ribbons.
#
# Honest framing: this is TRANSCRIPTIONAL CORRESPONDENCE, NOT lineage tracing.
# Different embryos per stage, no cells tracked.  The alluvial shows which
# populations are transcriptionally nearest stage-to-stage + the compositional
# shift; it does not claim a cell in HH25 became a cell in E06.
#
# Inputs (from the SHARED scratch dir, produced by the two export scripts):
#   multiome per-cell:  mm_<SAFE>_{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz,obs.tsv.gz}
#   multiome pseudobulk: mm_<SAFE>_pb.tsv.gz + mm_<SAFE>_pb_meta.tsv
#   old per-cell:       old_<stage>_{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz,obs.tsv.gz}
#   old pseudobulk:     old_<stage>_pb.tsv.gz + old_<stage>_pb_meta.tsv
#   (SAFE: "HH17/16"->"HH1716"; stages otherwise verbatim)
#
# Outputs (copied back with FLAT names).  Repo
# destination for the three figure-facing tables: figures/data/scrna_e6e9e11/.  The rest
# are not distributed:
#   xstage_auroc_matrix.tsv          full MetaNeighbor AUROC (harmonized celltype x celltype)
#   xstage_corr_ribbons.tsv          tidy adjacent-stage ribbons: from(stage,pop)->to(stage,pop),
#                                      auroc, spearman, n_from, n_to  [the alluvial input]
#   xstage_pb_spearman.tsv           full pseudobulk Spearman matrix (long)
#   xstage_composition.tsv           per (stage,population) cell counts (stratum heights)
#   xstage_gene_universe.txt         the frozen shared gene universe (one symbol/line)
#   xstage_correspondence_provenance.json
#
# Prereq (FLAG FOR USER TO VERIFY on the cluster): MetaNeighbor (Bioconductor)
# + SummarizedExperiment/SingleCellExperiment in /home/cl266/R/4.4.2.  If
# MetaNeighbor is NOT installed, set USE_METANEIGHBOR=FALSE below -> the script
# falls back to pseudobulk-Spearman-only ribbons (still a valid, transparent
# correspondence; the alluvial reads `spearman` for ribbon width and leaves
# auroc=NA).  pyMN (python MetaNeighbor) is the other fallback if R install is
# blocked — but the SCE path here is preferred.
#
# Execution model: run on the O2 cluster via SLURM (modest compute, but reads several
# mtx -> not a login-node job).  R env: module load gcc/14.2.0 R/4.4.2.
# Exploratory only.  Reads the scratch exports; writes ONLY the xstage_link/
# correspondence tables.
# =====================================================================
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressWarnings(suppressMessages({
  .libPaths(c(paste0(HOME_ROOT, "/R/4.4.2"), .libPaths()))
  library(Matrix); library(jsonlite)
}))
set.seed(1)
say <- function(...) cat(sprintf(...), "\n"); t0 <- Sys.time()

OUT <- Sys.getenv("XSTAGE_LINK", paste0(SCRATCH_USER, "/old_duck_review/xstage_link"))  # override for the cr9-CB rerun
USE_METANEIGHBOR <- TRUE          # set FALSE to force the Spearman-only fallback
N_HVG <- 2000                     # shared highly-variable genes for both measures

# Stage order along the developmental axis (left -> right in the alluvial).
STAGE_ORDER <- c("HH13", "HH17/16", "HH21", "HH25", "E06", "E09", "E11")
SAFE <- function(s) gsub("/", "", s)
DATASET <- function(s) if (startsWith(s, "HH")) "multiome" else "old_scrna"
PREFIX  <- function(s) if (startsWith(s, "HH")) sprintf("mm_%s", SAFE(s)) else sprintf("old_%s", s)

# ---------------------------------------------------------------------
# (0) Population-label HARMONIZATION across the 7 stages.
#     The raw vocabularies differ in surface form (old "Chondrogenic LPM" vs
#     multiome "Chondrogenic LPM"; old "ambiguous/mixed LPM" has no multiome
#     analogue; multiome has Early/Cardiopharyngeal LPM absent late).  We map
#     every raw label to a SINGLE canonical stratum vocabulary.  Labels NOT in
#     this map are kept verbatim (so a new label surfaces as its own stratum
#     rather than being silently merged) and flagged in provenance.
#     Strata that exist only at some stages (e.g. Chondrogenic emerges late,
#     Early LPM only early) are EXPECTED — the alluvial must show emergence/loss,
#     not force a match.  Contamination ('excl.') was already dropped upstream.
# ---------------------------------------------------------------------
HARMONIZE <- c(
  # --- shared canonical names (identity maps, listed for completeness) ---
  "Ventral spl."         = "Ventral spl.",
  "Dorsal sub-epi."      = "Dorsal sub-epi.",
  # 2026-08-22 rename: at E09/E11 (~HH35/37) "Chondrogenic LPM" reads wrong —
  # the population is the differentiating syringeal chondrocyte (triad with
  # Syringeal fibroblast / Syringeal smooth muscle). Old-duck arm only; the retired
  # pre-seeded-relabel multiome label and the chick arm's "Chondrogenic LPM" are untouched.
  "Chondrogenic LPM"     = "Syringeal chondrocyte",
  "Chondrogenic"         = "Syringeal chondrocyte", # legacy short form -> canonical
  "Syringeal chondrocyte" = "Syringeal chondrocyte",
  "Syringeal fibroblast"    = "Syringeal fibroblast",     # old-duck compendium (E09/E11)
  "Syringeal smooth muscle" = "Syringeal smooth muscle",  # old-duck compendium (E09/E11)
  "Syringeal"            = "Syringeal",
  "Skeletogenic LPM"     = "Skeletogenic LPM", # canonical identity (current label on both platforms)
  "Somatic"              = "Skeletogenic LPM", # legacy fold: a stale pre-rename mm_* export still
                                               # carrying "Somatic" maps to the canonical stratum
  "BARX1+"               = "BARX1+",
  "Venous-pole spl."     = "Venous-pole spl.",
  "Mesothelium"          = "Mesothelium",
  "Early LPM"            = "Early LPM",
  "Cardiopharyngeal LPM" = "Cardiopharyngeal LPM",
  "Dorsal-biased LPM"    = "Dorsal-biased LPM", # HH17/16 multiome (pre-TES dorsal analogue)
  # --- old-data-only ---
  "ambiguous/mixed LPM"  = "Unresolved LPM"     # rename for a cleaner stratum label
)
harmonize <- function(x) ifelse(x %in% names(HARMONIZE), HARMONIZE[x], x)

# ---------------------------------------------------------------------
# (1) Load each stage's per-cell matrix + obs; reconcile the NKX6-1 alias so
#     the gene namespaces line up before we freeze the shared universe.
#     Multiome native: 'NKX6-1'.  Old cr6.zju1: 'LOC113843492'.  We rename the
#     old feature to 'NKX6-1' so the intersection keeps it (task §3).
# ---------------------------------------------------------------------
ALIAS_OLD_TO_MM <- c("LOC113843492" = "NKX6-1")   # extend if more aliases surface

load_stage <- function(s) {
  pre <- file.path(OUT, PREFIX(s))
  mtx <- sprintf("%s_matrix.mtx.gz", pre)
  if (!file.exists(mtx)) { say("  MISSING %s -> skip stage %s", mtx, s); return(NULL) }
  m   <- as(Matrix::readMM(gzfile(mtx)), "CsparseMatrix")        # features x cells
  fe  <- readLines(gzfile(sprintf("%s_features.tsv.gz", pre)))
  bc  <- readLines(gzfile(sprintf("%s_barcodes.tsv.gz", pre)))
  obs <- read.table(gzfile(sprintf("%s_obs.tsv.gz", pre)), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (DATASET(s) == "old_scrna")
    fe <- ifelse(fe %in% names(ALIAS_OLD_TO_MM), ALIAS_OLD_TO_MM[fe], fe)
  rownames(m) <- make.unique(fe); colnames(m) <- bc
  obs <- obs[match(colnames(m), if ("barcode" %in% names(obs)) obs$barcode else obs$cell_id), ]
  # multiome export writes the population column as `population`; the old export
  # writes it as `label` — coalesce so both sides resolve.
  popcol <- if ("population" %in% names(obs)) obs$population else obs$label
  pop <- harmonize(popcol)
  list(stage = s, mat = m, pop = pop, cell = colnames(m), dataset = DATASET(s))
}

say("== loading per-stage matrices ==")
S <- lapply(STAGE_ORDER, load_stage); S <- S[!sapply(S, is.null)]
present_stages <- sapply(S, `[[`, "stage")
say("stages present: %s", paste(present_stages, collapse = ", "))

# ---------------------------------------------------------------------
# (2) Freeze the shared gene universe IDENTICALLY for every stage: the
#     intersection of all per-stage feature sets (already alias-reconciled).
#     The multiome (crarc_zju1) and old (cr6.zju1) namespaces share ~24.5K
#     symbols; this intersection is the frozen universe used for both measures.
# ---------------------------------------------------------------------
universe <- Reduce(intersect, lapply(S, function(x) rownames(x$mat)))
say("shared gene universe (intersection over %d stages): %d genes", length(S), length(universe))
writeLines(universe, file.path(OUT, "xstage_gene_universe.txt"))

# ---------------------------------------------------------------------
# (3) Pick shared HVGs on the universe (mean-variance, pooled across stages,
#     computed on CP10K-log to be normalization-robust).  Frozen once, reused
#     by BOTH MetaNeighbor and the pseudobulk Spearman so the two measures see
#     the same feature space.
# ---------------------------------------------------------------------
cp10k_log <- function(m) {
  cs <- Matrix::colSums(m); cs[cs == 0] <- 1
  lm <- m %*% Matrix::Diagonal(x = 1e4 / cs)
  lm@x <- log1p(lm@x); lm
}
# pooled HVG by variance of log-CP10K on the universe genes
pooled_var <- numeric(length(universe)); names(pooled_var) <- universe
nseen <- 0L
for (x in S) {
  lm <- cp10k_log(x$mat[universe, , drop = FALSE])
  mu <- Matrix::rowMeans(lm)
  v  <- Matrix::rowMeans(lm^2) - mu^2
  pooled_var <- pooled_var + as.numeric(v); nseen <- nseen + 1L
}
pooled_var <- pooled_var / nseen
hvg <- names(sort(pooled_var, decreasing = TRUE))[seq_len(min(N_HVG, length(universe)))]
say("shared HVGs: %d (pooled log-CP10K variance)", length(hvg))

# ---------------------------------------------------------------------
# (4) Composition table (stratum heights = per (stage,population) cell counts).
# ---------------------------------------------------------------------
comp <- do.call(rbind, lapply(S, function(x)
  data.frame(stage = x$stage, dataset = x$dataset,
             population = names(table(x$pop)),
             n_cells = as.integer(table(x$pop)), stringsAsFactors = FALSE)))
comp$stage <- factor(comp$stage, levels = STAGE_ORDER)
comp <- comp[order(comp$stage), ]
write.table(comp, file.path(OUT, "xstage_composition.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
say("composition table: %d (stage,population) strata", nrow(comp))

# ---------------------------------------------------------------------
# (5a) MetaNeighbor AUROC (per-cell).  Build a combined SCE on the shared HVGs;
#      study_id = stage, cell_type = harmonized population.  MetaNeighborUS
#      (unsupervised) returns a celltype x celltype AUROC matrix labelled
#      "study|celltype"; we parse it back to (stage,population) pairs.
# ---------------------------------------------------------------------
auroc_long <- NULL
mn_ok <- FALSE
if (USE_METANEIGHBOR) {
  mn_ok <- requireNamespace("MetaNeighbor", quietly = TRUE) &&
           requireNamespace("SummarizedExperiment", quietly = TRUE)
  if (!mn_ok) say("WARNING: MetaNeighbor / SummarizedExperiment NOT installed -> Spearman-only fallback")
}
if (USE_METANEIGHBOR && mn_ok) {
  say("== MetaNeighborUS (study_id=stage, cell_type=population, %d HVGs) ==", length(hvg))
  suppressWarnings(suppressMessages({
    library(SummarizedExperiment); library(SingleCellExperiment); library(MetaNeighbor)
  }))
  # cbind all stages on the shared HVGs (RAW counts; MetaNeighbor rank-standardizes internally)
  mats <- lapply(S, function(x) x$mat[hvg, , drop = FALSE])
  big  <- do.call(cbind, mats)
  study <- unlist(lapply(S, function(x) rep(x$stage, ncol(x$mat))))
  ctype <- unlist(lapply(S, function(x) x$pop))
  # globally-unique cell names (stage-prefixed) to avoid collisions across datasets
  colnames(big) <- make.unique(paste(study, unlist(lapply(S, `[[`, "cell")), sep = ":"))
  se <- SummarizedExperiment(assays = list(counts = big),
                             colData = DataFrame(study_id = study, cell_type = ctype))
  vg <- rownames(big)  # already the frozen shared HVGs
  aurocs <- MetaNeighborUS(var_genes = vg, dat = se,
                           study_id = se$study_id, cell_type = se$cell_type,
                           fast_version = TRUE)
  write.table(data.frame(celltype = rownames(aurocs), aurocs, check.names = FALSE),
              file.path(OUT, "xstage_auroc_matrix.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  # parse "study|celltype" -> long (MetaNeighbor uses '|' as the separator)
  parse_id <- function(id) {
    parts <- strsplit(id, "\\|"); st <- sapply(parts, `[`, 1)
    po <- sapply(parts, function(p) paste(p[-1], collapse = "|")); data.frame(stage = st, population = po)
  }
  rn <- parse_id(rownames(aurocs)); cn <- parse_id(colnames(aurocs))
  am <- as.data.frame(as.table(as.matrix(aurocs)))
  colnames(am) <- c("from_id", "to_id", "auroc")
  am <- cbind(am,
              from_stage = rn$stage[match(am$from_id, rownames(aurocs))],
              from_pop   = rn$population[match(am$from_id, rownames(aurocs))],
              to_stage   = cn$stage[match(am$to_id, colnames(aurocs))],
              to_pop     = cn$population[match(am$to_id, colnames(aurocs))])
  auroc_long <- am
  say("MetaNeighbor AUROC matrix: %d x %d cell types", nrow(aurocs), ncol(aurocs))
}

# ---------------------------------------------------------------------
# (5b) Pseudobulk Spearman (transparent cross-check).  Load the per-stage pb
#      tables, restrict to the shared universe, CP10K-log, correlate over HVGs.
# ---------------------------------------------------------------------
say("== pseudobulk Spearman (shared HVGs) ==")
pb_list <- list(); pb_meta <- list()
for (s in present_stages) {
  pre <- file.path(OUT, PREFIX(s))
  pbf <- sprintf("%s_pb.tsv.gz", pre); mf <- sprintf("%s_pb_meta.tsv", pre)
  if (!file.exists(pbf)) { say("  MISSING pb %s", pbf); next }
  pb <- read.table(gzfile(pbf), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  mt <- read.table(mf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (DATASET(s) == "old_scrna")
    rownames(pb) <- ifelse(rownames(pb) %in% names(ALIAS_OLD_TO_MM),
                           ALIAS_OLD_TO_MM[rownames(pb)], rownames(pb))
  colnames(pb) <- paste(s, harmonize(colnames(pb)), sep = "|")
  mt$key <- paste(s, harmonize(mt$population), sep = "|"); mt$stage <- s
  pb_list[[s]] <- pb; pb_meta[[s]] <- mt
}
common_g <- Reduce(intersect, c(list(universe), lapply(pb_list, rownames)))
# unname() the list before cbind: a NAMED list makes do.call(cbind, .) PREFIX every
# column with the list name (-> "HH13.HH13|pop"), which then mis-parses the stage.
pbmat <- do.call(cbind, unname(lapply(pb_list, function(p) p[common_g, , drop = FALSE])))
# CP10K-log per pseudobulk column, then Spearman over HVGs
cs <- colSums(pbmat); cs[cs == 0] <- 1
pbn <- log1p(sweep(pbmat, 2, 1e4 / cs, `*`))
hvg_pb <- intersect(hvg, common_g)
sp_mat <- cor(as.matrix(pbn[hvg_pb, , drop = FALSE]), method = "spearman")
spl <- as.data.frame(as.table(sp_mat)); colnames(spl) <- c("from_id", "to_id", "spearman")
splitid <- function(v) { p <- strsplit(as.character(v), "\\|"); list(s = sapply(p, `[`, 1), pop = sapply(p, function(x) paste(x[-1], collapse = "|"))) }
fi <- splitid(spl$from_id); ti <- splitid(spl$to_id)
spl$from_stage <- fi$s; spl$from_pop <- fi$pop; spl$to_stage <- ti$s; spl$to_pop <- ti$pop
write.table(spl, file.path(OUT, "xstage_pb_spearman.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
say("pseudobulk Spearman: %d x %d pseudobulks over %d HVGs", ncol(sp_mat), ncol(sp_mat), length(hvg_pb))

# ---------------------------------------------------------------------
# (6) ADJACENT-STAGE ribbons — the alluvial input.  For each consecutive stage
#     pair (i -> i+1) along STAGE_ORDER, emit every (from_pop -> to_pop) edge
#     with its AUROC (if MetaNeighbor ran) and Spearman.  The plot decides the
#     width metric + a display threshold; we emit the full bipartite edge set so
#     emergence/loss (a population with no strong partner) is visible.
# ---------------------------------------------------------------------
ncell_of <- setNames(comp$n_cells, paste(comp$stage, comp$population, sep = "|"))
ribbons <- list()
for (i in seq_len(length(present_stages) - 1)) {
  a <- present_stages[i]; b <- present_stages[i + 1]
  pa <- unique(comp$population[comp$stage == a]); pb_ <- unique(comp$population[comp$stage == b])
  for (fp in pa) for (tp in pb_) {
    au <- NA_real_; spv <- NA_real_
    if (!is.null(auroc_long)) {
      hit <- auroc_long[auroc_long$from_stage == a & auroc_long$from_pop == fp &
                          auroc_long$to_stage == b & auroc_long$to_pop == tp, "auroc"]
      if (length(hit)) au <- hit[1]
    }
    hs <- spl[spl$from_stage == a & spl$from_pop == fp & spl$to_stage == b & spl$to_pop == tp, "spearman"]
    if (length(hs)) spv <- hs[1]
    ribbons[[length(ribbons) + 1]] <- data.frame(
      from_stage = a, from_pop = fp, to_stage = b, to_pop = tp,
      auroc = au, spearman = spv,
      n_from = as.integer(ncell_of[paste(a, fp, sep = "|")]),
      n_to   = as.integer(ncell_of[paste(b, tp, sep = "|")]),
      stringsAsFactors = FALSE)
  }
}
rib <- do.call(rbind, ribbons)
write.table(rib, file.path(OUT, "xstage_corr_ribbons.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
say("adjacent-stage ribbons: %d edges across %d transitions",
    nrow(rib), length(present_stages) - 1)

# ---------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------
raw_labels <- sort(unique(unlist(lapply(S, function(x) x$pop))))
unmapped <- setdiff(unlist(lapply(S, function(x) names(table(harmonize(x$pop))))),
                    unname(HARMONIZE))
writeLines(toJSON(list(
  step = "xstage_correspondence",
  stages_present = present_stages, stage_order = STAGE_ORDER,
  n_gene_universe = length(universe), n_hvg = length(hvg),
  metaneighbor_ran = isTRUE(mn_ok && USE_METANEIGHBOR),
  harmonized_vocabulary = sort(unique(unname(HARMONIZE))),
  raw_population_labels = raw_labels,
  date = as.character(Sys.time()),
  note = "Transcriptional correspondence, NOT lineage tracing."),
  pretty = TRUE, auto_unbox = TRUE),
  file.path(OUT, "xstage_correspondence_provenance.json"))
say("\n== xstage_correspondence DONE in %.1f min ==",
    as.numeric(difftime(Sys.time(), t0, "mins")))
cat("SBATCH_DONE_OK\n")
