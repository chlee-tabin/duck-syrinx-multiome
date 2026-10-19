#!/usr/bin/env Rscript
# g4_hh1716_motif_enrichment.R
#
# Phase 3 (per fig04_chromatin_exploratory_plan.md): motif enrichment in
# HH17/16 L>R vs R>L peak sets, with focus on SMAD-family (BMP/NODAL
# transducers) — the load-bearing falsification test for the
# BMP-via-chromatin model.
#
# Approach: motifmatchr scans every peak in the f2 union for matches to
# the JASPAR2020 vertebrate motifs we care about (SMAD-family,
# FOX-family, PITX-family, splanchnic-LPM TFs).  Then per-motif
# Fisher exact test:
#
#   L>R enrichment: (L>R peaks with motif M) vs (background peaks with M)
#   R>L enrichment: (R>L peaks with motif M) vs (background peaks with M)
#
# Background = the full union peak set (all peaks called in any pool).
#
# Outputs:
#   figures/data/g4_hh1716_motif_enrichment.tsv      per-motif p, OR, n_hit, n_total
#   figures/data/g4_hh1716_peak_motif_matrix.tsv     peak × motif binary (large)

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))

.libPaths(c("~/R/4.4.2", .libPaths()))

suppressPackageStartupMessages({
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Anasplatyrhynchos.ZJU1.0)
  library(GenomicRanges)
  library(rtracklayer)
})

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH <- Sys.getenv("SCRATCH_ROOT", SCRATCH_ROOT)
# DA_DIR and PEAKS_BED env-overridable for the ArchR re-run chain
# (pipeline/run_atac_hh1716_chain_archr.sh sets DA_DIR=da_archr and
# points PEAKS_BED at the ArchR union BED).
ATAC_DA <- Sys.getenv("DA_DIR", file.path(SCRATCH, "phase2/atac/da"))
DESEQ_PATH <- file.path(ATAC_DA, "g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv")
PEAKS_BED <- Sys.getenv("PEAKS_BED",
                        file.path(SCRATCH, "phase2/atac/peaks/union_peaks.bed"))

OUT_DIR <- "figures/data"
OUT_ENRICH <- file.path(OUT_DIR, "g4_hh1716_motif_enrichment.tsv")
OUT_MATRIX <- file.path(OUT_DIR, "g4_hh1716_peak_motif_matrix.tsv.gz")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Significance thresholds for L>R / R>L peak set membership
PADJ_CUTOFF <- as.numeric(Sys.getenv("ATAC_PADJ_CUTOFF", "0.10"))
LFC_CUTOFF <- as.numeric(Sys.getenv("ATAC_LFC_CUTOFF", "0"))

# Same motif library as g3_pitx2_motif_scan.R for consistency
TARGET_MOTIFS <- c(
  SMAD2_3_4_human = "MA0513.1",
  SMAD3           = "MA0795.1",
  SMAD4           = "MA1153.1",
  SMAD5           = "MA1557.1",
  Smad2_3_mouse   = "MA1622.1",
  FOXH1           = "MA0479.1",
  FOXA2           = "MA0047.3",
  Foxa2_mouse     = "MA0047.2",
  FOXC1           = "MA0032.2",
  FOXC2           = "MA0846.1",
  Foxf1_mouse     = "MA1606.1",
  PITX1           = "MA0682.2",
  PITX2           = "MA1547.1",
  PITX3           = "MA0714.1",
  NKX2_5          = "MA0063.2",
  NKX6_1          = "MA0674.1",
  TBX5            = "MA0807.1",
  GATA4           = "MA0482.2",
  ISL1            = "MA1608.1",
  HAND1_TCF3      = "MA0092.1",
  MEF2C           = "MA0497.1"
)

cat("=== g4_hh1716_motif_enrichment ===\n")
cat("DESeq2 :", DESEQ_PATH, "\n")
cat("Peaks  :", PEAKS_BED, "\n")
cat("padj cutoff:", PADJ_CUTOFF, "\n")
cat("LFC cutoff :", LFC_CUTOFF, "\n\n")

# ── Pre-flight ──────────────────────────────────────────────────────────
if (!file.exists(DESEQ_PATH)) {
  stop("DESeq2 result missing: ", DESEQ_PATH, " — run g2_lpm_hh1716_deseq2.py first")
}
if (!file.exists(PEAKS_BED)) {
  stop("Union peaks missing: ", PEAKS_BED, " — run f2 first")
}

# ── Load DESeq2 results + classify L>R / R>L peak sets ─────────────────
# Pandas writes an index with empty header; R's read.table interprets that
# as the first numeric column being unnamed — often returned as "V1" by
# read.table with check.names=FALSE.  Also handle "Unnamed: 0" (pandas
# read-back convention).  Fall back to first column.
deseq <- read.table(DESEQ_PATH, header = TRUE, sep = "\t", check.names = FALSE)
cat("DESeq2 columns:", paste(colnames(deseq), collapse = ", "), "\n")
# Resolve peak_id column: try named candidates, else assume first column is peak ids
candidates <- c("peak_id", "Unnamed: 0", "V1", "X")
peak_col <- NA_character_
for (c in candidates) {
  if (c %in% colnames(deseq)) { peak_col <- c; break }
}
if (is.na(peak_col)) peak_col <- colnames(deseq)[1]
cat("Using peak_col=", peak_col, "\n", sep = "")

deseq <- deseq[, c(peak_col, "log2FoldChange", "padj")]
colnames(deseq) <- c("peak_id", "log2FoldChange", "padj")
deseq$padj[is.na(deseq$padj)] <- 1

l_gt_r <- deseq$peak_id[deseq$padj < PADJ_CUTOFF & deseq$log2FoldChange > LFC_CUTOFF]
r_gt_l <- deseq$peak_id[deseq$padj < PADJ_CUTOFF & deseq$log2FoldChange < -LFC_CUTOFF]
cat("L>R peaks (padj<", PADJ_CUTOFF, " + log2FC>", LFC_CUTOFF, "):", length(l_gt_r), "\n", sep = "")
cat("R>L peaks (padj<", PADJ_CUTOFF, " + log2FC<-", LFC_CUTOFF, "):", length(r_gt_l), "\n", sep = "")

if (length(l_gt_r) + length(r_gt_l) == 0) {
  cat("WARNING: no significant L vs R peaks at padj<", PADJ_CUTOFF,
      ".  Wave-2 motif test is vacuous in this state.\n", sep = "")
  enrich_df <- data.frame(motif_label = character(),
                           direction = character(),
                           OR = double(), p = double(),
                           n_hit_set = integer(),
                           n_set = integer(),
                           n_hit_bg = integer(),
                           n_bg = integer(),
                           atac_status = character())
  write.table(enrich_df, OUT_ENRICH, sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat("Wrote empty enrichment table -> ", OUT_ENRICH, "\n")
  quit(save = "no", status = 0)
}

# ── Load union peaks as GRanges ─────────────────────────────────────────
peaks <- import(PEAKS_BED, format = "BED")
cat("Loaded", length(peaks), "union peaks from", PEAKS_BED, "\n")

# Build a label vector for joining DESeq2 peak_id -> peak index.
# Peak IDs come from SnapATAC2's internal naming which preserves the
# original BED 0-based-half-open coords (chrom:start_0based-end_0based).
# rtracklayer::import returns 1-based-inclusive start, so we have to
# subtract 1 from start when constructing labels.  Build BOTH 0-based
# and 1-based variants + "chrom-start-end" alt to cover all conventions.
peak_labels_0b <- paste0(as.character(seqnames(peaks)), ":",
                          start(peaks) - 1L, "-", end(peaks))
peak_labels_1b <- paste0(as.character(seqnames(peaks)), ":",
                          start(peaks), "-", end(peaks))
peak_idx_0b <- setNames(seq_along(peaks), peak_labels_0b)
peak_idx_1b <- setNames(seq_along(peaks), peak_labels_1b)

alt_labels_0b <- paste0(as.character(seqnames(peaks)), "-",
                         start(peaks) - 1L, "-", end(peaks))
alt_labels_1b <- paste0(as.character(seqnames(peaks)), "-",
                         start(peaks), "-", end(peaks))
peak_idx_alt_0b <- setNames(seq_along(peaks), alt_labels_0b)
peak_idx_alt_1b <- setNames(seq_along(peaks), alt_labels_1b)

peak_idx_combined <- c(peak_idx_0b, peak_idx_1b, peak_idx_alt_0b, peak_idx_alt_1b)

resolve <- function(ids) {
  out <- peak_idx_combined[ids]
  out[!is.na(out)]
}

l_gt_r_idx <- resolve(l_gt_r)
r_gt_l_idx <- resolve(r_gt_l)
cat("Resolved L>R peaks to union: ", length(l_gt_r_idx), " / ", length(l_gt_r), "\n", sep = "")
cat("Resolved R>L peaks to union: ", length(r_gt_l_idx), " / ", length(r_gt_l), "\n", sep = "")

# ── Load target motif PWMs (getMatrixByID — robust) ─────────────────────
target_pfm_list <- list()
for (label in names(TARGET_MOTIFS)) {
  mid <- TARGET_MOTIFS[[label]]
  pfm <- getMatrixByID(JASPAR2020, ID = mid)   # a missing motif is an error, never a skip
  if (is.null(pfm)) stop("motif ", label, " (", mid, ") not in JASPAR2020")
  target_pfm_list[[label]] <- pfm
}
if (length(target_pfm_list) == 0L) stop("No motifs resolved.")

target_pwms <- do.call(PWMatrixList,
                       lapply(target_pfm_list, toPWM, pseudocounts = 0.8))
names(target_pwms) <- names(target_pfm_list)
cat("Resolved", length(target_pwms), "PWMs.\n")

# ── motifmatchr scan on full peak set ───────────────────────────────────
cat("Scanning", length(peaks), "peaks × ", length(target_pwms), "motifs...\n")
match_obj <- matchMotifs(target_pwms, peaks,
                         genome = BSgenome.Anasplatyrhynchos.ZJU1.0,
                         p.cutoff = 5e-4,
                         out = "matches")
# Binary peak × motif matrix
peak_motif <- as.matrix(motifMatches(match_obj))
mode(peak_motif) <- "integer"

# Save the matrix (gzipped TSV)
con <- gzfile(OUT_MATRIX, "w")
write.table(cbind(peak_id = peak_labels_0b, peak_motif),
            con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat("Wrote peak × motif matrix -> ", OUT_MATRIX, "\n")

# ── Fisher enrichment per motif × direction ─────────────────────────────
n_total <- length(peaks)
l_set_size <- length(l_gt_r_idx)
r_set_size <- length(r_gt_l_idx)

enrich_rows <- list()
for (motif_label in colnames(peak_motif)) {
  hits <- peak_motif[, motif_label]
  for (direction in c("L_gt_R", "R_gt_L")) {
    set_idx <- if (direction == "L_gt_R") l_gt_r_idx else r_gt_l_idx
    n_set <- length(set_idx)
    if (n_set == 0L) next
    n_hit_set <- sum(hits[set_idx])
    n_hit_bg <- sum(hits) - n_hit_set
    n_bg <- n_total - n_set
    mat <- matrix(c(n_hit_set,         n_set - n_hit_set,
                    n_hit_bg,          n_bg - n_hit_bg),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("set", "bg"), c("hit", "no_hit")))
    test <- fisher.test(mat, alternative = "greater")
    enrich_rows[[length(enrich_rows) + 1L]] <- data.frame(
      motif_label = motif_label,
      motif_id    = TARGET_MOTIFS[[motif_label]],
      direction   = direction,
      n_hit_set   = n_hit_set,
      n_set       = n_set,
      n_hit_bg    = n_hit_bg,
      n_bg        = n_bg,
      OR          = unname(test$estimate),
      p           = unname(test$p.value),
      stringsAsFactors = FALSE
    )
  }
}

enrich_df <- do.call(rbind, enrich_rows)
# BH adjust within direction
for (d in unique(enrich_df$direction)) {
  mask <- enrich_df$direction == d
  enrich_df$padj[mask] <- p.adjust(enrich_df$p[mask], method = "BH")
}
enrich_df <- enrich_df[order(enrich_df$direction, enrich_df$p), ]
write.table(enrich_df, OUT_ENRICH, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("Wrote enrichment table -> ", OUT_ENRICH, "\n")

cat("\n=== Top motifs by direction ===\n")
for (d in unique(enrich_df$direction)) {
  cat("\n[", d, "]\n", sep = "")
  print(head(enrich_df[enrich_df$direction == d,
                       c("motif_label", "n_hit_set", "n_set",
                         "OR", "p", "padj")], 10))
}

cat("\nDONE.\n")
