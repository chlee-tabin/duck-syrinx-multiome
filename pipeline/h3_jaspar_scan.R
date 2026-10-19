#!/usr/bin/env Rscript
# h3_jaspar_scan.R — unbiased JASPAR scan of the PITX2 locus.
#
# Tests whether
# unbiased scanning of the full JASPAR vertebrate motif library recovers
# the curated PITX/SMAD2/3 finding from §5.1 — and whether any non-
# PITX/non-SMAD motif also enriches.
#
# Inputs (env-configurable):
#   G2C_TSV          $SCRATCH/figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv
#   MODE             strict | raw   (cutoff regime, see h2_gkmsvm_train.py)
#   OUT_DIR          $SCRATCH/phase2/atac/jaspar_scan
#   P_CUTOFF         5e-4
#
# Positive set: L-biased peaks (log2FC>0, binom_p<0.05 [+ padj<0.1 if strict]).
# Background set: bilateral peaks (|log2FC|<0.1 & binom_p>0.5), 10x size,
# random sample.
#
# Outputs:
#   $OUT_DIR/h3_enrichment_<MODE>.tsv  one row per PWM:
#     motif_id, motif_name, n_pos_with_hit, n_neg_with_hit, frac_pos,
#     frac_neg, fold_enrich, fisher_p, fisher_padj
#   $OUT_DIR/h3_pwm_meta_<MODE>.tsv   PWM IDs + names (for joining)

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
})

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH  <- Sys.getenv("SCRATCH_ROOT", SCRATCH_ROOT)
G2C_TSV  <- Sys.getenv("G2C_TSV",
    file.path(SCRATCH, "figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv"))
MODE     <- Sys.getenv("MODE", "raw")  # strict | raw
P_CUTOFF <- as.numeric(Sys.getenv("P_CUTOFF", "5e-4"))
NEG_RATIO <- as.integer(Sys.getenv("NEG_RATIO", "10"))
SEED     <- as.integer(Sys.getenv("SEED", "42"))
OUT_DIR  <- Sys.getenv("OUT_DIR",
    file.path(SCRATCH, "phase2/atac/jaspar_scan"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ENRICH <- file.path(OUT_DIR, sprintf("h3_enrichment_%s.tsv", MODE))
OUT_PWM    <- file.path(OUT_DIR, sprintf("h3_pwm_meta_%s.tsv", MODE))

set.seed(SEED)

cat("=== h3_jaspar_scan ===\n")
cat("  G2C_TSV:  ", G2C_TSV, "\n")
cat("  MODE:     ", MODE, "\n")
cat("  P_CUTOFF: ", P_CUTOFF, "\n")
cat("  NEG_RATIO:", NEG_RATIO, "x\n")
cat("  OUT_DIR:  ", OUT_DIR, "\n")

# ── Load g2c TSV + define pos/neg ───────────────────────────────────────
df <- read.delim(G2C_TSV, header = TRUE, stringsAsFactors = FALSE)
cat("Loaded ", nrow(df), " peaks\n", sep = "")

# Filter
if (MODE == "strict") {
    pos_mask <- df$log2FC > 0 & df$binom_p < 0.05 & df$binom_padj < 0.1
} else if (MODE == "raw") {
    pos_mask <- df$log2FC > 0 & df$binom_p < 0.05
} else {
    stop(sprintf("Unknown MODE=%s", MODE))
}
neg_mask <- abs(df$log2FC) < 0.1 & df$binom_p > 0.5

pos_df <- df[pos_mask, ]
neg_pool <- df[neg_mask, ]
cat("  positive set: ", nrow(pos_df), " peaks\n", sep = "")
cat("  bilateral pool: ", nrow(neg_pool), " peaks\n", sep = "")

n_neg <- min(NEG_RATIO * nrow(pos_df), nrow(neg_pool))
neg_df <- neg_pool[sample(nrow(neg_pool), n_neg), ]
cat("  negative set: ", nrow(neg_df), " peaks\n", sep = "")

# ── Build GRanges ────────────────────────────────────────────────────────
parse_peak_id <- function(pids) {
    parts <- do.call(rbind, strsplit(pids, "[:-]"))
    GRanges(seqnames = parts[, 1],
            ranges = IRanges(start = as.integer(parts[, 2]) + 1L,
                             end = as.integer(parts[, 3])))
}
pos_gr <- parse_peak_id(pos_df$peak_id)
neg_gr <- parse_peak_id(neg_df$peak_id)
all_gr <- c(pos_gr, neg_gr)
labels <- c(rep("pos", length(pos_gr)), rep("neg", length(neg_gr)))

# Drop ranges on chromosomes not in the BSgenome (e.g., scaffolds)
bsg <- BSgenome.Anasplatyrhynchos.ZJU1.0
valid_chroms <- seqnames(bsg)
keep <- as.character(seqnames(all_gr)) %in% valid_chroms
n_dropped <- sum(!keep)
if (n_dropped > 0) {
    cat("  dropping ", n_dropped, " peaks on chromosomes not in BSgenome\n", sep = "")
    all_gr <- all_gr[keep]
    labels <- labels[keep]
}
cat("  final: ", length(all_gr), " peaks (", sum(labels == "pos"), " pos / ",
    sum(labels == "neg"), " neg)\n", sep = "")

# ── Load JASPAR2020 vertebrate PWMs ─────────────────────────────────────
cat("Loading JASPAR2020 vertebrate motifs...\n")
opts <- list(collection = "CORE", tax_group = "vertebrates",
             matrixtype = "PWM")
pwms <- TFBSTools::getMatrixSet(JASPAR2020, opts)
cat("  loaded ", length(pwms), " PWMs\n", sep = "")

pwm_meta <- data.frame(
    motif_id = sapply(pwms, function(m) m@ID),
    motif_name = sapply(pwms, function(m) m@name),
    stringsAsFactors = FALSE
)
write.table(pwm_meta, OUT_PWM, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  wrote PWM meta to ", OUT_PWM, "\n", sep = "")

# ── matchMotifs ─────────────────────────────────────────────────────────
cat("Running matchMotifs (p.cutoff=", P_CUTOFF, ")...\n", sep = "")
t0 <- Sys.time()
hits <- matchMotifs(pwms, all_gr, genome = bsg, p.cutoff = P_CUTOFF,
                    out = "matches")
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("  matchMotifs done in ", sprintf("%.1f", elapsed), "s\n", sep = "")

hit_mat <- as.matrix(motifMatches(hits))  # peaks × PWMs, logical
cat("  hit matrix: ", nrow(hit_mat), " peaks × ", ncol(hit_mat), " PWMs\n", sep = "")

# ── Per-PWM enrichment ──────────────────────────────────────────────────
is_pos <- labels == "pos"
n_pos <- sum(is_pos)
n_neg <- sum(!is_pos)

cat("Computing Fisher per PWM...\n")
results <- data.frame(
    motif_id = pwm_meta$motif_id,
    motif_name = pwm_meta$motif_name,
    n_pos_with_hit = colSums(hit_mat[is_pos, , drop = FALSE]),
    n_neg_with_hit = colSums(hit_mat[!is_pos, , drop = FALSE]),
    stringsAsFactors = FALSE
)
results$frac_pos <- results$n_pos_with_hit / n_pos
results$frac_neg <- results$n_neg_with_hit / n_neg
results$fold_enrich <- (results$frac_pos + 1e-6) / (results$frac_neg + 1e-6)

# One-sided Fisher (enrichment in positives)
results$fisher_p <- mapply(function(a, b) {
    fisher.test(matrix(c(a, n_pos - a, b, n_neg - b),
                       nrow = 2, byrow = TRUE),
                alternative = "greater")$p.value
}, results$n_pos_with_hit, results$n_neg_with_hit)
results$fisher_padj <- p.adjust(results$fisher_p, method = "BH")

# Sort: smallest p-value first
results <- results[order(results$fisher_p), ]

write.table(results, OUT_ENRICH, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nTop 20 enriched PWMs:\n")
print(head(results, 20))
cat("\nWrote enrichment to ", OUT_ENRICH, "\n", sep = "")
cat("=== h3_jaspar_scan DONE ===\n")
