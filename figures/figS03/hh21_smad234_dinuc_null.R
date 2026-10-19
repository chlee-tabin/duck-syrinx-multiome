#!/usr/bin/env Rscript
# hh21_smad234_dinuc_null.R — k-MER-PRESERVING shuffle null for the HH21 promoter panel.
# SHUFFLE_K selects the order: 1=mono, 2=dinucleotide (default), 3=trinucleotide.
#
# THE QUESTION. SMAD2::SMAD3::SMAD4 (JASPAR MA0513.1) is the only motif that repeatedly surfaces in
# the L-up arm of the HH21 promoter scan (+0.204, BH q=0.055) without crossing. A tail analysis
# (hh21_trans_motif_seq.py `tail`) shows the excess lives entirely in motif-dense promoters, which
# is what promoter COMPOSITION looks like rather than a pathway program. This settles it: if
# dinucleotide composition alone reproduces the L-up coefficient, the signal is composition.
#
# WHY DINUCLEOTIDE (k=2) IS THE DEFAULT, NOT MONONUCLEOTIDE. The GLM already carries a scalar GC covariate, so a
# mononucleotide null adds nothing. Dinucleotide structure (CpG depletion, purine tracts, simple
# repeats) is exactly what a scalar GC term cannot absorb and exactly what inflates PWM hit counts
# in low-complexity promoters. universalmotif's k=2 euler method is Altschul-Erickson: it preserves
# each sequence's dinucleotide composition EXACTLY, not merely in expectation.
#
# SCOPE. Two motifs only: SMAD2_3_4 (the question) and PITX2 (anchor — its observed L-up effect is
# ~0, so its null must also be ~0; if PITX2's null is off-zero the whole null is miscalibrated).
# The bilateral reference is SUBSAMPLED (default 3000) because the SE is dominated by the n=96 L-up
# arm, and the full 10,320 would make the shuffle scan hours rather than minutes. The local step
# recomputes the OBSERVED coefficient on the SAME subsample, so observed and null are like-for-like.
#
# Out: figures/data/hh21_smad234_{mono,dinuc,trinuc}_null.tsv  (per SHUFFLE_K)
#      gene_id, group, gc, rna_baseMean, obs_<motif>, null_<motif>_1 .. _K
#
# Run: sbatch pipeline/run_hh21_smad234_dinuc_null.sh
.libPaths(c("~/R/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(motifmatchr); library(TFBSTools); library(JASPAR2020)
  library(BSgenome.Anasplatyrhynchos.ZJU1.0)
  library(Biostrings); library(GenomicRanges); library(universalmotif)
})

DATA <- Sys.getenv("DATA_DIR", "figures/data")
BED <- file.path(DATA, "hh21_seq1_promoters.bed")          # the definitive single-TSS promoters
MAP <- file.path(DATA, "hh21_seq_promoter_map.tsv")
KORD <- as.integer(Sys.getenv("SHUFFLE_K", "2"))   # 1 = mono, 2 = dinucleotide, 3 = trinucleotide
if (!grepl("^[0-9]+$", Sys.getenv("SHUFFLE_K", "2")) || is.na(KORD) || !(KORD %in% 1:3))
  stop(sprintf("SHUFFLE_K must be 1, 2 or 3 (got %s) — an out-of-range value would run the\n  full K-replicate loop and then write an NA-named file", Sys.getenv("SHUFFLE_K")))
OUT <- file.path(DATA, sprintf("hh21_smad234_%s_null.tsv", c("mono","dinuc","trinuc")[KORD]))
K <- as.integer(Sys.getenv("N_SHUFFLE", "100"))
N_BIL <- as.integer(Sys.getenv("N_BILATERAL", "3000"))
PCUT <- as.numeric(Sys.getenv("MOTIF_P_CUTOFF", "1e-3"))
SEED <- as.integer(Sys.getenv("SEED", "1"))
BG <- c(A = 0.293, C = 0.207, G = 0.207, T = 0.293)        # SAME fixed bg as the main scan
set.seed(SEED)

# 5th column is prep()'s provenance stamp (see hh21_trans_motif_scan.R); tolerate its absence so
# an older BED still runs, but carry it into the null table when present.
bed <- read.table(BED, sep = "\t", stringsAsFactors = FALSE)
if (ncol(bed) == 5) {
  colnames(bed) <- c("chrom", "start", "end", "gene_id", "prep_stamp")
} else if (ncol(bed) == 4) {
  colnames(bed) <- c("chrom", "start", "end", "gene_id")
  bed$prep_stamp <- NA_character_
} else {
  stop(sprintf("null: %s has %d columns; expected 4 or 5", BED, ncol(bed)))
}
map <- read.table(MAP, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
n_bed0 <- nrow(bed)
bed <- merge(bed, map[, c("gene_id", "group", "rna_baseMean")], by = "gene_id")
# inner join: a BED/map pair from different preps would silently shrink the reference set,
# and the post-subsample log line prints 3000 either way.
if (nrow(bed) != n_bed0)
  stop(sprintf("BED/map merge dropped %d of %d genes — the two came from different preps",
               n_bed0 - nrow(bed), n_bed0))
cat("null: pre-subsample groups — ", paste(names(table(bed$group)), table(bed$group),
    sep = "=", collapse = "  "), "\n", sep = "")

keep <- c(which(bed$group != "bilateral"),
          sample(which(bed$group == "bilateral"), min(N_BIL, sum(bed$group == "bilateral"))))
bed <- bed[sort(keep), ]
cat("null: k=", KORD, " shuffle; ", nrow(bed), " promoters (", sum(bed$group == "L-up"), " L-up, ",
    sum(bed$group == "R-up"), " R-up, ", sum(bed$group == "bilateral"), " bilateral); K=", K, "\n", sep = "")

gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end))
seq <- getSeq(BSgenome.Anasplatyrhynchos.ZJU1.0, gr)

TARGETS <- c(SMAD2_3_4 = "MA0513.1", PITX2 = "MA1547.1")
pfm <- lapply(TARGETS, function(id) getMatrixByID(JASPAR2020, ID = id))
pwms <- do.call(PWMatrixList, lapply(pfm, toPWM, pseudocounts = 0.8))
names(pwms) <- names(pfm)

count_hits <- function(s) {
  posl <- matchMotifs(pwms, s, bg = BG, out = "positions", p.cutoff = PCUT)
  m <- vapply(posl, function(x) as.integer(lengths(x)), integer(length(s)))
  colnames(m) <- names(posl); m
}

obs <- count_hits(seq)
cat("null: observed mean sites/promoter — ",
    paste(sprintf("%s=%.3f", colnames(obs), colMeans(obs)), collapse = "  "), "\n", sep = "")

# integrity check: the shuffle must preserve dinucleotide composition EXACTLY
chk <- seq[1:min(50, length(seq))]
chk_s <- shuffle_sequences(chk, k = KORD, method = "euler", rng.seed = SEED)
d0 <- oligonucleotideFrequency(chk, width = KORD)
d1 <- oligonucleotideFrequency(chk_s, width = KORD)
if (!identical(as.vector(d0), as.vector(d1)))
  stop(sprintf("shuffle did NOT preserve %d-mer composition — null is invalid", KORD))
cat("null: ", KORD, "-mer composition preserved exactly on a 50-sequence check\n", sep = "")

nullmat <- vector("list", length(TARGETS))
names(nullmat) <- names(TARGETS)
for (nm in names(TARGETS)) nullmat[[nm]] <- matrix(0L, length(seq), K)
for (k in seq_len(K)) {
  sh <- shuffle_sequences(seq, k = KORD, method = "euler", rng.seed = SEED + k)
  cnt <- count_hits(sh)
  for (nm in names(TARGETS)) nullmat[[nm]][, k] <- cnt[, nm]
  if (k %% 10 == 0) cat("null: replicate ", k, "/", K, "\n", sep = "")
}

out <- data.frame(gene_id = bed$gene_id, group = bed$group,
                  gc = round(as.numeric(letterFrequency(seq, "GC", as.prob = TRUE)), 4),
                  rna_baseMean = bed$rna_baseMean, stringsAsFactors = FALSE)
for (nm in names(TARGETS)) out[[paste0("obs_", nm)]] <- obs[, nm]
for (nm in names(TARGETS)) {
  colnames(nullmat[[nm]]) <- paste0("null_", nm, "_", seq_len(K))
  out <- cbind(out, nullmat[[nm]])
}
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("null: wrote ", OUT, " (", nrow(out), " promoters x ", K, " replicates)\n", sep = "")
for (nm in names(TARGETS))
  cat("null: ", nm, " observed mean ", sprintf("%.3f", mean(obs[, nm])),
      " vs shuffled mean ", sprintf("%.3f", mean(nullmat[[nm]])), "\n", sep = "")
cat("=== DONE hh21_smad234_dinuc_null ===\n")
