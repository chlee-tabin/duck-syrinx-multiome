#!/usr/bin/env Rscript
# cr9_syrinx_lr_deseq2.R — paired L-vs-R DESeq2 per (stage x sex) on the
# old-data syrinx-population pseudobulks (E09/E11), for the expanded S16 C
# forest. Genome-wide fit (all genes), BH padj; design ~ donor + side.
# LOCAL step (duckfig2 R): run from repo root after the count export lands:
#   conda run -n duckfig2 Rscript --vanilla pipeline/scrna_e6e9e11/cr9_syrinx_lr_deseq2.R
# Reads  figures/data/scrna_e6e9e11/cr9_syrinx_pseudobulk_{counts,meta}_{E09,E11}.tsv*
# Writes figures/data/scrna_e6e9e11/cr9_syrinx_lr_deseq2_{E09,E11}_{M,F}.tsv
suppressPackageStartupMessages(library(DESeq2))
say <- function(...) cat(sprintf(...), "\n")
DIR <- "figures/data/scrna_e6e9e11"
if (!dir.exists(DIR)) stop("run from repo root: ", DIR, " not found")

for (s in c("E09", "E11")) {
  counts <- read.delim(gzfile(file.path(DIR,
    sprintf("cr9_syrinx_pseudobulk_counts_%s.tsv.gz", s))), check.names = FALSE)
  meta <- read.delim(file.path(DIR,
    sprintf("cr9_syrinx_pseudobulk_meta_%s.tsv", s)),
    colClasses = "character")  # identifier columns stay character (an all-"F" sex column otherwise type-converts to logical)
  meta$n_cells <- as.integer(meta$n_cells)
  m <- as.matrix(counts[, -1]); rownames(m) <- counts$gene
  if (!identical(colnames(m), meta$column))
    stop(sprintf("[%s] counts/meta column mismatch", s))
  for (sx in c("M", "F")) {
    ms <- meta[meta$sex == sx, ]
    # paired = exactly one L and one R pseudobulk for the donor
    tab <- table(ms$donor_id, ms$side)
    paired <- rownames(tab)[tab[, "L"] == 1 & tab[, "R"] == 1]
    dropped <- setdiff(unique(ms$donor_id), paired)
    if (length(dropped))
      say("[%s %s] unpaired donors dropped: %s", s, sx,
          paste(dropped, collapse = ","))
    ms <- ms[ms$donor_id %in% paired, ]
    if (length(paired) < 2)
      stop(sprintf("[%s %s] only %d paired donor(s); need >= 2", s, sx,
                   length(paired)))
    sub <- m[, ms$column, drop = FALSE]
    cd <- data.frame(donor = factor(ms$donor_id),
                     side = factor(ms$side, levels = c("R", "L")))
    dds <- DESeqDataSetFromMatrix(sub, cd, ~ donor + side)
    dds <- DESeq(dds, quiet = TRUE)
    # alpha must match the downstream padj<0.05 threshold (h2 convention;
    # DESeq2 optimizes independent filtering against it)
    r <- results(dds, contrast = c("side", "L", "R"), alpha = 0.05)
    out <- data.frame(gene_id = rownames(r), as.data.frame(r),
                      n_paired_donors = length(paired), row.names = NULL)
    f <- file.path(DIR, sprintf("cr9_syrinx_lr_deseq2_%s_%s.tsv", s, sx))
    write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
    say("[%s %s] pairs=%d genes=%d padj<0.05=%d -> %s", s, sx, length(paired),
        nrow(out), sum(out$padj < 0.05, na.rm = TRUE), f)
  }
}
say("DONE")
