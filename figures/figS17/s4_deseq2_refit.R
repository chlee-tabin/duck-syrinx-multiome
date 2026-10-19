#!/usr/bin/env Rscript
# s4_deseq2_refit.R -- E-S1 sensitivity arm S4: local DESeq2 refit with a
# composition covariate on the committed h1 inference-pool pseudobulk counts.
#
# Pools exactly as pipeline/h2_deseq2_snrna.R's POOL scope does (sum h1 counts
# over the stage's inference-pool subclusters per genotype x library_m x side;
# negative counts clipped to 0 and rounded; principled cell floor
# max(20, round(0.05 x median pool n_cells)); genotypes with both sides only;
# zero-count genes dropped; poscounts size factors; parametric dispersion; Wald
# test on side_L_vs_R), then fits two designs on the same matrix:
#   parity      ~ genotype + side                       <- must reproduce the committed
#                                                          h2 POOL row (checked by the caller)
#   composition ~ genotype + syringeal_fraction + side  <- the S4 arm
# syringeal_fraction = Syringeal nuclei / pool nuclei for each half (h1 n_cells).
#
# Usage: Rscript s4_deseq2_refit.R <h1_counts.tsv.gz> <h1_metadata.tsv> <stage>
#                                  <pool populations, "|"-joined> <gene> <out.tsv>
# Called by stats_esr1_lr_hh25.py. Fails loud: no fallback estimator.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L)
  stop("usage: s4_deseq2_refit.R <counts> <metadata> <stage> <pops> <gene> <out.tsv>")
counts_path <- args[1]; meta_path <- args[2]; stage <- args[3]
pops <- strsplit(args[4], "|", fixed = TRUE)[[1]]; gene <- args[5]; out_path <- args[6]
suppressPackageStartupMessages(library(DESeq2))

counts <- as.matrix(read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE))
counts[counts < 0] <- 0
counts <- round(counts)
storage.mode(counts) <- "integer"
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(colnames(counts) %in% meta$sample_id)) stop("count columns missing from the h1 metadata")
if (!gene %in% rownames(counts)) stop(gene, " is not a row of ", counts_path)

sm <- meta[meta$stage_anchor == stage, , drop = FALSE]
if (!nrow(sm)) stop("no h1 samples at stage ", stage)
if (!setequal(unique(sm$subcluster), pops))
  stop("h1 subclusters at ", stage, " {", paste(sort(unique(sm$subcluster)), collapse = ", "),
       "} != the inference pool {", paste(sort(pops), collapse = ", "), "}")
if (!"Syringeal" %in% pops) stop("syringeal_fraction needs Syringeal in the pool")
if (any(!sm$side %in% c("L", "R"))) stop("h1 side outside {L, R} at ", stage)

key <- paste(sm$genotype, sm$library_m, sm$side, sep = "|")
pool <- do.call(rbind, lapply(split(seq_len(nrow(sm)), key), function(i) {
  data.frame(genotype = sm$genotype[i[1]], library_m = sm$library_m[i[1]],
             side = sm$side[i[1]], n_cells = sum(sm$n_cells[i]),
             syr_cells = sum(sm$n_cells[i][sm$subcluster[i] == "Syringeal"]),
             ids = paste(sm$sample_id[i], collapse = "###"), stringsAsFactors = FALSE)
}))
floor_cells <- max(20L, as.integer(round(median(pool$n_cells) * 0.05)))
pool <- pool[pool$n_cells >= floor_cells, , drop = FALSE]
tab <- table(pool$genotype, pool$side)
paired <- rownames(tab)[rowSums(tab > 0) >= 2]
pool <- pool[pool$genotype %in% paired, , drop = FALSE]
if (anyDuplicated(paste(pool$genotype, pool$side)))
  stop("more than one pool sample per genotype x side after flooring")
pool$sample_id <- paste(pool$genotype, pool$library_m, pool$side, "POOL", stage, sep = "|")
cts <- sapply(strsplit(pool$ids, "###", fixed = TRUE),
              function(ids) as.integer(rowSums(counts[, ids, drop = FALSE])))
dimnames(cts) <- list(rownames(counts), pool$sample_id)
cd <- data.frame(genotype = factor(pool$genotype),
                 side = factor(pool$side, levels = c("R", "L")),
                 syringeal_fraction = pool$syr_cells / pool$n_cells,
                 row.names = pool$sample_id)
cat(sprintf("[s4_deseq2_refit] %s: %d pool samples / %d paired genotypes (floor %d cells)\n",
            stage, nrow(cd), length(paired), floor_cells))

fit_one <- function(design, label) {
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = cd, design = design)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, fitType = "parametric", quiet = TRUE)
  res <- results(dds, name = "side_L_vs_R", alpha = 0.05)
  if (!gene %in% rownames(res)) stop(label, ": ", gene, " dropped by the zero-count filter")
  r <- res[gene, ]
  data.frame(model = label, gene = gene, baseMean = r$baseMean,
             log2FoldChange = r$log2FoldChange, lfcSE = r$lfcSE, pvalue = r$pvalue,
             padj = r$padj, n_genes_tested = nrow(dds),
             n_genes_padj_nonNA = sum(!is.na(res$padj)),
             n_samples = nrow(cd), n_paired_genotypes = length(paired),
             floor_cells = floor_cells, deseq2_version = as.character(packageVersion("DESeq2")))
}

out <- rbind(fit_one(~ genotype + side, "parity"),
             fit_one(~ genotype + syringeal_fraction + side, "composition"))
print(out)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
