#!/usr/bin/env Rscript
# cr9_screen_export.R
# Re-export the FULL manuscript S8/S9/S10 gene panels (WNT ligands/modulators,
# BMP-vs-TGFb SMADs, BMP ligands) PER-CELL from the cr9 substage objects, carrying
# stage / side / sex / embryo / clustering so the screen can be resolved by genuine
# D/V subpopulation x stage (x L/R) instead of the manuscript's "pooled mesoderm".
#
# @data on the substage objects = log1p(CP10K)  ->  CP10K = expm1(@data).
# We export BOTH log-normalized (for violins, matches manuscript "Expression Level")
# and linear CP10K (for the per-cell "expression boost" L/R comparison).
#
#   sbatch run_cr9_screen_export.sh   # writes old_duck_cr9/substage/cr9_screen_*_percell.tsv.gz
# Repo destination (copied back from the cluster): figures/data/scrna_e6e9e11/;
# figS08/figS09/figS10 read these via figures/figS08/cr9_screen_io.R.
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")

suppressPackageStartupMessages({library(Seurat)})
say <- function(...) cat(sprintf(...), "\n")

SUB <- paste0(SCRATCH_USER, "/old_duck_cr9/substage"); OUT <- SUB

WNT  <- c("WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6",
          "WNT7A","WNT7B","WNT9A","WNT9B","WNT11","WNT16","SFRP2","SHISA2","WIF1")   # Fig S8
SMAD <- c("SMAD1","SMAD5","SMAD9","SMAD2","SMAD3")                                    # Fig S9
BMP  <- c("BMP1","BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP15")           # Fig S10
EXTRA<- c("GREM1","GREM2","NOG","CHRD","BMPR1A","BMPR1B","BMPR2",
          "ESR1","ESR2","ESRRG","PITX2")   # antagonists/receptors + refs for context
GENES <- unique(c(WNT, SMAD, BMP, EXTRA))
MK <- c("donor_id","embryo_id","sex.called","sex_numeric","breed","LR.called",
        "library","nCount_RNA")

# Compendium population labels — joined HERE (on the cluster) against the SCRATCH
# run-matched clusterlabels, so the export is self-contained and run-consistent with
# these exact objects' RNA_snn_res.0.2 (no stale-leiden-ID mismatch). That run-matched
# file is committed for provenance as
# figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv, and is a DIFFERENT
# substage run from the earlier cr6-era substage_compendium_clusterlabels.tsv (no cr9_
# prefix; 10,160 vs 8,736 E09 cells; not distributed) — do NOT re-join these exports to it.
LAB <- read.delim(file.path(SUB, "substage_compendium_clusterlabels.tsv"),
                  stringsAsFactors = FALSE)

for (s in c("E06","E09","E11")) {
  f <- file.path(SUB, sprintf("substage_%s_obj.rds", s))
  obj <- readRDS(f); DefaultAssay(obj) <- "RNA"
  dat <- GetAssayData(obj, assay="RNA", layer="data")
  present <- intersect(GENES, rownames(dat)); missing <- setdiff(GENES, rownames(dat))
  say("[%s] cells=%d  present=%d/%d  MISSING={%s}",
      s, ncol(obj), length(present), length(GENES), paste(missing, collapse=","))
  logm  <- as.matrix(dat[present, , drop=FALSE])
  cp10k <- expm1(logm)                                  # @data=log1p(CP10K)
  md <- obj@meta.data
  mk       <- intersect(MK, colnames(md))
  res_cols <- intersect(c("RNA_snn_res.0.2","RNA_snn_res.0.3"), colnames(md))
  labs <- LAB[LAB$stage == s, ]
  pmap <- setNames(labs$population, as.character(labs$cluster))
  population <- pmap[as.character(md[["RNA_snn_res.0.2"]])]
  if (any(is.na(population)))
    stop(sprintf("[%s] unmapped res0.2 clusters vs clusterlabels: {%s}", s,
                 paste(sort(unique(md[["RNA_snn_res.0.2"]][is.na(population)])), collapse=",")))
  say("[%s] population: %s", s, paste(sprintf("%s=%d", names(table(population)), table(population)), collapse=" "))
  base <- data.frame(cell_id = colnames(obj), stage = s, population = unname(population),
                     md[, mk, drop=FALSE], md[, res_cols, drop=FALSE],
                     check.names = FALSE, stringsAsFactors = FALSE)
  logdf <- as.data.frame(t(logm));  colnames(logdf) <- paste0("log__",   colnames(logdf))
  cpdf  <- as.data.frame(t(cp10k)); colnames(cpdf)  <- paste0("cp10k__", colnames(cpdf))
  out <- cbind(base, logdf, cpdf)
  gz <- gzfile(file.path(OUT, sprintf("cr9_screen_%s_percell.tsv.gz", s)), "w")
  write.table(out, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
  say("[%s] wrote %d cells x %d cols -> cr9_screen_%s_percell.tsv.gz", s, nrow(out), ncol(out), s)
}
say("DONE")
