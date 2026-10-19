#!/usr/bin/env Rscript
# g3f_hh21_motif_scan.R — DUCK motif presence + INSTANCE COUNTS at the g4f floor-20 HH21 windows.
#
# Re-scope of g3d to the consistent (floor-20, same engine as panel A) HH21 ventral L-biased
# cis-block + R-desert Class III + a/b promoters + ns controls (WIN_BED = g4f_hh21_motif_windows.bed,
# built by make_g4f_hh21_windows.py).  Adds the per-window motif INSTANCE COUNT (matchMotifs
# out="positions") alongside binary presence, on the SAME fixed-bg JASPAR panel as g3c/g3d.
#
# GC guardrail unchanged: every motif tagged GC-immune/GC-rich; raw AGAC SBE (GC-immune) kept;
# counts at the 81%-GC c-promoter are interpreted only for GC-immune motifs.
#
# Output: figures/data/g3f_hh21_duck_motifs.tsv (long: window,class,dist,gc,motif,present,n_sites,
#   motif_gc,gc_class) + g3f_hh21_window_gc.tsv (per-window gc + raw SBE AGAC count).
#
# Run: sbatch pipeline/run_g3f_g3g_domain.sh
.libPaths(c("~/R/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(motifmatchr); library(TFBSTools); library(JASPAR2020)
  library(BSgenome.Anasplatyrhynchos.ZJU1.0)
  library(Biostrings); library(GenomicRanges)
})
OUT <- "figures/data"
WIN_BED <- Sys.getenv("WIN_BED", "")
if (!nzchar(WIN_BED)) stop("WIN_BED must be set (pair it with the matching OUT_TAG). Runners: pipeline/run_g3f_g3g_domain.sh (shipped; committed figures/data/g4f_domain_motif_windows.bed)")
OUT_TAG <- Sys.getenv("OUT_TAG", "hh21")   # "domain" for the multi-stage fig04 panel C re-scope
PCUT <- as.numeric(Sys.getenv("MOTIF_P_CUTOFF", "1e-3"))
TSS_C <- 3971385L
set.seed(1)
duck <- BSgenome.Anasplatyrhynchos.ZJU1.0
BG <- c(A = 0.293, C = 0.207, G = 0.207, T = 0.293)  # SAME fixed bg as g3c/g3d/g3e (comparable)

bed <- read.table(WIN_BED, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(bed)[1:4] <- c("chrom", "start", "end", "name")
parts <- do.call(rbind, strsplit(bed$name, "\\|")); bed$label <- parts[, 1]; bed$class <- parts[, 2]
bed$dist <- as.integer((bed$start + bed$end) / 2 - TSS_C)
gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end)); seq <- getSeq(duck, gr); names(seq) <- bed$label
cat("g3f: ", length(seq), " duck windows\n", sep = ""); print(table(bed$class))

gc_frac <- function(s) as.numeric(letterFrequency(s, "GC", as.prob = TRUE))
sbe_raw <- sapply(seq, function(s) countPattern("AGAC", s) + countPattern("GTCT", s))
win_gc <- data.frame(window = bed$label, class = bed$class, dist = bed$dist,
                     gc = round(sapply(seq, gc_frac), 3), sbe_AGAC_count = as.integer(sbe_raw))
write.table(win_gc, file.path(OUT, sprintf("g3f_%s_window_gc.tsv", OUT_TAG)), sep = "\t", quote = FALSE, row.names = FALSE)

TARGETS <- c(SMAD2_3_4 = "MA0513.1", SMAD3 = "MA0795.1", SMAD4 = "MA1153.1", SMAD5 = "MA1557.1",
             FOXH1 = "MA0479.1", NKX6_1 = "MA0674.1", Foxf1 = "MA1606.1", FOXC1 = "MA0032.2",
             GATA4 = "MA0482.2", ISL1 = "MA1608.1", TBX5 = "MA0807.1", PITX2 = "MA1547.1")
pfm <- list()
for (lab in names(TARGETS)) {
  m <- getMatrixByID(JASPAR2020, ID = TARGETS[[lab]])   # a missing motif is an error
  if (is.null(m)) stop("g3f: motif ", lab, " (", TARGETS[[lab]], ") not in JASPAR2020")
  pfm[[lab]] <- m
}
# other BMP R-SMADs SMAD1 / SMAD9 by name (JASPAR2020 may lack standalone matrices — they share
# the SBE with SMAD5; report definitively + scan if present)
for (nm in c("SMAD1", "SMAD9")) {
  s <- tryCatch(getMatrixSet(JASPAR2020, list(name = nm)), error = function(e) NULL)
  if (!is.null(s) && length(s) > 0) { pfm[[nm]] <- s[[1]]; cat("g3f: added", nm, "=", ID(s[[1]]), "\n") }
  else cat("g3f: JASPAR2020 has NO standalone", nm, "matrix (BMP arm represented by SMAD5 + SBE)\n")
}
pwms <- do.call(PWMatrixList, lapply(pfm, toPWM, pseudocounts = 0.8)); names(pwms) <- names(pfm)
motif_gc <- sapply(pfm, function(m) { M <- m@profileMatrix; sum(M[c("C", "G"), ]) / sum(M) })
gc_class <- ifelse(motif_gc >= 0.60, "GC-rich", "GC-immune")

# export the JASPAR2020 PFMs (same matrices used in the scan) for downstream logo rendering
pfm_long <- do.call(rbind, lapply(names(pfm), function(nm) {
  M <- as.matrix(pfm[[nm]]@profileMatrix)            # 4 x L, rows A,C,G,T
  data.frame(motif = nm, pos = rep(seq_len(ncol(M)), each = 4),
             base = rep(rownames(M), times = ncol(M)), count = as.numeric(M), stringsAsFactors = FALSE)
}))
write.table(pfm_long, file.path(OUT, "g3f_motif_pfms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# presence (binary)
mm <- as.matrix(motifMatches(matchMotifs(pwms, seq, bg = BG, out = "scores", p.cutoff = PCUT)))
rownames(mm) <- bed$label
# INSTANCE COUNTS (out="positions": GRangesList per motif; seqnames = window label or index)
posl <- matchMotifs(pwms, seq, bg = BG, out = "positions", p.cutoff = PCUT)
cnt <- matrix(0L, length(seq), length(posl), dimnames = list(bed$label, names(posl)))
for (mi in seq_along(posl)) {            # posl[[mi]] = IRangesList, one (poss. empty) IRanges per window
  el <- lengths(posl[[mi]]); nm <- names(el)
  if (!is.null(nm) && all(nm %in% bed$label)) cnt[nm, names(posl)[mi]] <- as.integer(el)
  else if (length(el) == nrow(cnt)) cnt[, names(posl)[mi]] <- as.integer(el)
}
out <- data.frame()
for (mi in seq_len(ncol(mm))) for (wi in seq_len(nrow(mm))) {
  mo <- colnames(mm)[mi]; w <- rownames(mm)[wi]
  out <- rbind(out, data.frame(window = w, class = bed$class[wi], dist = bed$dist[wi],
    gc = win_gc$gc[wi], motif = mo, present = as.integer(mm[wi, mi]),
    n_sites = as.integer(cnt[w, mo]), motif_gc = round(motif_gc[mo], 2),
    gc_class = gc_class[mo], stringsAsFactors = FALSE))
}
write.table(out, file.path(OUT, sprintf("g3f_%s_duck_motifs.tsv", OUT_TAG)), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote figures/data/g3f_", OUT_TAG, "_duck_motifs.tsv (", nrow(out), " rows)\n", sep = "")
cat("=== DONE g3f ===\n")
