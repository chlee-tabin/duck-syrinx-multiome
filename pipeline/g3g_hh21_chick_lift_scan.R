#!/usr/bin/env Rscript
# g3g_hh21_chick_lift_scan.R — CHICK lift + motif presence/COUNTS for the g4f floor-20 HH21 windows.
#
# Chick counterpart of g3f: Smith-Waterman-lift each PROXIMAL duck window (the L cis-block + a/b
# promoters + ns controls; the R-desert Class III peaks are 400+ kb distal and not in the proximal
# chick region, so they stay duck-only) onto chick (galGal6a), call Class I (conserved) / Class II
# (duck-specific) by % identity, and scan the conserved chick orthologs with presence + INSTANCE
# COUNTS on the SAME fixed-bg panel as g3f.
#
# Output: figures/data/g3g_hh21_chick_motifs.tsv (window,class,dist,conservation,pct_id,motif,
#   duck_present,chick_present,duck_n_sites,chick_n_sites) + g3g_hh21_lift.tsv (per-window lift).
.libPaths(c("~/R/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(motifmatchr); library(TFBSTools); library(JASPAR2020)
  library(BSgenome.Anasplatyrhynchos.ZJU1.0); library(Biostrings); library(GenomicRanges)
})
OUT <- "figures/data"
WIN_BED <- Sys.getenv("WIN_BED", "")
if (!nzchar(WIN_BED)) stop("WIN_BED must be set (pair it with the matching OUT_TAG). Runners: pipeline/run_g3f_g3g_domain.sh (shipped; committed figures/data/g4f_domain_motif_windows.bed)")
OUT_TAG <- Sys.getenv("OUT_TAG", "hh21")   # "domain" for the multi-stage fig04 panel C re-scope
CHICK_PROX_FA <- Sys.getenv("CHICK_PROX_FA", "chick_pitx2_proximal_hh21.fa")
CHICK_PROX_START <- as.integer(Sys.getenv("CHICK_PROX_START", "57782829"))
PCUT <- as.numeric(Sys.getenv("MOTIF_P_CUTOFF", "1e-3"))
ID_THRESH <- as.numeric(Sys.getenv("ID_THRESH", "0.55")); COV_THRESH <- as.numeric(Sys.getenv("COV_THRESH", "0.50"))
# proximal classes only (R-desert is distal -> duck-only, not lifted here)
KEEP <- strsplit(Sys.getenv("KEEP_CLASS",
  "Lcis_HH21_paired,Lcis_HH21_unpaired,ctrl_aprom,ctrl_bprom,ctrl_ns"), ",")[[1]]
# TSS anchor for the `dist` annotation column. Default = the PITX2c TSS (the
# original locus); overridable so the same lift+scan can be pointed at another
# locus (g3h: duck BMP7). Annotation only — it does not affect lift or scan.
TSS_C <- as.integer(Sys.getenv("TSS_ANCHOR", "3971385")); set.seed(1)
duck <- BSgenome.Anasplatyrhynchos.ZJU1.0
BG <- c(A = 0.293, C = 0.207, G = 0.207, T = 0.293)

bed <- read.table(WIN_BED, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(bed)[1:4] <- c("chrom", "start", "end", "name")
p <- do.call(rbind, strsplit(bed$name, "\\|")); bed$label <- p[, 1]; bed$class <- p[, 2]
bed <- bed[bed$class %in% KEEP, , drop = FALSE]
bed$dist <- as.integer((bed$start + bed$end) / 2 - TSS_C)
gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end)); duck_seq <- getSeq(duck, gr); names(duck_seq) <- bed$label
chick_prox <- readDNAStringSet(CHICK_PROX_FA)[[1]]
cat("g3g: ", length(duck_seq), " proximal windows; chick region ", length(chick_prox), " bp (start ", CHICK_PROX_START, ")\n", sep = "")

# motif panel (shared)
TARGETS <- c(SMAD2_3_4 = "MA0513.1", SMAD3 = "MA0795.1", SMAD4 = "MA1153.1", SMAD5 = "MA1557.1",
             FOXH1 = "MA0479.1", NKX6_1 = "MA0674.1", Foxf1 = "MA1606.1", FOXC1 = "MA0032.2",
             GATA4 = "MA0482.2", ISL1 = "MA1608.1", TBX5 = "MA0807.1", PITX2 = "MA1547.1")
pfm <- list(); for (lab in names(TARGETS)) {   # a missing motif is an error
  m <- getMatrixByID(JASPAR2020, ID = TARGETS[[lab]])
  if (is.null(m)) stop("g3g: motif ", lab, " (", TARGETS[[lab]], ") not in JASPAR2020")
  pfm[[lab]] <- m }
for (nm in c("SMAD1", "SMAD9")) {  # other BMP R-SMADs by name (share the SBE with SMAD5)
  s <- tryCatch(getMatrixSet(JASPAR2020, list(name = nm)), error = function(e) NULL)
  if (!is.null(s) && length(s) > 0) { pfm[[nm]] <- s[[1]]; cat("g3g: added", nm, "=", ID(s[[1]]), "\n") }
  else cat("g3g: JASPAR2020 has NO standalone", nm, "matrix\n")
}
pwms <- do.call(PWMatrixList, lapply(pfm, toPWM, pseudocounts = 0.8)); names(pwms) <- names(pfm)

scan_pc <- function(seqs) {  # -> list(present=matrix, count=matrix), rows=names(seqs)
  mm <- as.matrix(motifMatches(matchMotifs(pwms, seqs, bg = BG, out = "scores", p.cutoff = PCUT)))
  rownames(mm) <- names(seqs)
  posl <- matchMotifs(pwms, seqs, bg = BG, out = "positions", p.cutoff = PCUT)
  cnt <- matrix(0L, length(seqs), length(posl), dimnames = list(names(seqs), names(posl)))
  for (mi in seq_along(posl)) {           # IRangesList, one per sequence; lengths() = sites per seq
    el <- lengths(posl[[mi]]); nm <- names(el)
    if (!is.null(nm) && all(nm %in% names(seqs))) cnt[nm, names(posl)[mi]] <- as.integer(el)
    else if (length(el) == nrow(cnt)) cnt[, names(posl)[mi]] <- as.integer(el)
  }
  list(present = mm, count = cnt)
}

# Smith-Waterman lift
mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
lift <- data.frame(); chick_hits <- DNAStringSet()
for (i in seq_along(duck_seq)) {
  aln <- pairwiseAlignment(duck_seq[[i]], chick_prox, type = "local", substitutionMatrix = mat, gapOpening = 5, gapExtension = 2)
  pid <- pid(aln, type = "PID1") / 100; cov <- width(pattern(aln)) / length(duck_seq[[i]])
  s <- start(subject(aln)); e <- end(subject(aln)); cons <- (pid >= ID_THRESH) && (cov >= COV_THRESH)
  lift <- rbind(lift, data.frame(window = bed$label[i], class = bed$class[i], dist = bed$dist[i],
    pct_id = round(pid, 3), cov = round(cov, 3), chick_start = CHICK_PROX_START + s - 1L,
    chick_end = CHICK_PROX_START + e - 1L,
    conservation = ifelse(cons, "ClassI_conserved", "ClassII_duck_specific"), stringsAsFactors = FALSE))
  if (cons) { cs <- DNAStringSet(subseq(chick_prox, s, e)); names(cs) <- bed$label[i]; chick_hits <- c(chick_hits, cs) }
}
write.table(lift, file.path(OUT, sprintf("g3g_%s_lift.tsv", OUT_TAG)), sep = "\t", quote = FALSE, row.names = FALSE)
cat("lift:\n"); print(table(lift$conservation))

dk <- scan_pc(duck_seq)
ck <- if (length(chick_hits) > 0) scan_pc(chick_hits) else list(present = NULL, count = NULL)
# GC-immune Smad-binding-element (SBE) raw count AGAC/GTCT per window (duck + chick ortholog)
sbe_count <- function(s) countPattern("AGAC", s) + countPattern("GTCT", s)
duck_sbe <- setNames(vapply(duck_seq, sbe_count, integer(1)), names(duck_seq))
chick_sbe <- if (length(chick_hits) > 0) setNames(vapply(chick_hits, sbe_count, integer(1)), names(chick_hits)) else integer(0)
cons <- setNames(lift$conservation, lift$window)
out <- data.frame()
for (wi in seq_len(nrow(dk$present))) { w <- rownames(dk$present)[wi]
  has_ck <- !is.null(ck$present) && w %in% rownames(ck$present)
  for (mo in colnames(dk$present)) {
    out <- rbind(out, data.frame(window = w, class = bed$class[match(w, bed$label)],
      dist = bed$dist[match(w, bed$label)], conservation = cons[[w]], motif = mo,
      duck_present = as.integer(dk$present[w, mo]), chick_present = if (has_ck) as.integer(ck$present[w, mo]) else NA_integer_,
      duck_n_sites = as.integer(dk$count[w, mo]), chick_n_sites = if (has_ck) as.integer(ck$count[w, mo]) else NA_integer_,
      duck_sbe = as.integer(duck_sbe[w]), chick_sbe = if (has_ck) as.integer(chick_sbe[w]) else NA_integer_,
      stringsAsFactors = FALSE)) } }
write.table(out, file.path(OUT, sprintf("g3g_%s_chick_motifs.tsv", OUT_TAG)), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote figures/data/g3g_", OUT_TAG, "_chick_motifs.tsv (", nrow(out), " rows)\n=== DONE g3g ===\n", sep = "")
