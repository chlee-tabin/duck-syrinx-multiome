#!/usr/bin/env Rscript
# hh21_trans_motif_scan.R — per-promoter-peak motif INSTANCE COUNTS + GC for the HH21
# response-vs-bilateral promoter set (see figures/figS03/hh21_trans_motif.py).
#
# Same engine, same JASPAR panel, same fixed background and same GC guardrail as
# pipeline/g3f_hh21_motif_scan.R — so counts are comparable with the fig04 panel-C scans.
# Differs only in scope (thousands of promoter peaks, not a curated window list) and in
# writing WIDE output (one row per peak), because the long rbind loop in g3f is O(n^2).
#
# In:  $SCAN_BED    4-col BED (chrom, start, end, id), 0-based/half-open.
#              default figures/data/hh21_trans_promoters.bed (the ATAC-peak promoters)
# Out: $SCAN_OUT    <$SCAN_ID_COL>, width, gc, gc_bases, n_TAATCC_raw, n_<MOTIF> x12
#              default figures/data/hh21_trans_promoter_motifs.tsv
#      figures/data/hh21_trans_motif_panel.tsv       motif, jaspar_id, motif_gc, gc_class
.libPaths(c("~/R/4.4.2", .libPaths()))
suppressPackageStartupMessages({
  library(motifmatchr); library(TFBSTools); library(JASPAR2020)
  library(BSgenome.Anasplatyrhynchos.ZJU1.0)
  library(Biostrings); library(GenomicRanges)
})

DATA <- Sys.getenv("DATA_DIR", "figures/data")
# Generic over any 4-column BED (chrom, start, end, id) ON DUCK ZJU1.0 — the BSgenome and the
# PWM background BG below are duck-specific and must be changed TOGETHER for another genome
# (a swapped BSgenome with duck BG silently shifts every p.cutoff=1e-3 hit count by GC).
# Defaults = peak mode. Sequence mode MUST set all three env vars; a partially-applied
# prefix would otherwise fall back to peak mode and overwrite the peak scan in place.
BED <- Sys.getenv("SCAN_BED", file.path(DATA, "hh21_trans_promoters.bed"))
OUT_TSV <- Sys.getenv("SCAN_OUT", file.path(DATA, "hh21_trans_promoter_motifs.tsv"))
ID_COL <- Sys.getenv("SCAN_ID_COL", "peak_id")
PANEL_OUT <- Sys.getenv("SCAN_PANEL_OUT", file.path(DATA, "hh21_trans_motif_panel.tsv"))
# all-or-nothing: any one of the three set means sequence mode, so all three must be set
set3 <- c(SCAN_BED = Sys.getenv("SCAN_BED"), SCAN_OUT = Sys.getenv("SCAN_OUT"),
          SCAN_ID_COL = Sys.getenv("SCAN_ID_COL")) != ""
if (any(set3) && !all(set3))
  stop(sprintf("SCAN_BED/SCAN_OUT/SCAN_ID_COL must be set together (missing: %s) — a partial\n  prefix would silently rescan the peak BED over the peak output",
               paste(names(set3)[!set3], collapse = ", ")))
PCUT <- as.numeric(Sys.getenv("MOTIF_P_CUTOFF", "1e-3"))
set.seed(1)
BG <- c(A = 0.293, C = 0.207, G = 0.207, T = 0.293)  # SAME fixed bg as g3c/g3d/g3e/g3f

# 4-col = the v1 peak BED; 5-col = a v2/v3 promoter BED carrying prep()'s provenance stamp in
# column 5. The stamp travels WITH the coordinates so the scan output can be checked against the
# code that produced the windows — stamping only the map left a stale scan undetectable whenever
# the gene set and widths were unchanged.
bed <- read.table(BED, sep = "\t", stringsAsFactors = FALSE)
if (ncol(bed) == 5) {
  colnames(bed) <- c("chrom", "start", "end", "id", "prep_stamp")
} else if (ncol(bed) == 4) {
  colnames(bed) <- c("chrom", "start", "end", "id")
  bed$prep_stamp <- NA_character_          # v1 peak BED predates the stamp
} else {
  stop(sprintf("scan: %s has %d columns; expected 4 (peak BED) or 5 (stamped promoter BED)",
               BED, ncol(bed)))
}
if (ncol(bed) == 5 && length(unique(bed$prep_stamp)) != 1L)
  stop("scan: BED carries more than one prep_stamp — it was concatenated from two preps")
gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end))   # BED is 0-based/half-open
seq <- getSeq(BSgenome.Anasplatyrhynchos.ZJU1.0, gr)
cat("scan: ", length(seq), " intervals, widths ",
    paste(range(width(seq)), collapse = "-"), "\n", sep = "")

# SAME 12-motif panel as g3f_hh21_motif_scan.R
TARGETS <- c(SMAD2_3_4 = "MA0513.1", SMAD3 = "MA0795.1", SMAD4 = "MA1153.1", SMAD5 = "MA1557.1",
             FOXH1 = "MA0479.1", NKX6_1 = "MA0674.1", Foxf1 = "MA1606.1", FOXC1 = "MA0032.2",
             GATA4 = "MA0482.2", ISL1 = "MA1608.1", TBX5 = "MA0807.1", PITX2 = "MA1547.1")
pfm <- lapply(TARGETS, function(id) getMatrixByID(JASPAR2020, ID = id))
if (any(vapply(pfm, is.null, logical(1)))) stop("scan: a JASPAR2020 matrix failed to load")
pwms <- do.call(PWMatrixList, lapply(pfm, toPWM, pseudocounts = 0.8))
names(pwms) <- names(pfm)
motif_gc <- vapply(pfm, function(m) { M <- m@profileMatrix; sum(M[c("C", "G"), ]) / sum(M) }, 0)

write.table(data.frame(motif = names(TARGETS), jaspar_id = unname(TARGETS),
                       motif_gc = round(motif_gc, 3),
                       # descriptive bins only. A PWM under 60% GC is NOT thereby "composition-immune" —
# the k-mer shuffle null (hh21_smad234_dinuc_null.R) is the real composition control.
                       gc_class = ifelse(motif_gc >= 0.60, "GC-rich", "GC-poor")),
            PANEL_OUT, sep = "\t",
            quote = FALSE, row.names = FALSE)

posl <- matchMotifs(pwms, seq, bg = BG, out = "positions", p.cutoff = PCUT)
cnt <- vapply(posl, function(x) as.integer(lengths(x)), integer(length(seq)))
colnames(cnt) <- paste0("n_", names(posl))

out <- data.frame(id = bed$id,
                  width = as.integer(width(seq)),
                  gc = round(as.numeric(letterFrequency(seq, "GC", as.prob = TRUE)), 4),
                  # gc_bases lets a caller re-derive GC after merging variable-width intervals
                  gc_bases = as.integer(letterFrequency(seq, "GC")),
                  # raw PITX/bicoid core, both strands. A COMPARATOR, not a composition
                  # control: TAATCC is strongly AT-rich and itself composition-dependent.
                  n_TAATCC_raw = as.integer(vcountPattern("TAATCC", seq) +
                                            vcountPattern("GGATTA", seq)),
                  cnt, prep_stamp = bed$prep_stamp,
                  check.names = FALSE, stringsAsFactors = FALSE)
names(out)[1] <- ID_COL
write.table(out, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat("scan: wrote ", OUT_TSV, " (", nrow(out), " intervals)\n", sep = "")
cat("scan: mean PITX2 sites/peak = ", round(mean(out$n_PITX2), 3),
    " | mean raw TAATCC = ", round(mean(out$n_TAATCC_raw), 3), "\n", sep = "")
cat("=== DONE hh21_trans_motif_scan ===\n")
