#!/usr/bin/env Rscript
# g3b_pitx2_candidate_peaks_motif_scan.R
#
# Per-peak motif scan for the 8 candidate PITX2-locus peaks
# (Peak #1/#2/#3/#4/#5 + SMAD5+ composite + PITX2a promoter peak +
# ASE_region — syntenic to Welsh 2015 chick ASE in PITX2c last intron),
# as a follow-up to the g3 §5.7 relaxed-cutoff scan that profiled
# Peak #1 only (and within a 601-bp `promoter|PITX2c` BED window
# rather than the 501-bp Peak #1 window used here, and on a 21-motif
# catalog rather than the 27-motif catalog here).
#
# Question we are answering: which of the candidate peaks are PITX-family
# motif-dense vs motif-desert at 3 relaxed cutoffs (p<5e-4 / p<5e-3 /
# p<1e-2), and which carry the SMAD2/3 + SMAD5 + FOX + MEIS/PBX
# signatures invoked by the §2 BMP-priming / §3.4 PITX2c-motif-desert
# hypotheses.
#
# Inputs:
#   references/pitx2_candidate_peaks.bed    (8 named peaks; this script)
#   BSgenome.Anasplatyrhynchos.ZJU1.0
#   JASPAR2020 vertebrate PWMs               (TFBSTools / motifmatchr stack)
#
# Outputs:
#   figures/data/g3b_candidate_peaks_motif_scan_p<cutoff>.tsv          (per-hit)
#   figures/data/g3b_candidate_peaks_motif_scan_summary_p<cutoff>.tsv (peak × motif)
#   figures/data/g3b_candidate_peaks_motif_family_table.tsv           (composite — all 3 cutoffs)
#
# Env vars:
#   REPO_BED        (default: references/pitx2_candidate_peaks.bed)
#   MOTIF_P_CUTOFF  (default: 5e-4 — but this driver runs all 3)
#
# Run on O2:
#   sbatch -p priority --mem=16G -t 0-0:30 -c 2 \
#     --wrap='module load gcc/14.2.0 R/4.4.2 && Rscript pipeline/g3b_pitx2_candidate_peaks_motif_scan.R'

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

# ── Paths ─────────────────────────────────────────────────────────────
SCRATCH  <- Sys.getenv("SCRATCH_ROOT", SCRATCH_ROOT)
REPO_BED <- Sys.getenv("REPO_BED", "references/pitx2_candidate_peaks.bed")
BED_PATH <- if (file.exists(REPO_BED)) REPO_BED else file.path(SCRATCH, REPO_BED)

OUT_DIR  <- "figures/data"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Same curated motif set as g3 so per-peak counts are directly comparable
# to g3_pitx2_motif_scan_summary_p*.tsv.
TARGET_MOTIFS <- c(
  # BMP/NODAL transducers
  SMAD2_3_4_human = "MA0513.1",
  SMAD3           = "MA0795.1",
  SMAD4           = "MA1153.1",
  SMAD5           = "MA1557.1",
  Smad2_3_mouse   = "MA1622.1",
  # NODAL receptor + pan-axial mesoderm Fox-family
  FOXH1           = "MA0479.1",
  FOXA2           = "MA0047.3",
  Foxa2_mouse     = "MA0047.2",
  FOXC1           = "MA0032.2",
  FOXC2           = "MA0846.1",
  Foxf1_mouse     = "MA1606.1",
  # PITX-family (auto-loop + paralogs)
  PITX1           = "MA0682.2",
  PITX2           = "MA1547.1",
  PITX3           = "MA0714.1",
  # Splanchnic / cardiopharyngeal LPM master regulators
  NKX2_5          = "MA0063.2",
  NKX6_1          = "MA0674.1",
  TBX5            = "MA0807.1",
  GATA4           = "MA0482.2",
  ISL1            = "MA1608.1",
  HAND1_TCF3      = "MA0092.1",
  MEF2C           = "MA0497.1",
  # MEIS / PBX — homeobox cofactors invoked by §5.7 narrative (some
  # auto-regulators of PITX-family loci sit on MEIS-PBX heterodimers)
  MEIS1           = "MA0498.2",
  MEIS2           = "MA0774.1",
  MEIS3           = "MA0775.1",
  PBX1            = "MA0070.1",
  PBX2            = "MA1113.2",
  PBX3            = "MA1114.1"
)

# Motif families for the per-peak × per-family rollup at the end.
MOTIF_FAMILIES <- list(
  PITX     = c("PITX1", "PITX2", "PITX3"),
  SMAD2_3  = c("SMAD2_3_4_human", "SMAD3", "Smad2_3_mouse"),
  SMAD5    = c("SMAD5"),
  SMAD4    = c("SMAD4"),
  FOX      = c("FOXH1", "FOXA2", "Foxa2_mouse", "FOXC1", "FOXC2", "Foxf1_mouse"),
  TBX5     = c("TBX5"),
  GATA4    = c("GATA4"),
  ISL1     = c("ISL1"),
  NKX2_5   = c("NKX2_5"),
  NKX6_1   = c("NKX6_1"),
  HAND1    = c("HAND1_TCF3"),
  MEF2     = c("MEF2C"),
  MEIS     = c("MEIS1", "MEIS2", "MEIS3"),
  PBX      = c("PBX1", "PBX2", "PBX3")
)

# Cutoffs to run.  g3 runs were single-cutoff; here we sweep all 3 in
# one driver to keep the per-cutoff TSVs aligned.
P_CUTOFFS <- c("5e-4" = 5e-4, "5e-3" = 5e-3, "1e-2" = 1e-2)

cat("=== g3b_pitx2_candidate_peaks_motif_scan ===\n")
cat("BED:        ", BED_PATH, "\n")
cat("Genome:     BSgenome.Anasplatyrhynchos.ZJU1.0\n")
cat("Motif set:  JASPAR2020 (", length(TARGET_MOTIFS), " targets)\n", sep = "")
cat("p.cutoffs:  ", paste(names(P_CUTOFFS), collapse = " / "), "\n", sep = "")
cat("OUT_DIR:    ", OUT_DIR, "\n\n")

# ── Load BED + genome ──────────────────────────────────────────────────
genome <- BSgenome.Anasplatyrhynchos.ZJU1.0
bed    <- import(BED_PATH, format = "BED")
if (length(bed) == 0L) stop("BED file empty: ", BED_PATH)

cat("Loaded", length(bed), "regions:\n")
for (i in seq_along(bed)) {
  cat(sprintf("  %2d. %-25s %s:%d-%d (%d bp)\n",
              i, bed[i]$name, as.character(seqnames(bed[i])),
              start(bed[i]), end(bed[i]), width(bed[i])))
}
cat("\n")

# ── Load PWMs (one-by-one via getMatrixByID for robustness) ────────────
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
cat("Resolved", length(target_pwms), "PWMs.\n\n")

# ── Helper: run scan at a single cutoff and write outputs ──────────────
scan_at_cutoff <- function(cutoff_label, p_cutoff) {
  cat("=== Scanning at p<", cutoff_label, " (", p_cutoff, ") ===\n", sep = "")
  match_positions <- matchMotifs(target_pwms, bed,
                                  genome = genome,
                                  out = "positions",
                                  p.cutoff = p_cutoff)

  records <- vector("list", 0L)
  for (label in names(match_positions)) {
    hits <- match_positions[[label]]
    if (length(hits) == 0L) next
    mid <- TARGET_MOTIFS[[label]]
    ovl <- findOverlaps(hits, bed)
    if (length(ovl) == 0L) next
    for (k in seq_along(ovl)) {
      hit_i <- queryHits(ovl)[k]
      region_i <- subjectHits(ovl)[k]
      bed_row <- bed[region_i]
      hit <- hits[hit_i]
      records[[length(records) + 1L]] <- data.frame(
        motif_id     = mid,
        motif_label  = label,
        region       = bed_row$name,
        region_chrom = as.character(seqnames(bed_row)),
        region_start = start(bed_row),
        region_end   = end(bed_row),
        hit_start    = start(hit),
        hit_end      = end(hit),
        hit_strand   = as.character(strand(hit)),
        hit_score    = round(score(hit), 4),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(records) == 0L) {
    hits_df <- data.frame(
      motif_id = character(), motif_label = character(),
      region = character(), region_chrom = character(),
      region_start = integer(), region_end = integer(),
      hit_start = integer(), hit_end = integer(),
      hit_strand = character(), hit_score = double()
    )
  } else {
    hits_df <- do.call(rbind, records)
  }

  out_hits <- file.path(OUT_DIR,
                        sprintf("g3b_candidate_peaks_motif_scan_p%s.tsv", cutoff_label))
  write.table(hits_df, out_hits, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  ", nrow(hits_df), " hit records → ", out_hits, "\n", sep = "")

  # Per-peak × per-motif summary
  summary_df <- data.frame()
  for (region_name in unique(bed$name)) {
    bed_row <- bed[bed$name == region_name][1]
    width_bp <- width(bed_row)
    for (label in names(target_pwms)) {
      mid <- TARGET_MOTIFS[[label]]
      n_hits <- nrow(hits_df[hits_df$region == region_name &
                               hits_df$motif_label == label, ])
      summary_df <- rbind(summary_df, data.frame(
        region       = region_name,
        region_width = width_bp,
        motif_id     = mid,
        motif_label  = label,
        n_hits       = n_hits,
        hits_per_kb  = round(n_hits / (width_bp / 1000), 3),
        stringsAsFactors = FALSE
      ))
    }
  }
  out_sum <- file.path(OUT_DIR,
                       sprintf("g3b_candidate_peaks_motif_scan_summary_p%s.tsv",
                               cutoff_label))
  write.table(summary_df, out_sum, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  ", nrow(summary_df), " summary rows → ", out_sum, "\n", sep = "")

  list(hits = hits_df, summary = summary_df)
}

# ── Run at all cutoffs ─────────────────────────────────────────────────
all_results <- list()
for (lab in names(P_CUTOFFS)) {
  all_results[[lab]] <- scan_at_cutoff(lab, P_CUTOFFS[[lab]])
  cat("\n")
}

# ── Composite per-peak × per-motif-family table across all 3 cutoffs ──
# Output schema: one row per (region, family, cutoff) with n_hits and
# hits_per_kb, plus a wide pivot that gives the same "Peak | family @
# cutoff" shape requested in the task brief.

cat("=== Building per-peak × family rollup ===\n")
fam_long <- data.frame()
for (lab in names(P_CUTOFFS)) {
  summary_df <- all_results[[lab]]$summary
  for (region_name in unique(summary_df$region)) {
    width_bp <- summary_df$region_width[summary_df$region == region_name][1]
    for (fam in names(MOTIF_FAMILIES)) {
      labels_in_fam <- MOTIF_FAMILIES[[fam]]
      sub <- summary_df[summary_df$region == region_name &
                          summary_df$motif_label %in% labels_in_fam, ]
      n_hits <- sum(sub$n_hits)
      fam_long <- rbind(fam_long, data.frame(
        region       = region_name,
        region_width = width_bp,
        family       = fam,
        cutoff       = lab,
        n_hits       = n_hits,
        hits_per_kb  = round(n_hits / (width_bp / 1000), 3),
        stringsAsFactors = FALSE
      ))
    }
  }
}

out_family <- file.path(OUT_DIR, "g3b_candidate_peaks_motif_family_table.tsv")
write.table(fam_long, out_family, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", nrow(fam_long), "family-cutoff rows to", out_family, "\n\n")

# ── Pretty-print to log: wide table per family, cols = peak × cutoff ──
cat("=== Per-family n_hits (rows = peak, cols = cutoff) ===\n\n")
for (fam in names(MOTIF_FAMILIES)) {
  cat(sprintf("--- %s family (%s) ---\n", fam,
              paste(MOTIF_FAMILIES[[fam]], collapse = ", ")))
  sub <- fam_long[fam_long$family == fam, ]
  wide <- reshape(sub[, c("region", "cutoff", "n_hits")],
                  idvar = "region", timevar = "cutoff", direction = "wide")
  colnames(wide) <- sub("^n_hits\\.", "p<", colnames(wide))
  print(wide, row.names = FALSE)
  cat("\n")
}

cat("DONE.\n")
