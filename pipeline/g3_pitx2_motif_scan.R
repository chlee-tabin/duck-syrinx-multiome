#!/usr/bin/env Rscript
# g3_pitx2_motif_scan.R
#
# Static motif catalog at the PITX2 locus.  No L/R, no statistics —
# answers the question: which TF motifs sit at the PITX2 promoter
# peaks + neighborhood (±50 kb), with focus on the BMP/NODAL pathway
# anchors (SMAD-family + FOX-family) that drive the fig04 hypothesis?
#
# Inputs:
#   references/pitx2_promoter_regions.bed    (canonical PITX2 regions)
#   BSgenome.Anasplatyrhynchos.ZJU1.0         (custom forge from ArchR track)
#   JASPAR2020 vertebrate PWMs                (TFBSTools / motifmatchr stack)
#
# Outputs:
#   figures/data/g3_pitx2_motif_scan.tsv           (per-hit detail)
#   figures/data/g3_pitx2_motif_scan_summary.tsv   (region × motif count table)

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

# ── Inputs ─────────────────────────────────────────────────────────────
SCRATCH  <- Sys.getenv("SCRATCH_ROOT", SCRATCH_ROOT)
REPO_BED <- "references/pitx2_promoter_regions.bed"
BED_PATH <- if (file.exists(REPO_BED)) REPO_BED else file.path(SCRATCH, REPO_BED)

OUT_DIR  <- "figures/data"
OUT_HITS    <- file.path(OUT_DIR, "g3_pitx2_motif_scan.tsv")
OUT_SUMMARY <- file.path(OUT_DIR, "g3_pitx2_motif_scan_summary.tsv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Motif families of interest (BMP / NODAL / FOX / PITX2 self-regulation /
# splanchnic-LPM transcription factors).  IDs verified against JASPAR2020
# SQLite on 2026-05-11.  Note JASPAR2020 has no standalone SMAD1; the
# SMAD2::SMAD3::SMAD4 trimer captures the canonical SBE.  Smad2::Smad3
# is the mouse dimer model (separate motif from human).
TARGET_MOTIFS <- c(
  # BMP/NODAL transducers
  SMAD2_3_4_human = "MA0513.1",  # SMAD2::SMAD3::SMAD4 — canonical SBE
  SMAD3           = "MA0795.1",
  SMAD4           = "MA1153.1",  # co-SMAD
  SMAD5           = "MA1557.1",  # BMP-target
  Smad2_3_mouse   = "MA1622.1",  # Smad2::Smad3 dimer
  # NODAL receptor + pan-axial mesoderm Fox-family
  FOXH1           = "MA0479.1",  # NODAL receptor
  FOXA2           = "MA0047.3",
  Foxa2_mouse     = "MA0047.2",
  FOXC1           = "MA0032.2",
  FOXC2           = "MA0846.1",
  Foxf1_mouse     = "MA1606.1",  # FoxF1 — splanchnic LPM marker
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
  MEF2C           = "MA0497.1"
)

P_CUTOFF <- as.numeric(Sys.getenv("MOTIF_P_CUTOFF", "5e-4"))

cat("=== g3_pitx2_motif_scan ===\n")
cat("BED:        ", BED_PATH, "\n")
cat("Genome:     BSgenome.Anasplatyrhynchos.ZJU1.0\n")
cat("Motif set:  JASPAR2020 (", length(TARGET_MOTIFS), " targets)\n", sep = "")
cat("p.cutoff:   ", P_CUTOFF, "\n", sep = "")
cat("OUT_HITS:   ", OUT_HITS, "\n")
cat("OUT_SUMMARY:", OUT_SUMMARY, "\n\n")

# ── Load BED + genome ──────────────────────────────────────────────────
genome <- BSgenome.Anasplatyrhynchos.ZJU1.0
bed    <- import(BED_PATH, format = "BED")
if (length(bed) == 0L) stop("BED file empty: ", BED_PATH)

# rtracklayer keeps the BED's `name` column on the GRanges.
cat("Loaded", length(bed), "regions:\n")
for (i in seq_along(bed)) {
  cat(sprintf("  %2d. %-40s %s:%d-%d (%d bp)\n",
              i, bed[i]$name, as.character(seqnames(bed[i])),
              start(bed[i]), end(bed[i]), width(bed[i])))
}
cat("\n")

# ── Load PWMs (one-by-one via getMatrixByID for robustness) ────────────
# `getMatrixSet` with species/tax_group filter returned 0 matrices on
# o2's JASPAR2020 install; pulling individual IDs via getMatrixByID
# bypasses the filter and works reliably.
target_pfm_list <- list()
for (label in names(TARGET_MOTIFS)) {
  mid <- TARGET_MOTIFS[[label]]
  pfm <- getMatrixByID(JASPAR2020, ID = mid)   # a missing motif is an error, never a skip
  if (is.null(pfm)) stop("motif ", label, " (", mid, ") not in JASPAR2020")
  target_pfm_list[[label]] <- pfm
}

if (length(target_pfm_list) == 0L) stop("No motifs resolved.  Check IDs.")

# Convert PFMs → PWMs (uniform background)
target_pwms <- do.call(PWMatrixList,
                       lapply(target_pfm_list, toPWM, pseudocounts = 0.8))
# Name the PWMs by our short labels (FOXH1, SMAD3, ...).
names(target_pwms) <- names(target_pfm_list)

cat("Resolved", length(target_pwms), "PWMs.  Scanning...\n")

# ── Scan ───────────────────────────────────────────────────────────────
match_positions <- matchMotifs(target_pwms, bed,
                                genome = genome,
                                out = "positions",
                                p.cutoff = P_CUTOFF)

# match_positions is a list-like (length == n_motifs) where each element is
# a GRangesList (length == n_regions) of hit positions in that region.

id_to_label <- setNames(names(TARGET_MOTIFS), as.character(TARGET_MOTIFS))
# motifmatchr returns one GRanges per motif (length = n_motifs); each
# GRanges holds all hits across the union of input regions in GENOME
# coords.  Use findOverlaps(hits, bed) to map back to input regions.

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

write.table(hits_df, OUT_HITS, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", nrow(hits_df), "hit records to", OUT_HITS, "\n\n")

# ── Region × motif count summary ───────────────────────────────────────
summary_df <- data.frame()
for (region_name in unique(bed$name)) {
  bed_row <- bed[bed$name == region_name][1]
  width_bp <- width(bed_row)
  for (label in names(target_pwms)) {
    mid <- TARGET_MOTIFS[[label]]
    n_hits <- nrow(hits_df[hits_df$region == region_name &
                             hits_df$motif_label == label, ])
    # Hits per kb so wide regions don't dominate by size.
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

write.table(summary_df, OUT_SUMMARY, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote ", nrow(summary_df), " summary rows to ", OUT_SUMMARY, "\n", sep = "")

# Pretty-print to log
cat("\n=== Hits per kb (region × motif) ===\n")
wide <- reshape(summary_df[, c("region", "motif_label", "hits_per_kb")],
                idvar = "region", timevar = "motif_label", direction = "wide")
colnames(wide) <- sub("^hits_per_kb\\.", "", colnames(wide))
print(format(wide, digits = 3), row.names = FALSE)

cat("\nDONE.\n")
