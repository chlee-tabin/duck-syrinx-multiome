#!/usr/bin/env Rscript
# panel_C_pileup.R — fig03 panel C: focal PITX2-locus L/R chromatin pileup.
#
# COVERAGE = PER-EMBRYO-HALF (pseudoreplication-rectified, 2026-06-17).  Each
#    (stage x compartment x genotype x side) embryo-HALF contributes EQUALLY to the
#    averaged track (pooled_pc() = rowMeans over the per-genotype-half _percell
#    tracks), so high-cell-count embryos no longer dominate the track SHAPE.  This
#    is the coverage analog of the g4f embryo-half stats rectification.  The L/R
#    bias survives equal weighting (inference-pool + ventral locus L/R 1.24-1.46:
#    HH17/16 pool 1.24, HH21/HH25 ventral 1.45-1.46, vs background ~1.1-1.4; the
#    dorsal lanes show no significant locus signal but are UNDERPOWERED at this
#    promoter — see the SIG block).  Significance STARS read the (b') paired
#    within-population DESeq2.
#
# Compartment-resolved L/R chromatin accessibility at the PITX2 locus, the
# coverage-track companion to panel B's L/R peak counts.  Six "mirror" lanes
# (Left = +blue up, Right = -red down): HH13 pool (Wave-1), HH17/16 pool, then
# {ventral splanchnic, dorsal sub-epithelial} x {HH21, HH25}.  The arc:
#   * VENTRAL splanchnic (inference pool)  -> Left-biased PITX2c c-module
#   * DORSAL sub-epithelial               -> underpowered control (no detected
#     bias; the within-locus specificity claim rides on c vs a/b)
#
# Adapted from an earlier BMP-mechanism pileup (not distributed).  TWO changes
# make it self-contained and statistically consistent with panels A/B:
#   1. Significance stars come from the (b') PAIRED within-population laterality
#      DESeq2 (figures/data/bprime_locus_fdr.tsv —
#      ONE pseudobulk per embryo-half WITHIN one d2g population, ~ genotype +
#      side, complete pairs, floor 20, locus-restricted BH).
#      The retired earlier stars used per-subcluster pseudobulks (h7c/h7_unpaired),
#      which carry the pseudoreplication the (b') design was built to remove.
#   2. The locus peak set (ticks + background mask) is the pool-independent
#      ArchR-union FAM241A-extended ±500 kb set (figures/data/pitx2_locus_tad_peaks.tsv; NOT a measured TAD), not the
#      earlier h7c called set.  No h7c/h7 import.
# The pileup is a PER-EMBRYO-HALF visual (equal-weight average L vs R accessibility
# shape); the L/R CLAIM is the (b') paired within-population DESeq2 (the stars).
# Per-lane raw L/R depth factors are annotated, and each lane's background is
# matched to L/R = 1 so a concentrated c-module L-bias survives but uniform depth
# does not.
#
# Inputs  : figures/data/f4_compartment_perembryo_pitx2_coverage.tsv.gz  (per-base
#           per-cell coverage; <stage>_<compartment>_<genotype>_<side>[_percell] cols) =
#           the f4 export on the current lr_inference_multiome.tsv sides; a split-tube
#           half (17-5/17-6) sums its two aliquots, as in g1
#           figures/data/f4_compartment_perembryo_pitx2_cell_counts.tsv  (nuclei per column)
#           figures/data/f4_hh13_perembryo_pitx2_coverage.tsv.gz  (HH13 pooled lane,
#           per-genotype-half)
#           figures/data/bprime_locus_fdr.tsv
#                                        ((b') paired within-population L/R DESeq2 — stars)
#           figures/data/g4f_sample_meta_floor20.tsv  (coverage-cohort restriction)
#           figures/data/pitx2_locus_tad_peaks.tsv              (ArchR-union peak set)
#           figures/data/pitx2_locus_gc_content.tsv             (GC% / CpG-island track)
# Output  : figures/fig03/panel_C_pileup.png   (GC + mirror pileup + gene model + peaks + stars)

local_lib <- "~/R/4.4.2"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

# ── paths (resolve relative to this script so it runs from anywhere) ─────────
HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
REPO <- normalizePath(file.path(HERE, "..", ".."))
DATA <- file.path(REPO, "figures", "data")
COV  <- file.path(DATA, "f4_compartment_perembryo_pitx2_coverage.tsv.gz")
COV_N <- file.path(DATA, "f4_compartment_perembryo_pitx2_cell_counts.tsv")
for (.f in c(COV, COV_N)) if (!file.exists(.f)) stop("missing ", .f, " — rerun f4 (pipeline/run_f4_perembryo.sh)")
source(file.path(REPO, "figures", "_shared", "pitx2_gene_model.R"))
OUT_PILEUP <- file.path(HERE, "panel_C_pileup.png")

# ── palette (house L/R from the color-registry SSOT; minimal color per project convention) ──
source(file.path(REPO, "figures", "_shared", "colors.R"))
lr_cols <- get_palette("lateral", labels = c("L", "R"))
L_BLUE <- unname(lr_cols["L"]); R_RED <- unname(lr_cols["R"])   # registry: L #4393C3, R #D6604D
INK <- "#2C3E50"; MUT <- "#5A6B7B"; FAINT <- "#9AA6B2"
C_BAND  <- L_BLUE      # PITX2c module highlight (the L-biased arm) = registry L blue
AB_BAND <- "#7F8C99"   # a/b promoter highlight (no-consistent-L/R-bias arm; neutral grey)
BIN <- 100             # bp smoothing window (peaks are 501 bp; 100 keeps shape, cuts per-base noise)

# ── featured elements: a/b/c PROMOTERS (load-bearing) + +5.4 kb (putative) ──
# Labeling decision (2026-06-17): the PITX2a/b/c PROMOTERS are the load-bearing
# units and are FOREGROUNDED (bold labels; the c-promoter is the CpG-island
# motif-DESERT memory element = the L arm; SMAD5 — the confirmed cis BMP/SBE input —
# annotates the a-promoter; a/b are the within-locus control, no consistent L/R
# bias).  Of the
# three exploratory distal marks only +5.4 kb survives, as a single PUTATIVE
# enhancer (italic, "(putative)"); +1.5 kb (over the PITX2c gene body, weakest) and
# +11.3 kb (weak, dorsal-trend) are DROPPED.  No "ASE" label (the ASE element is a
# different peak and the designation is retired); no NKX6-1 on distal peaks.
# `kind` drives the foregrounding (prom = bold + solid thick tick; putative =
# italic + thin dashed tick).
FEAT <- data.frame(
    name   = c("PITX2a\nprom\n(SMAD5)", "PITX2b\nprom", "PITX2c\nprom", "+5.4 kb\n(putative)"),
    start  = c(3961167, 3963950, 3971182, 3976533),
    end    = c(3961668, 3964451, 3971683, 3977034),
    module = c("ab", "ab", "c", "c"),
    kind   = c("prom", "prom", "prom", "putative"),
    stringsAsFactors = FALSE
)
FEAT$mid     <- (FEAT$start + FEAT$end) / 2
FEAT$face    <- ifelse(FEAT$kind == "prom", "bold", "italic")
FEAT$tick_lw <- ifelse(FEAT$kind == "prom", 0.9, 0.5)
# the a/b promoter callouts sit 2.8 kb apart and collide at the narrow
# (tall-variant) render width — splay them horizontally off their ticks
# (a right-aligned, b left-aligned), genome-browser style
FEAT$hjust <- c(1, 0, 0.5, 0.5)
FEAT$tick_lty <- ifelse(FEAT$kind == "prom", "solid", "21")

# ── significance source per lane = (b') PAIRED within-population DESeq2 ──────
#   Adopted 2026-08-19: the starred/primary
#   per-peak locus statistic at HH21/HH25 is the (b') design — ONE pseudobulk per
#   embryo-half WITHIN one d2g population (~ genotype + side, complete pairs,
#   floor 20, locus-restricted BH).  Ventral lanes star the Ventral-splanchnic-
#   population fit (8 L / 8 R embryo-halves, 0 singletons, at both stages);
#   dorsal lanes read their own (b') dorsal fit but are UNDERPOWERED (4 complete
#   pairs, SE ~ 1, estimates as large as ventral — the D/V contrast at this
#   promoter is power-limited, not effect-limited) -> "underpowered" note, and
#   n.s. means no star.  The inclusive-unpaired fit on the same unit (~ side,
#   singletons included) is the DISCLOSED SENSITIVITY arm (PROVENANCE.md) —
#   identical data at HH21/HH25 ventral, double the residual df.  HH17/16 lanes
#   carry NO stars: the HH17/16 sides are chromatin-defined (d2n), so a PITX2c-
#   promoter star would restate the side-defining signal (panel A withholds its
#   P there for the same reason); the pooled-unpaired trough value stays
#   disclosed in PROVENANCE/Methods only.  HH13 stays descriptive (pre-TES
#   pooled lane).
SIG <- list(
    "HH13_pool"       = list(stage = NA,        comp = NA),
    "HH17/16_ventral" = list(stage = NA,        comp = NA),
    "HH17/16_dorsal"  = list(stage = NA,        comp = NA),
    "HH21_ventral"    = list(stage = "HH21",    comp = "ventral"),
    "HH21_dorsal"     = list(stage = "HH21",    comp = "dorsal"),
    "HH25_ventral"    = list(stage = "HH25",    comp = "ventral"),
    "HH25_dorsal"     = list(stage = "HH25",    comp = "dorsal")
)

star_lab <- function(p) ifelse(is.na(p), "",
                        ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# ── load (b') paired within-population locus DESeq2 + peak coords ───────────
# Floor 20 = the settled default tier (S7 sweep + head-to-head — HH21/HH25
# are floor-invariant).  bprime "ventral" = the Ventral splanchnic POPULATION
# (one pseudobulk per embryo-half within it), NOT the pooled inference pool the
# coverage lanes average — the unit change is the point of the S1/S2 adoption:
# no pseudoreplication (embryo term), no composition leakage (single population).
# This committed table is the source; moved exploratory/ -> figures/data/
# 2026-08-26 (figures must not depend on exploratory/; producer moved to
# scripts/compare_bprime.py — BPRIME_RAW/BPRIME_META env-gated — which writes here).
BPRIME_FILE <- file.path(REPO, "figures", "data", "bprime_locus_fdr.tsv")
bprime <- read.delim(BPRIME_FILE, check.names = FALSE)
bprime <- bprime[bprime$design == "paired" & bprime$floor == 20, ]
stopifnot(nrow(bprime) > 0)
bp_m <- regmatches(bprime$peak_id, regexec("^(.*):([0-9]+)-([0-9]+)$", bprime$peak_id))
bprime$start <- as.integer(vapply(bp_m, function(x) x[3], ""))
bprime$end   <- as.integer(vapply(bp_m, function(x) x[4], ""))

# fail-fast: every featured peak must resolve to a (b') locus peak. ArchR peak
# boundaries can shift a few bp on a rerun → the exact-start star join
# (sg$start == FEAT$start) would silently drop a star. Stop loudly instead of
# quietly weakening the figure.
.feat_missing <- setdiff(FEAT$start, bprime$start)
if (length(.feat_missing))
    stop("FEAT peak start(s) absent from ", basename(BPRIME_FILE), " (peak-set drift?): ",
         paste(.feat_missing, collapse = ", "), " — re-derive FEAT from pitx2_locus_tad_peaks.tsv")

# read a lane's (b') paired L/R result -> data.frame(start, end, log2FoldChange, padj_locus)
read_sig <- function(spec) {
    if (is.na(spec$stage)) return(NULL)
    d <- bprime[bprime$stage == spec$stage & bprime$compartment == spec$comp, ]
    if (nrow(d) == 0) return(NULL)
    d[, c("start", "end", "log2FoldChange", "padj_locus")]
}

# ── locus peak set (ArchR-union, FAM241A-ext. ±500 kb; ticks + background mask) ──
tad <- read.delim(file.path(DATA, "pitx2_locus_tad_peaks.tsv"), check.names = FALSE)
peaks_all <- tad[, c("start", "end")]

# Display window: focal 40 kb (gene + 4 featured elements: a/b/c promoters + the +5.4 kb putative peak).
source(file.path(HERE, "_focal_window.R"))   # PITX2_FOCAL_WIN — shared with panel B (single source of truth; HERE-relative so panel_C stays cwd-independent)
WIN <- PITX2_FOCAL_WIN
peaks_win <- peaks_all[peaks_all$start >= WIN[1] & peaks_all$end <= WIN[2], ]
n_in_win <- nrow(peaks_win)
peak_note <- sprintf("%d ArchR-union peaks in view", n_in_win)
cat(sprintf("Window %d-%d (%.1f kb); %d peaks in view\n", WIN[1], WIN[2], diff(WIN) / 1e3, n_in_win))

# ── GC content over the window (the CpG island the PITX2c promoter sits on) ──
# Sliding-window GC% from the canonical duck reference (scripts/gc_content_pitx2_window.py,
# run on the cluster → ferried). The PITX2c promoter is a high-GC CpG island — the
# chromatin-memory substrate of the model (a memory-bearing, methylation-sensitive
# promoter), drawn as the top context track.
GC_TSV <- file.path(DATA, "pitx2_locus_gc_content.tsv")
if (!file.exists(GC_TSV))
    stop("missing ", GC_TSV, " — run scripts/run_gc_content_pitx2.sh on the cluster + ferry back")
gc <- read.delim(GC_TSV, check.names = FALSE)
gc <- gc[gc$pos >= WIN[1] & gc$pos <= WIN[2], ]
gc$kb <- gc$pos / 1e3
cat(sprintf("GC track: %d windows, %.0f-%.0f%% (island peak %.0f%%)\n",
            nrow(gc), min(gc$gc_pct), max(gc$gc_pct), max(gc$gc_pct)))

# ── load coverage ───────────────────────────────────────────────────────────
cat("Loading", basename(COV), "\n")
df <- read.delim(COV, check.names = FALSE)
# PSEUDOREPLICATION RECTIFICATION (2026-06-17): coverage is now per-(stage x
# compartment x GENOTYPE x side) — columns <stage>_<compartment>_<genotype>_<side>
# [_percell].  pooled_pc() averages the per-genotype _percell tracks so each
# embryo-HALF is weighted equally (was: cell-pooled, so high-cell embryos
# dominated the track SHAPE).  Significance STARS are unaffected (rectified g4f
# floor-20).
gpc <- grep("^HH[0-9_]+_(ventral|dorsal)_.+_[LR]_percell$", names(df), value = TRUE)
stopifnot(length(gpc) >= 12)  # >= 1 genotype-half per (stage x compartment x side)
cat(sprintf("  %d per-genotype-half _percell tracks (stage x compartment x genotype x side)\n",
            length(gpc)))

# ── HH13 (Wave-1) pooled L/R coverage — SEPARATE file ───────────────────────
# Pre-TES -> no D/V compartment split -> a single pooled inference-pool L vs R
# track.  Pools all HH13 genotype-halves (8 L / 7 R, singleton_rank-inclusive);
# UNDERPOWERED -> descriptive coverage shape only, NO DESeq2 stars (only 1L/2R are
# *firm*, and there is no g4f HH13 stat).  Rename to the HH13_pool_<side> scheme
# pooled_pc expects (compartment = "pool"); merge onto df by the shared pos grid.
HH13_COV <- file.path(DATA, "f4_hh13_perembryo_pitx2_coverage.tsv.gz")
if (!file.exists(HH13_COV)) stop("missing ", HH13_COV)
cat("Loading", basename(HH13_COV), "(HH13 Wave-1 pooled lane, per-embryo-half)\n")
hh13 <- read.delim(HH13_COV, check.names = FALSE)
stopifnot(nrow(hh13) == nrow(df), all(hh13$pos == df$pos))
# HH13 is pre-TES (no D/V split): columns HH13_<genotype>_<side>[_percell].
# Average the per-genotype _percell across embryo-halves -> the pooled lane
# (8 L / 7 R singleton_rank-inclusive halves; descriptive, no stars).
for (side in c("L", "R")) {
    cols <- grep(sprintf("^HH13_.+_%s_percell$", side), names(hh13), value = TRUE)
    stopifnot(length(cols) >= 1)
    df[[sprintf("HH13_pool_%s_percell", side)]] <- rowMeans(as.data.frame(hh13[cols]), na.rm = TRUE)
}

# smooth each track over the FULL series (avoids NA bands at window edges),
# then window.  Size factors use RAW (unsmoothed) per-cell over background.
ma <- function(x, w) as.numeric(stats::filter(x, rep(1 / w, w), sides = 2))
in_win <- df$pos >= WIN[1] & df$pos <= WIN[2]

# background mask = window positions NOT inside any called peak
pk_mask <- rep(FALSE, nrow(df))
for (i in seq_len(nrow(peaks_all))) pk_mask <- pk_mask | (df$pos >= peaks_all$start[i] & df$pos <= peaks_all$end[i])
bg_mask <- in_win & !pk_mask
cat(sprintf("Background positions (peak-free, in window): %s of %s\n",
            format(sum(bg_mask), big.mark = ","), format(sum(in_win), big.mark = ",")))

# ── 7 lanes: HH13 POOLED (Wave-1, pre-TES; underpowered/descriptive), then
#    ventral|dorsal split at HH17/16, HH21, HH25 (HH17/16 dorsal = Dorsal-biased
#    LPM, the earliest dorsal compartment; matches fig04A's 7 columns).  Each
#    split lane's coverage is cohort-matched to its OWN g4f (stage,compartment) stars.
lane_specs <- list(
    list(key = "HH13_pool", label = "HH13\n(pool)",
         srcs = list(c("HH13", "pool"))),
    list(key = "HH17/16_ventral", label = "HH17/16\nVentral (inference pool) \u2020", srcs = list(c("HH17_16", "ventral"))),
    list(key = "HH17/16_dorsal",  label = "HH17/16\nDorsal-biased LPM",        srcs = list(c("HH17_16", "dorsal"))),
    list(key = "HH21_ventral", label = "HH21\nVentral splanchnic", srcs = list(c("HH21", "ventral"))),
    list(key = "HH21_dorsal",  label = "HH21\nDorsal sub-epi.",    srcs = list(c("HH21", "dorsal"))),
    list(key = "HH25_ventral", label = "HH25\nVentral splanchnic", srcs = list(c("HH25", "ventral"))),
    list(key = "HH25_dorsal",  label = "HH25\nDorsal sub-epi.",    srcs = list(c("HH25", "dorsal")))
)
lane_keys   <- vapply(lane_specs, function(s) s$key, "")
lane_levels <- vapply(lane_specs, function(s) s$label, "")

# ── Cohort match: restrict the pileup to the g4f floor-20 DESeq2 cohort
#    (the same embryo-half census panel B's counts + fig04A's columns use), so the
#    displayed coverage is on a defined cohort rather than the f4 export's looser
#    per-half eligibility (which admits halves, e.g. HH17/16 ventral 19L/18R, that
#    never entered any fit).  The STARS are the (b') paired
#    within-population fit — at HH21/HH25 ventral its cohort is the same 8 embryos
#    (8 L / 8 R halves) as the g4f cohort shown.  HH13 is pre-pooled + descriptive
#    (no stars) → NOT filtered.
smeta_floor20 <- read.delim(file.path(DATA, "g4f_sample_meta_floor20.tsv"), check.names = FALSE)
cohort_pc <- unique(sprintf("%s_%s_%s_%s_percell",
                            gsub("/", "_", smeta_floor20$d2g_stage),
                            smeta_floor20$compartment,
                            smeta_floor20$unified_genotype_id,
                            smeta_floor20$inferred_side))
# Every scoped census half must have its own coverage column, and that column must
# hold at least the census nuclei (f4 halves include sub-floor subclusters the
# floor-20 census leaves out, so equality is not expected; fewer nuclei would mean the
# coverage and the census select different nuclei, e.g. one aliquot of a split-tube
# half).  Intersecting silently would draw a lane on a cohort that matches neither.
.census_scoped <- grep("^HH(17_16|21|25)_(ventral|dorsal)_", cohort_pc, value = TRUE)
.no_cov <- setdiff(.census_scoped, names(df))
if (length(.no_cov) > 0)
    stop(sprintf(paste0("panel C: %d g4f floor-20 census half(s) have no coverage column ",
                        "(coverage and census are on different side calls; rerun f4 ",
                        "on the current lr_inference_multiome.tsv): %s"),
                 length(.no_cov), paste(.no_cov, collapse = ", ")))
.f4n <- read.delim(COV_N, check.names = FALSE)
.cen <- smeta_floor20 |>
    mutate(col = sprintf("%s_%s_%s_%s", gsub("/", "_", d2g_stage), compartment,
                         unified_genotype_id, inferred_side)) |>
    filter(grepl("^HH(17_16|21|25)_", col)) |>
    left_join(.f4n |> select(group, f4_n = n_cells), by = c(col = "group"))
.short <- .cen |> filter(is.na(f4_n) | f4_n < n_cells)
if (nrow(.short) > 0)
    stop("panel C: f4 column has fewer nuclei than the g4f census for: ",
         paste(sprintf("%s (f4 %s < census %d)", .short$col, .short$f4_n, .short$n_cells), collapse = ", "))
.scoped_pc <- grep("^HH(17_16|21|25)_(ventral|dorsal)_.+_[LR]_percell$", names(df), value = TRUE)
cat(sprintf("Cohort match: %d/%d non-HH13 _percell tracks kept (g4f floor-20 census); %d dropped to keep the coverage on a defined cohort\n",
            length(intersect(.scoped_pc, cohort_pc)), length(.scoped_pc), length(setdiff(.scoped_pc, cohort_pc))))

# PER-EMBRYO-HALF averaged per-cell series for a lane-side: gather every
# genotype-half _percell track across the lane's sources and average them with
# EQUAL WEIGHT — each embryo-half counts once, regardless of cell count (the
# pseudoreplication rectification).  Non-HH13 sources are restricted to the g4f
# floor-20 cohort (above) so coverage and stars share one embryo set.  Returns 0
# for a side with no qualifying halves (a fully missing lane is caught loudly by
# the bgL>0/bgR>0 guard in the lane loop below).
lane_half_cols <- function(srcs, side) {
    half_cols <- character(0)
    for (s in srcs) {
        if (s[2] == "pool") {                                       # HH13 pre-TES pooled lane
            pat  <- sprintf("^%s_.+_%s_percell$", s[1], side)       #   (descriptive, no stars;
            cols <- grep(pat, names(df), value = TRUE)              #    NOT cohort-filtered)
        } else {                                                    # compartment x genotype
            pat  <- sprintf("^%s_%s_.+_%s_percell$", s[1], s[2], side)
            cols <- intersect(grep(pat, names(df), value = TRUE), cohort_pc)  # g4f cohort
        }
        half_cols <- c(half_cols, cols)
    }
    unique(half_cols)
}
pooled_pc <- function(srcs, side) {
    half_cols <- lane_half_cols(srcs, side)
    if (length(half_cols) == 0) return(numeric(nrow(df)))
    rowMeans(as.data.frame(df[half_cols]), na.rm = TRUE)
}
# embryo-halves behind a lane-side (the HH13 pool lane is pre-averaged in df, so
# count its per-genotype-half columns in hh13) — stdout record only
n_halves <- function(srcs, side) {
    if (all(vapply(srcs, function(s) s[2] == "pool", TRUE)))
        return(length(grep(sprintf("^HH13_.+_%s_percell$", side), names(hh13), value = TRUE)))
    length(lane_half_cols(srcs, side))
}
# PITX2c promoter (featured peak) — for the stdout record only (not drawn)
CPROM <- FEAT[FEAT$name == "PITX2c\nprom", c("start", "end")]
stopifnot(nrow(CPROM) == 1)
in_cprom <- df$pos >= CPROM$start & df$pos <= CPROM$end

rows <- list(); depth_tab <- list()
for (spec in lane_specs) {
    ln <- spec$key
    pcL <- pooled_pc(spec$srcs, "L"); pcR <- pooled_pc(spec$srcs, "R")
    bgL <- sum(pcL[bg_mask], na.rm = TRUE); bgR <- sum(pcR[bg_mask], na.rm = TRUE)
    stopifnot(bgL > 0, bgR > 0)  # empty-side lane → guard against NaN/Inf background scale
    locL <- sum(pcL[in_win], na.rm = TRUE); locR <- sum(pcR[in_win], na.rm = TRUE)
    Kref <- (bgL + bgR) / 2; sL <- Kref / bgL; sR <- Kref / bgR  # background-matched scale
    cpL <- sum(pcL[in_cprom], na.rm = TRUE); cpR <- sum(pcR[in_cprom], na.rm = TRUE)
    depth_tab[[ln]] <- data.frame(lane = ln, label = spec$label, bg_LR = bgL / bgR, locus_LR = locL / locR,
                                  cprom_LR_norm = (cpL / bgL) / (cpR / bgR),
                                  n_L = n_halves(spec$srcs, "L"), n_R = n_halves(spec$srcs, "R"))
    smL <- ma(pcL, BIN)[in_win]; smR <- ma(pcR, BIN)[in_win]; pw <- df$pos[in_win]
    rows[[ln]] <- rbind(
        data.frame(lane = ln, label = spec$label, side = "L", pos = pw, cov = smL * sL),
        data.frame(lane = ln, label = spec$label, side = "R", pos = pw, cov = smR * sR))
}
long <- do.call(rbind, rows); depth_tab <- do.call(rbind, depth_tab)
long$cov_m <- ifelse(long$side == "R", -long$cov, long$cov)  # mirror R below 0
long$kb <- long$pos / 1e3
long$lane_f <- factor(long$label, levels = lane_levels)
depth_tab$lane_f <- factor(depth_tab$label, levels = lane_levels)

cat("\nPer-lane L/R depth factors (raw, background vs whole locus), PITX2c-promoter L/R after\n",
    "background (depth) matching, and embryo-halves per side:\n", sep = "")
print(transform(depth_tab, bg_LR = round(bg_LR, 3), locus_LR = round(locus_LR, 3),
                cprom_LR_norm = round(cprom_LR_norm, 3))[, c("lane", "bg_LR", "locus_LR",
                                                             "cprom_LR_norm", "n_L", "n_R")],
      row.names = FALSE)

# ── per-lane y anchors (free_y) + lane annotations ──────────────────────────
lane_summ <- long %>% group_by(lane, lane_f) %>%
    summarise(top = max(cov_m, na.rm = TRUE), bot = min(cov_m, na.rm = TRUE), .groups = "drop")
# anchors: leave generous headroom above L for stars + info text (avoid collisions)
ANCHOR_TOP <- 1.58   # headroom factor above each lane's L-max (for depth/star text)
anchors <- bind_rows(
    transform(lane_summ, kb = WIN[1] / 1e3, cov_m = top * ANCHOR_TOP),
    transform(lane_summ, kb = WIN[1] / 1e3, cov_m = bot * 1.10))

# y-axis breaks: label 0 + each lane's rounded positive max (a.u.). With free_y +
# the anchor headroom, pretty_breaks() drops the top label on the low-amplitude
# lanes (HH17/16 pool, dorsal); this labels the data max on EVERY lane instead.
lane_y_breaks <- function(lims) {
    data_hi <- max(lims, na.rm = TRUE) / ANCHOR_TOP   # undo the top-anchor inflation
    out <- sort(unique(c(0, signif(data_hi, 1))))
    out[out >= min(lims) & out <= max(lims)]
}

# ── significance stars at featured peaks (per lane, from the (b') paired fit) ──
stars <- list()
for (ln in lane_keys) {
    sg <- read_sig(SIG[[ln]])
    if (is.null(sg)) next
    for (i in seq_len(nrow(FEAT))) {
        hit <- sg[sg$start == FEAT$start[i], ]
        if (nrow(hit) == 0) next
        lab <- star_lab(hit$padj_locus[1])
        if (lab == "") next
        ls <- lane_summ[lane_summ$lane == ln, ]
        stars[[length(stars) + 1]] <- data.frame(
            lane_f = ls$lane_f, kb = FEAT$mid[i] / 1e3, cov_m = ls$top * 1.26, lab = lab)
    }
}
stars <- if (length(stars)) do.call(rbind, stars) else
    data.frame(lane_f = factor(character(), levels = lane_levels), kb = numeric(), cov_m = numeric(), lab = character())

# ── peak / featured-band geometry ───────────────────────────────────────────
peaks_win$kbmid <- (peaks_win$start + peaks_win$end) / 2 / 1e3
feat_band <- FEAT %>% transmute(xmin = start / 1e3, xmax = end / 1e3,
                                fill = ifelse(module == "c", C_BAND, AB_BAND))

# ════════════════════════════════════════════════════════════════════════════
# TRACK 1 — ArchR-union peaks (grey ticks) + 4 featured (colored ticks + labels)
# ════════════════════════════════════════════════════════════════════════════
p_peaks <- ggplot() +
    geom_rect(data = feat_band, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = fill),
              inherit.aes = FALSE, alpha = 0.85) +
    scale_fill_identity() +
    geom_segment(data = peaks_win, aes(x = kbmid, xend = kbmid, y = 0, yend = 0.62),
                 color = FAINT, linewidth = 0.35) +
    geom_segment(data = FEAT, aes(x = mid / 1e3, xend = mid / 1e3, y = 0, yend = 1.0),
                 color = ifelse(FEAT$module == "c", C_BAND, AB_BAND),
                 linewidth = FEAT$tick_lw, linetype = FEAT$tick_lty) +
    # promoters: bold dark labels (foregrounded); +5.4 kb: italic muted (putative)
    geom_text(data = FEAT, aes(x = mid / 1e3, y = 1.25, label = name),
              size = 3.6, lineheight = 0.82, vjust = 0, hjust = FEAT$hjust,
              fontface = FEAT$face, color = ifelse(FEAT$kind == "prom", INK, MUT)) +
    # upper-right whitespace, ABOVE the ticks (clears the foregrounded a/b/c ticks+labels
    # on the left and sits clear of the grey ArchR ticks below)
    # (PNAS font pass: at print size the note sits just above the grey ticks on
    # the right, clear of the featured-peak labels above)
    annotate("text", x = WIN[2] / 1e3, y = 0.68, hjust = 1, vjust = 0,
             label = peak_note, size = 3.4, color = MUT) +
    scale_x_continuous(limits = WIN / 1e3, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.95), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +   # right-aligned PITX2a label may breathe past the left edge
    theme_void() + theme(plot.margin = margin(2, 10, 0, 5), legend.position = "none")

# ════════════════════════════════════════════════════════════════════════════
# TRACK 2 — the 6 mirror pileup lanes
# ════════════════════════════════════════════════════════════════════════════
p_tracks <- ggplot(long, aes(kb, cov_m, fill = side)) +
    # featured-peak bands span all lanes (eye-guide down through the stack).
    geom_rect(data = subset(feat_band, fill == C_BAND), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = C_BAND, alpha = 0.10) +
    geom_rect(data = subset(feat_band, fill == AB_BAND), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = AB_BAND, alpha = 0.10) +
    geom_area(position = "identity", alpha = 0.92, linewidth = 0) +
    geom_hline(yintercept = 0, color = INK, linewidth = 0.3) +
    geom_blank(data = anchors, aes(kb, cov_m), inherit.aes = FALSE) +
    # (no in-lane text, 2026-10-14: depth L/R is printed to stdout + PROVENANCE;
    #  underpowered / side-defining status is in the legend; HH17/16 ventral carries a dagger)
    # DESeq2 significance stars at featured peaks
    geom_text(data = stars, aes(x = kb, y = cov_m, label = lab),
              inherit.aes = FALSE, size = 4.4, color = INK, vjust = 0, fontface = "bold") +
    # lane names as one-line strips ABOVE each lane (PNAS font pass): left-side
    # strips at print size ate half the track width
    facet_wrap(~ lane_f, ncol = 1, scales = "free_y", strip.position = "top",
               labeller = as_labeller(function(x) gsub("\n", " ", x))) +
    scale_fill_manual(values = c(L = L_BLUE, R = R_RED),
                      labels = c(L = "Left (blue)", R = "Right (red)"), name = NULL) +
    scale_x_continuous(limits = WIN / 1e3, expand = c(0, 0)) +
    scale_y_continuous(breaks = lane_y_breaks) +
    labs(x = NULL, y = "Accessibility (a.u.)") +
    # PNAS final font pass (2026-10): sized for the tall variant, which Fig 5
    # prints at x0.653 (17.8 cm) -> >= ~7 pt printed.
    theme_classic(base_size = 12) +
    theme(
        strip.placement = "outside",
        strip.text.x = element_text(hjust = 0, size = 11, face = "bold", color = INK,
                                    margin = margin(3, 0, 1, 0)),
        strip.background = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10.5, color = MUT),
        legend.position = "top", legend.justification = "right",
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 11),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 10, 0, 5))

# ════════════════════════════════════════════════════════════════════════════
# TRACK 3 — PITX2 gene model (3 isoforms; shared helper) + neighbour genes
# ════════════════════════════════════════════════════════════════════════════
# Neighbour ncRNA models (crarc_zju1_PITX2 GTF; docs/pitx2_locus_annotation.md §2).
# All uncharacterized LOCs (no symbol/ortholog); identities POSITIONAL (PITX2-AS =
# PITX2 natural antisense; PANCR-like = PITX2-adjacent ncRNA analog, cf. human PANCR).
# The +5.4 kb element lies within LOC106017192 / LOC113843539 — drawn to show the overlap.
NEIGH <- list(
    list(name = "LOC113843538\n(PITX2-AS)",   strand = "-", kind = "ncRNA",
         exons = list(c(3963664, 3963765), c(3963991, 3964114))),
    list(name = "LOC113843539\n(PITX2-AS)",   strand = "-", kind = "ncRNA",
         exons = list(c(3973586, 3973815), c(3976331, 3976474), c(3976582, 3976858))),
    list(name = "LOC106017192\n(PANCR-like)", strand = "+", kind = "ncRNA",
         exons = list(c(3976875, 3977088), c(4032001, 4032085), c(4032251, 4033262))),
    list(name = "LOC106016804", strand = "-", kind = "ncRNA",  exons = list(c(3891849, 3954422))),
    list(name = "LOC101795212", strand = "-", kind = "coding", exons = list(c(4003281, 4037227)))
)

p_gene <- ggplot() +
    geom_rect(data = feat_band, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -0.95, ymax = 4.2),
              fill = rep(feat_band$fill, 1), alpha = 0.10)
p_gene <- draw_pitx2_gene_model(p_gene, label_size = 3.8)
EXH <- 0.13
for (g in NEIGH) {
    glo <- min(vapply(g$exons, function(e) e[1], 0)); ghi <- max(vapply(g$exons, function(e) e[2], 0))
    if (ghi < WIN[1] || glo > WIN[2]) next
    col <- if (g$kind == "coding") INK else MUT
    il <- max(glo, WIN[1]); ir <- min(ghi, WIN[2])
    p_gene <- p_gene + annotate("segment", x = il / 1e3, xend = ir / 1e3, y = 0, yend = 0,
                                color = col, linewidth = 0.4)              # intron line
    for (e in g$exons) {                                                   # exon boxes (in-window)
        if (e[2] < WIN[1] || e[1] > WIN[2]) next
        p_gene <- p_gene + annotate("rect", xmin = max(e[1], WIN[1]) / 1e3, xmax = min(e[2], WIN[2]) / 1e3,
                                    ymin = -EXH, ymax = EXH, fill = col, color = NA)
    }
    tss <- if (g$strand == "+") glo else ghi                               # TSS strand arrow
    if (tss >= WIN[1] && tss <= WIN[2]) {
        adir <- if (g$strand == "+") 1 else -1
        p_gene <- p_gene +
            annotate("segment", x = tss / 1e3, xend = tss / 1e3, y = EXH, yend = EXH + 0.17,
                     color = col, linewidth = 0.4) +
            annotate("segment", x = tss / 1e3, xend = (tss + adir * 300) / 1e3, y = EXH + 0.17, yend = EXH + 0.17,
                     color = col, linewidth = 0.4, arrow = arrow(length = unit(1, "mm"), type = "closed"))
    }
    p_gene <- p_gene + annotate("text", x = (il + ir) / 2 / 1e3, y = -0.52, label = g$name,
                                size = 3.0, color = col, fontface = "italic", lineheight = 0.8)
}
p_gene <- p_gene +
    scale_x_continuous(name = "NC_051775.1 position (kb)", limits = WIN / 1e3, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.95, 4.4)) +
    coord_cartesian(clip = "off") +   # right-aligned PITX2a isoform label sits at the window edge (2026-09-07: "P" was clipped)
    theme_void(base_size = 12) +
    theme(axis.title.x = element_text(size = 12, color = INK, margin = margin(t = 2)),
          axis.text.x = element_text(size = 10.8, color = INK),
          axis.ticks.x = element_line(linewidth = 0.3), axis.ticks.length.x = unit(1.4, "mm"),
          axis.line.x = element_line(linewidth = 0.3),
          plot.margin = margin(0, 10, 4, 5))

# ════════════════════════════════════════════════════════════════════════════
# TOP TRACK — GC% (the CpG island the PITX2c promoter sits on)
# ════════════════════════════════════════════════════════════════════════════
gc_peak_kb <- gc$kb[which.max(gc$gc_pct)]
p_gc <- ggplot(gc, aes(kb, gc_pct)) +
    # featured-peak bands span the track (eye-guide down through the stack)
    geom_rect(data = subset(feat_band, fill == C_BAND), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = C_BAND, alpha = 0.10) +
    geom_rect(data = subset(feat_band, fill == AB_BAND), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = AB_BAND, alpha = 0.10) +
    geom_area(fill = INK, alpha = 0.22, linewidth = 0) +
    geom_line(color = INK, linewidth = 0.35) +
    # CpG-island GC criterion (Gardiner-Garden & Frommer: GC ≥ 50%)
    geom_hline(yintercept = 50, color = MUT, linetype = "dashed", linewidth = 0.3) +
    # PNAS final files: the threshold label sits in the clear band above the
    # trace at the right end, keyed by a short dashed stroke (was on the trace)
    annotate("segment", x = WIN[2] / 1e3 - 0.41 * diff(WIN) / 1e3,
             xend = WIN[2] / 1e3 - 0.355 * diff(WIN) / 1e3, y = 90, yend = 90,
             color = MUT, linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = WIN[2] / 1e3 - 0.01 * diff(WIN) / 1e3, y = 90,
             label = "CpG-island GC (50%)", hjust = 1, vjust = 0.5, size = 3.5,
             color = MUT) +
    annotate("text", x = gc_peak_kb, y = 99, label = "CpG island", hjust = 0.5, vjust = 1,
             size = 3.8, color = INK, fontface = "italic") +
    scale_x_continuous(limits = WIN / 1e3, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
    labs(x = NULL, y = "GC %") +
    theme_classic(base_size = 12) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10.5, color = MUT),
          # top room so the composite letter "C." clears the "100" tick
          plot.margin = margin(16, 10, 0, 5))

# ── stack: GC / peaks / lanes / gene model (no panel title — see LEGEND.md) ──
combined <- p_gc / p_peaks / p_tracks / p_gene +
    plot_layout(heights = c(1.7, 1.8, 12.4, 3.1))   # tracks: 7 lanes (was 11.2 for 6)

ggsave(OUT_PILEUP, combined, width = 10.0, height = 11.1, dpi = 300, bg = "white")   # height: 7 lanes (was 10.0 for 6)
cat("wrote", OUT_PILEUP, "\n")
# tall variant for the manuscript Fig 5 side-column layout (A/B/D left, C right,
# 2026-08-24): same content, narrow-tall canvas so the lanes stretch
# vertically without font shrink at print
OUT_TALL <- file.path(HERE, "panel_C_pileup_tall.png")
ggsave(OUT_TALL, combined, width = 5.5, height = 13.0, dpi = 300, bg = "white")
cat("wrote", OUT_TALL, "\n")

# ── featured-peak (b') paired summary to stdout (the record) ────────────────
cat("\nFeatured-peak (b') paired L/R (log2FC / padj_locus), floor 20:\n")
for (ln in lane_keys) {
    sg <- read_sig(SIG[[ln]]); if (is.null(sg)) next
    cat(sprintf("-- %s --\n", ln))
    for (i in seq_len(nrow(FEAT))) {
        hit <- sg[sg$start == FEAT$start[i], ]
        if (nrow(hit) == 0) next
        cat(sprintf("   %-12s log2FC %+.2f  padj_locus %.2g %s\n",
                    gsub("\n", " ", FEAT$name[i]), hit$log2FoldChange[1], hit$padj_locus[1],
                    star_lab(hit$padj_locus[1])))
    }
}
