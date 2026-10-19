#!/usr/bin/env Rscript
# archr_build_peaks.R — chick galGal6a ArchR peak calling (snATAC STEP 1b).
#
# Faithful port of the duck pipeline/archr_build_peaks.R to galGal6a / GRCg6a.
# Mirrors the duck grouping decision (confirmed adequate by the STEP-0 power census
# h1b_power_census.py): groupBy = peak_pool = (population x stage), pooled across A+B
# genotypes for peak DISCOVERY (peaks are genomic regions; group/side is irrelevant here),
# MIN_POOL_CELLS=500, reproducibility="2" (peak in >=2 pseudobulk replicates, ArchR-default
# Sample/library-level = production; genotype-level tested as a sensitivity check, focal result
# invariant — docs/archr_peak_reproducibility_decision.md, not distributed), 501 bp fixed
# width.  8 of 16 pools clear the floor (the big HH25 pools + HH13 Early LPM/Somatic +
# HH15-16 Somatic/Dorsal-biased); the union peak set covers the PITX2 locus (chr4 ~57.79 Mb)
# via the HH25 pools.  The L/R DA (downstream) then runs only on Group B firm-oriented cells.
#
# ALL duck gotchas ported (see the duck script for full pilot-job rationale):
#   * custom chromSizes GRanges from chroms.txt (bypasses filterChrGR underscore/standard
#     stripping); main_seqs = chr1-28,30-33,chrZ,chrW (34 chroms; NO chrM, NO *_random).
#   * genes derived by aggregating transcripts per gene_id (NA gene-chrom crash guard) +
#     assert every main_seq has >=1 gene.
#   * GeneScoreMatParams useGeneBoundaries=FALSE (genes near contig origin).
#   * nChunk=1 (avoid 5x per-chrom chunk overhead); addArchRLocking(TRUE) (H5 link race).
#   * effective genome size = non-N ACGT over main_seqs (MACS3 -g for non-standard genome).
#
# Outputs (OUT_ROOT = chick_demux/atac_archr/):
#   ArrowFiles/<LIB>.arrow, proj_chick/, proj_chick_callable/,
#   peaks/union_peaks.bed, peaks/per_pool_peaks.tsv, peaks/reproducibility_metadata.tsv
#
# Run: sbatch pipeline/chick/archr_build_peaks.sh   (prereq: archr_build_bsgenome.R done)

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(ArchR)
    library(BSgenome.Ggallus.crarc.galGal6a)
    library(GenomicRanges)
    library(rtracklayer)
    library(Rsamtools)
    library(data.table)
    library(optparse)
})

option_list <- list(
    make_option("--libraries", type = "character", default = NULL,
                help = "Comma-separated subset (pilot). Default: ARC01,ARC02,ARC09,ARC10"),
    make_option("--threads", type = "integer", default = 4),
    make_option("--out-root", type = "character",
                default = paste0(SCRATCH_USER, "/chick_demux/atac_archr")),
    make_option("--macs3-path", type = "character",
                default = paste0(HOME_ROOT, "/.conda/envs/snapatac2/bin/macs3")),
    make_option("--ref-dir", type = "character",
                default = paste0(SCRATCH_USER, "/chick_demux/ref")),
    make_option("--meta", type = "character",
                default = paste0(SCRATCH_USER, "/chick_demux/atac_archr/chick_archr_coupled_metadata.tsv"))
)
opt <- parse_args(OptionParser(option_list = option_list))

OUT_ROOT  <- opt[["out-root"]]
ARROW_DIR <- file.path(OUT_ROOT, "ArrowFiles")
PROJ_DIR  <- file.path(OUT_ROOT, "proj_chick")
PEAK_DIR  <- file.path(OUT_ROOT, "peaks")
META_TSV  <- opt$meta
REF_DIR   <- opt[["ref-dir"]]
GTF_GZ    <- file.path(REF_DIR, "genes.gtf.gz")
CHROMS_TXT <- file.path(REF_DIR, "chroms.txt")
dir.create(ARROW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PEAK_DIR,  recursive = TRUE, showWarnings = FALSE)

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── library -> fragments (chick_demux/atac_qc/<lib>/outs, staged by hQ0) ───────
FRAG_ROOT <- paste0(SCRATCH_USER, "/chick_demux/atac_qc")
ARC_LIBS <- list(
    ARC01 = file.path(FRAG_ROOT, "ARC01/outs/atac_fragments.tsv.gz"),
    ARC02 = file.path(FRAG_ROOT, "ARC02/outs/atac_fragments.tsv.gz"),
    ARC09 = file.path(FRAG_ROOT, "ARC09/outs/atac_fragments.tsv.gz"),
    ARC10 = file.path(FRAG_ROOT, "ARC10/outs/atac_fragments.tsv.gz")
)
if (!is.null(opt$libraries)) {
    sel <- strsplit(opt$libraries, ",", fixed = TRUE)[[1]]
    miss <- setdiff(sel, names(ARC_LIBS))
    if (length(miss)) stop("Unknown library(s): ", paste(miss, collapse = ", "))
    ARC_LIBS <- ARC_LIBS[sel]
    ts("pilot mode — libraries: ", paste(names(ARC_LIBS), collapse = ", "))
}

SKIP_POOLS <- character(0)   # no gut-confound pool to mask in chick (duck masked Pharyngeal endoderm)
MIN_POOL_CELLS <- 500L
MIN_TSS   <- 2.0
MIN_FRAGS <- 1000L

# ── RNG seed ──────────────────────────────────────────────────────────────────
# Unseeded until 2026-08-06; see the long note in pipeline/archr_build_peaks.R.
# addGroupCoverages() subsamples cells into the reproducibility="2" replicates, so
# an unseeded rebuild yields a different union -- measured on chick: 103,836 vs
# 112,799 peaks (-7.9%), PITX2 +/-500 kb 69 vs 86 (-19.8%), with the 2 PITX2c-promoter
# peaks surviving (per-embryo promoter dFRiP r = 1.000, 15/17 embryos unchanged) --
# docs/archr_peak_reproducibility_decision.md (not distributed).
# The shipped 112,799-peak chick set is the artifact of record; this seed is for
# FUTURE builds only. Reproducible at threads = 1 (mclapply workers are unseeded).
ARCHR_SEED <- 1L
set.seed(ARCHR_SEED)

# ── ArchR genome + gene annotations ───────────────────────────────────────────
ts("building genomeAnnotation + geneAnnotation from BSgenome + GTF")
ts("RNG seed = ", ARCHR_SEED, " (threads = ", opt$threads,
   if (opt$threads > 1) "; >1 thread: mclapply worker streams are unseeded" else "", ")")
addArchRThreads(threads = opt$threads, force = TRUE)
addArchRLocking(locking = TRUE)   # serialise H5 link inserts (intra-library chunk race)

bsg <- BSgenome.Ggallus.crarc.galGal6a

# main chroms = chroms.txt (chr1-28,30-33,chrZ,chrW); excludes chrM + ~428 *_random
# scaffolds.  Restricting ArchR to these keeps the Arrow build tractable and PITX2 is on
# chr4 — all fig biology is on main chroms.  chroms.txt is comma/space/newline-separated
# (the crarc reference packs the names on one line) — split on any of those, like the h1
# scaffold's txt.replace("\n",",").split(",").
chroms_raw <- paste(readLines(CHROMS_TXT), collapse = ",")
main_seqs <- trimws(unlist(strsplit(chroms_raw, "[,[:space:]]+")))
main_seqs <- main_seqs[nzchar(main_seqs)]
ts("parsed ", length(main_seqs), " names from chroms.txt (e.g. ",
   paste(head(main_seqs, 4), collapse = ", "), ")")
main_seqs <- intersect(main_seqs, seqnames(bsg))
ts("main_seqs (in BSgenome): ", length(main_seqs), " chroms (e.g. ",
   paste(head(main_seqs), collapse = ", "), ")")
stopifnot("chr4" %in% main_seqs)   # PITX2 locus chrom

seqlens <- seqlengths(bsg)
main_chromSizes <- GRanges(
    seqnames = main_seqs,
    ranges   = IRanges(start = 1L, end = seqlens[main_seqs]),
    strand   = "*",
    seqlengths = seqlens[main_seqs]
)
genomeAnnotation <- createGenomeAnnotation(
    genome = bsg, chromSizes = main_chromSizes,
    blacklist = GRanges(), filter = FALSE)

ts("parsing GTF: ", GTF_GZ)
gtf <- rtracklayer::import(GTF_GZ, format = "gtf")
gtf <- keepSeqlevels(gtf, intersect(seqlevels(gtf), main_seqs), pruning.mode = "coarse")
ts("  GTF restricted to ", length(seqlevels(gtf)), " main chroms — ", length(gtf), " features")
exons_gr <- gtf[which(gtf$type == "exon")]
tx_gr    <- gtf[which(gtf$type == "transcript")]
native_genes_gr <- gtf[which(gtf$type == "gene")]
if (!"gene_id" %in% colnames(mcols(exons_gr))) stop("gene_id missing from GTF exons")
if (!"gene_id" %in% colnames(mcols(tx_gr)))    stop("gene_id missing from GTF transcripts")

# Derive gene-level GRanges by aggregating transcripts per gene_id (NA gene-chrom guard).
ts("aggregating transcripts -> genes (one row per gene_id) ...")
tx_dt <- as.data.table(tx_gr)
genes_dt <- tx_dt[, .(seqnames = as.character(seqnames[1L]),
                      start = min(start), end = max(end),
                      strand = as.character(strand[1L])), by = gene_id]
ts("  ", nrow(genes_dt), " genes from ", length(tx_gr), " transcripts; chroms covered: ",
   length(unique(genes_dt$seqnames)), " (target ", length(main_seqs), ")")

sym_lookup <- character(0)
if (length(native_genes_gr) > 0) {
    cand <- intersect(c("gene_name", "gene"), colnames(mcols(native_genes_gr)))
    if (length(cand)) {
        ng_sym <- rep(NA_character_, length(native_genes_gr))
        for (cc in cand) {
            v <- mcols(native_genes_gr)[[cc]]
            ng_sym <- ifelse(is.na(ng_sym) | ng_sym == "", v, ng_sym)
        }
        sym_lookup <- setNames(ng_sym, native_genes_gr$gene_id)
    } else sym_lookup <- setNames(native_genes_gr$gene_id, native_genes_gr$gene_id)
}
resolved_sym <- if (length(sym_lookup)) sym_lookup[genes_dt$gene_id] else NA_character_
resolved_sym <- ifelse(is.na(resolved_sym) | resolved_sym == "", genes_dt$gene_id, resolved_sym)

genes_for_archr <- GRanges(
    seqnames = genes_dt$seqnames,
    ranges   = IRanges(start = genes_dt$start, end = genes_dt$end),
    strand   = genes_dt$strand, gene_id = genes_dt$gene_id,
    symbol   = as.character(resolved_sym))
missing_chroms <- setdiff(main_seqs, unique(as.character(seqnames(genes_for_archr))))
if (length(missing_chroms)) stop("genes_for_archr missing chroms: ",
                                 paste(missing_chroms, collapse = ", "),
                                 " — .fastFeatureCounts will crash. Investigate GTF.")
ts("genes_for_archr: ", length(genes_for_archr), " ranges; all main_seqs covered")

tss_gr <- GRanges(
    seqnames = seqnames(tx_gr),
    ranges = IRanges(start = ifelse(as.character(strand(tx_gr)) == "-",
                                    end(tx_gr), start(tx_gr)), width = 1L),
    strand = strand(tx_gr), tx_id = tx_gr$transcript_id,
    tx_name = tx_gr$transcript_id, gene_id = tx_gr$gene_id)

exon_cols <- intersect(c("gene_name", "gene"), colnames(mcols(exons_gr)))
if (length(exon_cols) == 0) {
    exon_symbol <- exons_gr$gene_id
} else {
    exon_symbol <- rep(NA_character_, length(exons_gr))
    for (cc in exon_cols) {
        v <- mcols(exons_gr)[[cc]]
        exon_symbol <- ifelse(is.na(exon_symbol) | exon_symbol == "", v, exon_symbol)
    }
    exon_symbol <- ifelse(is.na(exon_symbol) | exon_symbol == "", exons_gr$gene_id, exon_symbol)
    exon_symbol <- as.character(exon_symbol)
}
geneAnnotation <- createGeneAnnotation(
    TSS = tss_gr,
    exons = GRanges(seqnames = seqnames(exons_gr), ranges = ranges(exons_gr),
                    strand = strand(exons_gr), gene_id = exons_gr$gene_id, symbol = exon_symbol),
    genes = genes_for_archr, annoStyle = "GENCODE")
ts("  geneAnnotation genes=", length(geneAnnotation$genes), " TSS=", length(geneAnnotation$TSS))

# ── load metadata sidecar (barcode -> peak_pool) ──────────────────────────────
ts("loading metadata: ", META_TSV)
meta <- fread(META_TSV, sep = "\t", header = TRUE, na.strings = "")
ts("  rows: ", nrow(meta), "  peak_pool NA: ", sum(is.na(meta$peak_pool)))
# ArchR cell barcode = "<sampleName>#<cb_fragments>"; sidecar already provides archr_cb.
setkey(meta, archr_cb)

# ── 1. Arrow files ────────────────────────────────────────────────────────────
setwd(OUT_ROOT)
arrow_paths <- character(0)
need_build  <- character(0)
for (lib in names(ARC_LIBS)) {
    cand <- c(file.path(OUT_ROOT, paste0(lib, ".arrow")),
              file.path(ARROW_DIR, paste0(lib, ".arrow")))
    found <- cand[file.exists(cand) & vapply(cand, function(p)
        if (file.exists(p)) file.info(p)$size > 1e6 else FALSE, logical(1))]
    if (length(found)) { ts("[skip] Arrow exists — ", found[1]); arrow_paths[lib] <- found[1]; next }
    if (!file.exists(ARC_LIBS[[lib]])) stop("missing fragments for ", lib, ": ", ARC_LIBS[[lib]])
    need_build <- c(need_build, lib)
}
if (length(need_build) > 0) {
    ts("createArrowFiles (parallel across ", length(need_build), " libs, threads=", opt$threads, ")")
    input_files <- vapply(need_build, function(lib) ARC_LIBS[[lib]], character(1))
    names(input_files) <- need_build
    arrs <- createArrowFiles(
        inputFiles = input_files, sampleNames = need_build, outputNames = need_build,
        minTSS = MIN_TSS, minFrags = MIN_FRAGS,
        addTileMat = TRUE, addGeneScoreMat = TRUE,
        GeneScoreMatParams = list(useGeneBoundaries = FALSE),
        genomeAnnotation = genomeAnnotation, geneAnnotation = geneAnnotation,
        force = FALSE, threads = opt$threads, nChunk = 1L)
    for (i in seq_along(arrs)) arrow_paths[need_build[i]] <- arrs[i]
}
ts("Arrow files ready (", length(arrow_paths), " libs)")

# ── 2. ArchRProject + peak_pool attach ────────────────────────────────────────
ts("building ArchRProject")
if (dir.exists(PROJ_DIR)) {
    proj <- loadArchRProject(path = PROJ_DIR, showLogo = FALSE)
} else {
    proj <- ArchRProject(ArrowFiles = arrow_paths, outputDirectory = PROJ_DIR,
                         copyArrows = TRUE, genomeAnnotation = genomeAnnotation,
                         geneAnnotation = geneAnnotation, showLogo = FALSE)
}
cell_ids <- rownames(getCellColData(proj))
ts("  ArchRProject cells (post-QC): ", length(cell_ids))

attach_df <- meta[J(cell_ids), on = "archr_cb",
                  .(cb = archr_cb, peak_pool, unified_genotype_id, population, stage, side, group)]
stopifnot(identical(attach_df$cb, cell_ids))
for (col in c("peak_pool", "unified_genotype_id", "population", "stage", "side", "group")) {
    proj <- addCellColData(ArchRProj = proj, data = attach_df[[col]], cells = cell_ids,
                           name = col, force = TRUE)
}
cc <- getCellColData(proj)
ts("  cells with non-NA peak_pool: ", sum(!is.na(cc$peak_pool)),
   " (coupling rate ", round(100 * sum(!is.na(cc$peak_pool)) / length(cell_ids), 1), "%)")

# Save the FULL project (with metadata attached) so downstream steps (locus PeakMatrix on
# the sub-floor focal cells) can reload all coupled cells — not just the callable subset.
ts("saveArchRProject (full) — ", PROJ_DIR)
saveArchRProject(ArchRProj = proj, outputDirectory = PROJ_DIR, load = FALSE)

valid <- !is.na(cc$peak_pool) & !(cc$peak_pool %in% SKIP_POOLS)
pool_counts  <- table(cc$peak_pool[valid])
callable     <- names(pool_counts)[pool_counts >= MIN_POOL_CELLS]
subthreshold <- names(pool_counts)[pool_counts <  MIN_POOL_CELLS]
ts("  callable pools (>= ", MIN_POOL_CELLS, "): ", length(callable))
for (p in callable) ts("    [call] ", p, " (", pool_counts[[p]], ")")
for (p in subthreshold) ts("    [skip<floor] ", p, " (", pool_counts[[p]], ")")

keep_cells <- cell_ids[!is.na(cc$peak_pool) & cc$peak_pool %in% callable]
ts("  cells feeding peak calling: ", length(keep_cells))
proj_call <- subsetArchRProject(ArchRProj = proj, cells = keep_cells,
                                outputDirectory = file.path(OUT_ROOT, "proj_chick_callable"),
                                dropCells = FALSE, force = TRUE)

# ── 3. group coverages + reproducible peak set ────────────────────────────────
ts("addGroupCoverages(groupBy='peak_pool')")
proj_call <- addGroupCoverages(ArchRProj = proj_call, groupBy = "peak_pool",  # ArchR-default Sample/
                               minCells = MIN_POOL_CELLS, maxCells = 10000L,  # library-level repro = production
                               force = TRUE, threads = opt$threads)

ts("computing effective genome size (non-N ACGT over main_seqs) ...")
genome_size_effective <- sum(vapply(main_seqs, function(seq)
    sum(letterFrequency(bsg[[seq]], letters = "ACGT", as.prob = FALSE)), numeric(1)))
ts("  effective genome size: ", format(genome_size_effective, big.mark = ","), " bp")

ts("addReproduciblePeakSet(reproducibility='2', MACS3=", opt[["macs3-path"]], ")")
proj_call <- addReproduciblePeakSet(
    ArchRProj = proj_call, groupBy = "peak_pool",
    pathToMacs2 = opt[["macs3-path"]], genomeSize = genome_size_effective,
    reproducibility = "2", peaksPerCell = 500L, maxPeaks = 150000L,
    minCells = MIN_POOL_CELLS, extendSummits = 250L, cutOff = 0.1,
    force = TRUE, threads = opt$threads)
peak_set <- getPeakSet(proj_call)
ts("peak set size: ", length(peak_set), " peaks")

# ── 4. exports ────────────────────────────────────────────────────────────────
union_bed <- file.path(PEAK_DIR, "union_peaks.bed")
export(peak_set, con = union_bed, format = "BED")
ts("wrote union BED: ", union_bed)

pool_rds_dir <- file.path(OUT_ROOT, "proj_chick_callable", "PeakCalls")
pool_files <- list.files(pool_rds_dir, pattern = "-reproduciblePeaks\\.gr\\.rds$",
                         full.names = TRUE, recursive = TRUE)
if (length(pool_files) == 0) {
    pool_files <- list.files(file.path(PROJ_DIR, "PeakCalls"),
                             pattern = "-reproduciblePeaks\\.gr\\.rds$",
                             full.names = TRUE, recursive = TRUE)
}
archr_sanitise <- function(x) gsub("[^a-zA-Z0-9._]", ".", x)
sanitised_to_orig <- setNames(names(pool_counts), archr_sanitise(names(pool_counts)))
per_pool_rows <- list(); repro_rows <- list()
for (f in pool_files) {
    pool_san <- sub("-reproduciblePeaks\\.gr\\.rds$", "", basename(f))
    pool <- sanitised_to_orig[[pool_san]]
    if (is.null(pool) || is.na(pool)) { ts("WARNING: unmapped pool '", pool_san, "' — skip"); next }
    gr <- readRDS(f)
    repro_rows[[length(repro_rows) + 1]] <- data.frame(
        peak_pool = pool, n_cells = as.integer(pool_counts[[pool]]),
        # n_replicates = genotype cardinality (informational); the peak reproducibility filter
        # itself is library-level (ArchR-default sampleLabels="Sample"; see header + decision doc).
        n_replicates = length(unique(cc$unified_genotype_id[cc$peak_pool == pool])),
        n_peaks = length(gr), stringsAsFactors = FALSE)
    per_pool_rows[[pool]] <- data.frame(
        peak_pool = pool, chrom = as.character(seqnames(gr)),
        start = start(gr) - 1L, end = end(gr), stringsAsFactors = FALSE)
}
if (length(repro_rows)) {
    write.table(do.call(rbind, repro_rows), file.path(PEAK_DIR, "reproducibility_metadata.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(do.call(rbind, per_pool_rows), file.path(PEAK_DIR, "per_pool_peaks.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    ts("wrote reproducibility_metadata.tsv + per_pool_peaks.tsv")
} else warning("no per-pool peak RDS under ", pool_rds_dir)

ts("saveArchRProject — ", file.path(OUT_ROOT, "proj_chick_callable"))
saveArchRProject(ArchRProj = proj_call,
                 outputDirectory = file.path(OUT_ROOT, "proj_chick_callable"), load = FALSE)
ts("=== archr_build_peaks.R done ===")
