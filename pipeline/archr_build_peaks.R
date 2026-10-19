#!/usr/bin/env Rscript
# archr_build_peaks.R — ArchR track A peak calling (plan §3.2 / §3.3).
#
# Parallel to the SnapATAC2 f2_peak_calling.py track.  Produces a unified
# peak set via ArchR's addReproduciblePeakSet with:
#   * groupBy   = "peak_pool"           (d2g_subcluster × d2g_stage; §3.3)
#   * reproducibility = "2"             (peak in ≥2 pseudobulk replicates; ArchR-default
#                                         Sample/library-level = the PRODUCTION set.  A
#                                         genotype-level variant was tested as a sensitivity
#                                         check, focal PITX2c result invariant — see
#                                         docs/archr_peak_reproducibility_decision.md, not distributed)
#   * peakMethod = "Macs2"              (ArchR's MACS wrapper; we point
#                                         pathToMacs2 at MACS3 binary in
#                                         the snapatac2 conda env)
#   * fixedWidth = 501 bp (extendSummits = 250)
#
# Granularity: peak calls per peak_pool; the 22 ≥500-cell pools call
# here, the 7 sub-500 LPM pools (plan §3.4) are handed to f2b's gkmQC
# rescue in the SnapATAC2 track — not an ArchR concern.
#
# Outputs (scratch phase2/atac_archr/):
#   ArrowFiles/<LIB>.arrow              — per-library Arrow files
#   proj_fig03/                         — ArchRProject (with obs, coverages, peak set)
#   peaks/union_peaks.bed               — merged 501 bp peak set
#   peaks/per_pool_peaks.tsv            — long format (peak_pool, chrom, start, end)
#   peaks/reproducibility_metadata.tsv  — (peak_pool, n_cells, n_replicates, n_peaks)
#
# Run on an o2 compute node via pipeline/run_archr_peaks.sh.  Pilot mode
# (--libraries ARC07[,ARC03-2]) exercises the full pipeline on 1-2
# libraries; default runs the full 12-library cohort.

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
REFERENCE_ROOT <- .duck_env("DUCK_REFERENCE_ROOT", "/n/groups/tabin/chlee/singlecell/reference")
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(ArchR)
    library(BSgenome.Anasplatyrhynchos.ZJU1.0)
    library(GenomicRanges)
    library(rtracklayer)
    library(Rsamtools)
    library(data.table)
    library(optparse)
})

# ── option parsing ─────────────────────────────────────────────────────
option_list <- list(
    make_option("--libraries", type = "character", default = NULL,
                help = "Comma-separated subset of libraries (pilot mode).  Default: all 12"),
    make_option("--threads", type = "integer", default = 4,
                help = "Thread count passed to ArchR addArchRThreads"),
    make_option("--out-root", type = "character",
                default = paste0(SCRATCH_ROOT, "/phase2/atac_archr"),
                help = "Output root directory (no trailing slash)"),
    make_option("--macs3-path", type = "character",
                default = paste0(HOME_ROOT, "/.conda/envs/snapatac2/bin/macs3"),
                help = "Path to MACS3 binary (ArchR takes MACS2-API but runs any MACS)"),
    make_option("--scratch", type = "character",
                default = SCRATCH_ROOT,
                help = "FC_08656 scratch root")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── paths ──────────────────────────────────────────────────────────────
SCRATCH  <- opt$scratch
OUT_ROOT <- opt[["out-root"]]
ARROW_DIR <- file.path(OUT_ROOT, "ArrowFiles")
PROJ_DIR  <- file.path(OUT_ROOT, "proj_fig03")
PEAK_DIR  <- file.path(OUT_ROOT, "peaks")
META_TSV  <- file.path(SCRATCH, "phase2/atac/step_f1_coupled_atac_metadata.tsv")
GTF_GZ    <- paste0(REFERENCE_ROOT, "/crarc_zju1_PITX2/genes/genes.gtf.gz")

dir.create(ARROW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PEAK_DIR,  recursive = TRUE, showWarnings = FALSE)

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── library → fragments (mirrors pipeline/f1_coupled_atac.py) ─────────
# Single source of truth: these paths must track the Python dict in
# f1_coupled_atac.main().  Any deviation breaks cross-track comparison.
ARC_LIBS <- list(
    `PILOTARC3`  = file.path(SCRATCH, "PILOTARC3/outs/atac_fragments.tsv.gz"),
    `PILOTARC4A` = file.path(SCRATCH, "PILOTARC4A/outs/atac_fragments.tsv.gz"),
    `PILOTARC4B` = file.path(SCRATCH, "PILOTARC4B/outs/atac_fragments.tsv.gz"),
    `ARC05`      = file.path(SCRATCH, "ARC05/outs/atac_fragments.tsv.gz"),
    `ARC06`      = file.path(SCRATCH, "ARC06/outs/atac_fragments.tsv.gz"),
    `ARC07`      = file.path(SCRATCH, "ARC07/outs/atac_fragments.tsv.gz"),
    `ARC08`      = file.path(SCRATCH, "ARC08/outs/atac_fragments.tsv.gz"),
    `ARC03-1A`   = file.path(SCRATCH, "ARC03-1A/outs/atac_fragments.tsv.gz"),
    `ARC03-1B`   = file.path(SCRATCH, "ARC03-1B/outs/atac_fragments.tsv.gz"),
    `ARC03-2`    = file.path(SCRATCH, "ARC03-2/outs/atac_fragments.tsv.gz"),
    `ARC04`      = file.path(SCRATCH, "ARC04/outs/atac_fragments.tsv.gz"),
    `ARC0S-duck` = file.path(SCRATCH, "ARC0S-duck/outs/atac_fragments.tsv.gz")
)

if (!is.null(opt$libraries)) {
    sel <- strsplit(opt$libraries, ",", fixed = TRUE)[[1]]
    miss <- setdiff(sel, names(ARC_LIBS))
    if (length(miss)) stop("Unknown library(s): ", paste(miss, collapse = ", "))
    ARC_LIBS <- ARC_LIBS[sel]
    ts("pilot mode — libraries: ", paste(names(ARC_LIBS), collapse = ", "))
}

# SKIP_POOLS mirrors f2_peak_calling.py SKIP_POOLS (plan §3.3 gut L/R
# confound).  Any pool listed here is masked from peak calling even if
# it exceeds the 500-cell floor.
SKIP_POOLS <- c("HH21|Pharyngeal endoderm", "HH25|Pharyngeal endoderm")
MIN_POOL_CELLS <- 500L  # plan §3.3 / §3.4 floor for MACS3 (sub-500 → gkmQC)

# ── QC thresholds (match f1_coupled_atac.py MIN_TSS / MIN_FRAGS) ──────
MIN_TSS   <- 2.0
MIN_FRAGS <- 1000L

# ── RNG seed ──────────────────────────────────────────────────────────
# The build was unseeded until 2026-08-06.  The RNG consumer that matters is
# addGroupCoverages(), which subsamples cells into the reproducibility="2"
# pseudobulk replicates; MACS3 and the iterative-overlap merge are deterministic
# given those coverages.  So an unseeded run yields a DIFFERENT union: a measured
# duck rebuild gave 277,604 vs 283,818 peaks (-2.2%), PITX2 +/-500 kb 310 vs 329
# (-5.8%), with the focal PITX2c-promoter peaks surviving (+/-4-7 bp summit shift)
# -- docs/archr_peak_reproducibility_decision.md (not distributed).
#
# THE SHIPPED BED IS THE ARTIFACT OF RECORD.  figures/data/archr_union_peaks.bed.gz
# (283,818 peaks) is what every published DA / motif / PITX2-locus number is
# computed against, and it is the ONLY surviving copy (absent from scratch and
# from standby atac_backup/atac_archr/peaks/).  This seed makes FUTURE builds
# reproducible; it does NOT reproduce the shipped set, and re-running would move
# every peak-level denominator in the paper.  Do not rebuild to "restore" it.
#
# Caveat: ArchR parallelises via mclapply, and worker RNG streams are not derived
# from this seed, so set.seed() alone guarantees reproducibility only at
# threads = 1.  Pin --threads 1 for a build intended to be byte-reproducible.
ARCHR_SEED <- 1L
set.seed(ARCHR_SEED)

# ── ArchR genome + gene annotations ───────────────────────────────────
ts("building genomeAnnotation + geneAnnotation from BSgenome + GTF")
ts("RNG seed = ", ARCHR_SEED, " (threads = ", opt$threads,
   if (opt$threads > 1) "; >1 thread: mclapply worker streams are unseeded" else "", ")")
addArchRThreads(threads = opt$threads, force = TRUE)
# H5 locking ENABLED.  Pilot 37866500 (2026-04-24) FAILED at tabix-to-Arrow
# link creation with `H5Lcreate_external ... HDF5. Links. Unable to
# initialize object` after 88% of tabix read.  Root cause: even with one
# library at a time, ArchR's mclapply spawns N=threads workers, all
# writing per-chrom chunks into the SAME master Arrow file.  When two
# workers finish near-simultaneously, the H5 link-table insertion races
# and corrupts the table.  The previous "safe because libraries process
# sequentially" reasoning ignored intra-library chunk concurrency.
# Locking serialises the link inserts; tabix reads can still parallelise
# at the OS level above the lock.
addArchRLocking(locking = TRUE)

bsg <- BSgenome.Anasplatyrhynchos.ZJU1.0

# ArchR's createArrowFiles chunks the genome per seqname (~5 chunks per
# chromosome).  The ZJU1.0 assembly has 34 chromosomes (NC_*) + 722
# tiny scaffolds (NW_*) — the latter produce 722 single-chunk tabix
# reads that scan the fragments file end-to-end each time, blowing the
# Arrow-build wall-time from ~30 min to 10+ hours.  Fig03 biology
# (PITX2 locus, genome-wide LPM DA) is entirely on chromosomes, so
# restrict ArchR's genome to NC_* only.  That gives 34 chroms × ~5
# chunks = ~170 chunks vs. 3780 today.  (Cancelled pilot 37827149
# confirmed: chunk 1 of 3780 took ~3 min; scaling implies ~190h wall.)
all_NC <- seqnames(bsg)[grep("^NC_", seqnames(bsg))]
# Exclude NC_009684.1 (Anas platyrhynchos mitochondrion, ~16.6 kb) on
# top of the NW_* scaffold drop.  Two reasons:
#   1. The NCBI Gnomon GTF annotates mtDNA genes with circular-genome
#      wrap-around coordinates: NC_009684.1:-1144-4833:+ (negative
#      start) and NC_009684.1:50-6088:+ (start near origin).  ArchR's
#      `.addGeneScoreMat` rejects these with "gene boundaries -
#      minError" (pilot 37878133, after the genes-from-transcripts fix
#      pushed past the Promoter/TileMatrix steps).  Clamping would
#      distort wrap-around genes; excluding the chromosome is cleaner.
#   2. mtDNA peaks are not relevant to fig03 biology (LPM mesenchyme,
#      PITX2 locus).  ATAC fragments on mtDNA are largely
#      mitochondrial-contamination noise.
main_seqs <- setdiff(all_NC, "NC_009684.1")
ts("filtering BSgenome to autosomes/sex only: ", length(main_seqs), " of ",
   length(seqnames(bsg)), " seqs (dropping ",
   length(seqnames(bsg)) - length(all_NC), " NW_* scaffolds + ",
   length(all_NC) - length(main_seqs), " mt-DNA contig)")

# Build our own chromSizes GRanges restricted to main chroms.  ArchR's
# `filter = TRUE` path uses filterChrGR, which has `underscore = TRUE`
# and `standard = TRUE` defaults — both strip NCBI-style names
# (NC_051772.1 has an underscore; not in UCSC "standard" chromosome
# list).  Passing our own chromSizes object bypasses that filter.
seqlens <- seqlengths(bsg)
main_chromSizes <- GRanges(
    seqnames = main_seqs,
    ranges   = IRanges(start = 1L, end = seqlens[main_seqs]),
    strand   = "*",
    seqlengths = seqlens[main_seqs]
)

# ``blacklist`` is optional — duck has no published ENCODE-style
# blacklist.  Pass empty GRanges rather than NULL to satisfy the
# ArchR validator.
blacklist <- GRanges()
genomeAnnotation <- createGenomeAnnotation(
    genome = bsg,
    chromSizes = main_chromSizes,
    blacklist = blacklist,
    filter = FALSE           # we already supplied chromSizes
)

# Gene annotation: parse the crarc GTF.  rtracklayer::import on a
# gzipped GTF is the canonical entry point; the result is a GRanges.
ts("parsing GTF: ", GTF_GZ)
gtf <- rtracklayer::import(GTF_GZ, format = "gtf")
# Scope to same chromosomes as genomeAnnotation (main chroms only;
# drops ~722 scaffolds' worth of genes).  Keeps seqlevels in sync,
# which ArchR enforces downstream.
gtf <- keepSeqlevels(gtf, intersect(seqlevels(gtf), main_seqs),
                     pruning.mode = "coarse")
ts("  GTF restricted to ", length(seqlevels(gtf)), " main chroms — ",
   length(gtf), " features retained")
exons_gr <- gtf[which(gtf$type == "exon")]
tx_gr    <- gtf[which(gtf$type == "transcript")]
native_genes_gr <- gtf[which(gtf$type == "gene")]

# Harmonise required fields — ArchR expects specific column names on
# each GRanges.
if (!"gene_id" %in% colnames(mcols(exons_gr))) stop("gene_id missing from GTF 'exon' features")
if (!"gene_id" %in% colnames(mcols(tx_gr)))    stop("gene_id missing from GTF 'transcript' features")

# CRITICAL: derive gene-level GRanges by aggregating transcripts per
# gene_id, NOT from `gtf$type == "gene"` rows.  The NCBI Gnomon GTF for
# ZJU1.0 has zero gene-type entries on NC_051772.1 (chr1), despite
# 8,186 transcripts there.  ArchR's `.fastFeatureCounts` does
# `split(Promoter, seqnames(Promoter))[chrArrow]` and a missing chrom
# in featureList yields an NA-named NULL element (R's named-list
# subsetting semantics) — the next line passes `chr=NA` to
# `.getFragsFromArrow` and crashes with "Chromosome NA not in
# ArrowFile" (pilots 37870636, 37874522).
#
# Native gene records still provide nicer symbol attributes; use them
# as a side lookup, but the gene boundaries themselves come from
# transcript spans.
suppressPackageStartupMessages(library(data.table))
ts("aggregating transcripts → genes (one row per gene_id) ...")
tx_dt <- as.data.table(tx_gr)
genes_dt <- tx_dt[, .(seqnames = as.character(seqnames[1L]),
                      start    = min(start),
                      end      = max(end),
                      strand   = as.character(strand[1L])),
                  by = gene_id]
ts("  ", nrow(genes_dt), " genes built from ", length(tx_gr), " transcripts")
ts("  chroms covered : ",
   length(unique(genes_dt$seqnames)),
   " (target = ", length(main_seqs), ")")

# Symbol lookup from native gene records when available.
sym_lookup <- character(0)
if (length(native_genes_gr) > 0) {
    cand_cols <- intersect(c("gene_name", "gene"),
                           colnames(mcols(native_genes_gr)))
    if (length(cand_cols)) {
        ng_sym <- rep(NA_character_, length(native_genes_gr))
        for (cc in cand_cols) {
            v <- mcols(native_genes_gr)[[cc]]
            ng_sym <- ifelse(is.na(ng_sym) | ng_sym == "", v, ng_sym)
        }
        sym_lookup <- setNames(ng_sym, native_genes_gr$gene_id)
    } else {
        sym_lookup <- setNames(native_genes_gr$gene_id, native_genes_gr$gene_id)
    }
}

resolved_sym <- if (length(sym_lookup)) sym_lookup[genes_dt$gene_id] else NA_character_
resolved_sym <- ifelse(is.na(resolved_sym) | resolved_sym == "",
                       genes_dt$gene_id, resolved_sym)

genes_for_archr <- GRanges(
    seqnames = genes_dt$seqnames,
    ranges   = IRanges(start = genes_dt$start, end = genes_dt$end),
    strand   = genes_dt$strand,
    gene_id  = genes_dt$gene_id,
    symbol   = as.character(resolved_sym)
)
# Sanity: every chrom in main_seqs MUST have at least one gene now,
# otherwise the same NA-name bug recurs.
covered <- unique(as.character(seqnames(genes_for_archr)))
missing_chroms <- setdiff(main_seqs, covered)
if (length(missing_chroms)) {
    stop("genes_for_archr is missing chroms expected in the Arrow file: ",
         paste(missing_chroms, collapse = ", "),
         " — `.fastFeatureCounts` will crash on these.  Investigate the GTF.")
}
ts("genes_for_archr: ", length(genes_for_archr), " ranges across ",
   length(covered), " chroms (all main_seqs covered)")

# Transcript-start TSS GRanges: one range per transcript at its TSS.
tss_gr <- GRanges(
    seqnames = seqnames(tx_gr),
    ranges   = IRanges(start = ifelse(as.character(strand(tx_gr)) == "-",
                                      end(tx_gr), start(tx_gr)),
                       width = 1L),
    strand   = strand(tx_gr),
    tx_id    = tx_gr$transcript_id,
    tx_name  = tx_gr$transcript_id,
    gene_id  = tx_gr$gene_id
)

# Resolve a display symbol for each exon the same way we did for genes.
exon_cols <- intersect(c("gene_name", "gene"), colnames(mcols(exons_gr)))
if (length(exon_cols) == 0) {
    exon_symbol <- exons_gr$gene_id
} else {
    exon_symbol <- rep(NA_character_, length(exons_gr))
    for (cc in exon_cols) {
        v <- mcols(exons_gr)[[cc]]
        exon_symbol <- ifelse(is.na(exon_symbol) | exon_symbol == "", v, exon_symbol)
    }
    exon_symbol <- ifelse(is.na(exon_symbol) | exon_symbol == "",
                          exons_gr$gene_id, exon_symbol)
    exon_symbol <- as.character(exon_symbol)
}

geneAnnotation <- createGeneAnnotation(
    TSS   = tss_gr,
    exons = GRanges(seqnames = seqnames(exons_gr),
                    ranges   = ranges(exons_gr),
                    strand   = strand(exons_gr),
                    gene_id  = exons_gr$gene_id,
                    symbol   = exon_symbol),
    genes = genes_for_archr,
    annoStyle = "GENCODE"    # harmless — ArchR just uses it as a tag
)

ts("  genomeAnnotation seqlevels: ", length(seqlevels(genomeAnnotation$chromSizes)))
ts("  geneAnnotation gene count : ", length(geneAnnotation$genes))
ts("  geneAnnotation TSS count  : ", length(geneAnnotation$TSS))

# ── load f1 metadata (barcode → peak_pool lookup) ─────────────────────
ts("loading f1 metadata: ", META_TSV)
meta <- fread(META_TSV, sep = "\t", header = TRUE, na.strings = "")
# The header has an empty leading field because pandas writes the
# index column unnamed; rename it.
if (colnames(meta)[1] == "V1" || colnames(meta)[1] == "") {
    setnames(meta, 1L, "cb_fragments")
}
ts("  metadata rows: ", nrow(meta))
ts("  peak_pool NA : ", sum(is.na(meta$peak_pool)))

# ArchR per-cell barcode format: "<sampleName>#<barcode>".  The
# fragments file uses "<cb16>-1", matching meta$cb_fragments.
meta[, archr_cb := paste0(arc_library, "#", cb_fragments)]
setkey(meta, archr_cb)

# ── 1. Arrow files (one per library) ──────────────────────────────────
setwd(OUT_ROOT)
# ArchR lazily checks a global genomeAnnotation in a handful of helpers.
# We pass our custom duck genomeAnnotation + geneAnnotation explicitly
# to createArrowFiles / ArchRProject / addGroupCoverages, so no
# addArchRGenome() call is needed.  Avoiding the global also sidesteps
# a BSgenome.Hsapiens.UCSC.hg38 (~860 MiB) hard dependency.

arrow_paths <- character(0)
need_build <- character(0)
for (lib in names(ARC_LIBS)) {
    # ArchR writes Arrow files to cwd (= OUT_ROOT after setwd above), not
    # to the ARROW_DIR we created.  Check both locations so subsequent
    # runs can resume from a previously-built Arrow file.
    cand <- c(file.path(OUT_ROOT,    paste0(lib, ".arrow")),
              file.path(ARROW_DIR,   paste0(lib, ".arrow")))
    found <- cand[file.exists(cand) & vapply(cand, function(p)
        if (file.exists(p)) file.info(p)$size > 1e6 else FALSE, logical(1))]
    if (length(found)) {
        ts("[skip] Arrow exists — ", found[1])
        arrow_paths[lib] <- found[1]
        next
    }
    if (!file.exists(ARC_LIBS[[lib]])) {
        stop("missing fragments for ", lib, ": ", ARC_LIBS[[lib]])
    }
    need_build <- c(need_build, lib)
}

# Build all uncached Arrows in a SINGLE createArrowFiles call so that
# ArchR's mclapply (gated on threads = opt$threads) parallelises across
# libraries.  Job 37897987 confirmed: calling createArrowFiles once
# per library inside an outer R loop = strictly sequential, threads
# wasted.  The mclapply parallelism only kicks in when multiple
# inputFiles are passed in a single call.  Each library is still
# single-threaded internally (locking=TRUE), so peak per-worker RAM
# stays at ~43 GB observed in pilot 37895184.
if (length(need_build) > 0) {
    ts("createArrowFiles (parallel across ", length(need_build),
       " libs at threads=", opt$threads, "): ",
       paste(need_build, collapse = ", "))
    input_files <- vapply(need_build, function(lib) ARC_LIBS[[lib]],
                          character(1))
    names(input_files) <- need_build
    arrs <- createArrowFiles(
        inputFiles  = input_files,
        sampleNames = need_build,
        outputNames = need_build,
        minTSS      = MIN_TSS,
        minFrags    = MIN_FRAGS,
        addTileMat  = TRUE,
        addGeneScoreMat = TRUE,
        # ArchR's `useGeneBoundaries=TRUE` (default) fires a false-
        # positive validation for any gene with start < min(extendUpstream)
        # = 1000.  After clamping `s = pmax(1, s)`, ArchR then checks
        # `pminGene - pReverseMin < start(extendedGeneRegion)` —
        # always TRUE when geneStart < 1001 (one or two such genes per
        # chrom, near the contig origin).  Pilot 37882011 hit this on
        # NC_051783.1 after the mtDNA-exclusion fix unblocked everything
        # else.  Setting useGeneBoundaries=FALSE skips the validation;
        # gene-score uses simple fixed-window extension.  Acceptable
        # since gene-score is downstream metadata, not the peak-calling
        # output.
        GeneScoreMatParams = list(useGeneBoundaries = FALSE),
        genomeAnnotation = genomeAnnotation,
        geneAnnotation   = geneAnnotation,
        force       = FALSE,
        threads     = opt$threads,
        # nChunk=5 default creates 170 genome-wide chunks (34 chroms × 5),
        # each of which scans the chrom's fragments subset.  Pilot 37828381
        # was on track for ~11h wall.  nChunk=1 → 34 chunks (one per
        # chrom) cuts per-chunk setup overhead ~5× while keeping each
        # chunk's memory ≤ ~1 GiB (even for chr1 with 32M fragments).
        nChunk      = 1L
    )
    # createArrowFiles returns the Arrow paths in the same order as
    # inputs.  Stitch back into arrow_paths keyed by lib name.
    for (i in seq_along(arrs)) {
        arrow_paths[need_build[i]] <- arrs[i]
    }
}

ts("Arrow files ready (", length(arrow_paths), " libs):")
for (lib in names(arrow_paths)) ts("  ", lib, " -> ", arrow_paths[[lib]])

# ── 2. ArchRProject + peak_pool attach ───────────────────────────────
ts("building ArchRProject")
if (dir.exists(PROJ_DIR)) {
    ts("  existing project dir — loading")
    proj <- loadArchRProject(path = PROJ_DIR, showLogo = FALSE)
} else {
    proj <- ArchRProject(
        ArrowFiles = arrow_paths,
        outputDirectory = PROJ_DIR,
        copyArrows = TRUE,
        genomeAnnotation = genomeAnnotation,
        geneAnnotation   = geneAnnotation,
        showLogo = FALSE
    )
}

cell_ids <- rownames(getCellColData(proj))
ts("  ArchRProject cells (post-QC): ", length(cell_ids))

# Attach peak_pool + unified_genotype_id to cellColData.
attach_df <- meta[J(cell_ids), on = "archr_cb",
                  .(cb = archr_cb, peak_pool, unified_genotype_id,
                    d2g_subcluster, d2g_stage)]
stopifnot(identical(attach_df$cb, cell_ids))

# ArchR requires columns to be named vectors aligned to cells()
proj <- addCellColData(ArchRProj = proj, data = attach_df$peak_pool,
                       cells = cell_ids, name = "peak_pool", force = TRUE)
proj <- addCellColData(ArchRProj = proj, data = attach_df$unified_genotype_id,
                       cells = cell_ids, name = "unified_genotype_id", force = TRUE)
proj <- addCellColData(ArchRProj = proj, data = attach_df$d2g_subcluster,
                       cells = cell_ids, name = "d2g_subcluster", force = TRUE)
proj <- addCellColData(ArchRProj = proj, data = attach_df$d2g_stage,
                       cells = cell_ids, name = "d2g_stage", force = TRUE)

cc <- getCellColData(proj)
ts("  cells with non-NA peak_pool: ", sum(!is.na(cc$peak_pool)))

# Subset to cells that have both a peak_pool and pass the SKIP_POOLS
# mask; peak calling runs only over these.
valid <- !is.na(cc$peak_pool) & !(cc$peak_pool %in% SKIP_POOLS)
pool_counts <- table(cc$peak_pool[valid])
callable <- names(pool_counts)[pool_counts >= MIN_POOL_CELLS]
subthreshold <- names(pool_counts)[pool_counts < MIN_POOL_CELLS]

ts("  callable pools (>= ", MIN_POOL_CELLS, " cells): ", length(callable))
ts("  sub-threshold pools (→ gkmQC rescue in f2b): ", length(subthreshold))
for (p in subthreshold) ts("    ", p, " (", pool_counts[[p]], ")")

# Restrict project for peak calling to callable pools.
keep_cells <- cell_ids[!is.na(cc$peak_pool) & cc$peak_pool %in% callable]
ts("  cells feeding peak calling: ", length(keep_cells))
proj_call <- subsetArchRProject(
    ArchRProj = proj,
    cells = keep_cells,
    outputDirectory = file.path(OUT_ROOT, "proj_fig03_callable"),
    dropCells = FALSE, force = TRUE
)

# ── 3. group coverages + reproducible peak set ────────────────────────
ts("addGroupCoverages(groupBy = 'peak_pool')")
# reproducibility="2" replicates are ArchR-default Sample/library-level — the PRODUCTION set.  A
# genotype-level variant (sampleLabels="unified_genotype_id") was tested as a sensitivity check; the
# focal PITX2c-promoter result is invariant in both species (chick ΔFRiP r=1.000; duck promoter peaks
# survive ±4-7 bp), so the production library-level set is retained.  See
# docs/archr_peak_reproducibility_decision.md *(not distributed)*.
proj_call <- addGroupCoverages(
    ArchRProj = proj_call,
    groupBy = "peak_pool",
    minCells = MIN_POOL_CELLS,
    maxCells = 10000L,      # cap pseudobulk size; ArchR default 1500 is too small
    force = TRUE,
    threads = opt$threads
)

# Effective genome size for MACS3 -g.  ArchR's addReproduciblePeakSet
# auto-detects this for hg/mm/etc but errors out for non-standard
# genomes ("Non-standard genome detected. Argument genomeSize is
# required!" — pilot 37887058).  Compute precisely as the count of
# non-N bases across main_seqs (ACGT in BSgenome).
ts("computing effective genome size (non-N bases across main_seqs) ...")
genome_size_effective <- sum(vapply(main_seqs, function(seq) {
    sum(letterFrequency(bsg[[seq]], letters = "ACGT", as.prob = FALSE))
}, numeric(1)))
ts("  effective genome size: ", format(genome_size_effective, big.mark = ","),
   " bp (", round(genome_size_effective / sum(seqlens[main_seqs]) * 100, 1),
   "% of assembly)")

ts("addReproduciblePeakSet(reproducibility = '2', MACS3 binary = ", opt[["macs3-path"]], ")")
proj_call <- addReproduciblePeakSet(
    ArchRProj = proj_call,
    groupBy = "peak_pool",
    pathToMacs2 = opt[["macs3-path"]],  # ArchR calls it MACS2 but shell-runs whatever you hand it
    genomeSize = genome_size_effective, # required for non-standard genomes
    reproducibility = "2",              # peak in >=2 library/Sample pseudobulks (ArchR-default sampleLabels="Sample" = PRODUCTION; genotype-level tested as a sensitivity check, focal result invariant — docs/archr_peak_reproducibility_decision.md, not distributed)
    peaksPerCell = 500L,
    maxPeaks = 150000L,
    minCells = MIN_POOL_CELLS,
    extendSummits = 250L,               # 501 bp fixed width (plan §3.2)
    cutOff = 0.1,                       # MACS3 q-value
    force = TRUE,
    threads = opt$threads
)

peak_set <- getPeakSet(proj_call)
ts("peak set size: ", length(peak_set), " peaks")

# ── 4. exports ────────────────────────────────────────────────────────
union_bed <- file.path(PEAK_DIR, "union_peaks.bed")
ts("writing union BED: ", union_bed)
export(peak_set, con = union_bed, format = "BED")

# Per-pool peak lists.  ArchR stores these inside the project under
# proj_fig03_callable/PeakCalls/peak_pool/<pool>-reproduciblePeaks.gr.rds
# (the `peak_pool/` subdirectory comes from groupBy="peak_pool" — confirmed
# by pilot 37895184 .out log).  recursive=TRUE picks them up regardless
# of which subdir ArchR uses.
pool_rds_dir <- file.path(OUT_ROOT, "proj_fig03_callable", "PeakCalls")
pool_files <- list.files(pool_rds_dir,
                         pattern = "-reproduciblePeaks\\.gr\\.rds$",
                         full.names = TRUE, recursive = TRUE)
if (length(pool_files) == 0) {
    # Older ArchR versions stored these one level up; check there too.
    alt_dir <- file.path(PROJ_DIR, "PeakCalls")
    pool_files <- list.files(alt_dir,
                             pattern = "-reproduciblePeaks\\.gr\\.rds$",
                             full.names = TRUE, recursive = TRUE)
}
# ArchR sanitises pool names in filenames: any character that's not
# [a-zA-Z0-9._] becomes ".".  Build sanitised → original lookup so
# downstream pool_counts / cc lookups (keyed on original names) work.
# Job 37930531 (full cohort) hit "subscript out of bounds" here because
# the loop tried `pool_counts[["HH17.16.Cardiopharyngeal.LPM"]]` while
# pool_counts had `"HH17/16|Cardiopharyngeal LPM"` as the key.
archr_sanitise <- function(x) gsub("[^a-zA-Z0-9._]", ".", x)
sanitised_to_orig <- setNames(names(pool_counts),
                              archr_sanitise(names(pool_counts)))

per_pool_rows <- list()
repro_rows <- list()
for (f in pool_files) {
    pool_san <- sub("-reproduciblePeaks\\.gr\\.rds$", "", basename(f))
    pool <- sanitised_to_orig[[pool_san]]
    if (is.null(pool) || is.na(pool)) {
        ts("WARNING: cannot map sanitised pool name '", pool_san, "' back to original; skipping")
        next
    }
    gr <- readRDS(f)
    n_pk <- length(gr)
    n_cells <- as.integer(pool_counts[[pool]])
    # n_rep = genotype cardinality (informational, exported as n_replicates); the peak
    # reproducibility filter itself is library-level (ArchR-default sampleLabels="Sample").
    n_rep <- length(unique(cc$unified_genotype_id[cc$peak_pool == pool]))
    repro_rows[[length(repro_rows) + 1]] <- data.frame(
        peak_pool = pool,
        n_cells = n_cells,
        n_replicates = n_rep,
        n_peaks = n_pk,
        stringsAsFactors = FALSE
    )
    per_pool_rows[[pool]] <- data.frame(
        peak_pool = pool,
        chrom = as.character(seqnames(gr)),
        start = start(gr) - 1L,      # 0-based BED
        end   = end(gr),
        stringsAsFactors = FALSE
    )
}
if (length(repro_rows)) {
    repro_df <- do.call(rbind, repro_rows)
    write.table(repro_df, file.path(PEAK_DIR, "reproducibility_metadata.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    per_pool_df <- do.call(rbind, per_pool_rows)
    write.table(per_pool_df, file.path(PEAK_DIR, "per_pool_peaks.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    ts("wrote reproducibility_metadata.tsv (", nrow(repro_df), " pools)")
    ts("wrote per_pool_peaks.tsv (", nrow(per_pool_df), " rows)")
} else {
    warning("no per-pool peak RDS files found under ", pool_rds_dir)
}

# ── 5. save project so subsequent steps (DA, coverage) can resume ─────
# proj_call was written into .../proj_fig03_callable by subsetArchRProject;
# save back to the SAME dir (not PROJ_DIR, which holds the pre-subset
# full-cohort project).  Downstream (pseudobulk DA, locus pileups) loads
# from proj_fig03_callable.
callable_dir <- file.path(OUT_ROOT, "proj_fig03_callable")
ts("saveArchRProject — ", callable_dir)
saveArchRProject(ArchRProj = proj_call, outputDirectory = callable_dir, load = FALSE)

ts("=== archr_build_peaks.R done ===")
