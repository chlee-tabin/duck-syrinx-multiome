#!/usr/bin/env Rscript
# archr_pitx2_locus_pseudobulk.R — chick STEP 2/3 data foundation.
#
# Extract the PITX2 +-500 kb PeakMatrix from the built chick ArchR arrows and aggregate to
# pseudobulks in the duck g3_pitx2_locus.py schema (so the local DESeq2 reuses that code path).
#
# The build (archr_build_peaks.R) only saved proj_chick_callable (4 >=500-cell pools / 5,065
# cells) — NOT the full proj_chick — so the sub-floor focal cells (HH25 Ventral spl. + Dorsal
# sub-epi. control) live only in the Arrow files.  We therefore REBUILD the full coupled-cell
# project from the Arrow files (re-attaching the sidecar metadata), add the union peak set from
# the saved proj_chick_callable, addPeakMatrix, and extract.  Peak CALLING (4 pools) only gated
# discovery; the L/R DA tests ANY coupled cell against the union peaks.
#
# NOTE: the genomeAnnotation/geneAnnotation preamble is copied verbatim from
# archr_build_peaks.R (proven) — kept inline (not refactored into a sourced helper) to avoid
# touching the frozen, already-run build script.  If the build is re-run (e.g. minFrags relax),
# factor both into a shared _chick_archr_annotation.R then.
#
# Outputs (chick_demux/atac_archr/locus/): locus_peak_coords.tsv, locus_pseudobulk_counts.tsv,
#   locus_pseudobulk_metadata.tsv, locus_percell_counts.tsv.gz
#
# Run: sbatch pipeline/chick/archr_pitx2_locus_pseudobulk.sh

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(ArchR)
    library(BSgenome.Ggallus.crarc.galGal6a)
    library(GenomicRanges)
    library(rtracklayer)
    library(data.table)
    library(SummarizedExperiment)
})

# OUT_ROOT defaults to the production build; override with LOCUS_OUT_ROOT to re-extract from an
# isolated alternative build (e.g. the genotype-reproducibility atac_archr_genorepro out-root).
# META_TSV defaults to OUT_ROOT/<sidecar> but the genorepro out-root does NOT
# carry the sidecar (it is a build *input*), so set LOCUS_META_TSV to the original when overriding.
OUT_ROOT  <- Sys.getenv("LOCUS_OUT_ROOT", paste0(SCRATCH_USER, "/chick_demux/atac_archr"))
PROJ_CALL <- file.path(OUT_ROOT, "proj_chick_callable")    # saved; carries the union peak set
PROJ_LOCUS <- file.path(OUT_ROOT, "proj_chick_locus")      # rebuilt-from-arrows project
PROJ_LOCUS_SUB <- file.path(OUT_ROOT, "proj_chick_locus_sub")  # coupled-cell subset (must differ)
LOCUS_DIR <- file.path(OUT_ROOT, "locus")
META_TSV  <- Sys.getenv("LOCUS_META_TSV", file.path(OUT_ROOT, "chick_archr_coupled_metadata.tsv"))
REF_DIR   <- paste0(SCRATCH_USER, "/chick_demux/ref")
GTF_GZ    <- file.path(REF_DIR, "genes.gtf.gz")
CHROMS_TXT <- file.path(REF_DIR, "chroms.txt")
dir.create(LOCUS_DIR, recursive = TRUE, showWarnings = FALSE)

ARC_LIBS <- c("ARC01", "ARC02", "ARC09", "ARC10")
# Arrows may live at OUT_ROOT/<lib>.arrow OR OUT_ROOT/ArrowFiles/<lib>.arrow depending on the build;
# archr_build_peaks.R searches BOTH (its `cand` list), so mirror that here — else a fresh build that
# wrote arrows into ArrowFiles/ trips the stopifnot below before any locus extraction.
ARROW_DIR <- file.path(OUT_ROOT, "ArrowFiles")
ARROWS <- vapply(ARC_LIBS, function(lib) {
    cand <- c(file.path(OUT_ROOT, paste0(lib, ".arrow")),
              file.path(ARROW_DIR, paste0(lib, ".arrow")))
    hit <- cand[file.exists(cand)]
    if (!length(hit)) stop("no arrow for ", lib, " at: ", paste(cand, collapse = " | "))
    hit[1]
}, character(1))

# PITX2 +-500 kb on galGal6a chr4 (references/chick_pitx2_promoter_regions.bed locus_500kb).
LOCUS_CHROM <- "chr4"; LOCUS_START <- 57285031L; LOCUS_END <- 58298939L
PROMC_START <- 57794528L; PROMC_END <- 57795129L; GENE_TSS <- 57785032L

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
addArchRThreads(threads = 4, force = TRUE)
addArchRLocking(locking = TRUE)

# ── genomeAnnotation + geneAnnotation (verbatim from archr_build_peaks.R) ──────
bsg <- BSgenome.Ggallus.crarc.galGal6a
chroms_raw <- paste(readLines(CHROMS_TXT), collapse = ",")
main_seqs <- trimws(unlist(strsplit(chroms_raw, "[,[:space:]]+")))
main_seqs <- intersect(main_seqs[nzchar(main_seqs)], seqnames(bsg))
stopifnot("chr4" %in% main_seqs)
seqlens <- seqlengths(bsg)
main_chromSizes <- GRanges(seqnames = main_seqs,
                           ranges = IRanges(start = 1L, end = seqlens[main_seqs]),
                           strand = "*", seqlengths = seqlens[main_seqs])
genomeAnnotation <- createGenomeAnnotation(genome = bsg, chromSizes = main_chromSizes,
                                           blacklist = GRanges(), filter = FALSE)
ts("parsing GTF: ", GTF_GZ)
gtf <- rtracklayer::import(GTF_GZ, format = "gtf")
gtf <- keepSeqlevels(gtf, intersect(seqlevels(gtf), main_seqs), pruning.mode = "coarse")
exons_gr <- gtf[which(gtf$type == "exon")]
tx_gr    <- gtf[which(gtf$type == "transcript")]
native_genes_gr <- gtf[which(gtf$type == "gene")]
tx_dt <- as.data.table(tx_gr)
genes_dt <- tx_dt[, .(seqnames = as.character(seqnames[1L]), start = min(start),
                      end = max(end), strand = as.character(strand[1L])), by = gene_id]
sym_lookup <- character(0)
if (length(native_genes_gr) > 0) {
    cand <- intersect(c("gene_name", "gene"), colnames(mcols(native_genes_gr)))
    if (length(cand)) {
        ng_sym <- rep(NA_character_, length(native_genes_gr))
        for (cc in cand) { v <- mcols(native_genes_gr)[[cc]]; ng_sym <- ifelse(is.na(ng_sym) | ng_sym == "", v, ng_sym) }
        sym_lookup <- setNames(ng_sym, native_genes_gr$gene_id)
    } else sym_lookup <- setNames(native_genes_gr$gene_id, native_genes_gr$gene_id)
}
resolved_sym <- if (length(sym_lookup)) sym_lookup[genes_dt$gene_id] else NA_character_
resolved_sym <- ifelse(is.na(resolved_sym) | resolved_sym == "", genes_dt$gene_id, resolved_sym)
genes_for_archr <- GRanges(seqnames = genes_dt$seqnames,
                           ranges = IRanges(start = genes_dt$start, end = genes_dt$end),
                           strand = genes_dt$strand, gene_id = genes_dt$gene_id,
                           symbol = as.character(resolved_sym))
miss <- setdiff(main_seqs, unique(as.character(seqnames(genes_for_archr))))
if (length(miss)) stop("genes_for_archr missing chroms: ", paste(miss, collapse = ", "))
tss_gr <- GRanges(seqnames = seqnames(tx_gr),
                  ranges = IRanges(start = ifelse(as.character(strand(tx_gr)) == "-", end(tx_gr), start(tx_gr)), width = 1L),
                  strand = strand(tx_gr), tx_id = tx_gr$transcript_id,
                  tx_name = tx_gr$transcript_id, gene_id = tx_gr$gene_id)
exon_cols <- intersect(c("gene_name", "gene"), colnames(mcols(exons_gr)))
if (length(exon_cols) == 0) { exon_symbol <- exons_gr$gene_id } else {
    exon_symbol <- rep(NA_character_, length(exons_gr))
    for (cc in exon_cols) { v <- mcols(exons_gr)[[cc]]; exon_symbol <- ifelse(is.na(exon_symbol) | exon_symbol == "", v, exon_symbol) }
    exon_symbol <- as.character(ifelse(is.na(exon_symbol) | exon_symbol == "", exons_gr$gene_id, exon_symbol))
}
geneAnnotation <- createGeneAnnotation(
    TSS = tss_gr,
    exons = GRanges(seqnames = seqnames(exons_gr), ranges = ranges(exons_gr),
                    strand = strand(exons_gr), gene_id = exons_gr$gene_id, symbol = exon_symbol),
    genes = genes_for_archr, annoStyle = "GENCODE")

# ── rebuild project from arrows + re-attach sidecar ───────────────────────────
stopifnot(all(file.exists(ARROWS)))
ts("ArchRProject from arrows: ", paste(ARC_LIBS, collapse = ", "))
proj <- ArchRProject(ArrowFiles = ARROWS, outputDirectory = PROJ_LOCUS, copyArrows = FALSE,
                     genomeAnnotation = genomeAnnotation, geneAnnotation = geneAnnotation,
                     showLogo = FALSE)
cell_ids <- rownames(getCellColData(proj))
ts("  project cells: ", length(cell_ids))
meta <- fread(META_TSV, sep = "\t", header = TRUE, na.strings = "")
setkey(meta, archr_cb)
att <- meta[J(cell_ids), on = "archr_cb",
            .(cb = archr_cb, peak_pool, unified_genotype_id, population, stage, side, group)]
stopifnot(identical(att$cb, cell_ids))
for (col in c("peak_pool", "unified_genotype_id", "population", "stage", "side", "group"))
    proj <- addCellColData(ArchRProj = proj, data = att[[col]], cells = cell_ids, name = col, force = TRUE)
cc <- getCellColData(proj)
keep <- cell_ids[!is.na(cc$peak_pool)]
ts("  coupled cells (peak_pool not NA): ", length(keep))
proj <- subsetArchRProject(ArchRProj = proj, cells = keep, outputDirectory = PROJ_LOCUS_SUB,
                           dropCells = FALSE, force = TRUE)

# ── union peak set (from the saved callable project) + PeakMatrix ─────────────
ts("loading union peak set from ", PROJ_CALL)
pk <- getPeakSet(loadArchRProject(path = PROJ_CALL, showLogo = FALSE))
ts("  union peaks: ", length(pk))
proj <- addPeakSet(ArchRProj = proj, peakSet = pk, force = TRUE)
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)
se <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix")
gr <- rowRanges(se)
locus <- which(as.character(seqnames(gr)) == LOCUS_CHROM & start(gr) < LOCUS_END & end(gr) > LOCUS_START)
ts("  locus peaks (chr4 +-500kb): ", length(locus), " of ", length(gr))
if (length(locus) == 0L) stop("no peaks in PITX2 +-500 kb window")

se_l <- se[locus, ]; gl <- rowRanges(se_l)
peak_id <- paste0(as.character(seqnames(gl)), ":", start(gl), "-", end(gl))
counts <- as.matrix(assay(se_l)); rownames(counts) <- peak_id

coords <- data.frame(peak_id = peak_id, chrom = as.character(seqnames(gl)),
                     start = start(gl), end = end(gl),
                     in_promoter_c = (start(gl) < PROMC_END & end(gl) > PROMC_START),
                     in_gene_tss = (start(gl) <= GENE_TSS & end(gl) > GENE_TSS),
                     stringsAsFactors = FALSE)
write.table(coords, file.path(LOCUS_DIR, "locus_peak_coords.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
ts("  wrote locus_peak_coords.tsv (", nrow(coords), " peaks; ", sum(coords$in_promoter_c),
   " in PITX2c promoter, ", sum(coords$in_gene_tss), " at gene TSS)")

ccL <- as.data.frame(getCellColData(proj))[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(ccL), colnames(counts)))
percell <- data.frame(cell = colnames(counts), stage = ccL$stage, population = ccL$population,
                      side = ccL$side, group = ccL$group, genotype = ccL$unified_genotype_id,
                      nFrags = ccL$nFrags, t(counts), check.names = FALSE)
fwrite(percell, file.path(LOCUS_DIR, "locus_percell_counts.tsv.gz"), sep = "\t", compress = "gzip")
ts("  wrote locus_percell_counts.tsv.gz (", nrow(percell), " cells x ", nrow(coords), " peaks)")

side_ok <- ccL$side %in% c("L", "R")
ts("  cells with firm L/R side: ", sum(side_ok), " of ", nrow(ccL))
sample_id <- paste(ccL$stage, ccL$population, ccL$side, ccL$unified_genotype_id, sep = "|")
sample_id[!side_ok] <- NA
uniq <- sort(unique(sample_id[!is.na(sample_id)]))
ts("  pseudobulks (stage|population|side|genotype): ", length(uniq))
pb <- matrix(0L, nrow = length(uniq), ncol = nrow(counts), dimnames = list(uniq, peak_id))
meta_rows <- list()
for (sid in uniq) {
    cells <- which(sample_id == sid)
    pb[sid, ] <- rowSums(counts[, cells, drop = FALSE])
    parts <- strsplit(sid, "|", fixed = TRUE)[[1]]
    meta_rows[[sid]] <- data.frame(sample_id = sid, stage = parts[1], population = parts[2],
                                   side = parts[3], genotype = parts[4],
                                   group = ccL$group[cells[1]], n_cells = length(cells),
                                   stringsAsFactors = FALSE)
}
write.table(data.frame(sample_id = rownames(pb), pb, check.names = FALSE),
            file.path(LOCUS_DIR, "locus_pseudobulk_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
meta_df <- do.call(rbind, meta_rows)
write.table(meta_df, file.path(LOCUS_DIR, "locus_pseudobulk_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
ts("  wrote locus_pseudobulk_counts.tsv (", nrow(pb), " pseudobulks x ", nrow(coords), " peaks) + metadata")

ts("\n=== pseudobulk design preview (firm-oriented L/R only) ===")
prev <- aggregate(cbind(n_cells, n_geno = 1) ~ stage + population + side, data = meta_df, FUN = sum)
print(prev[order(prev$stage, prev$population, prev$side), ], row.names = FALSE)
ts("\n=== archr_pitx2_locus_pseudobulk.R done ===")
