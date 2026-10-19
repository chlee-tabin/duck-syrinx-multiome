#!/usr/bin/env Rscript
# archr_build_bsgenome.R — forge + install BSgenome.Ggallus.crarc.galGal6a
#
# Chick analog of pipeline/archr_build_bsgenome.R (duck ZJU1.0).  Forges a
# BSgenome from the EXACT galGal6a FASTA the chick ATAC fragments were aligned
# to (cellranger-arc crarc_galgal6a, UCSC chrom names) — so peak/locus coords
# match the fragments byte-for-byte.  We do NOT use the stock
# BSgenome.Ggallus.UCSC.galGal6: a silent main-chrom size mismatch vs the crarc
# reference would corrupt the PITX2-locus analysis (the worst kind of bug), and
# forging is the duck-faithful, coord-guaranteed path.
#
# Staged FASTA (chick_demux/ref/genome.fa): 1.08 Gb, 464 sequences = 34 main
# chroms (chr1-28, chr30-33, chrZ, chrW) + chrM + ~428 *_random scaffolds.  NO
# EGFP (that contig is only in the Group B EGFP reference; main chroms identical,
# so one BSgenome serves both groups).
#
# Run on an o2 compute node (~8 GiB, ~30 min):  sbatch pipeline/chick/archr_build_bsgenome.sh
# Idempotent: smoke-test skip if the package already loads.

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(Biostrings)
    library(BSgenomeForge)
    library(GenomeInfoDb)
    library(S4Vectors)
})

REF_FASTA <- paste0(SCRATCH_USER, "/chick_demux/ref/genome.fa")
REF_FAI   <- paste0(REF_FASTA, ".fai")

BUILD_DIR <- paste0(SCRATCH_USER, "/chick_demux/atac_archr/bsgenome_build")
SEQ_DIR   <- file.path(BUILD_DIR, "single_sequences")
PKG_DIR   <- file.path(BUILD_DIR, "pkg_src")
TAR_DIR   <- file.path(BUILD_DIR, "tarballs")
dir.create(SEQ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PKG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAR_DIR, recursive = TRUE, showWarnings = FALSE)

PKG_NAME    <- "BSgenome.Ggallus.crarc.galGal6a"
PKG_VERSION <- "1.0.0"
PKG_TITLE   <- "Full genome sequences for Gallus gallus (galGal6a, crarc)"
PKG_DESC    <- paste(
    "Full genome sequences for the chicken (Gallus gallus) galGal6a / GRCg6a as",
    "packaged in the cellranger-arc crarc_galgal6a reference (UCSC chrom names),",
    "stored as a Biostrings object.  Built in-house for the duck-syrinx-multiome",
    "chick-control snATAC track (2026-06-13).  Coord-identical to the chick ATAC",
    "fragments (atac_fragments.tsv.gz)."
)
PKG_SEED_BASENAME <- "bsgenome_ggallus_crarc_galgal6a-seed.dcf"

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# Early idempotency skip.
if (suppressPackageStartupMessages(requireNamespace(PKG_NAME, quietly = TRUE))) {
    ts(PKG_NAME, " already installed — smoke test + exit")
    suppressPackageStartupMessages(library(PKG_NAME, character.only = TRUE))
    bsg <- getBSgenome(PKG_NAME)
    ts("  seqlevels N=", length(seqlevels(bsg)), "  e.g. ",
       paste(head(seqlevels(bsg)), collapse = ", "))
    ts("  chr4 length=", seqlengths(bsg)[["chr4"]], " (expect 91315245)")
    quit(save = "no", status = 0)
}

# ── 1. split FASTA per seqname ─────────────────────────────────────────
ts("reading FAI (", REF_FAI, ")")
fai <- read.table(REF_FAI, header = FALSE, sep = "\t",
                  col.names = c("seqname", "length", "offset", "linebases", "linewidth"),
                  stringsAsFactors = FALSE)
ts("  ", nrow(fai), " sequences in FAI")

existing <- length(list.files(SEQ_DIR, pattern = "\\.fa$"))
if (existing >= nrow(fai)) {
    ts("SEQ_DIR already has ", existing, " per-seqname FASTAs — skipping split")
} else {
    ts("reading full FASTA into memory (~1-2 min)")
    dna <- readDNAStringSet(REF_FASTA, format = "fasta")
    dna_names <- sub(" .*", "", names(dna))
    names(dna) <- dna_names
    miss <- setdiff(fai$seqname, dna_names)
    if (length(miss)) stop("FASTA missing seqnames present in FAI: ", paste(head(miss), collapse = ", "))
    ts("writing ", length(dna), " per-seqname FASTAs to ", SEQ_DIR)
    for (nm in dna_names) {
        out <- file.path(SEQ_DIR, paste0(nm, ".fa"))
        if (!file.exists(out)) writeXStringSet(dna[nm], filepath = out, format = "fasta", width = 80L)
    }
    rm(dna); gc(verbose = FALSE)
    ts("split complete: ", length(list.files(SEQ_DIR, pattern = "\\.fa$")), " files")
}

# ── 2. seed DCF ─────────────────────────────────────────────────────────
seed_path <- file.path(BUILD_DIR, PKG_SEED_BASENAME)
# galGal6a mitochondrion is named chrM in the crarc UCSC-style reference.
circ <- intersect(fai$seqname, c("chrM"))
circ_literal <- if (length(circ)) {
    sprintf("c(%s)", paste(sprintf('"%s"', circ), collapse = ", "))
} else "character(0)"

seed_lines <- c(
    sprintf("Package: %s", PKG_NAME),
    sprintf("Title: %s", PKG_TITLE),
    sprintf("Description: %s", PKG_DESC),
    sprintf("Version: %s", PKG_VERSION),
    "organism: Gallus gallus",
    "common_name: Chicken",
    "genome: galGal6a",
    "provider: UCSC/cellranger-arc",
    "release_date: 2018",
    "source_url: https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/",
    "organism_biocview: Gallus_gallus",
    "BSgenomeObjname: Ggallus",
    sprintf("seqnames: c(%s)", paste(sprintf('"%s"', fai$seqname), collapse = ", ")),
    sprintf("circ_seqs: %s", circ_literal),
    sprintf("seqs_srcdir: %s", SEQ_DIR),
    "seqfiles_suffix: .fa"
)
writeLines(seed_lines, seed_path)
ts("wrote seed DCF: ", seed_path, " (circ_seqs = ", circ_literal, ")")

# ── 3. forge ─────────────────────────────────────────────────────────────
pkg_src <- file.path(PKG_DIR, PKG_NAME)
if (dir.exists(pkg_src)) {
    ts("removing stale package dir: ", pkg_src)
    unlink(pkg_src, recursive = TRUE, force = TRUE)
}
ts("forgeBSgenomeDataPkg — builds package skeleton + 2bit")
BSgenomeForge::forgeBSgenomeDataPkg(x = seed_path, destdir = PKG_DIR,
                                    replace = TRUE, verbose = TRUE)

# ── 4. R CMD build ───────────────────────────────────────────────────────
ts("R CMD build")
build_out <- system(sprintf("cd %s && R CMD build %s 2>&1",
                            shQuote(TAR_DIR), shQuote(pkg_src)), intern = TRUE)
writeLines(tail(build_out, 10))
tarball <- list.files(TAR_DIR, pattern = paste0("^", PKG_NAME, ".*\\.tar\\.gz$"), full.names = TRUE)
if (length(tarball) == 0) stop("No tarball produced by R CMD build")
tarball <- tarball[order(file.info(tarball)$mtime, decreasing = TRUE)][1]
ts("tarball: ", tarball)

# ── 5. R CMD INSTALL into ~/R/4.4.2 ──────────────────────────────────────
ts("R CMD INSTALL --library=~/R/4.4.2 ", tarball)
install_out <- system(sprintf("R CMD INSTALL --library=%s %s 2>&1",
                              shQuote(path.expand("~/R/4.4.2")), shQuote(tarball)),
                      intern = TRUE)
writeLines(tail(install_out, 10))

# ── 6. smoke test ────────────────────────────────────────────────────────
ts("smoke test — library(", PKG_NAME, ")")
suppressPackageStartupMessages(library(PKG_NAME, character.only = TRUE))
bsg <- getBSgenome(PKG_NAME)
ts("  seqlevels N=", length(seqlevels(bsg)))
ts("  chr4 length=", seqlengths(bsg)[["chr4"]], " (expect 91315245)")
stopifnot(seqlengths(bsg)[["chr4"]] == 91315245L)
ts("=== archr_build_bsgenome.R done ===")
