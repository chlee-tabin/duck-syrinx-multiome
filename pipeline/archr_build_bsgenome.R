#!/usr/bin/env Rscript
# archr_build_bsgenome.R — forge + install BSgenome.Anasplatyrhynchos.ZJU1.0
#
# Canonical duck FASTA: /n/groups/tabin/chlee/singlecell/reference/
#   crarc_zju1_PITX2/fasta/genome.fa (persistent, byte-identical to NCBI
#   GCF_015476345.1_ZJU1.0_genomic.fna).  756 sequences (see .fai).
#
# Strategy (BSgenomeForge canonical workflow):
#   1. Split genome.fa into per-seqname fasta files under a temp dir on
#      scratch (home filesystem is quota-bound — sources are one-shot
#      and large).
#   2. Write a seed DCF file describing the package and pointing at the
#      per-seqname sources.
#   3. BSgenomeForge::forgeBSgenomeDataPkg(seed) to build the package
#      skeleton (outputs a directory, not a tarball).
#   4. R CMD build to make the source tarball.
#   5. R CMD INSTALL into ~/R/4.4.2/ so ``getBSgenome`` picks it up.
#
# Run on an o2 compute node (8 GiB mem, 30 min).  Login partition kills
# the R process during the FASTA split.
#
# The build output directory can be purged after install.  The installed
# package carries a .2bit copy of the genome (~370 MiB for a 1.1 Gb
# genome); plan for ~450 MiB permanent footprint in ~/R/4.4.2/.

## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
REFERENCE_ROOT <- .duck_env("DUCK_REFERENCE_ROOT", "/n/groups/tabin/chlee/singlecell/reference")
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))

suppressPackageStartupMessages({
    .libPaths(c("~/R/4.4.2", .libPaths()))
    library(Biostrings)
    library(BSgenomeForge)
})

# ── paths ──────────────────────────────────────────────────────────────
REF_FASTA <- paste0(REFERENCE_ROOT, "/crarc_zju1_PITX2/fasta/genome.fa")
REF_FAI   <- paste0(REF_FASTA, ".fai")

SCRATCH   <- SCRATCH_ROOT
BUILD_DIR <- file.path(SCRATCH, "phase2/atac_archr/bsgenome_build")
SEQ_DIR   <- file.path(BUILD_DIR, "single_sequences")       # per-seqname fastas
PKG_DIR   <- file.path(BUILD_DIR, "pkg_src")                # forge output
TAR_DIR   <- file.path(BUILD_DIR, "tarballs")               # R CMD build output

dir.create(SEQ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PKG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAR_DIR, recursive = TRUE, showWarnings = FALSE)

# ── package metadata ───────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(GenomeInfoDb)
    library(S4Vectors)
})

PKG_NAME      <- "BSgenome.Anasplatyrhynchos.ZJU1.0"
PKG_VERSION   <- "1.0.0"
PKG_TITLE     <- "Full genome sequences for Anas platyrhynchos (ZJU1.0)"
PKG_DESC      <- paste(
    "Full genome sequences for the duck (Anas platyrhynchos) as provided",
    "by NCBI assembly GCF_015476345.1 (ZJU1.0), stored as a Biostrings",
    "object.  Built in-house for the duck-syrinx-multiome project (fig03",
    "ArchR track, 2026-04-24)."
)
PKG_SEED_BASENAME <- "bsgenome_anasplatyrhynchos_zju1.0-seed.dcf"

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ── 1. split FASTA per seqname ─────────────────────────────────────────
# Skip if already populated with the full 756 sequences.
existing <- length(list.files(SEQ_DIR, pattern = "\\.fa$"))
ts("reading FAI (", REF_FAI, ")")
fai <- read.table(REF_FAI, header = FALSE, sep = "\t",
                  col.names = c("seqname", "length", "offset", "linebases", "linewidth"),
                  stringsAsFactors = FALSE)
ts("  ", nrow(fai), " sequences in FAI")

if (existing >= nrow(fai)) {
    ts("SEQ_DIR already has ", existing, " per-seqname FASTAs — skipping split")
} else {
    ts("SEQ_DIR currently has ", existing, " — splitting full genome")
    # Use BSgenomeForge helper if available; otherwise hand-roll.
    # As of BSgenomeForge 1.x, ``fastaTo2bit`` handles the conversion but
    # we still need per-seqname fastas for the seed's ``seqs_srcdir``
    # convention.  Biostrings::readDNAStringSet on a 1.2 GB FASTA needs
    # ~4 GiB RAM; we do it in one pass.
    ts("reading full FASTA into memory (this takes ~1-2 min)")
    dna <- readDNAStringSet(REF_FASTA, format = "fasta")
    # Biostrings strips everything after the first whitespace.  Sanity-
    # check against the FAI seqnames.
    dna_names <- names(dna)
    dna_names <- sub(" .*", "", dna_names)
    names(dna) <- dna_names
    miss <- setdiff(fai$seqname, dna_names)
    extra <- setdiff(dna_names, fai$seqname)
    if (length(miss)) stop("FASTA missing seqnames present in FAI: ", paste(head(miss), collapse=", "))
    if (length(extra)) warning("FASTA has extra seqnames not in FAI: ", paste(head(extra), collapse=", "))
    ts("writing ", length(dna), " per-seqname FASTAs to ", SEQ_DIR)
    for (nm in dna_names) {
        out <- file.path(SEQ_DIR, paste0(nm, ".fa"))
        if (!file.exists(out)) writeXStringSet(dna[nm], filepath = out, format = "fasta", width = 80L)
    }
    rm(dna); gc(verbose = FALSE)
    ts("split complete: ", length(list.files(SEQ_DIR, pattern="\\.fa$")), " files")
}

# ── 2. seed DCF file ───────────────────────────────────────────────────
# The seed describes the package to forgeBSgenomeDataPkg.  Keep the DCF
# under BUILD_DIR so the seed itself isn't shipped in the package.
seed_path <- file.path(BUILD_DIR, PKG_SEED_BASENAME)

# Circular sequences — duck genome includes mitochondrion NC_040069.1.
circ <- intersect(fai$seqname, c("NC_040069.1"))

# BSgenomeForge::forgeBSgenomeDataPkg calls Seqinfo(genome=<genome>) if
# either `seqnames` OR `circ_seqs` is NA.  ZJU1.0 isn't in GenomeInfoDb's
# registry → that call aborts with ".make_Seqinfo_from_genome" error.
# So ALWAYS write both fields — even `circ_seqs: character(0)` keeps
# BSgenomeForge off the registry path.
circ_literal <- if (length(circ)) {
    sprintf("c(%s)", paste(sprintf('"%s"', circ), collapse = ", "))
} else {
    "character(0)"
}

seed_lines <- c(
    sprintf("Package: %s", PKG_NAME),
    sprintf("Title: %s", PKG_TITLE),
    sprintf("Description: %s", PKG_DESC),
    sprintf("Version: %s", PKG_VERSION),
    "organism: Anas platyrhynchos",
    "common_name: Pekin duck",
    "genome: ZJU1.0",
    "provider: NCBI",
    "release_date: 2020",
    "source_url: https://www.ncbi.nlm.nih.gov/assembly/GCF_015476345.1/",
    sprintf("organism_biocview: Anas_platyrhynchos"),
    "BSgenomeObjname: Anasplatyrhynchos",
    sprintf("seqnames: c(%s)", paste(sprintf('"%s"', fai$seqname), collapse = ", ")),
    sprintf("circ_seqs: %s", circ_literal),
    sprintf("seqs_srcdir: %s", SEQ_DIR),
    "seqfiles_suffix: .fa"
)
writeLines(seed_lines, seed_path)
ts("wrote seed DCF: ", seed_path, " (circ_seqs = ", circ_literal, ")")

# ── 3. forgeBSgenomeDataPkg ────────────────────────────────────────────
# Writes a package skeleton under destdir/<PKG_NAME>/ — including the
# single 2bit file.  If the package dir already exists we rebuild to
# pick up any seed edits (cheap relative to the initial split).
pkg_src <- file.path(PKG_DIR, PKG_NAME)
if (dir.exists(pkg_src)) {
    ts("removing stale package dir: ", pkg_src)
    unlink(pkg_src, recursive = TRUE, force = TRUE)
}

ts("forgeBSgenomeDataPkg — builds package skeleton + 2bit file")
BSgenomeForge::forgeBSgenomeDataPkg(
    x = seed_path,
    destdir = PKG_DIR,
    replace = TRUE,
    verbose = TRUE
)

# ── 4. R CMD build ─────────────────────────────────────────────────────
# From inside PKG_DIR so the tarball lands next to the source.
ts("R CMD build — builds source tarball")
build_cmd <- sprintf(
    "cd %s && R CMD build %s 2>&1",
    shQuote(TAR_DIR), shQuote(pkg_src)
)
build_out <- system(build_cmd, intern = TRUE)
ts("  last 10 lines of build output:")
writeLines(tail(build_out, 10))

tarball <- list.files(TAR_DIR, pattern = paste0("^", PKG_NAME, ".*\\.tar\\.gz$"), full.names = TRUE)
if (length(tarball) == 0) stop("No tarball produced by R CMD build")
tarball <- tarball[order(file.info(tarball)$mtime, decreasing = TRUE)][1]
ts("tarball: ", tarball)

# ── 5. R CMD INSTALL into ~/R/4.4.2/ ───────────────────────────────────
ts("R CMD INSTALL --library=~/R/4.4.2 ", tarball)
install_cmd <- sprintf(
    "R CMD INSTALL --library=%s %s 2>&1",
    shQuote(path.expand("~/R/4.4.2")), shQuote(tarball)
)
install_out <- system(install_cmd, intern = TRUE)
ts("  last 10 lines of install output:")
writeLines(tail(install_out, 10))

# ── 6. smoke test ──────────────────────────────────────────────────────
ts("smoke test — library(BSgenome.Anasplatyrhynchos.ZJU1.0)")
suppressPackageStartupMessages(library(BSgenome.Anasplatyrhynchos.ZJU1.0))
bsg <- getBSgenome(PKG_NAME)
ts("  seqlevels head: ", paste(head(seqlevels(bsg)), collapse = ", "))
ts("  seqinfo N    : ", length(seqlevels(bsg)))

ts("=== archr_build_bsgenome.R done ===")
