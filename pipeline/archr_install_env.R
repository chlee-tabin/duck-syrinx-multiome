#!/usr/bin/env Rscript
# archr_install_env.R — install ArchR + BSgenomeForge dependencies into ~/R/4.4.2/
#
# Run on an o2 compute node (interactive or sbatch).  NOT the login node —
# compile time for Bioconductor packages is ~20-40 min and login-partition
# quotas kill long-running R.
#
# Strategy:
#   1. Ensure BiocManager + remotes are present.
#   2. Install ArchR (GitHub: GreenleafLab/ArchR) — pins Bioconductor release.
#   3. Install BSgenomeForge + BSgenome (Bioconductor) for custom duck BSgenome.
#   4. Install peripheral ArchR deps that commonly fail (Cairo, chromVAR,
#      motifmatchr).  Not strictly required for peak calling.
#
# Size budget: existing lib = 938 MiB (2026-04-24); target ≤ 2.5 GiB total.
# HOME quota headroom at the time of writing was 23 GiB (quota-v2).
#
# Idempotent: every `install.packages` / `BiocManager::install` call checks
# whether the package already loads successfully; skips if so.

.libPaths(c("~/R/4.4.2", .libPaths()))

options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    Ncpus = max(1L, parallel::detectCores() - 1L),
    timeout = 600
)

cat("--- archr_install_env.R ---\n")
cat("R:       ", R.version$major, ".", R.version$minor, "\n", sep = "")
cat("libPaths:\n")
for (p in .libPaths()) cat("  ", p, "\n")
cat("Ncpus:   ", getOption("Ncpus"), "\n\n")

# ── helper: install only if not already importable ───────────────────────
# For installer="github", pass the repo spec as `source` (e.g.
# "GreenleafLab/ArchR@release_1.0.3") — `pkg` is used only for the
# pre/post-install import-check and message.  For "cran" / "bioc", the
# source defaults to pkg.
need <- function(pkg, installer, source = pkg, ...) {
    ok <- suppressPackageStartupMessages(
        requireNamespace(pkg, quietly = TRUE, lib.loc = .libPaths())
    )
    if (ok) {
        cat(sprintf("[skip] %s already present\n", pkg))
        return(invisible(TRUE))
    }
    cat(sprintf("[install] %s via %s (src=%s)\n", pkg, installer, source))
    t0 <- Sys.time()
    args <- list(...)
    res <- switch(
        installer,
        cran   = do.call(install.packages, c(list(source), args)),
        bioc   = do.call(BiocManager::install, c(list(source, update = FALSE, ask = FALSE), args)),
        github = do.call(remotes::install_github, c(list(source, upgrade = "never"), args)),
        stop(sprintf("Unknown installer %q", installer))
    )
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("[done]  %s (%.1fs)\n", pkg, dt))
    ok2 <- suppressPackageStartupMessages(
        requireNamespace(pkg, quietly = TRUE, lib.loc = .libPaths())
    )
    if (!ok2) stop(sprintf("FAILED: %s is not importable after install", pkg))
    invisible(res)
}

# ── 1. bootstrap installers ──────────────────────────────────────────────
need("BiocManager", "cran")
need("remotes",     "cran")
# Lightweight CLI parser used by archr_build_peaks.R.
need("optparse",    "cran")
# Data.table is an ArchR runtime dep; ensure it is present before ArchR install.
need("data.table",  "cran")

# Pin Bioc release; 3.20 pairs with R 4.4.2.
bioc_ver <- BiocManager::version()
cat(sprintf("Bioconductor version: %s\n", bioc_ver))

# ── 2. BSgenomeForge chain (custom duck genome) ──────────────────────────
# BSgenome (provides getBSgenome, forgeBSgenomeDataPkg)
need("BSgenome",       "bioc")
# BSgenomeForge — the modern forge API (Bioc 3.17+)
need("BSgenomeForge",  "bioc")

# ── 3. ArchR hard deps from Bioc (declared in ArchR DESCRIPTION).
# Without these, `install_github("GreenleafLab/ArchR@release_1.0.3")`
# fails with:
#   ERROR: dependencies 'chromVAR', 'ComplexHeatmap', 'motifmatchr',
#          'rhdf5', 'sparseMatrixStats' are not available for package 'ArchR'
# Install them explicitly here so the ArchR tarball build succeeds.
#
# Known O2 caveat (cluster software environment):
#   The R `curl` package is pinned at 5.2.1 because newer curl requires
#   libcurl >= 7.80, which o2 doesn't have.  We rely on `curl 5.2.1`
#   already being present.  Do NOT let this install chain upgrade it —
#   all our Bioc installs pass `update = FALSE`.
need("ComplexHeatmap",    "bioc")
need("rhdf5",             "bioc")
need("sparseMatrixStats", "bioc")
# Chain of chromVAR transitive deps.  Caveats:
#   * TFMPvalue >= 1.0.0 requires R >= 4.5.0 (o2 has R/4.4.2).  Install
#     from the CRAN archive via `remotes::install_version` (0.0.9 — the
#     last R-4.x-compatible release).
#   * pwalign / TFBSTools are Bioconductor, install via BiocManager.
# Order matters: TFMPvalue → TFBSTools → DT → chromVAR.
if (!requireNamespace("TFMPvalue", quietly = TRUE)) {
    cat("[install] TFMPvalue 0.0.9 via remotes::install_version\n")
    remotes::install_version("TFMPvalue", version = "0.0.9",
                             repos = "https://cloud.r-project.org",
                             upgrade = "never")
    if (!requireNamespace("TFMPvalue", quietly = TRUE)) {
        stop("FAILED: TFMPvalue is not importable after install")
    }
    cat("[done]  TFMPvalue (0.0.9)\n")
} else {
    cat("[skip] TFMPvalue already present\n")
}
need("pwalign",           "bioc")
need("TFBSTools",         "bioc")
need("DT",                "cran")
need("chromVAR",          "bioc")
need("motifmatchr",       "bioc")
# chromVARmotifs lives on GitHub (Greenleaf lab).  Needed transitively by
# ArchR; install_github sets upgrade="never" via need() so it won't
# clobber curl 5.2.1.
need("chromVARmotifs", "github", source = "GreenleafLab/chromVARmotifs")

# ── 4. ArchR core ───────────────────────────────────────────────────────
# ArchR lives on GitHub; we pin to tag release_1.0.3 for reproducibility
# (latest stable as of 2026-04).
#
# WORKAROUND — o2 curl pin:
#   ArchR 1.0.3 lists `devtools` in Imports.  devtools has a cascading
#   dep on `curl (>= 6)`, which needs libcurl >= 7.80 — and o2 has
#   7.76.1 (cluster system library).  devtools is used by ArchR only
#   for its .onAttach() banner and for `ArchR::installExtraPackages()`
#   (which we deliberately skip).  No peak-calling or Arrow code path
#   touches devtools.
#
# Strategy: download the ArchR source tarball, patch DESCRIPTION to
# drop `devtools` from Imports, then install from the patched tarball.
install_archr_without_devtools <- function() {
    if (requireNamespace("ArchR", quietly = TRUE)) {
        cat("[skip] ArchR already present\n")
        return(invisible(TRUE))
    }
    cat("[install] ArchR via patched tarball (drop devtools from Imports)\n")
    t0 <- Sys.time()
    tmp <- tempfile(fileext = ".tar.gz")
    # Use `curl -L` rather than R's download.file: the tag-archive URL
    # 302-redirects to codeload, and download.file on R 4.4.2 (default
    # method "libcurl") sometimes fails to follow the redirect cleanly.
    # release_1.0.3 is a branch (not a tag — the tag name is v1.0.3 but
    # we track the release branch head, which is what Greenleaf lab
    # treats as stable).  archive/refs/heads/<branch>.tar.gz is the
    # canonical GitHub URL form for a branch tarball.
    dl_cmd <- sprintf(
        "curl -L -s -o %s https://github.com/GreenleafLab/ArchR/archive/refs/heads/release_1.0.3.tar.gz",
        shQuote(tmp)
    )
    dl_rc <- system(dl_cmd)
    if (dl_rc != 0L || file.info(tmp)$size < 1e5) {
        stop(sprintf("FAILED curl download of ArchR tarball (rc=%d, size=%s)",
                     dl_rc, file.info(tmp)$size))
    }
    work <- tempfile("archr_src_")
    dir.create(work, recursive = TRUE)
    untar(tmp, exdir = work)
    # GitHub tarball extracts to ArchR-release_1.0.3/; move into ArchR/ so
    # R CMD INSTALL picks the right package name.
    extracted <- list.files(work, full.names = TRUE)
    stopifnot(length(extracted) == 1L)
    pkg_src <- file.path(work, "ArchR")
    file.rename(extracted[[1]], pkg_src)
    desc_path <- file.path(pkg_src, "DESCRIPTION")
    if (!file.exists(desc_path)) stop("ArchR DESCRIPTION missing after extract: ", desc_path)
    desc <- readLines(desc_path)
    # DESCRIPTION keeps a bulleted multi-line Imports block; we only
    # need to drop the line that is literally "    devtools,".
    before <- length(desc)
    desc <- desc[!grepl("^\\s*devtools\\s*,?\\s*$", desc)]
    cat(sprintf("   patched DESCRIPTION: %d → %d lines (removed devtools entry)\n",
                before, length(desc)))
    writeLines(desc, desc_path)

    # ArchR's .onAttach() iterates ArchRDependency and calls
    # packageVersion(pkg) on each member before attempting to load it.
    # If "devtools" is in the vector and isn't installed, .onAttach
    # aborts with `error: there is no package called 'devtools'`, which
    # R CMD INSTALL treats as a load failure → removes the package.
    # Drop "devtools" from the R-source vector too.
    global_defaults <- file.path(pkg_src, "R", "GlobalDefaults.R")
    if (file.exists(global_defaults)) {
        gdl <- readLines(global_defaults)
        gdl2 <- gsub('^(\\s*)"devtools"\\s*,?\\s*$', '', gdl, perl = TRUE)
        dropped <- sum(gdl != gdl2)
        if (dropped > 0) {
            writeLines(gdl2, global_defaults)
            cat(sprintf("   patched GlobalDefaults.R: removed %d 'devtools' entry\n", dropped))
        }
    } else {
        warning("GlobalDefaults.R not found in tarball — .onAttach may still fail")
    }
    # Install.
    inst_cmd <- sprintf("R CMD INSTALL --library=%s %s 2>&1",
                        shQuote(path.expand("~/R/4.4.2")), shQuote(pkg_src))
    out <- system(inst_cmd, intern = TRUE)
    cat("   last 8 lines of R CMD INSTALL:\n")
    writeLines(paste0("   ", tail(out, 8)))
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (!requireNamespace("ArchR", quietly = TRUE)) {
        stop("FAILED: ArchR is not importable after patched install")
    }
    cat(sprintf("[done]  ArchR patched (%.1fs)\n", dt))
    invisible(TRUE)
}

install_archr_without_devtools()

# NOTE: ArchR::installExtraPackages() is deliberately NOT called.
# That helper pulls `devtools`, which depends on `curl 6.x` — and
# curl 6.x needs libcurl >= 7.80 which o2 doesn't ship.  The peak-
# calling path used by archr_build_peaks.R does not need Cairo,
# presto, or Seurat (all "extras" added by installExtraPackages).
# If a downstream step needs them later, install them manually in
# an env that sidesteps the curl pin.

# ── 4. motif/chromatin peripherals (not strictly required for peak call
#      but used by forest plot, manhattan, chromVAR downstream).  These
#      can fail without blocking the peak-calling track — guarded.
try_need <- function(pkg, installer, ...) {
    res <- tryCatch(
        need(pkg, installer, ...),
        error = function(e) {
            message(sprintf("[skip/fail] %s: %s", pkg, conditionMessage(e)))
            NULL
        }
    )
    invisible(res)
}

try_need("chromVAR",     "bioc")
try_need("motifmatchr",  "bioc")
try_need("JASPAR2020",   "bioc")

# ── 5. sanity load of core packages ─────────────────────────────────────
cat("\n--- final load smoke test ---\n")
pkgs_to_load <- c("ArchR", "BSgenome", "BSgenomeForge",
                  "GenomicRanges", "Biostrings", "Rsamtools", "rtracklayer")
for (p in pkgs_to_load) {
    ok <- suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE))
    cat(sprintf("  %-16s %s\n", p, if (ok) "OK" else "MISSING"))
}

# Print install footprint for the next run to monitor.
lib_size <- system("du -sh ~/R/4.4.2 2>/dev/null | cut -f1", intern = TRUE)
if (length(lib_size) > 0) cat(sprintf("\n~/R/4.4.2 size: %s\n", lib_size))

cat("\n=== archr_install_env.R done ===\n")
