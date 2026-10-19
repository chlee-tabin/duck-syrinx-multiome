#!/usr/bin/env Rscript
# archr_pool_composition.R — what each ArchR peak-calling pool IS in the current
# annotation vocabulary.
#
# The 17 stage × cell-group pools that called the production peak set carry the
# labels in force at peak-calling time ("Somatic", "Splanchnic", "Chondrogenic
# LPM", an HH25 "Syringeal" that is a since-retired subcluster).  Those names are
# baked into proj_fig03_callable's cellColData$peak_pool and into
# peaks/reproducibility_metadata.tsv, which figS03 panel B1 (ms Fig S10) plots.
# This script measures — per pool, per nucleus — the CURRENT population of each
# pooled nucleus, so B1 can be labelled by what the pools contain rather than by
# a hand-written name crosswalk (figures/figS03/PROVENANCE.md carried one, but it
# could only bound the mixed pools, not measure them).
#
# Method (mirrors pipeline/chick/duck_pitx2_locus_percell_extract.py's SSOT
# re-derivation): cellColData(peak_pool) ⟕ f1 sidecar (arc_library, cb16 →
# d2g_stage, d2g_cluster, cell_type) ⟕ references/d2g_subcluster_annotations.tsv
# (stage, cluster → population).  Nuclei that carry a d2g (stage, leiden) pair
# absent from the SSOT are an ERROR (stale sidecar), not a silent drop.  Nuclei
# with no d2g assignment at all (the paraxial / epithelial pools that predate the
# LPM selection) are summarised by their global-atlas cell_type instead.
#
# Writes (under $SCRATCH/phase2/atac_archr/peaks/):
#   pool_composition.tsv         long: peak_pool, stage, pool_n_cells, label_kind,
#                                label, n_cells, frac        (→ figures/data/)
#   pool_cell_membership.tsv.gz  per nucleus: archr_cb, peak_pool, baked_label,
#                                d2g_stage, d2g_cluster, current_population, cell_type
# Compute: loads one ArchRProject + two TSVs; < 5 min, 1 core.
suppressPackageStartupMessages({
    library(ArchR)
    library(data.table)
})
ts <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))

OUT_ROOT <- file.path(SCRATCH_ROOT, "phase2/atac_archr")
PROJ_DIR <- file.path(OUT_ROOT, "proj_fig03_callable")
# The 2026-06 purge left proj_fig03_callable/ a skeleton; Save-ArchR-Project.rds was
# restored 2026-09-08 from STANDBY_FC/phase2_backup_20260512/atac_archr/ (written
# 2026-04-26 04:57, after peak_pool was attached and before peak calling).  The arrow
# files are gone, so read the object with readRDS — cellColData is all this needs.
PROJ_RDS <- file.path(PROJ_DIR, "Save-ArchR-Project.rds")
PEAK_DIR <- file.path(OUT_ROOT, "peaks")
F1_META  <- file.path(SCRATCH_ROOT, "phase2/atac/step_f1_coupled_atac_metadata.tsv")
SSOT     <- file.path(SCRATCH_ROOT, "references/d2g_subcluster_annotations.tsv")
# reproducibility_metadata.tsv was ferried into the repo; scratch peaks/ keeps only the BED.
REPRO    <- Filter(file.exists, c(file.path(PEAK_DIR, "reproducibility_metadata.tsv"),
                                  file.path(SCRATCH_ROOT, "figures/data/archr_reproducibility_metadata.tsv")))[1]
if (is.na(REPRO)) stop("reproducibility_metadata.tsv found neither under peaks/ nor figures/data/")
for (p in c(PROJ_RDS, F1_META, SSOT)) if (!file.exists(p)) stop("missing input: ", p)
ts("reproducibility metadata: ", REPRO)

ts("reading ArchRProject object ", PROJ_RDS)
proj <- readRDS(PROJ_RDS)
cc <- tryCatch(getCellColData(proj), error = function(e) {
    ts("  getCellColData() failed (", conditionMessage(e), ") — reading the cellColData slot directly")
    proj@cellColData })
cc <- as.data.table(as.data.frame(cc), keep.rownames = "archr_cb")
need <- c("peak_pool", "d2g_subcluster", "d2g_stage")
if (!all(need %in% names(cc))) stop("cellColData lacks: ", paste(setdiff(need, names(cc)), collapse = ", "))
cc <- cc[, .(archr_cb, peak_pool, baked_label = d2g_subcluster, pool_stage = d2g_stage)]
cc <- cc[!is.na(peak_pool)]
ts("  nuclei with a peak_pool: ", nrow(cc), " in ", uniqueN(cc$peak_pool), " pools")

# Every pool in reproducibility_metadata.tsv (the 17 that yielded a reproducible
# peak set — what B1 plots) must be in the project; the project also holds the
# callable pools that produced none (e.g. HH13|Early LPM), which are reported but
# carry no B1 bar.
repro <- fread(REPRO)
missing_pools <- setdiff(repro$peak_pool, unique(cc$peak_pool))
if (length(missing_pools))
    stop("pools in reproducibility_metadata.tsv absent from the project: ", paste(missing_pools, collapse = "; "))
extra <- setdiff(unique(cc$peak_pool), repro$peak_pool)
ts("  callable pools with no reproducible peak set (in project, not in B1): ", length(extra),
   if (length(extra)) paste0(" — ", paste(extra, collapse = "; ")) else "")
pool_n <- cc[, .N, by = peak_pool]
bad <- merge(pool_n, repro[, .(peak_pool, n_cells)], by = "peak_pool")[N != n_cells]
if (nrow(bad)) { cat("WARNING: pool cell counts differ from reproducibility_metadata.tsv:\n");
                 print(bad) } else ts("  pool membership matches reproducibility_metadata.tsv (all n_cells equal)")

ts("loading f1 sidecar ", F1_META)
meta <- fread(F1_META, sep = "\t", header = TRUE, na.strings = "")
if (names(meta)[1] %in% c("V1", "")) setnames(meta, 1L, "cb_fragments")
for (col in c("arc_library", "cb_fragments", "d2g_stage", "d2g_cluster", "cell_type"))
    if (!col %in% names(meta)) stop("f1 sidecar lacks column ", col)
# Same key construction as archr_build_peaks.R, so the join reproduces the
# original attach exactly.
meta[, archr_cb := paste0(arc_library, "#", cb_fragments)]
if (anyDuplicated(meta$archr_cb)) stop("f1 sidecar: duplicate archr_cb keys")
meta[, d2g_cluster := as.integer(round(as.numeric(d2g_cluster)))]
x <- merge(cc, meta[, .(archr_cb, d2g_stage, d2g_cluster, cell_type)], by = "archr_cb", all.x = TRUE)
ts("  pooled nuclei missing from the sidecar: ", sum(is.na(x$cell_type)))
x[is.na(cell_type), cell_type := "absent from current f1 sidecar"]

ssot <- fread(SSOT)
ssot[, cluster := as.integer(cluster)]
x <- merge(x, ssot[, .(d2g_stage = stage, d2g_cluster = cluster, current_population = population)],
           by = c("d2g_stage", "d2g_cluster"), all.x = TRUE)
# Fail loud: a (stage, leiden) that exists in the sidecar but not in the SSOT means
# the sidecar's d2g run and the annotation TSV disagree (stale era).
un <- x[!is.na(d2g_cluster) & is.na(current_population), unique(paste(d2g_stage, d2g_cluster))]
if (length(un)) stop("sidecar (stage, leiden) pairs with no SSOT entry: ", paste(un, collapse = "; "),
                     " — the f1 sidecar and references/d2g_subcluster_annotations.tsv are from different d2g runs")
st_mis <- x[!is.na(d2g_stage) & d2g_stage != pool_stage, .N]
ts("  nuclei whose current d2g_stage differs from the pool's stage: ", st_mis)

x[, label_kind := ifelse(!is.na(current_population), "d2g_population", "global_cell_type_outside_d2g")]
x[, label := ifelse(!is.na(current_population), current_population, cell_type)]
comp <- x[, .(n_cells = .N), by = .(peak_pool, label_kind, label)]
comp <- merge(comp, pool_n[, .(peak_pool, pool_n_cells = N)], by = "peak_pool")
comp[, frac := n_cells / pool_n_cells]
comp[, stage := sub("\\|.*$", "", peak_pool)]
setorder(comp, stage, peak_pool, -frac)
setcolorder(comp, c("peak_pool", "stage", "pool_n_cells", "label_kind", "label", "n_cells", "frac"))

out1 <- file.path(PEAK_DIR, "pool_composition.tsv")
fwrite(comp, out1, sep = "\t")
out2 <- file.path(PEAK_DIR, "pool_cell_membership.tsv.gz")
fwrite(x[, .(archr_cb, peak_pool, baked_label, d2g_stage, d2g_cluster, current_population, cell_type)],
       out2, sep = "\t", compress = "gzip")
ts("wrote ", out1, " (", nrow(comp), " rows) and ", out2, " (", nrow(x), " nuclei)")

# Human-readable summary in the log.
cat("\n=== pool composition (current vocabulary) ===\n")
for (p in c(repro$peak_pool, extra)) {
    d <- comp[peak_pool == p]
    cat(sprintf("%-32s n=%5d  %s\n", p, d$pool_n_cells[1],
                paste(sprintf("%s %.0f%%", d$label, 100 * d$frac), collapse = " | ")))
}
ts("=== archr_pool_composition.R done ===")
