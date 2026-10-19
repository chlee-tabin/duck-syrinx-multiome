#!/usr/bin/env Rscript
# =====================================================================
# substage_compendium_annotate.R — align old-duck E06/E09/E11 substage
# subcluster annotation to the COMPENDIUM method used for the multiome d2g
# (scripts/annotate_d2g_via_compendium.py) and chick (pipeline/chick/
# cAM_compendium_annotate.py), replacing the abandoned threshold classifier
# (annotate_and_plot.py::classify(), the PRRX1>=0.65 & FOXF1<0.25 priority rules).
#
# SAME ALGORITHM as the duck/chick scorers — U-method UEGs (Stein 2026:
# U = P_in - max P_out) + Wilcoxon DE + SIGNED lineage-panel scoring
# (up-panel minus down-panel, AddModuleScore) + margin-of-victory — implemented
# in R/Seurat to match the old-duck subclustering environment
# (substage_subcluster_v2.R).  The lineage compendium REUSES the LPM panels and
# ADDS the post-differentiation syringeal derivatives (cartilage / fibroblast /
# smooth muscle) absent at the earlier HH13-25 stages (a "tissue-appropriate
# panel" approach).
#
# Why: the threshold classifier mislabelled syringeal FIBROBLAST as "Somatic"
# (PRRX1-by-exclusion rule, never checked TBX5/TBX4) and dropped syringeal
# SMOOTH MUSCLE (ACTA2 0.79 / MYH11) to "ambiguous/mixed LPM" (no SM rule) —
# the exact priority-artifact + silent-drop failures the compendium method fixes.
#
# Inputs (scratch): substage/substage_E{06,09,11}_obj.rds  (Seurat; RNA_snn_res.0.2)
# Outputs (scratch substage/):
#   substage_compendium_clusterlabels.tsv      per (stage,cluster) label + evidence
#   substage_E*_percell_compendium.tsv.gz      per-cell (barcode,library,cluster,population)
# Runs as an O2 SLURM job.
# =====================================================================
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressWarnings(suppressMessages({
  .libPaths(c(paste0(HOME_ROOT, "/R/4.4.2"), .libPaths()))
  library(Seurat); library(Matrix)
}))
set.seed(1)
say <- function(...) cat(sprintf(...), "\n")

STAGES  <- c("E06", "E09", "E11")
RES_COL <- "RNA_snn_res.0.2"                 # canonical substage-v2 resolution
# SUBSTAGE_DIR env var selects the substage object dir (default = cr6); the cr9-CB
# rebase runs with SUBSTAGE_DIR=/n/scratch/.../old_duck_cr9/substage (run_cr9_compendium.sh).
OBJDIR  <- Sys.getenv("SUBSTAGE_DIR", paste0(SCRATCH_USER, "/old_duck_review/substage"))
OUT     <- OBJDIR
NKX61   <- "LOC113843492"                     # NKX6-1 in cr6.zju1/Gnomon

# ── Signed lineage panels (up, down, label) ──────────────────────────────
# Reused LPM panels (old-duck-appropriate splanchnic = FOXF1/GATA-family/ZFPM2,
# matching the substage_dotplot definition) + 3 differentiation panels.  Genes
# absent from a stage's object are dropped with a warning (panel still scores on
# the present genes), exactly like the d2g _score_panel.
PANELS <- list(
  Dorsal_subepi        = list(up=c(NKX61,"WNT2B","FOXD1","MSX1"),               down=c("HAND2","GATA4"),                       label="Dorsal sub-epi."),
  # PRRX1 removed from Ventral down: it is high in PRRX1+ mesenchyme generally and
  # wrongly tanked the Ventral score for E06 cl0, handing it to Dorsal by default.
  Ventral_spl          = list(up=c("FOXF1","GATA4","GATA5","GATA6","ZFPM2","HAND2"), down=c("BARX1"),                          label="Ventral spl."),
  Early_LPM            = list(up=c("TBX5","FOXF1","ISL1"),                       down=c(),                                      label="Early LPM"),
  # Skeletogenic LPM (was "Somatic"): the body-wall up-panel (TBX5/TBX4/PITX1) was
  # a misnomer — the identity vet (oldduck_somatic_identity_probe.R) showed the real
  # population is the pharyngeal-cranial SKELETOGENIC progenitor, defined by
  # EYA1/EYA2/EBF2 (PITX1~0, posterior-HOX-neg → NOT trunk somatopleure).  UP is now
  # IDENTICAL to the multiome Skeletogenic LPM panel (cross-platform consistency).
  # DOWN keeps SOX9/FNDC1 (+ HAND2) so the later, differentiated old-duck stages
  # (E09 cartilage SOX9/MATN4-high; E11 fibroblast) lose to Chondrogenic/Syringeal-
  # fibroblast — the stage-split then emerges from the data (E06=Skeletogenic LPM,
  # E09=Chondrogenic, E11=Syringeal fibroblast; Chondrogenic is already E09/E11-gated).
  Skeletogenic         = list(up=c("PRRX1","TWIST1","EBF2","EYA1","EYA2"),        down=c("GATA4","FOXF1","ZFPM2","HAND2","SOX9","FNDC1"), label="Skeletogenic LPM"),
  Epicardial           = list(up=c("WT1","TCF21"),                              down=c("FOXF1","GATA4","PRRX1"),               label="Mesothelium"),
  # NB: the multiome Splanchnic_vis ("Venous-pole spl.") panel is DROPPED for
  # old-duck — it has no venous pole / proepicardium in a syrinx dissection, and
  # its GATA4/ZFPM2-up + low-down signature spuriously captured E09 differentiating
  # cells (retained splanchnic memory) at 91%.
  # ── differentiation (E09/E11) ── RUNX2/WNT5A added so committing (not-yet-mature)
  # E09 cartilage still wins Chondrogenic over residual splanchnic.
  # label renamed 2026-08-22: "Chondrogenic LPM" -> "Syringeal chondrocyte"
  # (stage-appropriate at E09/E11 ~ HH35/37; triad with fibroblast / smooth muscle).
  # Panel KEY stays "Chondrogenic" (code identifier only).
  Chondrogenic         = list(up=c("SOX9","COL2A1","ACAN","CNMD","COL9A1","MATN1","RUNX2","WNT5A"), down=c(),                   label="Syringeal chondrocyte"),
  Syringeal_fibroblast = list(up=c("FNDC1","SCX","TNMD","LUM","POSTN","DCN","COL1A1","COL3A1"),
                                                                                 down=c("SOX9","COL2A1","ACAN","GATA4","FOXF1","TBX5"), label="Syringeal fibroblast"),
  Syringeal_smoothmuscle = list(up=c("ACTA2","MYH11","TAGLN","CNN1","DES","PDGFRB"), down=c(),                                  label="Syringeal smooth muscle")
)
STAGE_GATED <- list(
  Early_LPM              = c("E06"),
  # Skeletogenic LPM is the EARLY pharyngeal-skeletal progenitor state (E06 = the
  # cross-platform bridge from multiome HH25 Skeletogenic LPM).  By E09 it has
  # differentiated → Chondrogenic LPM / Syringeal fibroblast (identity-probe: EYA1
  # 52%@E06 → 5%@E09).  Gate to E06 so the PRRX1/TWIST1-driven (pan-mesenchyme)
  # spurious wins on E09/E11 fibroblast clusters are suppressed and the trajectory
  # reads cleanly.  Mirrors the Early_LPM gate.
  Skeletogenic           = c("E06"),
  Chondrogenic           = c("E09","E11"),
  Syringeal_fibroblast   = c("E09","E11"),
  Syringeal_smoothmuscle = c("E09","E11")
)
MIN_MARGIN_STRONG <- 0.10
MIN_MARGIN_MOD    <- 0.05

# Marker-specific lineages: each defined by a pathognomonic marker.  Instead of a
# hardcoded fraction floor (NKX6-1>=0.5 etc.), these panels may win ONLY if their
# diagnostic marker is among the cluster's top data-driven UEGs (U-method).  This
# both gates the "win-on-absence" artifact AND reclaims a strongly-marked cluster,
# is rank-based (robust to ambient detection), and removes every per-marker number
# — the sole tunable is UEG_TOPK, a conventional top-markers count.
MARKER_DIAGNOSTIC <- list(
  Dorsal_subepi          = c("LOC113843492", "WNT2B"),   # NKX6-1 / WNT2B
  Epicardial             = c("WT1", "TCF21"),
  Syringeal_smoothmuscle = c("ACTA2", "MYH11")
)
UEG_TOPK <- 30

# ── helpers ──────────────────────────────────────────────────────────────
score_panel <- function(obj, genes, key) {
  # signed AddModuleScore component: control-matched mean of present genes,
  # 0 if none present.  Returns per-cell numeric vector.
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0) return(setNames(rep(0, ncol(obj)), colnames(obj)))
  obj <- AddModuleScore(obj, features = list(present), name = key, ctrl = 50, seed = 1)
  v <- obj[[paste0(key, "1")]][, 1]
  setNames(v, colnames(obj))
}

uegs_for_stage <- function(det, clusters) {
  # det: genes x cells 0/1 detection (dgCMatrix); clusters: factor per cell.
  # U = P_in - max P_out per (gene,cluster).  Returns list cluster -> data.frame.
  cl_levels <- levels(clusters)
  # per-cluster detection fraction matrix: genes x clusters
  frac <- sapply(cl_levels, function(cl) {
    idx <- which(clusters == cl)
    Matrix::rowMeans(det[, idx, drop = FALSE])
  })
  colnames(frac) <- cl_levels
  out <- list()
  for (cl in cl_levels) {
    p_in  <- frac[, cl]
    p_out <- if (length(cl_levels) > 1) apply(frac[, setdiff(cl_levels, cl), drop = FALSE], 1, max) else 0
    u <- p_in - p_out
    ord <- order(u, decreasing = TRUE)
    out[[cl]] <- data.frame(gene = rownames(frac)[ord], u = u[ord],
                            p_in = p_in[ord], stringsAsFactors = FALSE)
  }
  out
}

compose_label <- function(stage, panel_scores, markers, ueg_genes) {
  # panel_scores: named numeric (panel-key -> mean signed score for this cluster).
  # markers: named numeric (gene -> detection fraction in this cluster).
  g <- function(x) if (x %in% names(markers)) as.numeric(markers[[x]]) else 0
  # (0) contamination pre-check (mirrors the d2g/cAM upstream mesenchyme filter):
  # epithelial / endothelial / immune / neural-crest clusters are not LPM lineage.
  if (g("EPCAM") >= 0.30 || g("CDH5") >= 0.20 || g("PECAM1") >= 0.20 ||
      g("PTPRC") >= 0.20 || g("SOX10") >= 0.20 || g("FOXD3") >= 0.20)
    return(list(label = "Contamination", top = "contam_precheck", top_score = NA_real_,
                runner = NA, margin = NA_real_, confidence = "strong"))
  # (0b) cycling QC pre-check: ccreg (linear S/G2M-score regression) is applied
  # upstream but does not fully dissolve proliferation structure in the late,
  # highly-mitotic E09/E11 tissue, so small cell-cycle-driven clusters persist —
  # the same residual the multiome quarantines via its apoptotic QC class. UEG-based
  # (consistent with the marker-diagnostic rule): if the cluster's top UEGs are
  # dominated by cell-cycle genes, quarantine it rather than force a lineage.
  CC_GENES <- c("MKI67","TOP2A","PCNA","CCNB2","BUB1B","CDC20","UBE2C","PTTG1",
                "NUSAP1","KIF15","TTK","SGO1","CENPF","CENPA","AURKB","CDK1","CKS2",
                "HMGB2","CENPE","KIF23","CCNB1","CDK1","TYMS","RRM2")
  if (sum(head(ueg_genes, 6) %in% CC_GENES) >= 2)
    return(list(label = "Cycling", top = "cycling_precheck", top_score = NA_real_,
                runner = NA, margin = NA_real_, confidence = "strong"))

  s <- panel_scores
  for (p in names(STAGE_GATED)) if (p %in% names(s) && !(stage %in% STAGE_GATED[[p]])) s[p] <- -Inf
  # (1) UEG corroboration (replaces the fixed fraction floors): a marker-specific
  # panel (Dorsal=NKX6-1, Mesothelium=WT1, Smooth muscle=ACTA2/MYH11) wins on
  # near-zero down-markers alone (the "win-on-absence" artifact). Require its
  # diagnostic marker to be among the cluster's top data-driven UEGs, else gate it.
  for (panel in names(MARKER_DIAGNOSTIC)) {
    if (panel %in% names(s) && !any(MARKER_DIAGNOSTIC[[panel]] %in% ueg_genes))
      s[[panel]] <- -Inf
  }
  s <- s[is.finite(s)]
  ranked <- sort(s, decreasing = TRUE)
  top <- names(ranked)[1]; top_s <- ranked[1]
  runner <- if (length(ranked) > 1) names(ranked)[2] else NA
  runner_s <- if (length(ranked) > 1) ranked[2] else 0
  margin <- top_s - runner_s
  conf <- if (margin >= MIN_MARGIN_STRONG) "strong" else if (margin >= MIN_MARGIN_MOD) "moderate" else "weak"
  list(label = PANELS[[top]]$label, top = top, top_score = top_s,
       runner = runner, margin = margin, confidence = conf)
}

# ── main ─────────────────────────────────────────────────────────────────
all_cluster_rows <- list()
for (stage in STAGES) {
  objf <- file.path(OBJDIR, sprintf("substage_%s_obj.rds", stage))
  if (!file.exists(objf)) { say("MISSING %s -> skip", objf); next }
  say("== %s ==  loading %s", stage, objf)
  obj <- readRDS(objf)
  DefaultAssay(obj) <- "RNA"
  if (!RES_COL %in% colnames(obj[[]])) stop(sprintf("%s lacks %s", stage, RES_COL))
  clusters <- factor(as.character(obj[[RES_COL]][, 1]))
  obj$.cl <- clusters
  say("  %d cells, %d clusters (%s)", ncol(obj), nlevels(clusters), RES_COL)

  # (1) signed lineage scores per cell -> mean per cluster
  panel_cell <- list()
  for (pk in names(PANELS)) {
    up   <- score_panel(obj, PANELS[[pk]]$up,   paste0("up_", pk))
    down <- score_panel(obj, PANELS[[pk]]$down, paste0("dn_", pk))
    panel_cell[[pk]] <- up - down
    pres_up <- intersect(PANELS[[pk]]$up, rownames(obj))
    say("    panel %-22s up=%d/%d present", pk, length(pres_up), length(PANELS[[pk]]$up))
  }
  panel_mat <- do.call(cbind, panel_cell)                # cells x panels
  colnames(panel_mat) <- names(PANELS)
  cl_means <- t(sapply(levels(clusters), function(cl)
    colMeans(panel_mat[clusters == cl, , drop = FALSE])))   # clusters x panels

  # (2) UEGs (U-method) per cluster — Seurat v4 'slot' / v5 'layer' robust
  det <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
  det <- as(det > 0, "CsparseMatrix")
  uegs <- uegs_for_stage(det, clusters)
  # per-cluster detection fractions for the contamination pre-check + dorsal floor
  # single source of truth: the gate is derived from CONTAM_GATE, not a hand-copied duplicate
  # (a second list silently rots if a contamination marker is added/removed in one place only).
  CONTAM_GATE <- c("EPCAM","CDH5","PECAM1","PTPRC","SOX10","FOXD3")          # contamination markers
  DIAGNOSTIC  <- c("LOC113843492","WT1","TBX5","FOXF1","SOX9","ACTA2")        # diagnostic markers
  CONTAM <- c(CONTAM_GATE, DIAGNOSTIC)
  # the contamination pre-check treats a gene absent from `markers` as 0% detected, which
  # silently disables that gene's contamination clause — warn loudly so an inert gate is visible
  # rather than passing every cluster through it (gene-symbol drift is plausible on cr6->cr9).
  missing_contam <- setdiff(CONTAM_GATE, rownames(det))
  if (length(missing_contam))
    say("    WARNING: contamination-gate genes absent from namespace (their clause is inert): %s",
        paste(missing_contam, collapse = ", "))
  need <- intersect(CONTAM, rownames(det))
  cl_frac <- sapply(levels(clusters), function(cl)
    Matrix::rowMeans(det[need, clusters == cl, drop = FALSE]))
  rownames(cl_frac) <- need

  # (3) Wilcoxon top markers per cluster (evidence only)
  Idents(obj) <- clusters
  wil <- tryCatch(
    FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                   verbose = FALSE),
    error = function(e) { say("    FindAllMarkers failed: %s", conditionMessage(e)); NULL })

  # (4) compose per-cluster label + evidence
  for (cl in levels(clusters)) {
    ps <- setNames(as.numeric(cl_means[cl, ]), colnames(cl_means))
    markers <- setNames(as.numeric(cl_frac[, cl]), rownames(cl_frac))
    dg <- function(x) { i <- match(x, names(markers)); if (is.na(i)) NA_real_ else markers[i] }
    diag <- sprintf("NKX6-1=%.2f WT1=%.2f TBX5=%.2f FOXF1=%.2f SOX9=%.2f ACTA2=%.2f",
                    dg("LOC113843492"), dg("WT1"), dg("TBX5"), dg("FOXF1"), dg("SOX9"), dg("ACTA2"))
    ueg_genes <- head(uegs[[cl]]$gene[uegs[[cl]]$u > 0], UEG_TOPK)   # data-driven top UEGs
    ann <- compose_label(stage, ps, markers, ueg_genes)
    top_uegs <- paste(head(uegs[[cl]]$gene, 6), collapse = ",")
    top_de <- if (!is.null(wil) && cl %in% wil$cluster)
      paste(head(wil$gene[wil$cluster == cl], 6), collapse = ",") else ""
    n_cl <- sum(clusters == cl)
    all_cluster_rows[[length(all_cluster_rows) + 1]] <- data.frame(
      stage = stage, cluster = cl, n_cells = n_cl,
      population = ann$label, top_panel = ann$top,
      top_score = round(ann$top_score, 4), runner = ann$runner,
      margin = round(ann$margin, 4), confidence = ann$confidence,
      diag = diag, top_uegs = top_uegs, top_de = top_de, stringsAsFactors = FALSE)
    say("  cl%-3s n=%5d  %-24s [%s m=%.3f]  %s  UEG: %s", cl, n_cl,
        ann$label, ann$confidence, ann$margin, diag, top_uegs)
  }

  # per-cell population output for the xstage re-export
  st_rows <- do.call(rbind, Filter(function(d) d$stage[1] == stage, all_cluster_rows))
  pop_map <- setNames(st_rows$population, as.character(st_rows$cluster))
  md <- obj[[]]
  bc <- if ("barcode" %in% colnames(md)) md$barcode else rownames(md)
  lib <- if ("library" %in% colnames(md)) md$library else NA
  pc <- data.frame(barcode = bc, library = lib,
                   cluster = as.character(clusters),
                   population = unname(pop_map[as.character(clusters)]),
                   stage = stage, stringsAsFactors = FALSE)
  pcf <- file.path(OUT, sprintf("substage_%s_percell_compendium.tsv.gz", stage))
  gz <- gzfile(pcf, "w"); write.table(pc, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  say("  wrote %s (%d cells, %d populations)", pcf, nrow(pc), length(unique(pc$population)))
  rm(obj); gc()
}

out <- do.call(rbind, all_cluster_rows)
# Repo destination (copied back from the cluster):
# figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv (note the cr9_ prefix)
outf <- file.path(OUT, "substage_compendium_clusterlabels.tsv")
write.table(out, outf, sep = "\t", quote = FALSE, row.names = FALSE)
say("\nwrote %s (%d clusters)", outf, nrow(out))
say("\n== population composition (compendium) ==")
print(as.data.frame(table(out$stage, out$population)))
cat("SBATCH_DONE_OK\n")
