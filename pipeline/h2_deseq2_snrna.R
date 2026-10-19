#!/usr/bin/env Rscript
# Phase H Step 2: paired DESeq2 L vs R per (stage × subcluster_scope) on
# the snRNA gene-level pseudobulks produced by h1.
#
# Counterpart of g2 (ATAC) for the snRNA-seq DEG track.  Same paired
# design (~ genotype + side, contrast L vs R) and the same poscounts
# size factor estimator family for sparse single-cell pseudobulk.
#
# Two scope variants per stage:
#   1. POOL              all inference-pool subclusters per genotype × library × side × stage
#   2. PER_SUBCLUSTER    one DESeq2 run per (stage, subcluster) bucket
# (Scope 3 SPLANCHNIC_RELAXED — HH13-only MIN=15 exploratory — retired
#  2026-05-23: redundant under --min-cells-mode=principled (HH13 floor → 20),
#  hardcoded the now-stale "Splanchnic" subcluster name (renamed to
#  "Venous-pole spl." on 2026-05-04), and had no callers in pipeline/ or
#  scripts/.  Consumer block in generate_volcano.R was deleted.)
#
# Stage rollout: --stages HH13 (default) for the feasibility prototype;
#                --stages HH13,HH17/16,HH21,HH25 for the full sweep later.
#
# Inputs (figures/data/):
#   h1_pseudobulk_snrna_counts.tsv.gz
#   h1_pseudobulk_snrna_metadata.tsv
#
# Outputs (figures/data/):
#   h2_deg_<scope>_<stage_safe>.tsv          (per scope × stage, gene-sorted by padj)
#     Schema (14 cols):
#       gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
#       scope, stage, design, n_samples, n_L_samples, n_R_samples,
#       n_paired_genotypes
#     `design` ∈ {"paired", "unpaired", "unpaired+sex",
#                 "unpaired (fallback from unpaired+sex)"}; constant
#     within a file, recording the formula actually fit (downstream
#     consumers should filter on this column rather than infer from
#     filename — see build_share.py::assert_design_provenance).
#   h2_deg_summary.tsv                        (one row per scope × stage with
#                                              n_samples, n_paired_genotypes,
#                                              n_DEGs at multiple padj cutoffs)
#
# Run on O2:
#   sbatch pipeline/run_pipeline.sh --from h2 --to h2
# (runs every invocation: paired × {LPM, excluded, dorsal} + unpaired+sex × {LPM, excluded}.
#  DESeq2 supplied by ~/R/4.4.2/ — see pipeline/install_o2_r_packages.R.)
#
# Or locally (only if pseudobulk TSVs synced):
#   conda activate duckfig    # has bioconductor-deseq2
#   Rscript pipeline/h2_deseq2_snrna.R --stages HH13
#
# Flags:
#   --stages         comma-separated stage list (default: HH13)
#   --min-cells      fixed-mode floor for POOL + PER_SUBCLUSTER (default: 30)
#   --min-cells-mode "fixed" (use --min-cells globally) or "principled".
#                    Principled = max(20, 0.05 × median(n_cells)) per scope,
#                    computed independently for POOL and each PER_SUBCLUSTER
#                    bucket.  Recovers small-pseudobulk subclusters that
#                    fixed=30 silently excludes (HH13/HH17/HH16/HH25
#                    Splanchnic; HH17/16 Cardiopharyngeal LPM) while keeping
#                    a sensible 20-cell absolute floor.  Default: "principled"
#                    (legacy "fixed" mode kept only for parity diagnostics).
#   --pseudobulk-dir override input directory (default: figures/data/)
#   --output-dir     override output directory (default: figures/data/)
.libPaths(c("~/R/4.4.2", .libPaths()))

# Print warnings as they occur (default warn=0 buffers them until script
# exit, which makes the per-bucket fallback warning at line ~245 invisible
# in real time on long SLURM logs).  See M1.
options(warn = 1)

suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
    library(dplyr)
    library(tidyr)
})

# ── Args ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
    idx <- match(paste0("--", name), args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[idx + 1]
}
get_flag <- function(name) {
    paste0("--", name) %in% args
}

stages <- strsplit(get_arg("stages", "HH13"), ",")[[1]]
min_cells <- as.integer(get_arg("min-cells", "30"))
# Principled is the canonical mode (default).  Use --min-cells-mode=fixed
# to fall back to the global --min-cells floor (legacy, only for parity
# diagnostics).
min_cells_mode <- get_arg("min-cells-mode", "principled")
if (!min_cells_mode %in% c("fixed", "principled")) {
    stop(sprintf("--min-cells-mode must be 'fixed' or 'principled', got '%s'",
                 min_cells_mode))
}
unpaired <- get_flag("unpaired")
control_sex <- get_flag("control-sex")
if (control_sex && !unpaired) {
    stop("--control-sex requires --unpaired (sex is collinear with genotype in the paired design and would produce a rank-deficient model).",
         call. = FALSE)
}
e1b_path <- get_arg("e1b-tsv", "figures/data/e1b_stage_assignments.tsv")
in_dir <- get_arg("pseudobulk-dir", "figures/data")
out_dir <- get_arg("output-dir", "figures/data")
# Optional: override h1-output filenames + h2-output prefix.  Defaults
# preserve the LPM-pipeline behavior; the excluded and dorsal
# sister-pipelines override these — see pipeline/run_pipeline.sh's H2
# block for the per-scope invocations.
counts_name <- get_arg("counts-name", "h1_pseudobulk_snrna_counts.tsv.gz")
meta_name   <- get_arg("meta-name",   "h1_pseudobulk_snrna_metadata.tsv")
out_prefix  <- get_arg("out-prefix",  "h2_deg")
# Guard against silent paired-output overwrite when the user runs the
# unpaired model without overriding --out-prefix.
if (unpaired && !grepl("unpaired", out_prefix)) {
    stop(sprintf(
        "Refusing to write unpaired model output under prefix '%s' — would overwrite paired files. Pass --out-prefix h2_deg_unpaired (or h2_deg_unpaired_excluded for the variant runner).",
        out_prefix
    ))
}

cat("=== Phase H Step 2: DESeq2 ===\n")
cat(sprintf("  stages:             %s\n", paste(stages, collapse = ", ")))
cat(sprintf("  min_cells_mode:     %s\n", min_cells_mode))
cat(sprintf("  min_cells (fixed):  %d\n", min_cells))
cat(sprintf("  unpaired:           %s\n", unpaired))
cat(sprintf("  control_sex:        %s\n", control_sex))
cat(sprintf("  in_dir:             %s\n", in_dir))
cat(sprintf("  out_dir:            %s\n", out_dir))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load h1 pseudobulks ──────────────────────────────────────────────────
counts_path <- file.path(in_dir, counts_name)
meta_path   <- file.path(in_dir, meta_name)
if (!file.exists(counts_path)) stop(sprintf("Missing: %s", counts_path))
if (!file.exists(meta_path))   stop(sprintf("Missing: %s", meta_path))

cat(sprintf("Loading %s ...\n", counts_path))
counts <- read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
# CellBender-corrected master h5ad can produce small negative values
# in .X (ambient subtraction).  DESeq2 requires non-negative integer
# counts.  Round-floor: max(0, round(x)).  This is the standard fix
# for CellBender-corrected pseudobulks.
n_neg <- sum(counts < 0)
if (n_neg > 0) {
    cat(sprintf("  Clipping %d negative cells (CellBender residuals) to 0\n",
                n_neg))
    counts[counts < 0] <- 0
}
counts <- round(counts)
storage.mode(counts) <- "integer"
cat(sprintf("  %d genes × %d samples\n", nrow(counts), ncol(counts)))

cat(sprintf("Loading %s ...\n", meta_path))
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  %d sample rows\n", nrow(meta)))

# Sanity: every column in counts must exist as sample_id in meta
stopifnot(all(colnames(counts) %in% meta$sample_id))
meta <- meta[match(colnames(counts), meta$sample_id), ]
stopifnot(all(colnames(counts) == meta$sample_id))

# Load FIXED e1b sex map: genotype_id -> {M, F}.  Used when
# control_sex=TRUE to add sex as a fixed-effect covariate.  Sex is
# perfectly collinear with genotype in the paired design (each
# genotype = one embryo = one sex), so paired model can't add it
# without singularity.  In unpaired model, ~sex + side disentangles
# Z-chromosome dosage from L/R asymmetry (CFC1 is Z-linked — see
# docs/cfc1_nodal_lr_review.md §4 hypothesis D, not distributed).
sex_map <- list()
if (file.exists(e1b_path)) {
    e1b <- read.delim(e1b_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
    sex_map <- setNames(e1b$sex, as.character(e1b$genotype_id))
    cat(sprintf("Loaded sex map: %d genotypes (M=%d, F=%d, ?=%d)\n",
                length(sex_map),
                sum(sex_map == "M", na.rm = TRUE),
                sum(sex_map == "F", na.rm = TRUE),
                sum(is.na(sex_map) | sex_map == "?")))
} else if (control_sex) {
    stop(sprintf("ERROR: --control-sex requires e1b sidecar at %s",
                 e1b_path))
}

# Fail-fast on unresolved sex.  '?' encodes upstream sex-inference failure
# (e.g., FC_08656/20260329_step07c_pseudobulk_pca_subsets.py:164 default).
# Without this guard a bucket with M + '?' would be counted as 2 sex levels
# and fed to `~ sex + side`, treating '?' as an artificial covariate stratum
# .  Surfacing the data-integrity issue upstream
# is preferred over masking it in the design choice.
if (control_sex) {
    unresolved <- names(sex_map)[is.na(sex_map) | sex_map == "?"]
    if (length(unresolved) > 0) {
        stop(sprintf(
            "--control-sex requires confident M/F sex for every genotype, but %d genotype(s) in %s have unresolved sex (NA or '?'):\n%s\nResolve upstream (e.g., FC_08656/20260329_step07c_pseudobulk_pca_subsets.py) before rerunning h2.",
            length(unresolved), e1b_path,
            paste(sprintf("  - %s", unresolved), collapse = "\n")
        ), call. = FALSE)
    }
}
meta$sex <- sex_map[as.character(meta$genotype)]

# ── Helper: compute MIN_CELLS floor under principled rule ───────────────
#
# Principled rule: max(20, 0.05 × median(n_cells)) computed over the
# pseudobulk buckets passing into a single DESeq2 run.  20 is a hard floor
# (DESeq2 size-factor estimation degrades sharply below 20-cell pseudobulks
# regardless of cohort scale); 5% of median scales the floor up for very
# large pseudobulk distributions (e.g., HH25 Ventral spl. median=784 →
# floor=39).  Returns the fixed value when mode == "fixed".
compute_floor <- function(scope_meta, mode, fixed_value) {
    if (mode == "fixed") return(fixed_value)
    med <- median(scope_meta$n_cells)
    max(20L, as.integer(round(med * 0.05)))
}

# ── Helper: filter samples by stage / scope and run paired DESeq2 ────────
#
# Returns a tibble with one row per gene of (gene_id, log2FoldChange, lfcSE,
# pvalue, padj, baseMean, scope, stage, n_samples, n_paired_genotypes).
# Returns NULL if there are too few paired genotypes for a paired model.
run_deseq2_paired <- function(scope_name, stage, subset_meta, floor) {
    keep <- subset_meta$n_cells >= floor
    sub  <- subset_meta[keep, , drop = FALSE]
    if (nrow(sub) < 4) {
        cat(sprintf("  [%s × %s] Skip: only %d samples pass MIN=%d\n",
                    scope_name, stage, nrow(sub), floor))
        return(NULL)
    }

    if (unpaired) {
        # Unpaired model: keep all samples (singletons + pairs alike).
        # Need >=2 of each side for any contrast.
        n_L <- sum(sub$side == "L")
        n_R <- sum(sub$side == "R")
        if (n_L < 2 || n_R < 2) {
            cat(sprintf(
                "  [%s × %s] Skip (unpaired): n_L=%d, n_R=%d (need >=2 each side)\n",
                scope_name, stage, n_L, n_R))
            return(NULL)
        }
        n_paired <- NA  # not relevant in unpaired mode
        cat(sprintf("  [%s × %s] %d samples (unpaired; n_L=%d, n_R=%d); ",
                    scope_name, stage, nrow(sub), n_L, n_R))
    } else {
        # Paired model: drop genotypes with only one side.
        geno_sides <- table(sub$genotype, sub$side)
        paired_genos <- rownames(geno_sides)[
            rowSums(geno_sides > 0) >= 2
        ]
        if (length(paired_genos) < 2) {
            cat(sprintf("  [%s × %s] Skip: only %d paired genotypes (need >=2)\n",
                        scope_name, stage, length(paired_genos)))
            return(NULL)
        }
        sub <- sub[sub$genotype %in% paired_genos, , drop = FALSE]
        n_paired <- length(paired_genos)
        cat(sprintf("  [%s × %s] %d samples / %d paired genotypes; ",
                    scope_name, stage, nrow(sub), n_paired))
    }
    cts <- counts[, sub$sample_id, drop = FALSE]

    cd <- data.frame(
        sample = sub$sample_id,
        genotype = factor(sub$genotype),
        sex = factor(sub$sex),
        side = factor(sub$side, levels = c("R", "L")),
        n_cells = sub$n_cells,
        row.names = sub$sample_id
    )

    # Choose design + matching label in one if/else if/else chain so the
    # two stay structurally aligned.
    # NA-aware sex-level count: setdiff drops NA so a bucket with all-male
    # samples plus some unknown-sex samples doesn't slip past the >=2
    # threshold and feed an NA factor level into DESeq2.  '?' is caught
    # upstream by the fail-fast at sex-map load; the
    # NA guard here covers the residual case of a sample whose genotype
    # isn't present in sex_map.
    n_sex_levels <- length(setdiff(unique(cd$sex), NA))
    if (unpaired && control_sex && n_sex_levels >= 2) {
        design_formula <- ~ sex + side
        design_label <- "unpaired+sex"
    } else if (unpaired && control_sex) {
        stop(sprintf(
            "[%s × %s] only %d non-NA sex level(s) present; ~ sex + side is not estimable. Not falling back to ~ side silently.",
            scope_name, stage, n_sex_levels
        ), call. = FALSE)
    } else if (unpaired) {
        design_formula <- ~ side
        design_label <- "unpaired"
    } else {
        design_formula <- ~ genotype + side
        design_label <- "paired"
    }
    dds <- DESeqDataSetFromMatrix(
        countData = cts,
        colData = cd,
        design = design_formula
    )

    # Drop genes with zero counts across all samples — DESeq2 chokes on
    # them with poscounts (no positive count to estimate from).
    dds <- dds[rowSums(counts(dds)) > 0, ]
    cat(sprintf("%d genes after zero-row filter\n", nrow(dds)))

    # poscounts size factors handle sparse single-cell pseudobulks where
    # MOR (median-of-ratios) returns zeros from rows that are all-zero in
    # at least one sample.  Mirrors g2 ATAC convention.
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds, fitType = "parametric", quiet = TRUE)
    res <- results(dds, name = "side_L_vs_R", alpha = 0.05)

    # Convert to tibble; drop NA padj rows (independent filtering rejects)
    n_L_samples <- sum(sub$side == "L")
    n_R_samples <- sum(sub$side == "R")
    res_df <- as.data.frame(res) |>
        tibble::rownames_to_column("gene_id") |>
        dplyr::arrange(padj) |>
        dplyr::mutate(
            scope = scope_name,
            stage = stage,
            design = design_label,
            n_samples = nrow(sub),
            n_L_samples = n_L_samples,
            n_R_samples = n_R_samples,
            n_paired_genotypes = n_paired
        )
    list(per_gene = res_df,
         summary = data.frame(
             scope = scope_name, stage = stage,
             design = design_label,
             n_samples = nrow(sub),
             n_L_samples = n_L_samples,
             n_R_samples = n_R_samples,
             n_paired_genotypes = n_paired,
             n_genes_tested = nrow(dds),
             n_DEG_padj_05 = sum(res_df$padj < 0.05, na.rm = TRUE),
             n_DEG_padj_01 = sum(res_df$padj < 0.01, na.rm = TRUE),
             n_DEG_padj_001 = sum(res_df$padj < 0.001, na.rm = TRUE)
         ))
}

# ── Iterate stages × scopes ──────────────────────────────────────────────
stage_safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
all_summaries <- list()

for (stage in stages) {
    cat(sprintf("\n=== Stage %s ===\n", stage))
    stage_meta <- meta[meta$stage_anchor == stage, , drop = FALSE]
    if (nrow(stage_meta) == 0) {
        cat("  No samples for this stage; skipping.\n")
        next
    }
    cat(sprintf("  Stage cohort: %d samples across %d genotypes / "
                ,
                nrow(stage_meta), length(unique(stage_meta$genotype))))
    cat(sprintf("%d subclusters\n",
                length(unique(stage_meta$subcluster))))

    # ── Scope 1: POOL — one pseudobulk per genotype × library × side ────
    # Our h1 emits at the (genotype × library × side × subcluster × stage)
    # granularity already — to "pool", we sum across subclusters within
    # the same (genotype, library, side, stage) bucket.  Aggregate counts
    # column-by-column by re-keying.
    cat("  -- scope: POOL (sum across inference-pool subclusters) --\n")
    # Use "###" as the component-list separator (does NOT appear in any
    # sample_id field).  Earlier "+" separator collided with subcluster
    # names that contain a "+" (e.g., "BARX1+") at HH17/16 / HH21.
    pool_meta <- stage_meta |>
        dplyr::group_by(genotype, library_m, side, stage_anchor) |>
        dplyr::summarise(
            sample_id_list = paste(sort(sample_id), collapse = "###"),
            n_cells   = sum(n_cells),
            firmness_norm = paste(unique(firmness_norm), collapse = ";"),
            subcluster = "POOL",
            sex = dplyr::first(sex),
            n_stage_mismatch_cells = sum(n_stage_mismatch_cells),
            .groups = "drop"
        )
    # Synthesize a SINGLE pool sample_id used as both the column name in
    # pool_counts and the row identifier in DESeq2 colData.
    pool_meta$sample_id <- paste(pool_meta$genotype, pool_meta$library_m,
                                 pool_meta$side, "POOL", stage,
                                 sep = "|")
    pool_meta$stage_anchor <- stage
    # Build pooled counts matrix by summing component sample columns
    pool_counts <- matrix(0L, nrow = nrow(counts), ncol = nrow(pool_meta),
                          dimnames = list(rownames(counts), pool_meta$sample_id))
    for (i in seq_len(nrow(pool_meta))) {
        component_ids <- strsplit(pool_meta$sample_id_list[i], "###")[[1]]
        pool_counts[, i] <- as.integer(rowSums(counts[, component_ids, drop = FALSE]))
    }
    # Splice POOL counts into the global counts matrix temporarily —
    # via a small helper closure.
    counts_orig <- counts
    counts <- cbind(counts, pool_counts)

    pool_floor <- compute_floor(pool_meta, min_cells_mode, min_cells)
    if (min_cells_mode == "principled") {
        cat(sprintf("  POOL principled floor: %d (median n_cells = %d)\n",
                    pool_floor, as.integer(median(pool_meta$n_cells))))
    }
    res_pool <- run_deseq2_paired("POOL", stage, pool_meta, pool_floor)
    if (!is.null(res_pool)) {
        out_path <- file.path(out_dir,
                              sprintf("%s_POOL_%s.tsv", out_prefix, stage_safe(stage)))
        write.table(res_pool$per_gene, out_path, sep = "\t",
                    quote = FALSE, row.names = FALSE)
        cat(sprintf("  Wrote %s (%d genes)\n", out_path, nrow(res_pool$per_gene)))
        all_summaries[[length(all_summaries) + 1]] <- res_pool$summary
    }

    counts <- counts_orig

    # ── Scope 2: PER_SUBCLUSTER ──────────────────────────────────────────
    cat("  -- scope: PER_SUBCLUSTER --\n")
    for (sc_name in unique(stage_meta$subcluster)) {
        sc_meta <- stage_meta[stage_meta$subcluster == sc_name, , drop = FALSE]
        scope_label <- paste0("SUB:", sc_name)
        sc_floor <- compute_floor(sc_meta, min_cells_mode, min_cells)
        if (min_cells_mode == "principled") {
            cat(sprintf("  [%s] principled floor: %d (median n_cells = %d)\n",
                        sc_name, sc_floor,
                        as.integer(median(sc_meta$n_cells))))
        }
        res_sc <- run_deseq2_paired(scope_label, stage, sc_meta, sc_floor)
        if (!is.null(res_sc)) {
            sc_safe <- gsub("[^A-Za-z0-9]", "_", sc_name)
            out_path <- file.path(out_dir,
                sprintf("%s_SUB_%s_%s.tsv", out_prefix, sc_safe, stage_safe(stage)))
            write.table(res_sc$per_gene, out_path, sep = "\t",
                        quote = FALSE, row.names = FALSE)
            cat(sprintf("  Wrote %s (%d genes)\n",
                        out_path, nrow(res_sc$per_gene)))
            all_summaries[[length(all_summaries) + 1]] <- res_sc$summary
        }
    }

    # (Scope 3 SPLANCHNIC_RELAXED retired 2026-05-23 — see file-header docstring.)
}

# ── Write summary table ──────────────────────────────────────────────────
if (length(all_summaries) > 0) {
    summary_df <- do.call(rbind, all_summaries)
    summary_path <- file.path(out_dir, sprintf("%s_summary.tsv", out_prefix))
    write.table(summary_df, summary_path, sep = "\t",
                quote = FALSE, row.names = FALSE)
    cat(sprintf("\nWrote summary: %s\n", summary_path))
    print(summary_df)
} else {
    cat("\nNo DESeq2 runs produced output.\n")
}

cat("\n=== Done ===\n")
