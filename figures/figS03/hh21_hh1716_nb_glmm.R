#!/usr/bin/env Rscript

# Negative-binomial mixed-model arm test for the HH21 positive control and the
# paired HH17/16 sensitivity cohort.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("lme4", "glmmTMB", "emmeans")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Could not resolve this script's path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
counts_path <- file.path(repo_root, "figures", "data", "hh21_hh1716_response_counts.tsv")
halfmeta_path <- file.path(repo_root, "figures", "data", "hh21_hh1716_halfmeta.tsv")
output_path <- file.path(repo_root, "figures", "data", "hh21_hh1716_nb_glmm.tsv")

for (path in c(counts_path, halfmeta_path)) {
  if (!file.exists(path)) {
    stop("Required input does not exist: ", path)
  }
}

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

assert_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
}

parse_bool <- function(x, label) {
  if (is.logical(x)) {
    if (anyNA(x)) stop(label, " contains NA.")
    return(x)
  }
  y <- tolower(trimws(as.character(x)))
  if (anyNA(y) || any(!y %in% c("true", "false"))) {
    bad <- unique(as.character(x)[is.na(y) | !y %in% c("true", "false")])
    stop(label, " contains invalid boolean value(s): ", paste(bad, collapse = ", "))
  }
  y == "true"
}

key_string <- function(x, columns) {
  do.call(paste, c(x[columns], sep = "\r"))
}

counts <- read_tsv(counts_path)
halfmeta <- read_tsv(halfmeta_path)
# Content fingerprint of BOTH inputs this model reads, carried into the output so a consumer can
# prove the fit and the scores drawn beside it came from one module-test run. Counts alone is not
# enough: halfmeta supplies every offset and the paired flags, and a stale same-tier halfmeta
# would otherwise pass unnoticed. Byte-concatenation must match the
# producer's hashlib.sha256(counts_bytes + halfmeta_bytes).
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Missing required R package: digest (needed for the module-input fingerprint)")
}
.read_bytes <- function(path) readBin(path, "raw", file.info(path)$size)
# "sha:" prefix: a bare 12-hex digest can look numeric (13038807e193 parses as 1.3e+200) and
# read.delim would coerce the column to double, refusing a perfectly fresh set. Must match
# figures/figS03/_module_inputs.SHA_PREFIX.
MODULE_INPUTS_SHA <- paste0("sha:", substr(
  digest::digest(c(.read_bytes(counts_path), .read_bytes(halfmeta_path)),
                 algo = "sha256", serialize = FALSE),
  1L, 12L))

count_columns <- c("stage", "floor", "embryo", "side", "peak_id", "arm", "is_ref",
                   "pitx2_locus", "count")
meta_columns <- c(
  "stage", "floor", "embryo", "side", "size_factor", "total_fragments", "n_cells", "paired"
)
assert_columns(counts, count_columns, "Response-count table")
assert_columns(halfmeta, meta_columns, "Half-metadata table")
# Response-arm peak universe, derived from the counts file itself (July run: 212; 2026-08
# seeded relabel: 91; after the PITX2-locus hold-out: the independent responders only).
# Derived AFTER assert_columns, so a missing `arm` column raises the named column error
# rather than silently evaluating to 0.
EXPECTED_PEAKS <- length(unique(counts$peak_id[counts$arm %in% c("L-up", "R-up")]))
EXPECTED_ARM <- vapply(c("R-up", "L-up"),
                       function(a) length(unique(counts$peak_id[counts$arm == a])), integer(1))
counts <- counts[count_columns]
halfmeta <- halfmeta[meta_columns]

counts$is_ref <- parse_bool(counts$is_ref, "counts$is_ref")
counts$pitx2_locus <- parse_bool(counts$pitx2_locus, "counts$pitx2_locus")
# The pseudobulk tier is the contested parameter, so the fit must record which one it
# modelled — otherwise a panel cannot assert that its dots and its subtitle share a tier.
FLOOR <- unique(c(counts$floor, halfmeta$floor))
if (length(FLOOR) != 1L) {
  stop("counts/halfmeta mix pseudobulk tiers: ", paste(FLOOR, collapse = ", "))
}
halfmeta$paired <- parse_bool(halfmeta$paired, "halfmeta$paired")

if (anyNA(counts) || anyNA(halfmeta)) {
  stop("Input table contains NA values.")
}
if (any(!is.finite(counts$count)) || any(counts$count < 0) ||
    any(counts$count != round(counts$count))) {
  stop("counts$count must contain non-negative finite integers.")
}
if (any(!is.finite(halfmeta$size_factor)) || any(halfmeta$size_factor <= 0)) {
  stop("halfmeta$size_factor must be finite and strictly positive.")
}
if (any(!is.finite(halfmeta$total_fragments)) || any(halfmeta$total_fragments <= 0)) {
  stop("halfmeta$total_fragments must be finite and strictly positive.")
}
if (any(!is.finite(halfmeta$n_cells)) || any(halfmeta$n_cells <= 0)) {
  stop("halfmeta$n_cells must be finite and strictly positive.")
}

if (!setequal(unique(counts$stage), c("HH21", "HH17/16")) ||
    !setequal(unique(halfmeta$stage), c("HH21", "HH17/16"))) {
  stop("Both inputs must contain exactly stages HH21 and HH17/16.")
}
if (any(!counts$side %in% c("R", "L")) || any(!halfmeta$side %in% c("R", "L"))) {
  stop("side must contain only R and L.")
}
if (any(duplicated(key_string(counts, c("stage", "embryo", "side", "peak_id"))))) {
  stop("Response-count table has duplicate (stage, embryo, side, peak_id) rows.")
}
if (any(duplicated(key_string(halfmeta, c("stage", "embryo", "side"))))) {
  stop("Half-metadata table has duplicate (stage, embryo, side) rows.")
}

# Arm contract (module test, 2026-08-18): the response arms are the INDEPENDENT responder
# promoters; PITX2's own elements are split into the c promoter (PITX2c_ref), the a/b
# bilateral control (PITX2ab_ctrl) and its remaining locus elements (PITX2_locus_other).
# `is_ref` stays GENE-level, so it is TRUE for both PITX2 promoter arms.
PITX2_ARMS <- c("PITX2c_ref", "PITX2ab_ctrl")
RESPONSE_ARMS <- c("R-up", "L-up")
if (any(counts$is_ref != (counts$arm %in% PITX2_ARMS))) {
  stop("is_ref and the PITX2 promoter arms disagree.")
}
if (any(!counts$arm %in% c(RESPONSE_ARMS, PITX2_ARMS, "PITX2_locus_other"))) {
  stop("Unexpected arm value in response-count table.")
}
arm_peaks_per_stage <- function(arm) {
  vapply(c("HH21", "HH17/16"),
         function(stage) length(unique(counts$peak_id[counts$stage == stage & counts$arm == arm])),
         FUN.VALUE = integer(1))
}
ref_counts <- arm_peaks_per_stage("PITX2c_ref")
ab_counts <- arm_peaks_per_stage("PITX2ab_ctrl")
if (any(ref_counts != 1L) || any(ab_counts != 2L)) {
  stop(
    "Expected 1 PITX2c reference peak and 2 PITX2a/b control peaks per stage; observed c=",
    paste(names(ref_counts), ref_counts, sep = "=", collapse = ", "), " ab=",
    paste(names(ab_counts), ab_counts, sep = "=", collapse = ", ")
  )
}

joined <- merge(
  counts,
  halfmeta,
  by = c("stage", "embryo", "side"),
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)
if (nrow(joined) != nrow(counts)) {
  stop("Count-to-halfmeta join changed the row count.")
}
if (anyNA(joined$size_factor) || anyNA(joined$total_fragments) ||
    anyNA(joined$n_cells) || anyNA(joined$paired)) {
  bad <- unique(joined[is.na(joined$size_factor), c("stage", "embryo", "side")])
  stop(
    "Count-to-halfmeta join failed for at least one row: ",
    paste(apply(bad, 1L, paste, collapse = "/"), collapse = ", ")
  )
}

# Remove EVERY PITX2-locus element (the c promoter, the a/b controls and the +5.4 kb element)
# before any sampling or model fitting: the model tests whether the OTHER response genes'
# promoters carry the directionality, so PITX2's own cis-elements cannot sit inside an arm.
# The contrast for what the held-out element contributed is the previous committed run, which
# had it inside the L-up arm: there the HH21 embryo interaction slope was identifiable
# (+0.19 log2, p 0.019); with it out, the slope goes singular in every engine x offset and the
# directionality p is suppressed. Recorded in figS03/PROVENANCE.md, not re-run here.
model_input <- joined[joined$arm %in% RESPONSE_ARMS, ]
# Test what the message claims: `is_ref` alone covers only PITX2's own promoters, so the +5.4 kb
# element could re-enter a response arm and every guard here would still pass — and that element
# is precisely what made the HH21 slope identifiable before it was held out (a wrong headline).
if (any(model_input$is_ref) || any(model_input$pitx2_locus) ||
    any(!model_input$arm %in% RESPONSE_ARMS)) {
  stop("PITX2-locus peaks were not completely excluded.")
}
# The three checks above all test flags MINTED BY ONE PRODUCER from one window literal, so a
# mislabelled flag passes them all — and re-admitting the +5.4 kb element is exactly what turns
# the suppressed directionality back into a significant p=0.019 headline. Check the COORDINATES
# independently, against the same BED f5/f6/g3/the scan read.
locus_bed <- file.path(repo_root, "references", "pitx2_promoter_regions.bed")
if (!file.exists(locus_bed)) {
  stop("Cannot verify the PITX2-locus hold-out: ", locus_bed, " is missing.")
}
bed <- read.delim(locus_bed, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
locus_row <- bed[bed[[4]] == "locus_50kb|PITX2", , drop = FALSE]
if (nrow(locus_row) != 1L) {
  stop("Expected exactly one 'locus_50kb|PITX2' row in ", locus_bed, "; found ", nrow(locus_row))
}
locus_chrom <- locus_row[[1]][1]; locus_lo <- locus_row[[2]][1]; locus_hi <- locus_row[[3]][1]
pk <- regmatches(model_input$peak_id,
                 regexec("^(.+):([0-9]+)-([0-9]+)$", model_input$peak_id))
if (any(vapply(pk, length, integer(1)) != 4L)) {
  stop("Un-parseable peak_id(s) in the model input; cannot verify the PITX2-locus hold-out.")
}
pk_chrom <- vapply(pk, `[`, character(1), 2L)
pk_start <- as.numeric(vapply(pk, `[`, character(1), 3L))
pk_end <- as.numeric(vapply(pk, `[`, character(1), 4L))
inside <- pk_chrom == locus_chrom & pk_start < locus_hi & pk_end > locus_lo
if (any(inside)) {
  stop("Peak(s) inside the PITX2 locus window reached a response arm despite their flags: ",
       paste(unique(model_input$peak_id[inside]), collapse = ", "),
       " (window ", locus_chrom, ":", locus_lo, "-", locus_hi, " from ", basename(locus_bed), ")")
}

prepare_stage <- function(stage) {
  x <- model_input[model_input$stage == stage, ]
  if (stage == "HH17/16") {
    x <- x[x$paired, ]
  } else if (stage == "HH21") {
    if (any(!x$paired)) stop("HH21 contains an unpaired half.")
    x <- x[x$paired, ]
  } else {
    stop("Unexpected stage: ", stage)
  }

  # Guardrails are CONSISTENCY checks against the model input's own halfmeta / response set,
  # not literals: the July run had 8 / 4 paired embryos and 212 promoter peaks; the 2026-08
  # seeded-relabel run at the floor-20 tier has 8 / 8 and 90 (91 before the PITX2-locus
#   hold-out). A literal here would (and did)
  # stop the script on the very reruns it exists for.
  expected_embryos <- length(unique(halfmeta$embryo[halfmeta$stage == stage & halfmeta$paired]))
  if (expected_embryos < 3L) stop(stage, ": fewer than 3 paired embryos in halfmeta (", expected_embryos, ")")
  halves <- unique(x[c("stage", "embryo", "side")])
  n_embryo <- length(unique(x$embryo))
  n_peak <- length(unique(x$peak_id))
  if (n_embryo != expected_embryos) {
    stop(stage, ": expected ", expected_embryos, " paired embryos; observed ", n_embryo)
  }
  if (nrow(halves) != 2L * expected_embryos) {
    stop(stage, ": every sampled embryo must contribute exactly two halves.")
  }
  side_counts <- table(halves$embryo, halves$side)
  if (!all(c("R", "L") %in% colnames(side_counts)) || any(side_counts[, c("R", "L")] != 1L)) {
    stop(stage, ": every sampled embryo must have exactly one R and one L half.")
  }
  if (n_peak != EXPECTED_PEAKS) {
    stop(stage, ": expected ", EXPECTED_PEAKS, " non-reference response peaks (from the counts' response set); observed ", n_peak)
  }
  arm_peak_counts <- vapply(
    c("R-up", "L-up"),
    function(arm) length(unique(x$peak_id[x$arm == arm])),
    FUN.VALUE = integer(1)
  )
  if (!identical(unname(arm_peak_counts), unname(EXPECTED_ARM))) {
    stop(
      stage, ": expected R-up=", EXPECTED_ARM[["R-up"]], " and L-up=", EXPECTED_ARM[["L-up"]],
      " peaks (the counts' response set); observed ",
      paste(names(arm_peak_counts), arm_peak_counts, sep = "=", collapse = ", ")
    )
  }
  if (nrow(x) != nrow(halves) * n_peak) {
    stop(stage, ": expected one count per sampled half and response peak.")
  }

  x$side <- factor(x$side, levels = c("R", "L"))
  x$arm <- factor(x$arm, levels = c("R-up", "L-up"))
  x$embryo <- factor(x$embryo)
  x$peak_id <- factor(x$peak_id)
  # Numeric interaction column (= the armL-up:sideL design column, 1 only on the
  # Left half of an L-up peak). This is the embryo random SLOPE we care about:
  # the estimand is the side x arm INTERACTION (L-up-minus-R-up directionality),
  # so to keep the embryo the inferential unit FOR THAT ESTIMAND the embryo
  # random effect must carry the interaction slope, (1 + int_LL | embryo) — NOT
  # a common-side slope (1 + side | embryo), which lets the interaction borrow
  # false precision across the ~200 promoters and is anti-conservative for the
  # directionality (verified: se 0.025 vs the honest 0.044; p 2e-8 vs 2e-4).
  # Treatment coding is used because it fits cleanly (non-singular) in BOTH
  # engines. A sum-coded, reference-invariant slope is a verified sensitivity
  # (the directionality inference does not depend on the coding) — verified once
  # during development (coding-invariance), not repeated here.
  x$int_LL <- as.integer(x$arm == "L-up") * as.integer(x$side == "L")
  if (anyNA(x$side) || anyNA(x$arm)) {
    stop(stage, ": factor coding introduced NA values.")
  }
  if (!identical(levels(x$side), c("R", "L")) ||
      !identical(levels(x$arm), c("R-up", "L-up"))) {
    stop(stage, ": side/arm reference-level guardrail failed.")
  }
  x
}

stage_data <- list(
  HH21 = prepare_stage("HH21"),
  `HH17/16` = prepare_stage("HH17/16")
)

cat("Validated model cohorts (PITX2 reference peaks excluded):\n")
cohort_table <- do.call(
  rbind,
  lapply(names(stage_data), function(stage) {
    x <- stage_data[[stage]]
    data.frame(
      stage = stage,
      n_obs = nrow(x),
      n_embryo = length(unique(x$embryo)),
      n_halves = nrow(unique(x[c("embryo", "side")])),
      n_peak = length(unique(x$peak_id)),
      check.names = FALSE
    )
  })
)
print(cohort_table, row.names = FALSE)

# Target random structure: the embryo carries the INTERACTION slope
# (1 + int_LL | embryo) so the directionality (side x arm) — the estimand — is
# inferred against between-embryo variance in the interaction (embryo = the
# inferential unit for it). (1 + side | peak_id) gives per-promoter baseline +
# per-promoter L/R heterogeneity (promoter dependence preserved). If the
# interaction slope is unidentifiable (a singular fit), the ladder falls back to
# (1 | embryo), whose directionality SE is anti-conservative — flagged by
# ranef_structure so the fallback is transparent and the label-swap floor is the
# honest anchor there.
# An embryo random-SLOPE variance needs a reasonable number of groups to be
# reliably estimated; with only ~4 the estimate is at the boundary and
# engine-dependent. Require >= 6 embryos before ATTEMPTING the interaction slope.
# On the floor-20 tier both stages now clear that gate (8 pairs each), so the
# a-priori skip no longer fires for either — a non-identifiable trough fit is
# recorded by the ladder (converged = FALSE) rather than pre-empted.
# See fit_with_fallback().
MIN_EMBRYO_FOR_SLOPE <- 6L

ranef_ladder <- list(
  list(
    label = "(1 + int_LL | embryo) + (1 + side | peak_id)",
    formula = "(1 + int_LL | embryo) + (1 + side | peak_id)",
    has_random_slope = TRUE
  ),
  list(
    label = "(1 | embryo) + (1 + side | peak_id)",
    formula = "(1 | embryo) + (1 + side | peak_id)",
    has_random_slope = TRUE
  ),
  list(
    label = "(1 | embryo) + (1 | peak_id)",
    formula = "(1 | embryo) + (1 | peak_id)",
    has_random_slope = FALSE
  )
)

tmb_is_singular <- function(fit, tol = 1e-4) {
  variance_components <- glmmTMB::VarCorr(fit)$cond
  if (length(variance_components) == 0L) return(FALSE)
  any(vapply(variance_components, function(component) {
    covariance <- as.matrix(component)
    if (any(!is.finite(covariance))) return(TRUE)
    eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(eigenvalues)) || any(eigenvalues < -sqrt(.Machine$double.eps))) {
      return(TRUE)
    }
    min(sqrt(pmax(eigenvalues, 0))) < tol
  }, FUN.VALUE = logical(1)))
}

format_messages <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("none")
  paste(gsub("[\r\n]+", " ", as.character(x)), collapse = " | ")
}

fit_once <- function(engine, formula, data, stage, offset_type, step_number) {
  fit_warnings <- character()
  fit_error <- NULL
  fit <- withCallingHandlers(
    tryCatch(
      {
        if (engine == "glmer.nb") {
          set.seed(1)
          lme4::glmer.nb(
            formula,
            data = data,
            control = lme4::glmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 2e5)
            )
          )
        } else if (engine == "glmmTMB") {
          glmmTMB::glmmTMB(
            formula,
            data = data,
            family = glmmTMB::nbinom2,
            control = glmmTMB::glmmTMBControl(
              optCtrl = list(iter.max = 20000L, eval.max = 20000L)
            )
          )
        } else {
          stop("Unknown engine: ", engine)
        }
      },
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      cat(
        sprintf(
          "  WARNING [%s | %s | %s | step %d]: %s\n",
          stage, engine, offset_type, step_number, conditionMessage(w)
        )
      )
      invokeRestart("muffleWarning")
    }
  )

  if (is.null(fit)) {
    cat(
      sprintf(
        "  RESULT  [%s | %s | %s | step %d]: convergence_code=ERROR; singular=NA; message=%s\n",
        stage, engine, offset_type, step_number, fit_error
      )
    )
    return(list(
      fit = NULL,
      converged = FALSE,
      singular = NA,
      convergence_code = "ERROR",
      warnings = fit_warnings,
      error = fit_error
    ))
  }

  warning_is_bad <- any(grepl(
    paste(
      "failed to converge", "false convergence", "non-positive-definite",
      "degenerate Hessian", "iteration limit",
      "unable to evaluate scaled gradient",
      "estimate truncated at zero",   # MASS::theta.ml NB dispersion hit its 0 boundary
      sep = "|"),
    fit_warnings,
    ignore.case = TRUE
  ))
  if (engine == "glmer.nb") {
    opt_code <- fit@optinfo$conv$opt
    if (is.null(opt_code) || length(opt_code) == 0L) opt_code <- 0L
    conv_messages <- fit@optinfo$conv$lme4$messages
    non_singularity_messages <- conv_messages[
      !grepl("boundary \\(singular\\) fit", conv_messages, ignore.case = TRUE)
    ]
    converged <- all(opt_code == 0L) &&
      length(non_singularity_messages) == 0L && !warning_is_bad
    singular <- lme4::isSingular(fit, tol = 1e-4)
    convergence_code <- paste0(
      "opt=", paste(opt_code, collapse = ","),
      "; lme4=", format_messages(conv_messages)
    )
  } else {
    opt_code <- fit$fit$convergence
    pd_hessian <- isTRUE(fit$sdr$pdHess)
    converged <- length(opt_code) == 1L && opt_code == 0L && pd_hessian && !warning_is_bad
    singular <- tmb_is_singular(fit, tol = 1e-4)
    convergence_code <- paste0(
      "opt=", paste(opt_code, collapse = ","),
      "; pdHess=", pd_hessian,
      "; message=", format_messages(fit$fit$message)
    )
  }
  cat(
    sprintf(
      "  RESULT  [%s | %s | %s | step %d]: convergence_code=%s; converged=%s; singular=%s\n",
      stage, engine, offset_type, step_number, convergence_code, converged, singular
    )
  )
  list(
    fit = fit,
    converged = converged,
    singular = singular,
    convergence_code = convergence_code,
    warnings = fit_warnings,
    error = NULL
  )
}

fit_with_fallback <- function(data, stage, engine, offset_type) {
  if (!offset_type %in% c("size_factor", "total_fragments")) {
    stop("Unknown offset type: ", offset_type)
  }
  data$offset_value <- data[[offset_type]]
  if (any(!is.finite(data$offset_value)) || any(data$offset_value <= 0)) {
    stop(stage, " ", offset_type, ": offset must be finite and strictly positive.")
  }

  # Only ATTEMPT the embryo interaction slope when there are enough embryos to
  # identify its variance. At n=4 the slope is not reliably identifiable and the
  # two engines disagree at the optimizer boundary (glmer.nb reports singular ->
  # fallback; glmmTMB sits on a non-singular boundary and keeps it) — which would
  # make the directionality inference engine-DEPENDENT. Gating the ladder START on
  # n_embryo (not on each optimizer's singularity verdict) forces BOTH engines to
  # the same starting point, so the model directionality p is suppressed
  # CONSISTENTLY whenever the cohort is too small. (On the floor-20 tier both
  # stages have 8 pairs, so this gate is currently inert; it fires only if a
  # future tier/cohort drops a stage below MIN_EMBRYO_FOR_SLOPE.)
  n_embryo_stage <- length(unique(data$embryo))
  first_step <- if (n_embryo_stage >= MIN_EMBRYO_FOR_SLOPE) 1L else 2L
  cat(sprintf("\n[%s | %s | %s] n_embryo=%d; starting at ladder step %d (%s).\n",
              stage, engine, offset_type, n_embryo_stage, first_step,
              if (first_step == 1L) "target interaction slope"
              else "interaction slope skipped a priori (too few embryos)"))

  for (step_number in seq.int(first_step, length(ranef_ladder))) {
    structure <- ranef_ladder[[step_number]]
    formula <- as.formula(paste0(
      "count ~ arm * side + ", structure$formula,
      " + offset(log(offset_value))"
    ))
    cat(
      sprintf(
        "  FIT     [%s | %s | %s | step %d]: %s\n",
        stage, engine, offset_type, step_number, structure$label
      )
    )
    attempt <- fit_once(engine, formula, data, stage, offset_type, step_number)

    acceptable <- !is.null(attempt$fit) && isTRUE(attempt$converged) &&
      (!structure$has_random_slope || !isTRUE(attempt$singular))
    if (acceptable) {
      attempt$ranef_structure <- structure$label
      attempt$step_number <- step_number
      return(attempt)
    }

    if (step_number < length(ranef_ladder)) {
      reason <- if (is.null(attempt$fit)) {
        "fit error"
      } else if (!isTRUE(attempt$converged)) {
        "non-convergence"
      } else {
        "singular random-slope fit"
      }
      cat(
        sprintf(
          "  FALLBACK [%s | %s | %s]: %s; advancing to step %d.\n",
          stage, engine, offset_type, reason, step_number + 1L
        )
      )
    } else {
      # Ladder exhausted. Since 2026-08-17 this is RECORDED (converged = FALSE, plus a
      # "(ladder exhausted; NOT converged)" suffix on `ranef_structure` — there is no separate
      # flag column) rather than fatal, so a non-converging arm (the HH17/16 trough on the
      # seeded relabel) does not black-hole the whole table; the final guardrail below
      # requires convergence for the HH21 positive control only. A fit that errored outright
      # (fit NULL) still stops.
      if (is.null(attempt$fit)) {
        stop(stage, " ", engine, " ", offset_type,
             ": fallback ladder exhausted with a fit ERROR (not just non-convergence): ", attempt$error)
      }
      cat(sprintf("  UNCONVERGED [%s | %s | %s]: ladder exhausted; recording the last attempt with converged=FALSE (last code: %s)\n",
                  stage, engine, offset_type, attempt$convergence_code))
      attempt$ranef_structure <- paste0(structure$label, " (ladder exhausted; NOT converged)")
      attempt$step_number <- step_number
      attempt$converged <- FALSE
      return(attempt)
    }
  }
  stop(stage, " ", engine, " ", offset_type, ": internal fallback-ladder error.")
}

# emmeans contrast -> data.frame of point estimates on the link scale. We consume
# only `estimate` + `SE` (make_row builds the CI/p itself from est +- crit*se with
# the small-sample t reference) plus the orientation-label columns the callers
# assert on. emmeans' own CI/p/df are dropped: we don't request `infer=c(TRUE,TRUE)`,
# so no CI columns are produced; any default p/z columns are simply ignored.
summary_frame <- function(contrast_grid, context) {
  ans <- as.data.frame(summary(contrast_grid, adjust = "none", type = "link"))
  if (!all(c("estimate", "SE") %in% names(ans))) {
    stop(context, ": unexpected emmeans summary schema: ", paste(names(ans), collapse = ", "))
  }
  ans
}

one_sided_p <- function(estimate, p_two_sided, alternative) {
  if (!is.finite(estimate) || !is.finite(p_two_sided) || p_two_sided < 0 || p_two_sided > 1) {
    stop("Cannot derive a one-sided p-value from non-finite/invalid emmeans output.")
  }
  if (alternative == "greater") {
    p <- if (estimate >= 0) p_two_sided / 2 else 1 - p_two_sided / 2
  } else if (alternative == "less") {
    p <- if (estimate <= 0) p_two_sided / 2 else 1 - p_two_sided / 2
  } else {
    stop("Unknown one-sided alternative: ", alternative)
  }
  min(1, max(0, p))
}

fixed_effects <- function(fit, engine) {
  if (engine == "glmer.nb") lme4::fixef(fit) else glmmTMB::fixef(fit)$cond
}

extract_effects <- function(fit_result, data, stage, engine, offset_type) {
  fit <- fit_result$fit
  side_means <- emmeans::emmeans(fit, ~ side | arm, type = "link")
  side_contrasts <- emmeans::contrast(side_means, method = "revpairwise", adjust = "none")
  side_summary <- summary_frame(
    side_contrasts,
    paste(stage, engine, offset_type, "per-arm")
  )
  if (nrow(side_summary) != 2L ||
      !identical(as.character(side_summary$arm), c("R-up", "L-up")) ||
      any(as.character(side_summary$contrast) != "L - R")) {
    stop(stage, " ", engine, " ", offset_type, ": emmeans per-arm contrast orientation failed.")
  }

  joint_means <- emmeans::emmeans(fit, ~ side * arm, type = "link")
  direction_contrast <- emmeans::contrast(
    joint_means,
    interaction = c("revpairwise", "revpairwise"),
    adjust = "none"
  )
  direction_summary <- summary_frame(
    direction_contrast,
    paste(stage, engine, offset_type, "directionality")
  )
  if (nrow(direction_summary) != 1L ||
      as.character(direction_summary$side_revpairwise) != "L - R" ||
      as.character(direction_summary$arm_revpairwise) != "(L-up) - (R-up)") {
    stop(stage, " ", engine, " ", offset_type, ": emmeans interaction orientation failed.")
  }

  beta <- fixed_effects(fit, engine)
  required_beta <- c("sideL", "armL-up:sideL")
  if (!all(required_beta %in% names(beta))) {
    stop(
      stage, " ", engine, " ", offset_type,
      ": required fixed-effect coefficients missing: ",
      paste(setdiff(required_beta, names(beta)), collapse = ", ")
    )
  }
  expected <- c(
    `R-up` = unname(beta[["sideL"]]),
    `L-up` = unname(beta[["sideL"]] + beta[["armL-up:sideL"]]),
    directionality = unname(beta[["armL-up:sideL"]])
  )
  observed <- c(
    `R-up` = side_summary$estimate[side_summary$arm == "R-up"],
    `L-up` = side_summary$estimate[side_summary$arm == "L-up"],
    directionality = direction_summary$estimate
  )
  if (any(abs(observed - expected) > 1e-7 * (1 + abs(expected)))) {
    stop(stage, " ", engine, " ", offset_type, ": emmeans/factor-coding coefficient check failed.")
  }

  n_obs <- nrow(data)
  n_embryo <- length(unique(data$embryo))
  n_peak <- length(unique(data$peak_id))

  # Per-TERM small-sample / inferential-unit calibration (addresses the two
  # review findings). The embryo INTERACTION slope int_LL makes the embryo the
  # inferential unit for exactly ONE quantity — the side x arm DIRECTIONALITY —
  # and only when it is actually in the model:
  #   - directionality WITH int_LL retained: an embryo-level estimate -> t(
  #     n_embryo - 1) reference for BOTH the CI and the p (the asymptotic z
  #     over-states at 8 embryos). df carried in `df_used`. The run log prints
  #     the exact label-swap p beside it (read from the committed stats TSV).
  #   - directionality WITHOUT int_LL (the fallback, whenever the slope is not
  #     identifiable): the embryo INTERCEPT cancels from the contrast, so
  #     (1|embryo) does NOT make the embryo the unit -> model p promoter-driven
  #     -> SUPPRESS (NA), defer to the label-swap floor.
  #   - per-arm (L-up / R-up) side effects: int_LL is 0 for every R-up
  #     observation and the embryo intercept cancels from each per-arm L-R
  #     contrast, so NO per-arm p is embryo-level under this RE (calibrating them
  #     needs per-arm embryo slopes, unidentifiable at n=8). ALWAYS suppress the
  #     per-arm model p (NA) and defer per-arm inference to the label-swap arms
  #     (hh21_hh1716_module_labelswap_stats.tsv). Estimate + an asymptotic CI kept for
  #     effect-size context; `embryo_calibrated=FALSE` flags them.
  # From the ladder STEP, not by substring-matching `ranef_structure` — that label is a display
  # string the exhaustion path appends to, so a reworded suffix would silently flip this bit.
  int_LL_retained <- identical(fit_result$step_number, 1L)
  make_row <- function(term, estimate, se, alternative) {
    calibrated <- term == "directionality" && int_LL_retained
    df_used <- if (calibrated) n_embryo - 1L else NA_integer_
    crit <- if (calibrated) stats::qt(0.975, df = df_used) else stats::qnorm(0.975)
    if (calibrated) {
      p_two <- 2 * stats::pt(-abs(estimate / se), df = df_used)
      p1 <- one_sided_p(estimate, p_two, alternative)
    } else {
      p1 <- NA_real_   # not embryo-calibrated under this RE -> defer to label-swap
    }
    data.frame(
      stage = stage,
      floor = FLOOR,       # the pseudobulk tier this fit modelled (the contested parameter)
      module_inputs_sha = MODULE_INPUTS_SHA,   # which module-test run produced counts + halfmeta
      engine = engine,
      offset_type = offset_type,
      term = term,
      estimate_log2 = estimate / log(2),
      ci_lo_log2 = (estimate - crit * se) / log(2),
      ci_hi_log2 = (estimate + crit * se) / log(2),
      p_onesided = p1,
      est_natlog = estimate,
      se_natlog = se,
      df_used = df_used,
      embryo_calibrated = calibrated,
      ranef_structure = fit_result$ranef_structure,
      converged = isTRUE(fit_result$converged),
      singular = isTRUE(fit_result$singular),
      n_obs = n_obs,
      n_embryo = n_embryo,
      n_peak = n_peak,
      check.names = FALSE
    )
  }

  rows <- list()
  for (term in c("L-up", "R-up")) {
    i <- match(term, as.character(side_summary$arm))
    rows[[term]] <- make_row(
      term,
      side_summary$estimate[i],
      side_summary$SE[i],
      if (term == "L-up") "greater" else "less"
    )
  }
  rows$directionality <- make_row(
    "directionality",
    direction_summary$estimate,
    direction_summary$SE,
    "greater"
  )
  do.call(rbind, rows)
}

set.seed(1)
result_blocks <- list()
result_index <- 0L
for (stage in c("HH21", "HH17/16")) {
  data <- stage_data[[stage]]
  for (engine in c("glmer.nb", "glmmTMB")) {
    for (offset_type in c("size_factor", "total_fragments")) {
      fit_result <- fit_with_fallback(data, stage, engine, offset_type)
      result_index <- result_index + 1L
      result_blocks[[result_index]] <- extract_effects(
        fit_result, data, stage, engine, offset_type
      )
    }
  }
}

output <- do.call(rbind, result_blocks)
rownames(output) <- NULL
# `sign_unstable` marks a directionality whose two offset variants straddle zero (see the
# offset-stability guard below). It is a REPORTED state, not a failure: consumers must print
# "sign not determined" rather than quoting a direction.
output$sign_unstable <- FALSE
output_columns <- c(
  "stage", "floor", "module_inputs_sha", "engine", "offset_type", "term", "estimate_log2", "ci_lo_log2",
  "ci_hi_log2", "p_onesided", "est_natlog", "se_natlog", "df_used",
  "embryo_calibrated", "ranef_structure", "converged", "singular", "n_obs",
  "n_embryo", "n_peak", "sign_unstable"
)
if (!identical(names(output), output_columns)) {
  stop("Internal output-schema mismatch.")
}
# HH21 directionality must be embryo-calibrated (interaction slope retained) at
# both engines/offsets; otherwise the headline p is not honest.
hh21_dir <- output[output$stage == "HH21" & output$term == "directionality", ]
# nrow check FIRST: all(logical(0)) is TRUE, so an empty subset would pass the
# calibration assert vacuously (4 = 2 engines x 2 offsets).
if (nrow(hh21_dir) != 4L) {
  stop("Expected 4 HH21 directionality rows (2 engines x 2 offsets); got ", nrow(hh21_dir), ".")
}
if (!all(hh21_dir$embryo_calibrated)) {
  # RECORDED, not fatal (2026-08-18). The protection this guard existed for is enforced
  # upstream: when the interaction slope is not retained, make_row() sets p_onesided to NA,
  # so an uncalibrated p can never be printed as the headline. Crashing here additionally
  # discarded the whole table — including the effect sizes and the trough contrast — the
  # first time the HH21 slope went singular (it did on the PITX2-locus-excluded module,
  # where the +5.4 kb element no longer supplies between-embryo spread).
  cat("\n!! HH21 directionality is NOT embryo-calibrated in ",
      sum(!hh21_dir$embryo_calibrated), " of ", nrow(hh21_dir), " engine x offset fits: ",
      paste(unique(hh21_dir$ranef_structure[!hh21_dir$embryo_calibrated]), collapse = " | "),
      "\n   The model p is suppressed (NA) for those rows; the label-swap permutation is the",
      " anchor for the HH21 directionality.\n", sep = "")
}
.n_paired <- function(st) length(unique(halfmeta$embryo[halfmeta$stage == st & halfmeta$paired]))
# One condition per message: "Final output guardrail failed." told the reader nothing about
# which invariant broke.
if (nrow(output) != 24L) {
  stop("Expected 24 output rows (2 stages x 2 engines x 2 offsets x 3 terms); got ", nrow(output), ".")
}
if (any(!output$converged[output$stage == "HH21"])) {
  stop("An HH21 fit did not converge — the positive control must converge.")
}
if (any(output$n_peak != EXPECTED_PEAKS)) {
  stop("A fit used ", paste(unique(output$n_peak), collapse = "/"), " promoter peaks; expected ",
       EXPECTED_PEAKS, " (the response-arm universe in the counts file).")
}
for (st in c("HH21", "HH17/16")) {
  if (!all(output$n_embryo[output$stage == st] == .n_paired(st))) {
    stop(st, ": fitted n_embryo disagrees with the paired-embryo count in halfmeta (",
         .n_paired(st), ").")
  }
}

directionality <- output[output$term == "directionality", ]
for (stage in c("HH21", "HH17/16")) {
  for (engine in c("glmer.nb", "glmmTMB")) {
    z <- directionality[directionality$stage == stage & directionality$engine == engine, ]
    if (nrow(z) != 2L) stop("Missing offset sensitivity result for ", stage, " ", engine)
    offset_difference <- abs(diff(z$estimate_log2))
    same_sign <- length(unique(sign(z$estimate_log2))) == 1L
    if (!same_sign || offset_difference > 0.1) {
      # Offset stability is only an assertion about estimates the run STANDS BEHIND. A pair
      # already recorded as non-converged is reported, not fatal — otherwise a wobbly
      # non-identifiable arm black-holes the whole table, which is the failure mode the
      # ladder-exhaustion handling exists to remove.
      msg <- paste0(stage, " ", engine,
                    ": directionality is not stable across offsets (absolute log2 difference=",
                    signif(offset_difference, 4), ").")
      # SIGN STRADDLING ZERO WITHIN THE TOLERANCE is likewise reported, not fatal (2026-08-22).
      # The assertion this guard protects is "the run stands behind this direction". When the two
      # offset variants differ by less than the tolerance AND both confidence intervals contain
      # zero, the run stands behind no direction at all — the sign is decided by the offset
      # convention because the estimate IS zero. Failing there black-holes the whole table
      # (including the stable HH21 positive control) over an arm nobody may quote a sign for.
      # Encountered when the +/-1 kb promoter rescue (2026-08-22) added 11 responder
      # promoters to the module arms and moved the HH17/16 trough directionality from +0.040
      # to +/-0.008: null before, null after, but now close enough to zero to straddle it.
      # A sign flip with a LARGE gap, or one where a CI excludes zero, is still fatal.
      within_noise <- !same_sign && offset_difference <= 0.1 &&
        all(z$ci_lo_log2 <= 0 & z$ci_hi_log2 >= 0)
      if (all(z$converged) && !within_noise) stop(msg)
      if (within_noise) {
        output$sign_unstable[output$stage == stage & output$engine == engine &
                               output$term == "directionality"] <- TRUE
        cat("  note (estimates straddle zero within tolerance; recorded sign_unstable=TRUE, ",
            "no direction is asserted): ", msg, "\n", sep = "")
      } else {
        cat("  note (non-converged fit, reported not fatal): ", msg, "\n", sep = "")
      }
    }
  }
}

write.table(
  output,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)
if (!file.exists(output_path) || is.na(file.info(output_path)$size) ||
    file.info(output_path)$size <= 0L) {
  stop("Output TSV was not written or is empty: ", output_path)
}

# Committed within-embryo label-swap directionality p, READ from the artifact the panel writes
# (hardcoding it printed July values next to a fresh fit).  Absent file -> NA, never a literal.
labelswap_path <- file.path(repo_root, "figures", "data", "hh21_hh1716_module_labelswap_stats.tsv")
reference <- data.frame(stage = c("HH21", "HH17/16"), label_swap_p = NA_real_,
                        stringsAsFactors = FALSE)
if (file.exists(labelswap_path)) {
  ls_stats <- read_tsv(labelswap_path)
  # The label-swap table is written by the PANEL, which needs this file as input — so on the
  # regeneration path it is always from the previous render. Compare its input fingerprint
  # before quoting its p beside a fresh fit.
  ls_sha <- unique(as.character(ls_stats$module_inputs_sha))
  if (length(ls_sha) != 1L || !identical(ls_sha, MODULE_INPUTS_SHA)) {
    cat("  note: label-swap table is from module-test run(s) ", paste(ls_sha, collapse = ", "),
        " but this fit read ", MODULE_INPUTS_SHA,
        " — reporting NA rather than pairing a fresh estimate with a stale permutation p.\n",
        sep = "")
  } else {
    for (i in seq_len(nrow(reference))) {
      hit <- ls_stats$p[ls_stats$stage == reference$stage[i] & ls_stats$module == "directionality"]
      if (length(hit) == 1L) reference$label_swap_p[i] <- hit
    }
  }
} else {
  cat("  note: ", labelswap_path, " absent — label-swap column reported as NA",
      " (re-render panel_hh21_hh1716_module.py to produce it).\n", sep = "")
}

cat("\n===== PER-STAGE DIRECTIONALITY COMPARISON =====\n")
for (stage in c("HH21", "HH17/16")) {
  d <- directionality[directionality$stage == stage, ]
  ref <- reference[reference$stage == stage, ]
  comparison <- data.frame(
    engine = d$engine,
    offset = d$offset_type,
    GLMM_est_log2 = sprintf("%+.4f", d$estimate_log2),
    GLMM_95pct_CI = sprintf("[%+.4f, %+.4f]", d$ci_lo_log2, d$ci_hi_log2),
    ref_df = ifelse(is.na(d$df_used), "z", paste0("t", d$df_used)),
    GLMM_p_one_sided = ifelse(is.na(d$p_onesided), "NA (deferred)",
                              format.pval(d$p_onesided, digits = 4, eps = 1e-99)),
    label_swap_p = format.pval(rep(ref$label_swap_p, nrow(d)), digits = 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat("\n", stage, if (stage == "HH21") " (primary; firm sides)\n" else " (sensitivity; design-constrained sides)\n", sep = "")
  print(comparison, row.names = FALSE, right = FALSE)
  if (any(is.na(d$p_onesided))) {
    n_emb <- unique(d$n_embryo)
    reason <- if (length(n_emb) == 1L && !is.na(n_emb) && n_emb < MIN_EMBRYO_FOR_SLOPE) {
      paste0("too few embryos for the interaction slope (n = ", n_emb, " < ",
             MIN_EMBRYO_FOR_SLOPE, ") -> skipped a priori")
    } else {
      paste0("no embryo-calibrated slope at n = ", paste(n_emb, collapse = "/"),
             ": ", paste(unique(d$ranef_structure[is.na(d$p_onesided)]), collapse = " | "))
    }
    cat("  note: NA p =", reason,
        "; model directionality inference deferred to the label-swap floor.\n")
  }
}

cat("\n===== PER-ARM SIDE EFFECTS (L - R within arm; effect size only) =====\n")
cat("  per-arm model p is NOT embryo-calibrated under this RE (see make_row);",
    "\n  reported as NA -> per-arm inference defers to the label-swap arms",
    "(see hh21_hh1716_module_labelswap_stats.tsv, written by panel_hh21_hh1716_module.py).\n")
for (stage in c("HH21", "HH17/16")) {
  a <- output[output$stage == stage & output$term %in% c("L-up", "R-up"), ]
  arm_table <- data.frame(
    engine = a$engine,
    offset = a$offset_type,
    arm = a$term,
    estimate_log2 = sprintf("%+.4f", a$estimate_log2),
    CI_95pct = sprintf("[%+.4f, %+.4f]", a$ci_lo_log2, a$ci_hi_log2),
    p_one_sided = ifelse(is.na(a$p_onesided), "NA (label-swap)",
                         format.pval(a$p_onesided, digits = 4, eps = 1e-99)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat("\n", stage, "\n", sep = "")
  print(arm_table, row.names = FALSE, right = FALSE)
}

cat("\n===== OFFSET AND ENGINE AGREEMENT =====\n")
for (stage in c("HH21", "HH17/16")) {
  for (engine in c("glmer.nb", "glmmTMB")) {
    z <- directionality[directionality$stage == stage & directionality$engine == engine, ]
    cat(sprintf(
      "%s | %s | offset absolute difference = %.6f log2\n",
      stage, engine, abs(diff(z$estimate_log2))
    ))
  }
  for (offset_type in c("size_factor", "total_fragments")) {
    z <- directionality[
      directionality$stage == stage & directionality$offset_type == offset_type,
    ]
    z <- z[match(c("glmer.nb", "glmmTMB"), z$engine), ]
    agreement <- round(z$estimate_log2[1], 2) == round(z$estimate_log2[2], 2)
    cat(sprintf(
      "%s | %s | engines agree at 2 decimals = %s (absolute difference = %.6f log2)%s\n",
      stage, offset_type, agreement, abs(diff(z$estimate_log2)),
      if (!agreement) "; use glmer.nb as source of truth" else ""
    ))
  }
}

cat("\nWrote ", nrow(output), " rows to ", output_path, "\n", sep = "")
cat("\n===== SESSION INFO =====\n")
print(sessionInfo())
