# Shared loader for d2g stage-restricted subcluster annotations (R side).
#
# Single source of truth = references/d2g_subcluster_annotations.tsv
# (owned by pipeline/d2g_validate_subclusters.py). R figure scripts use
# this helper instead of hardcoding cluster->population maps.
#
# Leiden cluster IDs are not seed-stable. Do not hardcode them.

.subcluster_map_cache <- NULL

.find_repo_root <- function() {
  d <- getwd()
  marker <- c("references", "d2g_subcluster_annotations.tsv")
  while (d != "/" && !file.exists(do.call(file.path, c(list(d), as.list(marker))))) {
    d <- dirname(d)
  }
  if (d == "/") stop("Could not locate references/d2g_subcluster_annotations.tsv from ", getwd())
  d
}

.load_subcluster_annotations <- function() {
  if (is.null(.subcluster_map_cache)) {
    path <- file.path(.find_repo_root(), "references", "d2g_subcluster_annotations.tsv")
    if (!file.exists(path)) {
      stop(sprintf(
        "%s missing -- run pipeline/d2g_validate_subclusters.py after any d2g rerun.",
        path
      ))
    }
    .subcluster_map_cache <<- read.table(
      path, sep = "\t", header = TRUE,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  .subcluster_map_cache
}

get_leiden_col <- function(stage) {
  ann <- .load_subcluster_annotations()
  cols <- unique(ann$leiden_col[ann$stage == stage])
  if (length(cols) != 1) {
    stop(sprintf("stage=%s resolved to %d leiden cols: %s",
                 stage, length(cols), paste(cols, collapse = ",")))
  }
  cols[1]
}

get_domain_map <- function(stage, laterality_only = FALSE) {
  ann <- .load_subcluster_annotations()
  sub <- ann[ann$stage == stage, , drop = FALSE]
  if (nrow(sub) == 0) stop(sprintf("stage=%s not in annotations table", stage))
  if (laterality_only) {
    # Derive from POOLS so a future inference/bilateral edit flows here automatically.
    keep <- c(.POOLS[[stage]]$inference, .POOLS[[stage]]$bilateral)
    sub <- sub[sub$population %in% keep, , drop = FALSE]
  }
  setNames(sub$population, as.character(sub$cluster))
}

get_population_clusters <- function(stage, population) {
  ann <- .load_subcluster_annotations()
  sub <- ann[ann$stage == stage & ann$population == population, , drop = FALSE]
  sort(as.integer(sub$cluster))
}

# ── Laterality pool classification (mirror of subcluster_map.py POOLS) ──
# Per-stage pool membership for L/R inference panels. Keep in sync with
# figures/_shared/subcluster_map.py.
#
# 2026-05-04 update (docs/lpm_mesenchyme_subcluster_biology.md):
#   - "Splanchnic" → "Venous-pole spl." rename + move from inference → excluded
#     (off-axis dissection-edge introgression at the cardiac venous pole)
#   - "Epicardial mes." → "Mesothelium" rename (still excluded)
# 2026-06-18 update (annotation vet):
#   - "Somatic" → "Skeletogenic LPM" rename (FOXF1-low, HOX-negative,
#     TBX1+/EYA1+ pharyngeal-cranial skeletal mesenchyme; L/R-symmetric; still
#     excluded).  + "Contamination" added to excluded (ribosomal/erythroid
#     technical junk reclassified out by the compendium contamination_check).
.POOLS <- list(
  "HH13" = list(
    inference = c("Early LPM", "Cardiopharyngeal LPM"),
    bilateral = character(0),
    # 2026-07-02 (annotation vet): HH13 skeletogenic label softened to "PRRX1+/FOXF1-low
    # mes." (identity not earned pre-HH17) — STILL EXCLUDED (FOXF1-low). HH17/16-
    # HH25 keep "Skeletogenic LPM". Pool membership unchanged; only the name differs.
    excluded  = c("PRRX1+/FOXF1-low mes.", "Sclerotome", "Dermomyotome",
                  "Mesothelium", "Venous-pole spl.", "Contamination",
                  "ambiguous")
  ),
  "HH17/16" = list(
    # 2026-08-09 (seeded d2g): the reproducible res-0.3 partition
    # emits NO "Early LPM" at HH17/16 (old cl1 merged into Skeletogenic /
    # Venous-pole). Mirrors subcluster_map.py POOLS.
    inference = c("Cardiopharyngeal LPM", "BARX1+"),
    bilateral = c("Dorsal-biased LPM"),
    # 2026-07-02 vet: cl4 = PRRX1+/FOXF1-low fibroblast (branchiomeric gate);
    # cl2/cl9 remain genuine Skeletogenic LPM. Both excluded.
    excluded  = c("Skeletogenic LPM", "PRRX1+/FOXF1-low mes.",
                  "Sclerotome", "Dermomyotome",
                  "Mesothelium", "Venous-pole spl.", "Contamination",
                  "ambiguous")
  ),
  "HH21" = list(
    inference = c("Ventral spl.", "Syringeal", "BARX1+"),
    bilateral = c("Dorsal sub-epi."),
    excluded  = c("Skeletogenic LPM", "Sclerotome", "Dermomyotome",
                  "Mesothelium", "Venous-pole spl.", "Contamination",
                  "ambiguous")
  ),
  "HH25" = list(
    # 2026-07-02 vet: HH25 syrinx precursor (IRX1+) is cl3, relabeled
    # "Chondrogenic LPM" -> "Syringeal"; old "Syringeal" cl5 (no IRX1) -> Ventral
    # spl. "Chondrogenic LPM" no longer emitted at HH25 -> dropped from inference.
    inference = c("Ventral spl.", "Syringeal"),
    bilateral = c("Dorsal sub-epi."),
    excluded  = c("Skeletogenic LPM", "Sclerotome", "Dermomyotome",
                  "Mesothelium", "Venous-pole spl.", "Contamination",
                  "ambiguous")
  )
)

# %||% is Python's `or` for NULL fallback; defined here since not all
# R versions provide it (rlang exposes it but we don't take the dep).
# Defined before first use in get_pool_map.
`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Publication display names (mirror of subcluster_map.py) ──────────────────
# DISPLAY-ONLY, purely additive: the SSOT population string never changes, and
# no code may filter, join, or compare on a value from here. "Syringeal" names
# a marker- and dissection-INFERRED identity, not a spatially validated one — the
# only molecular in-situ at this locus is PITX2 HCR. The chondrogenic program
# (WNT5A/RUNX2/DPP10, PTCHD4) is HH25-ONLY; the HH21 cluster with the same label is
# IRX1-defined with no cartilage gene, so the name is earned at HH25 and carried back
# to HH21 by lineage continuity. Rationale:
# docs/lpm_mesenchyme_subcluster_biology.md "Naming honesty".
.POPULATION_DISPLAY_NAMES <- list(
  # SSOT label = c(long form: legends & prose, short form: panel/axis text)
  "Syringeal" = c("Syringeal chondrocyte precursor", "Syringeal chondro.")
)

# A typo'd key would silently never fire — pin it to the .POOLS vocabulary at
# source() time (.POOLS is a literal, so this costs no file I/O).
local({
  pool_labels <- unique(unlist(.POOLS, use.names = FALSE))
  bad <- setdiff(names(.POPULATION_DISPLAY_NAMES), pool_labels)
  if (length(bad) > 0) {
    stop(sprintf(".POPULATION_DISPLAY_NAMES keys absent from .POOLS: %s",
                 paste(bad, collapse = ", ")))
  }
  # Validate the VALUE shape too: a scalar passes a key-only check and then makes
  # display_population(short = TRUE) return NA_character_, which renders as a facet
  # literally labelled "NA".  Mirrors the Python assert.
  ok <- vapply(.POPULATION_DISPLAY_NAMES,
               function(v) is.character(v) && length(v) == 2L && all(nzchar(v)),
               logical(1))
  if (!all(ok)) {
    stop(sprintf(".POPULATION_DISPLAY_NAMES values must be c(long, short) non-empty character; malformed: %s",
                 paste(names(.POPULATION_DISPLAY_NAMES)[!ok], collapse = ", ")))
  }
})

#' Reader-facing label for an SSOT population name (vectorised).
#'
#' Call ONLY where a label is rendered for a human (facet strips, axis text,
#' legends). NEVER on a value used to filter, join, group, or compare AGAINST
#' SSOT DATA. One deliberate exception: when the value being matched is itself a
#' rendered label (e.g. the `group` column of figures/data/d2g_dotplot_spatial*.tsv),
#' the key must be derived here rather than hardcoded.
#'
#' `unhedged` is the caller's own label for populations carrying no hedge (a
#' figure often shortens names to fit its axis); pass it so a hedged name
#' overrides the figure-local shortening while every other entry keeps it.
#' Recycled against `population`. It is NOT an error fallback — most
#' populations legitimately have no hedge; the source-time check above catches
#' the real failure mode, a key matching no population.
display_population <- function(population, short = FALSE, unhedged = NULL) {
  idx <- if (short) 2L else 1L
  population <- as.character(population)
  if (!is.null(unhedged) && !length(unhedged) %in% c(1L, length(population))) {
    stop(sprintf(paste0("display_population(): `unhedged` has length %d, expected 1 or %d ",
                        "— silent recycling would mislabel rows"),
                 length(unhedged), length(population)))
  }
  unhedged <- if (is.null(unhedged)) population else
    rep_len(as.character(unhedged), length(population))
  vapply(seq_along(population), function(i) {
    hedged <- .POPULATION_DISPLAY_NAMES[[population[i]]]
    if (is.null(hedged)) unhedged[i] else hedged[idx]
  }, character(1), USE.NAMES = FALSE)
}

# Mirror of figures/_shared/subcluster_map.py EXCLUDED_OPTIONAL_PER_STAGE.
# Per-stage allow-list of excluded labels that may be absent without raising
# (biological no-shows: Mesothelium has 0 clusters at HH13; Venous-pole spl.
# absent at HH25).
.EXCLUDED_OPTIONAL_PER_STAGE <- list(
  "HH13"    = c("Mesothelium"),
  "HH17/16" = character(0),
  "HH21"    = character(0),
  # 2026-07-02 vet: Mesothelium may be absent at HH25 (cl8->Ventral, cl10->Venous-
  # pole after the respiratory-gate fix); Venous-pole present (cl10).
  "HH25"    = c("Venous-pole spl.", "Mesothelium")
)

# Mirror of subcluster_map.py _ALWAYS_OPTIONAL — paraxial (filtered upstream
# by pipeline/_d2g_lpm_selection.py) + the residual ambiguous category +
# technical-junk Contamination (present only at early stages).
.ALWAYS_OPTIONAL <- c("Sclerotome", "Dermomyotome", "ambiguous", "Contamination")

get_pool_map <- function(stage, kind = "inference") {
  if (!stage %in% names(.POOLS)) {
    stop(sprintf("stage=%s not in POOLS", stage))
  }
  if (!kind %in% c("inference", "bilateral", "excluded")) {
    stop(sprintf("kind=%s not in {inference, bilateral, excluded}", kind))
  }
  wanted <- .POOLS[[stage]][[kind]]
  ann <- .load_subcluster_annotations()
  sub <- ann[ann$stage == stage & ann$population %in% wanted, , drop = FALSE]

  # Validate every wanted label is present (catches typos + post-rename drift)
  # for all three pool kinds.  excluded kind permits biologically-absent
  # labels per .EXCLUDED_OPTIONAL_PER_STAGE | .ALWAYS_OPTIONAL.
  if (length(wanted) > 0) {
    present <- unique(ann$population[ann$stage == stage])
    missing <- setdiff(wanted, present)
    if (kind == "excluded") {
      optional <- union(
        .EXCLUDED_OPTIONAL_PER_STAGE[[stage]] %||% character(0),
        .ALWAYS_OPTIONAL
      )
      missing <- setdiff(missing, optional)
    }
    if (length(missing) > 0) {
      stop(sprintf(
        "POOLS[%s][%s] expects {%s} but d2g_subcluster_annotations.tsv has none -- regenerate via pipeline/d2g_validate_subclusters.py, update POOLS, or add to .EXCLUDED_OPTIONAL_PER_STAGE if biologically-absent.",
        stage, kind, paste(missing, collapse = ", ")
      ))
    }
  }
  setNames(sub$population, as.character(sub$cluster))
}

get_pool_clusters <- function(stage, kind = "inference") {
  sort(as.integer(names(get_pool_map(stage, kind))))
}
