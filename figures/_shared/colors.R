# Shared color registry for all R figures (ggplot2).
#
# Single source of truth = references/color_registry.tsv.
#
# Usage:
#   source(".../figures/_shared/colors.R")
#   get_color("stage", "HH25")                       # single hex
#   get_palette("stage")                             # named vector, all stages
#   get_palette("celltype", labels = c("LPM mesenchyme", "Sclerotome"))
#   scale_fill_manual(values = get_palette("lateral", labels = c("L", "R")))
#
# Unknown labels stop() with an error -- never silently fall back.

.color_registry_cache <- NULL

.find_repo_root_colors <- function() {
  d <- getwd()
  marker <- c("references", "color_registry.tsv")
  while (d != "/" && !file.exists(do.call(file.path, c(list(d), as.list(marker))))) {
    d <- dirname(d)
  }
  if (d == "/") {
    stop("Could not locate references/color_registry.tsv from ", getwd())
  }
  d
}

.load_color_registry <- function() {
  if (is.null(.color_registry_cache)) {
    path <- file.path(.find_repo_root_colors(), "references", "color_registry.tsv")
    if (!file.exists(path)) {
      stop(sprintf(
        "%s missing -- this is the single source of truth for figure colors.",
        path
      ))
    }
    df <- read.table(
      path, sep = "\t", header = TRUE,
      stringsAsFactors = FALSE, check.names = FALSE,
      quote = "", comment.char = ""
    )
    for (col in c("level", "label", "lineage", "hex")) {
      df[[col]] <- trimws(df[[col]])
    }
    # Sanity: every hex must be #RRGGBB
    bad <- !grepl("^#[0-9A-Fa-f]{6}$", df$hex)
    if (any(bad)) {
      stop(sprintf(
        "color_registry.tsv contains malformed hex values:\n%s",
        paste(capture.output(print(df[bad, , drop = FALSE])), collapse = "\n")
      ))
    }
    dup <- duplicated(df[, c("level", "label")])
    if (any(dup)) {
      stop("color_registry.tsv has duplicate (level, label) rows")
    }
    .color_registry_cache <<- df
  }
  .color_registry_cache
}

get_color <- function(level, label) {
  df <- .load_color_registry()
  match_rows <- df[df$level == level & df$label == label, , drop = FALSE]
  if (nrow(match_rows) == 0) {
    stop(sprintf(
      "color_registry: no entry for level='%s', label='%s'. Add it to references/color_registry.tsv.",
      level, label
    ))
  }
  match_rows$hex[1]
}

get_palette <- function(level, labels = NULL, lineage = NULL) {
  df <- .load_color_registry()
  sub <- df[df$level == level, , drop = FALSE]
  if (!is.null(lineage)) {
    sub <- sub[sub$lineage == lineage, , drop = FALSE]
  }
  registry <- setNames(sub$hex, sub$label)

  if (is.null(labels)) {
    return(registry)
  }

  missing <- setdiff(labels, names(registry))
  if (length(missing) > 0) {
    stop(sprintf(
      "color_registry: labels not found for level='%s'%s: %s. Add to references/color_registry.tsv.",
      level,
      if (!is.null(lineage)) sprintf(" lineage='%s'", lineage) else "",
      paste(missing, collapse = ", ")
    ))
  }
  registry[labels]
}

get_lineage <- function(level, label) {
  df <- .load_color_registry()
  match_rows <- df[df$level == level & df$label == label, , drop = FALSE]
  if (nrow(match_rows) == 0) {
    stop(sprintf(
      "color_registry: no entry for level='%s', label='%s'.",
      level, label
    ))
  }
  match_rows$lineage[1]
}

list_labels <- function(level, lineage = NULL) {
  df <- .load_color_registry()
  sub <- df[df$level == level, , drop = FALSE]
  if (!is.null(lineage)) {
    sub <- sub[sub$lineage == lineage, , drop = FALSE]
  }
  sub$label
}
