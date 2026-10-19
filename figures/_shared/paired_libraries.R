# Which genotypes are drawn as a half-vs-half point (fig02 panel C / Fig 4D).
#
# A genotype seen in two libraries is one embryo's LEFT half in one library and
# its RIGHT half in the other ONLY when the design says so.  The library pair
# alone does not decide it: lane-loading-record tube group 17-B(2) (embryos 17-5/17-6) is ONE left-half tube
# split across GEX04 and GEX0S-duck (a technical duplicate of the same tissue),
# while tube 17-A is a real L/R pair in the SAME two libraries.  So:
#
#   1. The d2l sidecar (figures/data/d2l_lr_inference.tsv -- raw d2l on purpose:
#      fig02 is the pre-fig03 progressive-disclosure surface) marks a same-side
#      duplicate genotype with firmness "same_side_dup" on BOTH library rows,
#      with one shared side.  Such a genotype is EXCLUDED (logged): its two
#      libraries hold the same half.  The side itself comes from the design;
#      d2l's inferred_side / pair_constraint_ok are NOT used to decide pairing.
#   2. Every other two-library genotype is drawn only if its library pair is an
#      opposite-side pair in design/embryo_lr_design.tsv (some embryo has L in
#      one of the two libraries and R in the other, after the tech-rep
#      collapse).  Anything else raises.
#
# Library A / B = lexical order of the two collapsed library names (the
# lower-numbered library is A), identical to the historical hardcoded list.
#
# Requires library_naming.R (collapse_lib) to be sourced first.

D2L_FIRMNESS <- c("direct", "placement", "singleton_rank", "weak", "same_side_dup")
SAME_SIDE_DUP <- "same_side_dup"

.pair_key <- function(libs) paste(sort(libs, method = "radix"), collapse = "|")

# Opposite-side and same-side library pairs from the design (tech reps collapsed).
design_library_pairs <- function(design_path) {
  if (!file.exists(design_path)) stop("missing design table: ", design_path)
  d <- read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in c("embryo", "library", "side")) {
    if (!col %in% names(d)) stop(design_path, " lacks column ", col)
  }
  bad_side <- setdiff(unique(d$side), c("L", "R"))
  if (length(bad_side)) stop("design side values outside {L,R}: ",
                             paste(bad_side, collapse = ", "))
  d$lib_m <- collapse_lib(d$library)
  d <- unique(d[, c("embryo", "lib_m", "side")])
  opposite <- character(0)
  same <- character(0)
  for (emb in unique(d$embryo)) {
    e <- d[d$embryo == emb, ]
    if (any(duplicated(e$lib_m)))
      stop("design embryo ", emb, " has both sides in one library_m: ",
           paste(e$lib_m[duplicated(e$lib_m)], collapse = ", "))
    if (nrow(e) == 1) next
    if (nrow(e) != 2)
      stop("design embryo ", emb, " spans ", nrow(e), " library_m (",
           paste(e$lib_m, collapse = ", "), "); only two-library embryos are modelled")
    key <- .pair_key(e$lib_m)
    if (e$side[1] != e$side[2]) opposite <- c(opposite, key) else same <- c(same, key)
  }
  # same_n: design same-side embryos per library pair (cross-checked vs d2l)
  list(opposite = sort(unique(opposite)), same = sort(unique(same)),
       same_n = table(same))
}

# Same-side duplicate genotypes from the raw d2l sidecar.
d2l_same_side_duplicates <- function(d2l_path) {
  if (!file.exists(d2l_path)) stop("missing d2l sidecar: ", d2l_path)
  lr <- read.delim(d2l_path, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in c("genotype", "library_m", "firmness", "inferred_side")) {
    if (!col %in% names(lr)) stop(d2l_path, " lacks column ", col)
  }
  unknown <- setdiff(unique(lr$firmness), D2L_FIRMNESS)
  if (length(unknown)) stop("d2l firmness value(s) outside the known alphabet {",
                            paste(D2L_FIRMNESS, collapse = ", "), "}: ",
                            paste(unknown, collapse = ", "))
  dup_genos <- unique(lr$genotype[lr$firmness == SAME_SIDE_DUP])
  out <- lapply(dup_genos, function(g) {
    r <- lr[lr$genotype == g, ]
    if (!all(r$firmness == SAME_SIDE_DUP) || nrow(r) != 2 ||
        length(unique(r$library_m)) != 2 || length(unique(r$inferred_side)) != 1)
      stop("d2l same-side duplicate ", g, " is malformed (need exactly two ",
           SAME_SIDE_DUP, " rows, two libraries, one side): ",
           paste(sprintf("%s/%s/%s", r$library_m, r$firmness, r$inferred_side),
                 collapse = "; "))
    data.frame(genotype = g, pair_key = .pair_key(r$library_m),
               side = r$inferred_side[1], stringsAsFactors = FALSE)
  })
  if (length(out) == 0)
    return(data.frame(genotype = character(0), pair_key = character(0),
                      side = character(0), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

# The d2l duplicates must be exactly the design's same-side embryos: one
# duplicate genotype per same-side design embryo, in the same library pair.
check_duplicates_against_design <- function(dups, design_pairs) {
  have <- table(dups$pair_key)
  want <- design_pairs$same_n
  keys <- union(names(have), names(want))
  for (k in keys) {
    h <- if (k %in% names(have)) as.integer(have[[k]]) else 0L
    w <- if (k %in% names(want)) as.integer(want[[k]]) else 0L
    if (h != w)
      stop(sprintf("same-side duplicates in %s: d2l sidecar has %d genotype(s), design has %d embryo(s)",
                   k, h, w))
  }
  invisible(TRUE)
}

# genotype -> (lib_A, lib_B) for half-vs-half genotypes.
# geno_libs: data.frame(genotype, libs = list-column of collapsed library_m).
# stage_label: only for messages.
resolve_half_pairs <- function(geno_libs, design_pairs, dups, stage_label) {
  rows <- list()
  for (i in seq_len(nrow(geno_libs))) {
    g <- geno_libs$genotype[i]
    libs <- sort(unique(geno_libs$libs[[i]]), method = "radix")
    if (length(libs) == 1) next
    if (length(libs) > 2)
      stop(sprintf("%s: genotype %s spans %d libraries (%s); a half-vs-half pair needs exactly two",
                   stage_label, g, length(libs), paste(libs, collapse = ", ")))
    key <- .pair_key(libs)
    if (g %in% dups$genotype) {
      dk <- dups$pair_key[dups$genotype == g]
      if (dk != key)
        stop(sprintf("%s: same-side duplicate %s is in %s here but %s in the d2l sidecar",
                     stage_label, g, key, dk))
      message(sprintf("    excluded (same-side duplicate, one %s half in both %s): %s",
                      dups$side[dups$genotype == g], key, g))
      next
    }
    if (!key %in% design_pairs$opposite)
      stop(sprintf("%s: genotype %s spans %s, which no design embryo holds as an L/R pair%s",
                   stage_label, g, key,
                   if (key %in% design_pairs$same) " (design: same-side libraries only)" else ""))
    rows[[length(rows) + 1]] <- data.frame(genotype = g, lib_A = libs[1], lib_B = libs[2],
                                           stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}
