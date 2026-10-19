#!/usr/bin/env Rscript
# s8_glmm.R -- E-S1 sensitivity arm S8 (docs/m5_esr1_prolif_stats_plan.md (not distributed) §3).
#
# Depth-adjusted per-cell ESR1 detection: binomial GLMM with a log-sequencing-
# depth offset and an embryo random intercept.  Three models:
#   side       detected ~ side + offset(log n_counts) + (1 | embryo)  <- the plan's spec
#   side_pop   ... + population                                       <- GLMM analogue of the
#                                                                        composition standardisation
#   side_slope ... + (1 + side | embryo)                              <- side effect tested at
#                                                                        the EMBRYO level
# `side` and `side_pop` put the side contrast at the CELL level: with ~16 k cells
# and 8 embryos their Wald p is anticonservative by orders of magnitude (the real
# replication unit is the embryo).  `side_slope` lets the side effect vary by
# embryo, which is the comparable test.  All three are SENSITIVITY arms and
# non-claim-bearing; the primary is the exact paired test on embryo-level
# standardised fractions.
#
# Usage: Rscript s8_glmm.R <in.tsv> <out.tsv>
# Input columns: genotype, inferred_side, population, detected (0/1), n_counts.
# Called by stats_esr1_lr_hh25.py, which supplies exactly the cells the
# standardised estimate is built from.  Fails loud: no fallback estimator.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: s8_glmm.R <in.tsv> <out.tsv>")
suppressPackageStartupMessages(library(glmmTMB))

d <- read.delim(args[1], stringsAsFactors = FALSE)
need <- c("genotype", "inferred_side", "population", "detected", "n_counts")
miss <- setdiff(need, names(d))
if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))
if (any(!d$inferred_side %in% c("L", "R"))) stop("inferred_side outside {L, R}")
if (any(d$n_counts <= 0)) stop("non-positive n_counts -- log-depth offset undefined")

# Right is the reference level, so the side coefficient reads "L vs R".
d$side <- factor(d$inferred_side, levels = c("R", "L"))
d$embryo <- factor(d$genotype)
d$population <- factor(d$population)
d$logdepth <- log(d$n_counts)
cat(sprintf("[s8_glmm] %d cells, %d embryos, %d populations\n",
            nrow(d), nlevels(d$embryo), nlevels(d$population)))

fit_one <- function(form, label) {
  m <- glmmTMB(form, family = binomial, data = d)
  if (!m$fit$convergence %in% c(0L)) stop(label, ": glmmTMB did not converge")
  co <- summary(m)$coefficients$cond
  if (!"sideL" %in% rownames(co)) stop(label, ": no sideL coefficient in the fit")
  est <- co["sideL", "Estimate"]; se <- co["sideL", "Std. Error"]
  data.frame(model = label, estimate = est,
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
             p = co["sideL", "Pr(>|z|)"], n_cells = nrow(d))
}

out <- rbind(
  fit_one(detected ~ side + offset(logdepth) + (1 | embryo), "side"),
  fit_one(detected ~ side + population + offset(logdepth) + (1 | embryo), "side_pop"),
  fit_one(detected ~ side + offset(logdepth) + (1 + side | embryo), "side_slope"))
print(out)
write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
