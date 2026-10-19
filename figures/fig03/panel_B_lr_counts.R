#!/usr/bin/env Rscript
# fig03 panel B — L-vs-R counts, DIVERGENT, three rows × two compartments.
#   Rows  : snRNA (# L/R DE GENES)  /  ATAC focal 40 kb (the panel-C window)  /
#           ATAC ±500 kb FAM241A-extended (# L/R-biased PEAKS) — a nested zoom:
#           transcriptome -> focal locus -> wider locus.
#   Cols  : Inference pool (ventral laterality pool)  /  Dorsal (bilateral control)
#
# The segway: with the sides fixed by panel A, the snRNA is nearly silent at HH17/16
# (inference pool: 2 DEG, both Right-biased — AHR, CHRM2 — on the corrected-design
# cohort of 13 paired embryos; 0 on the pre-correction 15) — yet in the SAME multiome
# cells the PITX2-locus chromatin is already L-biased (L peaks ≫ R), and stays
# L-biased through HH25, in the ventral
# pool but NOT in the dorsal control (0 L / 0 R locus peaks at every stage; the dorsal
# snRNA arm is 0 at HH17/16 and small-and-Right at HH21/HH25 — see fig03/LEGEND.md (B)).  The spatial
# substructure of those peaks — focal-L (PITX2c core + PANCR) vs the late upstream-R
# desert near FAM241A/Playrr — is panel C.  (Genome-wide ATAC L/R is depth-confounded/
# null — h4 perm emp_p 0.385 (pre-fix tail-biased estimator; see pipeline/_lr_binom.py —
# corrected values in figures/data/h4_summary_*.tsv: HH17/16 count-level emp_p
# L=0.91 / R=0.61, still null; HH25 exceeds the null both directions, read as
# pseudobulk composition structure) — so the honest ATAC signal is
# FOCAL; described in text.  Per-subcluster DEG overlap is figS03 UpSet.)
#
# L-biased UP (blue), R-biased DOWN (red).
#
# Data:
#   RNA per-gene paired DESeq2 (~genotype+side, padj<0.05; +log2FC = Left):
#     figures/data/h2_deg_POOL_*.tsv      (inference pool)
#     figures/data/h2_dorsal_POOL_*.tsv   (dorsal bilateral control)
#   ATAC per-peak L/R DESeq2 at two scales (locus-FDR<0.1):
#     figures/data/g4f_lr_long_floor20.tsv  (compartment ventral|dorsal, design=unpaired)
#     figures/data/pitx2_locus_tad_peaks.tsv (peak coords; focal = the panel-C 40 kb
#       window, wider = full 359-peak FAM241A-extended ±500 kb set — NOT a measured TAD)
#
#   conda run -n duckfig Rscript figures/fig03/panel_B_lr_counts.R
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

F    <- "figures/data"
HERE <- "figures/fig03"
STAGES    <- c("HH13", "HH17/16", "HH21", "HH25")
PADJ_RNA  <- 0.05
PADJ_ATAC <- 0.10
SCOPE_ORD <- c("Inference pool", "Dorsal (ctrl)")
# (wrapped to 3 lines: the row strips are horizontal at print size)
MOD_ORD   <- c("snRNA\nL vs R DEG",
               "ATAC · PITX2\nfocal (40 kb,\npanel C)",
               "ATAC · PITX2\n±500 kb\n(FAM241A-ext.)")
side_cols <- local({ source("figures/_shared/colors.R"); get_palette("lateral", labels = c("L", "R")) })

stage_safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)

# ── RNA: signed L/R DEG counts (inference pool + dorsal control) ─────────────
rna_count <- function(path) {
  d <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  s <- d[!is.na(d$padj) & d$padj < PADJ_RNA, ]
  list(L = sum(s$log2FoldChange > 0), R = sum(s$log2FoldChange < 0))
}
RNA_SRC <- tibble::tribble(
  ~scope,           ~prefix,
  "Inference pool", "h2_deg",       # POOL
  "Dorsal (ctrl)",  "h2_dorsal"     # POOL on dorsal pseudobulks
)
# A stage with no DEG table is NOT a measured zero, and rendering it as an empty bar makes those
# two cases identical to the eye. `if (!file.exists(p)) next` did exactly that: on the pre-relabel
# pseudobulks the HH25 dorsal RNA control was un-estimable, the column rendered blank, and both
# LEGEND.md and PROVENANCE.md read it back as "the Dorsal control is flat ... every stage". HH25 is
# where the ventral pool ramps to 185 L / 290 R, so the dorsal contrast at that stage is the
# panel's punchline and it was being published as a flat control while unmeasured.
#
# HH13 is the one legitimate absence — there is no dorsal compartment at HH13 by construction — so
# it is declared here rather than inferred from a missing file. Everything else must exist.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
# HH13: no dorsal compartment exists at that stage, by construction.
# HH25 WAS un-estimable on the pre-relabel pseudobulks (h2 skipped the stage: of 8 HH25 dorsal
#   genotypes, seven had zero left-side Dorsal sub-epi. cells above the 20-cell floor). On the
#   2026-08 seeded relabel it has SEVEN paired genotypes and h2_dorsal_POOL_HH25.tsv exists, so
#   the cell is rendered from data and this entry is dead. It is kept only as the record of why
#   a declared-absence list exists; if a future rerun loses the stage again, the declaration is
#   already here and the panel will not silently draw a blank bar as a flat control.
EXPECTED_ABSENT <- list("Dorsal (ctrl)" = "HH13")
rna_rows <- list()
for (i in seq_len(nrow(RNA_SRC))) {
  for (stage in STAGES) {
    p <- file.path(F, sprintf("%s_POOL_%s.tsv", RNA_SRC$prefix[i], stage_safe(stage)))
    sc <- RNA_SRC$scope[i]
    if (!file.exists(p)) {
      if (stage %in% (EXPECTED_ABSENT[[sc]] %||% character(0))) next
      stop(sprintf(paste0("panel B: %s is MISSING for scope '%s' at %s.\n",
                          "  This is an UNMEASURED cell, not a measured zero — rendering it as a\n",
                          "  blank bar would publish an absence as a flat control. Either run h2\n",
                          "  for that scope/stage, or add it to EXPECTED_ABSENT with a reason."),
                   basename(p), sc, stage))
    }
    cnt <- rna_count(p)
    rna_rows[[length(rna_rows) + 1]] <- tibble(
      modality = MOD_ORD[1], scope = sc, stage = stage, L = cnt$L, R = cnt$R)
  }
}
rna <- bind_rows(rna_rows)

# ── ATAC: L/R sig-peak counts per compartment, at TWO nested scales ──────────
#   focal = the panel-C 40 kb display window (PITX2 gene + a/b/c proms + +5.4 kb);
#   wider = the full 359-peak FAM241A-extended ±500 kb set (NOT a measured TAD).
atac_long <- suppressMessages(read_tsv(file.path(F, "g4f_lr_long_floor20.tsv"), show_col_types = FALSE))
if ("design" %in% names(atac_long)) atac_long <- atac_long %>% filter(design == "unpaired")
locus <- suppressMessages(read_tsv(file.path(F, "pitx2_locus_tad_peaks.tsv"), show_col_types = FALSE)) %>%
  select(peak_id, start, end)
atac_long <- atac_long %>% inner_join(locus, by = "peak_id")
COMP_TO_SCOPE <- c(ventral = "Inference pool", dorsal = "Dorsal (ctrl)")
source(file.path(HERE, "_focal_window.R"))   # PITX2_FOCAL_WIN — shared with panel C (single source of truth)
FOCAL <- PITX2_FOCAL_WIN   # the panel-C display window

atac_count <- function(df, mod_label) {
  tested <- df %>% distinct(compartment, stage)
  cnt <- df %>% filter(padj_locus < PADJ_ATAC) %>%
    group_by(compartment, stage) %>%
    summarise(L = sum(log2FoldChange > 0), R = sum(log2FoldChange < 0), .groups = "drop")
  tested %>% left_join(cnt, by = c("compartment", "stage")) %>%
    mutate(L = coalesce(L, 0L), R = coalesce(R, 0L),
           modality = mod_label, scope = COMP_TO_SCOPE[compartment]) %>%
    filter(!is.na(scope)) %>% select(modality, scope, stage, L, R)
}
atac <- bind_rows(
  # focal = peaks FULLY CONTAINED in the panel-C window (the predicate panel C
  # uses for its ticks/mask) so this count == the peaks panel C displays.
  atac_count(atac_long %>% filter(start >= FOCAL[1], end <= FOCAL[2]), MOD_ORD[2]),
  atac_count(atac_long, MOD_ORD[3])
)

# ── Combine + long ──────────────────────────────────────────────────────────
wide <- bind_rows(rna, atac) %>%
  filter(stage %in% STAGES)
# Every TESTED cell prints both side counts, a 0 included (grey), so an empty side
# reads as measured.  An untested cell (HH13 dorsal: no compartment, absent from
# `wide`) stays blank.
long <- wide %>%
  pivot_longer(c(L, R), names_to = "side", values_to = "n") %>%
  mutate(signed = ifelse(side == "L", n, -n),
         stage  = factor(stage, levels = STAGES),
         scope  = factor(scope, levels = SCOPE_ORD),
         modality = factor(modality, levels = MOD_ORD),
         side   = factor(side, levels = c("L", "R")))

p <- ggplot(long, aes(stage, signed, fill = side)) +
  geom_col(position = "identity", width = 0.72, colour = "grey15", linewidth = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_text(data = subset(long, n > 0),
            aes(label = n, vjust = ifelse(side == "L", -0.3, 1.2)),
            size = 3.7, colour = "grey20") +
  geom_text(data = subset(long, n == 0),
            aes(label = n, vjust = ifelse(side == "L", -0.3, 1.2)),
            size = 3.7, colour = "grey45") +
  facet_grid(modality ~ scope, scales = "free_y", switch = "y") +
  scale_fill_manual(values = side_cols, name = NULL,
                    labels = c(L = "Left-biased", R = "Right-biased")) +
  scale_y_continuous(labels = function(x) abs(x),
                     expand = expansion(mult = c(0.3, 0.3))) +
  labs(x = NULL, y = NULL) +
  # PNAS final font pass (2026-10): 300 dpi (compositor assumes 300) on A's
  # canvas width, so Fig 5 prints this at x0.653; text >= 7 pt printed.
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10.8, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10.8),
        strip.text.x = element_text(size = 11.5, face = "bold"),
        strip.text.y.left = element_text(size = 10.8, face = "bold", angle = 0,
                                         hjust = 1, lineheight = 0.88),
        strip.placement = "outside",
        legend.position = "top", legend.key.size = unit(0.36, "cm"),
        legend.text = element_text(size = 11),
        legend.margin = margin(0, 0, 0, 0),
        panel.spacing.y = unit(1.1, "lines"),
        plot.margin = margin(2, 6, 2, 2))

# (The hardcoded "silent gap" call-out over HH17/16 is gone: on the corrected-design
# cohort that cell holds 2 Right-biased DEGs, and the count label above says so.)

out <- file.path(HERE, "panel_B_lr_counts.png")
ggsave(out, p, width = 5.067, height = 3.9, dpi = 300)
cat("wrote", out, "\n")
