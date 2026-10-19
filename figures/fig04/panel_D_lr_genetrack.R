#!/usr/bin/env Rscript
# Panel D (ggplot2 re-render) — PITX2 Pitx2-TAD locus per-peak L/R log2FC,
# position-accurate, with THREE aligned tracks:
#   (1) gene models (squished; key genes bold, minor genes grey/small),
#   (2) promoter + notable-peak marks lane (a/b/c promoters; the +5.4 kb (putative)
#       distal L peak; the FAM241A-border R peak) — ggrepel + ggforce::geom_mark_rect,
#   (3) the per-peak L/R log2FC heatmap (4 testable compartment×stage columns).
#
# L/R convention (project standard): Left = blue (+log2FC), Right = red (−log2FC).
# Columns = the 4 L/R-testable contrasts (census lr_testable=1); HH13-ventral
# (1L/2R) and HH25-dorsal (1L/1R) are NOT testable → omitted (same grouping/pool
# as fig03: g4f "ventral"=inference POOL, "dorsal"=bilateral control).
#
# Gene-name conservatism: in the zju1 annotation the Kurpios-lab orthologs are
# LOC-prefixed; we label them LOCxxxxx\n(Playrr?) / (PANCR?) — candidate, not
# asserted. FAM241A uses its CORRECTED coding extent (NOT the 268 kb Gnomon
# over-extension XM_038178404.1).
#
#   conda run -n duckfig Rscript figures/fig04/panel_D_lr_genetrack.R
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(ggrepel); library(ggforce)
})

# fig04 panel A's engine (the filename is historical; it also renders other
# gene-track variants via G4F_OUTDIR).
HERE <- "figures/fig04"
F    <- "figures/data"
PADJ <- 0.1
MB   <- 1e6

# Output overrides (default = HERE = figures/fig04 + computed name).
# fig04/generate.py renders the floor-20 7-track variant straight into
# figures/fig04/panel_A_tad_genetrack.png via these (the TAD gene-track moved
# from fig03 to fig04), so the composite stays self-contained without forking
# this (canonical) renderer.
G4F_OUTDIR  <- Sys.getenv("G4F_OUTDIR", HERE)
G4F_OUTNAME <- Sys.getenv("G4F_OUTNAME", "")      # "" = computed from SUF/DESIGN/FILL

# ── data ────────────────────────────────────────────────────────────────────
# G4F_SUFFIX selects the floor variant: "" = canonical floor-50; "_floor20" =
# floor-20 sensitivity.  G4F_DESIGN selects the laterality replicate design over
# EMBRYO-HALVES (the rectified unit — one pooled pseudobulk per genotype×side, NOT
# per-subcluster): "unpaired" (~inferred_side, ALL halves; primary) or "paired"
# (~unified_genotype_id + inferred_side, both-sided genotypes only; confirmatory).
# Columns are DERIVED from the census testable flag (never hardcoded): unpaired
# uses `lr_testable` (>=2 L AND >=2 R halves); paired uses `paired_testable`
# (>=2 both-sided genotypes).  The axis N is EMBRYO-HALVES, never pseudobulks.
SUF    <- Sys.getenv("G4F_SUFFIX", "")
DESIGN <- Sys.getenv("G4F_DESIGN", "unpaired")    # "unpaired" (primary) | "paired"
stopifnot(DESIGN %in% c("unpaired", "paired"))
# G4F_FILL selects the per-peak EFFECT-SIZE encoded by tile colour:
#   "log2fc" = DESeq2 L/R log2FoldChange (fold-change; scale-free but noisy/explosive
#              for sparse-R peaks, over-weights tiny peaks).
#   "dfrip"  = ΔFRiP = mean(FRiP_L) − mean(FRiP_R) over embryo-halves, shown as Δ mean CPM
#              (= ΔFRiP × 1e6).  A depth-normalised DIFFERENCE: bounded + self-gating for
#              R-end sparsity (few insertions → both FRiP tiny → ΔFRiP≈0), so the locus-wide
#              picture is not dominated by sparse-R noise.  Magnitude-weighted (suppresses
#              genuine fold-changes at low-accessibility peaks) → COMPLEMENT to log2fc, not a
#              replacement.  Same DESeq2 padj_locus asterisks; same grey-for-absent.
#              (precedent: chick-snatac-archr cross_species_halfvshalf.py "diff" metric.)
FILL   <- Sys.getenv("G4F_FILL", "log2fc")        # "log2fc" | "dfrip"
stopifnot(FILL %in% c("log2fc", "dfrip"))
peaks <- read_tsv(file.path(F, "pitx2_locus_tad_peaks.tsv"), show_col_types = FALSE) %>%
  mutate(mid = (start + end) / 2)
lr  <- read_tsv(file.path(F, paste0("g4f_lr_long", SUF, ".tsv")), show_col_types = FALSE)
cen <- read_tsv(file.path(F, paste0("g4f_census", SUF, ".tsv")), show_col_types = FALSE)
if ("design" %in% names(lr)) lr <- lr %>% filter(design == DESIGN)

STAGE_ORD <- c("HH13", "HH17/16", "HH21", "HH25")
COMP_ORD  <- c("ventral", "dorsal")
TESTABLE_COL <- if (DESIGN == "paired") "paired_testable" else "lr_testable"
cols <- cen %>%
  filter(.data[[TESTABLE_COL]] == 1) %>%
  mutate(compartment = factor(compartment, levels = COMP_ORD),
         stage       = factor(stage,       levels = STAGE_ORD)) %>%
  arrange(compartment, stage) %>%                 # ventral stages, then dorsal
  mutate(compartment = as.character(compartment), stage = as.character(stage),
         # axis N = EMBRYO-HALVES behind each column (unpaired: nL/nR halves;
         # paired: # within-embryo pairs) — the honest replication, not pseudobulks
         xlab_disp = if (DESIGN == "paired")
             sprintf("%s\n%s\n%d pairs", compartment, stage, n_paired)
           else
             sprintf("%s\n%s\n%dL/%dR", compartment, stage, n_L_half, n_R_half),
         x_pos = row_number())                     # continuous x so labels sit left of col 1
# Full (column × ALL locus peaks) grid, then LEFT-join the tested L/R log2FC.
# Peaks absent (pooled count < floor) in a column stay NA -> rendered GREY
# (na.value), visually distinct from white = tested-but-bilateral (log2FC≈0).
hm <- tidyr::crossing(cols %>% select(compartment, stage, x_pos, xlab_disp),
                      peaks %>% select(peak_id, start, end, mid, region)) %>%
  left_join(lr %>% select(peak_id, compartment, stage, log2FoldChange, padj_locus),
            by = c("peak_id", "compartment", "stage")) %>%
  mutate(sig    = !is.na(padj_locus) & padj_locus < PADJ,
         tested = !is.na(log2FoldChange),
         y_mb   = mid / MB)

# ── tile fill: log2FC (DESeq2) or ΔFRiP (= Δ mean CPM over embryo-halves) ──
if (FILL == "dfrip") {
  traj <- read_tsv(file.path(F, paste0("g4f_cpm_trajectory_long", SUF, ".tsv")),
                   show_col_types = FALSE)
  dcpm <- traj %>% filter(side %in% c("L", "R")) %>%
    select(peak_id, compartment, stage, side, mean_cpm) %>%
    tidyr::pivot_wider(names_from = side, values_from = mean_cpm) %>%
    mutate(dfrip_cpm = L - R)                 # Δ mean CPM (L−R) = ΔFRiP × 1e6
  hm <- hm %>%
    left_join(dcpm %>% select(peak_id, compartment, stage, dfrip_cpm),
              by = c("peak_id", "compartment", "stage")) %>%
    mutate(fillval = ifelse(tested, dfrip_cpm, NA_real_))   # grey = absent (floor-fail), same as log2fc
  FILL_NAME <- "ΔFRiP\n(per 1e6 frags,\n+ = Left)"
} else {
  hm <- hm %>% mutate(fillval = log2FoldChange)
  FILL_NAME <- "L/R log2FC\n(+ = Left)"
}

ylim   <- range(peaks$mid) / MB + c(-0.01, 0.01)
PITX2C <- 3971385 / MB   # PITX2c TSS — the focal point

# ── optional GC% / CpG-island strip (G4F_GC_TSV=<path>; fig04 sets it) ────────
# A thin VERTICAL viridis strip (sequential, DISTINCT from the L/R blue-red so GC is
# not misread as laterality) showing GC% across the TAD: the CpG islands at the PITX2
# promoter cluster (~88%, the focal L-bias substrate), FAM241A, and downstream gene
# promoters pop bright against the GC-poor desert. Off by default → exploratory renders
# of this engine are unchanged.
GC_TSV <- Sys.getenv("G4F_GC_TSV", "")
add_gc <- nzchar(GC_TSV) && file.exists(GC_TSV)
if (add_gc) {
  gc <- read_tsv(GC_TSV, show_col_types = FALSE) %>%
    filter(pos >= ylim[1] * MB, pos <= ylim[2] * MB) %>%
    mutate(y_mb = pos / MB)
  gc_h <- if (nrow(gc) > 1) stats::median(diff(sort(gc$pos))) / MB * 1.05 else 0.0002
  cat(sprintf("GC strip: %d windows, %.0f-%.0f%%\n", nrow(gc), min(gc$gc_pct), max(gc$gc_pct)))
}

# ── (1) gene track (CORRECTED FAM241A; conservative LOC naming) ──────────────
genes <- tribble(
  ~gene,                  ~start,    ~end,      ~strand, ~key,
  "FAM241A",              3387963,   3398391,   "-",     TRUE,    # corrected coding (XM_027456395.2); NOT the 3,655,942 Gnomon model
  "LOC106016802\n(Playrr?)", 3665204,3752463,  "-",     TRUE,    # candidate Playrr ortholog
  "LOC106016804",         3891849,   3954422,   "-",     FALSE,
  "PITX2",                3961387,   3976254,   "+",     TRUE,
  "LOC113843538",         3963664,   3964114,   "-",     FALSE,
  "LOC113843539",         3973586,   3976858,   "-",     FALSE,
  "LOC106017192\n(PANCR?)",3976875,  4033262,   "+",     TRUE,    # candidate PANCR ortholog
  "LOC101795212",         4003281,   4037227,   "-",     FALSE,
  "LOC119716870",         4041303,   4071427,   "+",     FALSE,
  "ELOVL6",               4105940,   4165489,   "+",     TRUE,
  "LOC113843470",         4159541,   4162591,   "-",     FALSE,
  "LOC113843471",         4168730,   4179954,   "-",     FALSE,
  "LOC119716871",         4180310,   4192324,   "-",     FALSE,
  "EGF",                  4194846,   4263663,   "-",     FALSE,
  "LRIT3",                4273391,   4285294,   "-",     FALSE,
  "RRH",                  4285433,   4300222,   "-",     FALSE,
  "GAR1",                 4302081,   4306236,   "-",     FALSE,
  "CFI",                  4306428,   4318961,   "+",     FALSE,
  "PLA2G12A",             4319835,   4325918,   "+",     FALSE,
  "CASP6",                4327886,   4337093,   "+",     FALSE,
  "MCUB",                 4337483,   4364812,   "-",     FALSE,
  "SEC24B",               4369737,   4415416,   "-",     FALSE,
  "COL25A1",              4451303,   4458808,   "+",     FALSE    # clipped to window max
) %>%
  filter(end >= ylim[1]*MB, start <= ylim[2]*MB) %>%
  mutate(s_mb = start/MB, e_mb = end/MB,
         ystt_mb = ifelse(strand == "+", s_mb, e_mb),
         yend_mb = ifelse(strand == "+", e_mb, s_mb))

# greedy lane packing (squish overlapping genes into separate lanes)
genes <- genes %>% arrange(start); genes$lane <- NA_integer_; lane_end <- numeric(0)
for (i in seq_len(nrow(genes))) {
  placed <- FALSE
  for (l in seq_along(lane_end))
    if (genes$start[i] > lane_end[l]) { genes$lane[i] <- l; lane_end[l] <- genes$end[i]; placed <- TRUE; break }
  if (!placed) { genes$lane[i] <- length(lane_end) + 1; lane_end <- c(lane_end, genes$end[i]) }
}
nlane <- max(genes$lane)

p_genes <- ggplot(genes) +
  geom_segment(aes(x = lane, xend = lane, y = ystt_mb, yend = yend_mb,
                   colour = key, linewidth = key),
               arrow = arrow(length = unit(0.045, "in"), type = "closed")) +
  scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey62"), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.5), guide = "none") +
  geom_text_repel(data = filter(genes, key),
                  aes(x = lane, y = (s_mb + e_mb)/2, label = gene),
                  direction = "y", hjust = 1, nudge_x = -1.6, segment.size = 0.3,
                  segment.colour = "grey55", size = 2.7, fontface = "bold",
                  lineheight = 0.82, min.segment.length = 0, box.padding = 0.3, max.overlaps = Inf) +
  geom_text_repel(data = filter(genes, !key),
                  aes(x = lane, y = (s_mb + e_mb)/2, label = gene),
                  direction = "y", hjust = 1, nudge_x = -1.6, segment.size = 0.2,
                  segment.colour = "grey78", size = 1.7, colour = "grey55",
                  min.segment.length = 0, box.padding = 0.12, max.overlaps = Inf) +
  scale_y_reverse(limits = rev(ylim), expand = expansion(0)) +
  scale_x_continuous(limits = c(-2.7, nlane + 0.25)) +
  labs(y = "NC_051775.1 position (Mb)", x = NULL, title = "Gene models (squished)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.title = element_text(size = 8, hjust = 0),
        axis.title.y = element_text(size = 8, margin = margin(r = 1)),
        plot.margin = margin(2, 1, 2, 2))

# ── promoter + notable-peak marks: dogleg leaders into the heatmap left edge ──
# No glyphs — a tight left gutter holds the labels; from each label an oblique
# segment turns HORIZONTAL to land on that feature's band row at the heatmap edge.
# label_y = where the (vertically de-crowded) label sits; pos/MB = the true row it points to.
X_LAB <- 0.18; X_KNEE <- 0.36; X_EDGE <- 0.54   # tile col 1 spans 0.54–1.46
marks <- tribble(
  ~pos,      ~cls,     ~lab,                  ~label_y,
  3961387,   "prom",   "PITX2a",              3.905,
  3964253,   "prom",   "PITX2b",              3.945,
  3971385,   "focal",  "PITX2c",              3.987,
  3976783,   "Lpeak",  "+5.4 kb (putative)",  4.028,
  3444117,   "Rpeak",  "FAM241A-border (R)",  3.444
) %>% mutate(ty = pos / MB)
# Focal region = PITX2 ±50 kb (the canonical d2i / locus-FDR proximal window,
# NC_051775.1:3,911,387–4,026,254) — drawn as a bracket on the heatmap's RIGHT edge.
FOCAL_LO <- 3911387 / MB; FOCAL_HI <- 4026254 / MB
X_BR <- max(cols$x_pos) + 0.72   # bracket x, right of the dorsal column

vmax <- max(abs(hm$fillval), na.rm = TRUE)
p_heat <- ggplot(hm, aes(x = x_pos, y = y_mb)) +
  geom_tile(aes(fill = fillval), width = 0.92, height = 0.0028) +
  geom_point(data = filter(hm, sig), shape = 8, size = 0.9, stroke = 0.3, colour = "black") +
  # focal (PITX2c) line — restricted to the band/heatmap region only (tile left→right edge),
  # i.e. between the peak-annotation gutter and the focal bracket, not across either.
  annotate("segment", x = min(cols$x_pos) - 0.46, xend = max(cols$x_pos) + 0.46,
           y = PITX2C, yend = PITX2C, colour = "#1b7837", linewidth = 0.4, linetype = "22") +
  # dogleg leaders: oblique (label → knee@target_y) then horizontal (knee → edge)
  geom_segment(data = marks, aes(x = X_LAB, y = label_y, xend = X_KNEE, yend = ty, colour = cls),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = marks, aes(x = X_KNEE, y = ty, xend = X_EDGE, yend = ty, colour = cls),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = marks, aes(x = X_LAB, y = label_y, label = lab, colour = cls),
            inherit.aes = FALSE, size = 2.4, fontface = "bold", hjust = 1) +
  # focal-region bracket on the RIGHT edge: PITX2 ±50 kb (vertical bar + end-caps + label)
  annotate("segment", x = X_BR, xend = X_BR, y = FOCAL_LO, yend = FOCAL_HI,
           colour = "#1b7837", linewidth = 0.5) +
  annotate("segment", x = X_BR, xend = X_BR - 0.13, y = FOCAL_LO, yend = FOCAL_LO,
           colour = "#1b7837", linewidth = 0.5) +
  annotate("segment", x = X_BR, xend = X_BR - 0.13, y = FOCAL_HI, yend = FOCAL_HI,
           colour = "#1b7837", linewidth = 0.5) +
  annotate("text", x = X_BR + 0.12, y = (FOCAL_LO + FOCAL_HI)/2, label = "focal\n±50 kb",
           colour = "#1b7837", size = 2.3, fontface = "bold", hjust = 0, lineheight = 0.85) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
                       limits = c(-vmax, vmax), name = FILL_NAME,
                       na.value = "grey80") +   # grey = peak absent (<floor) in this column
  scale_colour_manual(values = c(focal = "#1b7837", prom = "grey30",
                                 Lpeak = "#2166ac", Rpeak = "#b2182b"), guide = "none") +
  scale_x_continuous(breaks = cols$x_pos, labels = cols$xlab_disp,
                     limits = c(-1.05, X_BR + 0.85)) +
  scale_y_reverse(limits = rev(ylim), expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.x = element_text(size = 6.6),
        legend.key.width = unit(0.3, "cm"), plot.margin = margin(2, 2, 2, 0))

# Title/subtitle intentionally OMITTED from the figure — they live in LEGEND.md
# (production-grade panel to be folded into fig03).
# heatmap track + figure width scale with the number of testable columns so the
# tiles keep a constant aspect (4 cols -> the canonical 9.2in; 7 cols -> wider).
ncol_test <- nrow(cols)
heat_w <- 0.42 * ncol_test + 0.40            # 4 cols -> 2.08 (canonical)
fill_desc <- if (FILL == "dfrip") {
    "fill = ΔFRiP (Δ mean CPM, L−R; depth-normalised difference, sparsity-robust)"
  } else { "fill = DESeq2 L/R log2FC" }
cap <- paste0("Replicate unit = embryo-half (",
              if (DESIGN == "paired") "paired ~genotype+side, within-embryo L vs R"
                else "unpaired ~side, all halves",
              ").  ", fill_desc, ".  grey tile = peak absent (pooled count < floor) in that ",
              "column; white = tested but bilateral (≈0).  ✳ padj_locus < 0.1 (DESeq2).")
# ── optional GC% strip (vertical viridis heatmap lane), between genes and heatmap ──
GC_W <- 0.30
if (add_gc) {
  p_gc <- ggplot(gc, aes(x = 1, y = y_mb, fill = gc_pct)) +
    geom_tile(width = 1, height = gc_h) +
    scale_fill_viridis_c(name = "GC %", limits = c(10, 90), breaks = c(25, 50, 75),
                         oob = scales::squish,
                         guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    scale_y_reverse(limits = rev(ylim), expand = expansion(0)) +
    scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "GC %") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 8, hjust = 0.5),
          legend.position = "bottom", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5),
          legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.22, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(2, 1, 2, 1))
  body <- p_genes + p_gc + p_heat + plot_layout(widths = c(0.92, GC_W, heat_w))
  fig_w <- 3.07 * (0.92 + GC_W + heat_w)
} else {
  body <- p_genes + p_heat + plot_layout(widths = c(0.92, heat_w))
  fig_w <- 3.07 * (0.92 + heat_w)
}
p <- body +
  plot_annotation(caption = cap,
                  theme = theme(plot.caption = element_text(size = 6, hjust = 0,
                                                            colour = "grey35")))

DSUF <- if (DESIGN == "paired") "_paired" else ""   # unpaired keeps the canonical name
FSUF <- if (FILL == "dfrip") "_dfrip" else ""
outname <- if (nzchar(G4F_OUTNAME)) {
  G4F_OUTNAME
} else {
  paste0("panel_D_lr_genetrack", SUF, DSUF, FSUF, ".png")
}
out <- file.path(G4F_OUTDIR, outname)
ggsave(out, p, width = fig_w, height = 9.4, dpi = 200)
cat("wrote", out, "\n")
