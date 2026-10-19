#!/usr/bin/env Rscript
# cAD plot — chick PITX2-locus pileup, faithful duck fig02D adaptation (stage x population).
# Reads the cAD cluster outputs (depth + gene model + polyA priming) and renders the duck
# generate_pitx2_pileup_multistage.R layout: stacked stage x population coverage tracks with a
# two-tier y-axis (dorsal natural / non-dorsal amplified), isoform first-exon highlight bands,
# polyA internal-priming markers, and a 3-isoform gene model (UTR/CDS/promoters) from the GTF.
#
#   Rscript pipeline/chick/cAD_plot_pileup.R
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)})

# CAD_DATA overrides the input dir (e.g. a side directory holding a fresh cAD run, so it can be
# rendered + compared before it replaces the committed copy); default = the committed data.
D <- Sys.getenv("CAD_DATA"); if (!nzchar(D)) D <- "figures/data/chick_control"
depth <- read_tsv(file.path(D, "cAD_pitx2_pileup_depth.tsv.gz"), show_col_types = FALSE, name_repair = "minimal")
trk   <- read_tsv(file.path(D, "cAD_pitx2_pileup_tracks.tsv"),  show_col_types = FALSE)
gm    <- read_tsv(file.path(D, "cAD_pitx2_gene_model.tsv"),     show_col_types = FALSE)
prime <- read_tsv(file.path(D, "cAD_pitx2_polya_priming.tsv"),  show_col_types = FALSE)

# ── render parameters (env toggles: manuscript figS04 vs exploratory) ──────────
# Manuscript figS04 Panel E (default): NO antisense track, NO title (compositor
# letters it), 300 dpi. Exploratory render keeps the antisense CONTROL track that
# proves the 5' sense signal is genuine a/b transcription (antisense ~0 there) +
# the title. The antisense control's validation role is cited in LEGEND.md.
ANTISENSE  <- Sys.getenv("CAD_ANTISENSE", "0") == "1"
SHOW_TITLE <- Sys.getenv("CAD_TITLE", "0") == "1"
# default render path: docs/chick_control/ (not distributed; set CAD_OUT to write elsewhere)
OUT_PNG    <- Sys.getenv("CAD_OUT", "docs/chick_control/cAD_pitx2_pileup.png")  # (not distributed)
OUT_DPI    <- as.numeric(Sys.getenv("CAD_DPI", "300"))
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)   # comma-separated n=

SMOOTH <- 50
smooth <- function(y) as.numeric(stats::filter(y, rep(1/SMOOTH, SMOOTH), sides = 2))

# track order + compartment (fill) + dorsal flag (two-tier y-axis)
ord <- c("HH13|LPM mesenchyme", "HH15-16|Dorsal-biased", "HH15-16|Ventral (pooled)",
         "HH25|Dorsal sub-epi.", "HH25|Ventral spl.", "HH25|Syringeal-Cardio")
comp <- c("HH13|LPM mesenchyme"="All LPM", "HH15-16|Dorsal-biased"="Dorsal",
          "HH15-16|Ventral (pooled)"="Ventral", "HH25|Dorsal sub-epi."="Dorsal",
          "HH25|Ventral spl."="Ventral", "HH25|Syringeal-Cardio"="Syringeal")
fill_vals <- c("All LPM"="gray35", "Dorsal"="black", "Ventral"="gray45", "Syringeal"="gray55",
               "Antisense"="firebrick")
dorsal_comp <- c("Dorsal")
NON_DORSAL_AMP <- 3.0

# antisense control cell count (capture before filtering trk to the sense-track order)
n_anti <- { v <- trk$n_cells[trk$track == "ANTISENSE control"]; if (length(v)) v[1] else NA_real_ }

trk <- trk %>% filter(track %in% ord)
present <- ord[ord %in% trk$track]
# build long depth: per track, reconstruct the cAD column key = gsub('|','__',track) + '_percell'
long <- lapply(present, function(t) {
  col <- paste0(trk$safe_key[trk$track == t], "_percell")
  if (!col %in% names(depth)) return(NULL)
  n <- trk$n_cells[trk$track == t]
  data.frame(position = depth$position, depth = smooth(depth[[col]]),
             track = t, compartment = comp[[t]], n_cells = n)
}) %>% bind_rows()
# antisense (-f 16) CONTROL as its own bottom track: antisense ~0 at the sense 5' peak
# (exploratory render only; manuscript figS04 Panel E omits it — see LEGEND.md)
if (ANTISENSE && "antisense_percell" %in% names(depth)) {
  long <- bind_rows(long, data.frame(
    position = depth$position, depth = smooth(depth$antisense_percell),
    track = "ANTISENSE control", compartment = "Antisense",
    n_cells = ifelse(is.na(n_anti), 0L, as.integer(n_anti))))
  present <- c(present, "ANTISENSE control")
}
# track_label levels pulled from long itself (robust to the appended antisense row)
lvl <- long %>% distinct(track, n_cells)
lvl <- lvl[match(present, lvl$track), ]
long$track_label <- factor(paste0(gsub("\\|", "\n", long$track), " (n=", fmt_n(long$n_cells), ")"),
                           levels = paste0(gsub("\\|", "\n", lvl$track), " (n=", fmt_n(lvl$n_cells), ")"))

GENE_START <- min(gm$start[gm$feature == "gene"]); GENE_END <- max(gm$end[gm$feature == "gene"])
PLOT_START <- GENE_START - 1500; PLOT_END <- GENE_END + 1500
long <- long %>% filter(position >= PLOT_START, position <= PLOT_END)

# two-tier y ceilings
dmax <- max(long$depth[long$compartment %in% dorsal_comp], na.rm = TRUE)
amax <- dmax / NON_DORSAL_AMP
ylim <- long %>% distinct(track_label, compartment) %>%
  mutate(y_max = ifelse(compartment %in% dorsal_comp, dmax, amax), pos_kb = PLOT_START/1e3)
cat(sprintf("y-tiers: dorsal=%.4f amplified=%.4f (%.1fx)\n", dmax, amax, NON_DORSAL_AMP))

# isoform first exons (smallest-start exon per transcript) -> highlight bands + labels a/b/c
tx_first <- gm %>% filter(feature == "exon") %>% group_by(transcript_id) %>%
  summarise(s = min(start), e = start[which.min(start)] , ee = end[which.min(start)], .groups="drop") %>%
  arrange(s)
tx_first$iso <- c("a","b","c")[seq_len(nrow(tx_first))]
bands <- data.frame(xmin = tx_first$e/1e3, xmax = tx_first$ee/1e3)
# PITX2c TSS (downstream of a/b on + strand): any SENSE signal 5' of this line can only be a/b
c_tss <- suppressWarnings(min(tx_first$s[tx_first$iso == "c"]))
cat(sprintf("PITX2c TSS = %s; sense signal 5' of it is a/b-only\n", c_tss))

px <- function(p) p/1e3
p_tracks <- ggplot(long, aes(position/1e3, depth)) +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, fill = "gray92", color = NA) +
  { if (nrow(prime)) geom_point(data = prime, aes(x = position/1e3, y = Inf),
            inherit.aes = FALSE, shape = 25, size = 0.7, fill = "gray60", color = "gray60") } +
  geom_area(aes(fill = compartment), alpha = 0.9, linewidth = 0) +
  { if (is.finite(c_tss)) geom_vline(xintercept = c_tss/1e3, linetype = "dashed",
            color = "firebrick", linewidth = 0.35) } +
  scale_fill_manual(values = fill_vals, guide = "none") +
  geom_blank(data = ylim, aes(x = pos_kb, y = y_max)) +
  facet_grid(track_label ~ ., scales = "free_y", switch = "y") +
  scale_x_continuous(limits = c(px(PLOT_START), px(PLOT_END)), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), breaks = scales::pretty_breaks(2)) +
  labs(y = sprintf("UMIs per cell (%d bp smooth)", SMOOTH)) +
  theme_classic(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, lineheight = 1.05),
        strip.background = element_blank(), strip.placement = "outside",
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6),
        panel.spacing = unit(0.3, "lines"), plot.margin = margin(5,10,0,5))

# gene model: per transcript (y=3 a / 2 b / 1 c) draw UTR (light) + CDS (dark) + intron line + TSS arrow
ymap <- setNames(c(3,2,1), tx_first$transcript_id)
p_gene <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -0.5, ymax = 4.5),
            fill = "gray92", color = NA)
for (i in seq_len(nrow(tx_first))) {
  tid <- tx_first$transcript_id[i]; y <- ymap[[tid]]
  sub <- gm %>% filter(transcript_id == tid)
  ex  <- sub %>% filter(feature == "exon") %>% arrange(start)
  p_gene <- p_gene + annotate("segment", x = px(min(ex$start)), xend = px(max(ex$end)),
                              y = y, yend = y, linewidth = 0.35)
  utrs <- sub %>% filter(feature == "UTR"); cds <- sub %>% filter(feature == "CDS")
  if (nrow(utrs)) p_gene <- p_gene + annotate("rect", xmin = px(utrs$start), xmax = px(utrs$end),
                              ymin = y-0.18, ymax = y+0.18, fill = "#606060", color = "black", linewidth = 0.15)
  if (nrow(cds))  p_gene <- p_gene + annotate("rect", xmin = px(cds$start), xmax = px(cds$end),
                              ymin = y-0.30, ymax = y+0.30, fill = "#303030", color = "black", linewidth = 0.15)
  p_gene <- p_gene +
    annotate("text", x = px(min(ex$start)) - 0.3, y = y, label = paste0("PITX2", tx_first$iso[i]),
             size = 2.4, fontface = "italic", hjust = 1) +
    annotate("segment", x = px(min(ex$start)), xend = px(min(ex$start)), y = y+0.35, yend = y+0.62, linewidth = 0.3) +
    annotate("segment", x = px(min(ex$start)), xend = px(min(ex$start))+0.4, y = y+0.62, yend = y+0.62,
             linewidth = 0.3, arrow = arrow(length = unit(1.1,"mm"), type = "closed"))
}
if (nrow(prime)) p_gene <- p_gene + geom_point(data = prime, aes(x = position/1e3, y = 4.3),
            shape = 25, size = 0.7, fill = "gray60", color = "gray60")
if (is.finite(c_tss)) p_gene <- p_gene +
  geom_vline(xintercept = c_tss/1e3, linetype = "dashed", color = "firebrick", linewidth = 0.35) +
  annotate("text", x = c_tss/1e3, y = 4.75, label = "PITX2c TSS", color = "firebrick",
           size = 2.3, hjust = -0.04, fontface = "bold") +
  annotate("text", x = c_tss/1e3 - 0.12, y = -0.35, label = "5' sense = a/b only",
           color = "firebrick", size = 2.1, hjust = 1)
p_gene <- p_gene +
  scale_x_continuous(name = sprintf("%s position (kb)  |  bands = isoform first exons; v = polyA internal-priming; dashed = PITX2c TSS",
                                    depth$chrom[1]),
                     limits = c(px(PLOT_START), px(PLOT_END)), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.5, 4.8)) +
  theme_void(base_size = 9) +
  theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 7),
        axis.ticks.x = element_line(linewidth = 0.3), axis.ticks.length.x = unit(1.5,"mm"),
        axis.line.x = element_line(linewidth = 0.3), plot.margin = margin(0,10,5,5))

combined <- p_tracks / p_gene + plot_layout(heights = c(10, 2))
if (SHOW_TITLE) combined <- combined +
  plot_annotation(title = "Chick PITX2-locus snRNA pileup (Group B, stage × population; sense-strand, per-cell)",
                  subtitle = "duck fig02D analog — total PITX2 (isoform first-exons unresolved, cZ); dorsal natural / non-dorsal amplified",
                  theme = theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8)))
dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, combined, width = 10.5, height = 5.6, dpi = OUT_DPI, bg = "white")
cat("Saved", OUT_PNG, "\n")
