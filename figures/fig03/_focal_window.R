# Shared focal display/count window for fig03 panels B and C — the SINGLE SOURCE
# OF TRUTH so the two panels cannot drift apart.  Panel C draws the L/R pileup
# over this window (background mask + ticks); panel B's "ATAC focal" row counts
# the L/R-significant peaks fully contained in it (same predicate), so the focal
# count is exactly the peak set panel C displays.  ~40 kb on NC_051775.1 covering
# the PITX2 gene, its a/b/c promoters, and the +5.4 kb putative element.
#
# NB one parallel copy lives outside R and cannot source this: the GC-track
# producer scripts/gc_content_pitx2_window.py (GC_START/GC_END defaults) must be
# kept equal to these bounds — it writes pitx2_locus_gc_content.tsv for panel C.
PITX2_FOCAL_WIN <- c(3958000, 3998000)
