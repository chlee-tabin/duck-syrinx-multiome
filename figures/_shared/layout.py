"""Canonical panel widths + font sizing for every figure generator.

Every panel generator should author at a width chosen from this module
and apply ``apply_layout_fonts()`` before plotting. The compositor renders
panels at their native size (no per-column stretching), so if two panels
share a row, they should share a width class (or two complementary
classes — e.g. TWO_THIRD + ONE_THIRD).

Page width target: 10.5 inches — deliberately wider than Nature's 183 mm
double-column so panels do not have to cram. Figures can be scaled down
by the journal layout stage if needed; we optimize for readability first.

Use in matplotlib:
    from _shared.layout import W_HALF, H_STD, DPI, apply_layout_fonts
    apply_layout_fonts()
    fig = plt.figure(figsize=(W_HALF, H_STD), dpi=DPI)

Use in svgwrite: pass the width in pixels at 96 dpi, e.g.
    size=(f"{W_FULL * 96}px", f"{height * 96}px")

Use in R/ggsave:
    ggsave(path, plot, width=W_HALF, height=H_STD, units="in", dpi=300)

The compositor uses DPI=300 to recover each PNG's physical width in
inches. Any panel that uses dpi != 300 will appear at the wrong scale.
"""
import os

# ── Canonical widths (inches) ────────────────────────────────────────────
# Widened 2026-04-17 from the original 7.2" Nature-double-column target —
# per-panel cramming was forcing tick labels on top of each other. Page
# can be shrunk in final layout if a journal demands 183 mm.
W_FULL       = 10.5
W_TWO_THIRD  = 7.00
W_HALF       = 5.25
W_ONE_THIRD  = 3.50
W_QUARTER    = 2.625

# Convenience height classes — authors may pick any height, but these
# work well with the widths above.
H_SHORT = 2.5
H_STD   = 3.75
H_TALL  = 5.5
H_XTALL = 7.5

# Rendering discipline — every author must use DPI for raster panels.
DPI = 300

# Gap between panels in the composite (points). 12 pt ≈ 1/6 inch.
COMPOSITE_GAP_PT = 12.0


# ── Font sizing for the canonical widths ─────────────────────────────────
# After the 2026-04-16 layout migration dropped panel widths without
# touching matplotlib's default font size, text-to-canvas ratios blew up
# on every panel. Call this helper from every generator so tick labels,
# legends, and axis titles scale with the new canvas.
_LAYOUT_RCPARAMS = {
    "font.size":        7,     # base text
    "axes.titlesize":   8,     # panel subtitles (outer letters come from compositor)
    "axes.titleweight": "bold",
    "axes.labelsize":   7,     # axis labels
    "xtick.labelsize":  6,
    "ytick.labelsize":  6,
    "legend.fontsize":  6,
    "legend.title_fontsize": 6,
    "figure.titlesize": 8,
    "figure.dpi":       DPI,
    "savefig.dpi":      DPI,
    "savefig.bbox":     "tight",
}


def apply_layout_fonts():
    """Set rcParams sized for the canonical widths in this module.

    Call once at the top of every matplotlib panel generator, after any
    ``apply_pub_style()`` call (or instead of it). Idempotent — rcParams
    are process-global but reapplying is safe.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update(_LAYOUT_RCPARAMS)


# ── Panel letters (house convention, 2026-09-07) ─────────────────────────
# Capital + trailing period, italic, NOT bold, one face everywhere: Arial
# Italic, embedded as TrueType.  Rationale: `family="Helvetica",
# fontstyle="italic"` silently resolved to Helvetica REGULAR (matplotlib
# cannot pick the oblique face out of macOS's Helvetica.ttc), so every
# composite since 2026-08-26 carried upright letters, and the DejaVu Sans
# Oblique the panel producers baked in reads as bold next to Helvetica.
def panel_letter_font(size_pt: float):
    """FontProperties for a panel letter — fails loud if no italic face."""
    from matplotlib import font_manager as fm
    path = fm.findfont(fm.FontProperties(
        family=["Arial", "Liberation Sans", "DejaVu Sans"], style="italic"))
    base = os.path.basename(path).lower()
    if "italic" not in base and "oblique" not in base:
        raise RuntimeError(f"panel_letter_font: resolved a non-italic face {path!r}")
    return fm.FontProperties(fname=path, size=size_pt)


def panel_letter(ax, label: str, size_pt: float = 9.0) -> None:
    """Bake a panel letter into an axes as its left-aligned title.
    Bare capital letters get the trailing period ("A."); sub-labels such
    as "A1" are drawn as-is."""
    text = label + "." if (len(label) == 1 and "A" <= label <= "Z") else label
    ax.set_title(text, fontproperties=panel_letter_font(size_pt), loc="left")
