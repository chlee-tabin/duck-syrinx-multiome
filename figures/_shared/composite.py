"""Figure compositor — PDF-only, vector-preserving, globally consistent zoom.

Design rules (per 2026-04-16 session):
  1. Output is a single PDF per figure. No PNG composite. No per-panel PDF
     sidecars. SVG sources remain the only vector panel form on disk.
  2. Panels are never upscaled beyond their native size when the global
     scale factor would be >1 — instead we keep S=1 and let the page be as
     wide as the widest row at native size. If a caller forces a smaller
     target width, every panel scales down by the SAME factor so effective
     zoom is consistent across the whole figure.
  3. Missing panels are a hard error. No silent skipping.
  4. Labels are vector text drawn via matplotlib overlay merged with pypdf.
     A None entry in the labels grid means "this tile already has internal
     letters baked in" (sub-composite PNGs like panel_libqc) — compositor
     skips it.

The entry point is ``save_composite``. ``compose_pdf`` and
``_panel_to_pdf`` are the internals.
"""
from __future__ import annotations

import os
import tempfile
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from .layout import COMPOSITE_GAP_PT, panel_letter_font

try:
    import cairosvg  # type: ignore
except ImportError:  # pragma: no cover
    cairosvg = None

try:
    from pypdf import PdfReader, PdfWriter, Transformation  # type: ignore
except ImportError:  # pragma: no cover
    PdfReader = PdfWriter = Transformation = None  # type: ignore


# ── Pre-compose utility: stack multiple PNGs vertically into one ────────
def stack_vertical(paths: list[str], out_path: str,
                   gap_px: int = 30, bg: str = "white") -> str:
    """Stack PNG panels vertically into a single PNG. Useful when two
    sub-variants of a panel should appear as one logical tile (e.g. Fig 2
    panel C: PITX2+ fraction above mean CP10K). Returns `out_path`.
    """
    from PIL import Image as _Image
    imgs = []
    for p in paths:
        img = _Image.open(p)
        imgs.append(img.convert("RGB") if img.mode != "RGB" else img)
    max_w = max(img.width for img in imgs)
    total_h = sum(img.height for img in imgs) + gap_px * (len(imgs) - 1)
    canvas = _Image.new("RGB", (max_w, total_h), bg)
    y = 0
    for img in imgs:
        canvas.paste(img, ((max_w - img.width) // 2, y))
        y += img.height + gap_px
    canvas.save(out_path, dpi=(300, 300))
    return out_path


# ── Panel → single-page PDF (vector for SVG, embedded raster for PNG) ────
def _panel_to_pdf(path: str, tmpdir: str) -> str:
    base = os.path.splitext(os.path.basename(path))[0]
    pdf_path = os.path.join(tmpdir, f"{base}.pdf")
    if path.lower().endswith(".svg"):
        if cairosvg is None:
            raise RuntimeError(f"cairosvg required to embed {path}")
        cairosvg.svg2pdf(url=path, write_to=pdf_path)
    elif path.lower().endswith(".pdf"):
        return path  # already a PDF; use directly
    else:
        img = Image.open(path)
        if img.mode == "RGBA":
            bg = Image.new("RGB", img.size, "white")
            bg.paste(img, mask=img.split()[3])
            img = bg
        else:
            img = img.convert("RGB")
        # Save at 300 DPI — native page size = pixels/300 inches
        img.save(pdf_path, "PDF", resolution=300.0)
    return pdf_path


# ── Vector label overlay (matplotlib → PDF) ──────────────────────────────
def _draw_label_overlay(
    positions: list[tuple[float, float, str]],
    target_w_pt: float,
    total_h_pt: float,
    out_path: str,
    label_size_pt: float,
) -> None:
    """Render a transparent-background PDF containing just the panel labels.
    `positions` entries are (x_pt, y_pt_from_bottom, text) — va='top',
    ha='left'. Vector text (Arial Italic, TrueType-embedded) survives PDF
    merging. House letter convention: capital italic with a
    trailing period ("A."), not bold — callers still pass the bare letter.
    The face itself comes from `layout.panel_letter_font` (one source).
    """
    fig = plt.figure(figsize=(target_w_pt / 72.0, total_h_pt / 72.0), dpi=72)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, target_w_pt)
    ax.set_ylim(0, total_h_pt)
    ax.set_axis_off()
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)
    fp = panel_letter_font(label_size_pt)
    for x_pt, y_pt, text in positions:
        if not (len(text) == 1 and "A" <= text <= "Z"):
            raise ValueError(f"panel label must be a bare capital letter, got {text!r}")
        ax.text(x_pt, y_pt, text + ".", fontproperties=fp, ha="left", va="top")
    # fonttype 42 = embed the TrueType outlines (Illustrator keeps the face;
    # the Type 3 default is not editable text there).
    with plt.rc_context({"pdf.fonttype": 42}):
        fig.savefig(out_path, format="pdf", transparent=True,
                    bbox_inches=None, pad_inches=0)
    plt.close(fig)


# ── Core PDF composer ────────────────────────────────────────────────────
def compose_pdf(
    rows: list[list[str]],
    labels: Optional[list[list[Optional[str]]]],
    out_pdf: str,
    target_w_pt: Optional[float] = None,
    gap_pt: float = 12.0,
    label_size_pt: float = 10.0,
    label_pad_pt: float = 6.0,
    row_valign: Optional[list[str]] = None,
    label_gutter_pt: float = 0.0,
) -> None:
    """Compose a vector-preserving PDF figure.

    Args:
        rows: list of rows, each row a list of panel file paths. Rows may
            have different column counts.
        labels: list parallel to `rows`; each sublist parallel to the row.
            Entry is the outer letter to draw (e.g. "A") or None to skip
            (tile has internal letters baked in).
        out_pdf: destination .pdf path.
        target_w_pt: page width in PDF points. If None (default), page is
            sized to the widest row's native width (S=1 everywhere — no
            upscaling). If set, every panel scales by the SAME factor
            target_w_pt/widest_row_native_w so effective zoom is consistent
            across the whole figure.
        gap_pt: spacing between panels (horizontal) and rows (vertical).
        label_size_pt: font size for outer letter labels.
        label_pad_pt: label offset from top-left of the panel cell.
        row_valign: per-row vertical alignment for shorter panels within
            the row height.  Each entry is "center" (default — shorter
            panels sit vertically centered with equal whitespace above and
            below the content) or "top" (shorter panels flush with the
            letter label, whitespace falls at the bottom).  Length must
            match `rows` when provided; None means "center everywhere".
        label_gutter_pt: extra blank height reserved ABOVE every row that
            carries an outer letter, so the letter sits clear of panels
            whose PNG has no top-left margin (2026-09-07:
            figS16/figS17 letters overlapped the y-axis labels).
    """
    if PdfReader is None:
        raise RuntimeError(
            "pypdf is required for the PDF compositor. "
            "Activate the duckfig conda env or `pip install pypdf`.")

    # Hard-fail on missing panels — no silent skipping.
    missing = [p for row in rows for p in row if not os.path.exists(p)]
    if missing:
        joined = "\n  ".join(missing)
        raise FileNotFoundError(
            f"compose_pdf: {len(missing)} panel file(s) not found. "
            f"Cannot build {out_pdf}:\n  {joined}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Convert each panel to PDF; read native (width, height) in points.
        per_panel: list[list[tuple[str, float, float]]] = []
        for row in rows:
            row_entries: list[tuple[str, float, float]] = []
            for path in row:
                panel_pdf = _panel_to_pdf(path, tmpdir)
                reader = PdfReader(panel_pdf)
                page = reader.pages[0]
                w = float(page.mediabox.width)
                h = float(page.mediabox.height)
                row_entries.append((panel_pdf, w, h))
            per_panel.append(row_entries)

        # Widest row's native width drives the global scale factor.
        row_native_widths: list[float] = []
        for row in per_panel:
            if not row:
                row_native_widths.append(0.0)
                continue
            n = len(row)
            w = sum(panel_w for _, panel_w, _ in row) + gap_pt * (n - 1)
            row_native_widths.append(w)
        max_row_w = max(row_native_widths) if row_native_widths else 1.0
        if target_w_pt is None:
            target_w_pt = max_row_w
        S = target_w_pt / max_row_w if max_row_w > 0 else 1.0

        # Scaled row heights (vertical extent in output page), plus the
        # label gutter on rows that carry an outer letter.
        row_gutter_pt: list[float] = [
            label_gutter_pt if (labels is not None and i < len(labels)
                                and any(labels[i])) else 0.0
            for i in range(len(per_panel))]
        row_heights_pt: list[float] = []
        for row, gutter in zip(per_panel, row_gutter_pt):
            if not row:
                row_heights_pt.append(0.0)
                continue
            row_heights_pt.append(max(h for _, _, h in row) * S + gutter)

        n_nonempty = sum(1 for h in row_heights_pt if h > 0)
        total_h_pt = sum(row_heights_pt) + gap_pt * max(0, n_nonempty - 1)

        writer = PdfWriter()
        master_page = writer.add_blank_page(
            width=target_w_pt, height=total_h_pt)

        # Place panels (pypdf origin is bottom-left).
        label_positions: list[tuple[float, float, str]] = []
        y_top = total_h_pt
        for row_idx, (row, row_h) in enumerate(zip(per_panel, row_heights_pt)):
            if not row:
                continue
            y_top -= row_h  # y_top now = bottom of this row
            # Center-align the row within the page. Rows whose panels do
            # not sum to the page width (e.g. two half-width panels next
            # to a W_FULL panel on another row) were previously pinned
            # flush-left, producing asymmetric right-side whitespace.
            row_native_w = row_native_widths[row_idx]
            x_left = (target_w_pt - row_native_w * S) / 2
            # Per-row vertical alignment (center by default).  "top" pins
            # the content top-edge to the row top (directly under the
            # letter label); "center" leaves equal whitespace above and
            # below.  Labels always draw at the row top regardless of
            # content alignment — they anchor the panel cell, not the
            # content itself.
            valign = "center"
            if row_valign is not None and row_idx < len(row_valign):
                valign = row_valign[row_idx]

            for col_idx, (panel_pdf, w, h) in enumerate(row):
                reader = PdfReader(panel_pdf)
                panel_page = reader.pages[0]
                scaled_w = w * S
                scaled_h = h * S
                content_h = row_h - row_gutter_pt[row_idx]
                if valign == "top":
                    y_offset = y_top + (content_h - scaled_h)
                else:  # "center"
                    y_offset = y_top + (content_h - scaled_h) / 2
                transform = (Transformation()
                             .scale(S, S)
                             .translate(x_left, y_offset))
                master_page.merge_transformed_page(panel_page, transform)

                # Collect label for overlay
                lbl: Optional[str] = None
                if labels is not None and row_idx < len(labels):
                    row_labels = labels[row_idx]
                    if col_idx < len(row_labels):
                        lbl = row_labels[col_idx]
                if lbl:
                    label_y_pt = y_top + row_h - label_pad_pt
                    label_positions.append(
                        (x_left + label_pad_pt, label_y_pt, lbl))

                # Gaps scale with S, matching the row-width/centring maths
                # above (row_native_w includes native gaps). No-op at S == 1.
                x_left += scaled_w + gap_pt * S
            y_top -= gap_pt

        # Overlay vector labels
        if label_positions:
            overlay_path = os.path.join(tmpdir, "_labels.pdf")
            _draw_label_overlay(
                label_positions, target_w_pt, total_h_pt,
                overlay_path, label_size_pt)
            overlay_reader = PdfReader(overlay_path)
            master_page.merge_page(overlay_reader.pages[0])

        with open(out_pdf, "wb") as f:
            writer.write(f)

    # Summary line
    dims = f"{target_w_pt:.0f}x{total_h_pt:.0f} pt ({target_w_pt/72:.2f}x{total_h_pt/72:.2f} in)"
    n_panels = sum(len(row) for row in rows)
    n_labels = len(label_positions)
    print(f"  Saved: {out_pdf} [{dims}, {n_panels} panels, {n_labels} labels, scale={S:.3f}]")


# ── Public entry point ───────────────────────────────────────────────────
def save_composite(
    rows: list[list[str]],
    out_pdf: str,
    labels: Optional[list[list[Optional[str]]]] = None,
    target_w_pt: Optional[float] = None,
    gap_pt: float = 12.0,
    label_size_pt: float = 10.0,
    label_pad_pt: float = 6.0,
    row_valign: Optional[list[str]] = None,
    label_gutter_pt: float = 0.0,
) -> None:
    """Assemble rows of panel paths (PNG/SVG/PDF) into a composite PDF.

    This is the single entry point used by every fig*/generate.py. PNG
    composites are not produced — the PDF is the publication artifact and
    the reviewable deliverable.

    Missing panel files raise FileNotFoundError — no silent skipping.

    Set `target_w_pt=None` (default) for "native-size" output: the page is
    as wide as the widest row's native width, every panel at its native
    scale (S=1). Pass an explicit value (e.g. 540 for a 7.5-inch page) to
    force a specific page width; every panel scales by the SAME factor.

    `row_valign` per-row list of "center" (default) or "top".  Use "top"
    only when a short panel genuinely benefits from anchoring under its
    letter label (e.g. a companion schematic next to a wide data plot).
    Most rows look more balanced with "center".
    """
    compose_pdf(
        rows=rows,
        labels=labels,
        out_pdf=out_pdf,
        target_w_pt=target_w_pt,
        gap_pt=gap_pt,
        label_size_pt=label_size_pt,
        label_pad_pt=label_pad_pt,
        row_valign=row_valign,
        label_gutter_pt=label_gutter_pt,
    )
