"""Fig S8 (fig01_supp) tick-label helpers -- local to this figure on purpose,
so nothing here can change another figure's output."""


def vertical_xticklabels(ax, labels, **kw):
    """Set x tick labels at exactly 90 deg, each label's top end anchored on
    its tick (rotation_mode='anchor') so it sits directly under the tick."""
    ax.set_xticklabels(labels, rotation=90, ha="right", va="center",
                       rotation_mode="anchor", **kw)


def tick_label_bottom(ax, indices=None) -> float:
    """Lowest edge of ax's x tick labels, in axes-y coordinates (negative).

    Draws the figure once to measure the rendered labels, so annotations
    placed below the axis (group brackets) sit directly under the labels.
    With indices, measure only those labels (per-group brackets).
    Raises if the axis has no visible tick labels to measure.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = ax.get_xticklabels()
    if indices is not None:
        labels = [labels[i] for i in indices]
    boxes = [t.get_window_extent(renderer) for t in labels
             if t.get_visible() and t.get_text()]
    if not boxes:
        raise RuntimeError("tick_label_bottom: axis has no x tick labels")
    y_disp = min(b.y0 for b in boxes)
    return ax.transAxes.inverted().transform((0, y_disp))[1]
