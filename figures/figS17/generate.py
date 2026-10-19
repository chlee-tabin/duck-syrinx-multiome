"""Compose figS17 (E-S1 ESR1 left-asymmetry supplement) -> figures/figS17.pdf."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _shared.composite import save_composite  # noqa: E402


def main() -> None:
    rows = [
        # panel_C_depth_glmm CUT 2026-08-24: uncited technical control;
        # its numbers live in the legend's depth-control note. File + producer kept.
        [str(HERE / "panel_A_esr1_paired.png"),
         str(HERE / "panel_B_dorsal_control.png")],
    ]
    labels = [["A", "B"]]
    out = HERE.parent / "figS17.pdf"
    save_composite(rows, out_pdf=str(out), labels=labels, label_gutter_pt=14.0)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
