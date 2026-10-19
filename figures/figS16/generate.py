"""Compose figS16 (P-S1 proliferation-timing supplement) -> figures/figS16.pdf."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _shared.composite import save_composite  # noqa: E402


def main() -> None:
    rows = [
        [str(HERE / "panel_A_cycling_paired.png"),
         str(HERE / "panel_B_falsification_share.png")],
        [str(HERE / "panel_C_mitogen_program.png")],
    ]
    labels = [["A", "B"], ["C"]]
    out = HERE.parent / "figS16.pdf"
    save_composite(rows, out_pdf=str(out), labels=labels, label_gutter_pt=14.0)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
