#!/usr/bin/env python3
"""Validate references/same_side_duplicate_genotypes.tsv against e3.

d2l and d2n read the committed same-side duplicate reference instead of the
downstream e3 composition table (``pipeline/_same_side_duplicates.py``,
``load_duplicate_reference``).  This validator recomputes the tube-split
evidence from e3 (``identify_duplicate_genotypes``: composition chi-square per
candidate genotype, one-tube nuclei-ratio test) and raises on any disagreement
with the reference -- a different candidate set, a flipped ``is_duplicate``,
different nuclei counts, or a chi2 / p / TV value outside rtol 1e-6.

Runs after e3 (cluster runner ``FC_08656/run_regenerate_figure_data.sh``, not distributed;
against the freshly written e3) and in ``make check`` (against the committed e3).

Usage:
    python3 pipeline/validate_same_side_duplicates.py
        [--composition PATH] [--reference PATH] [--design PATH]
    python3 pipeline/validate_same_side_duplicates.py --write-reference
        (regenerate the reference from e3 after a deliberate e3 change;
         review the diff before committing it)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _same_side_duplicates import (  # noqa: E402
    COMPOSITION_PATH, DESIGN_PATH, DUPLICATE_REFERENCE, EVIDENCE_COLS,
    identify_duplicate_genotypes, load_collapsed_design, validate_reference)

HEADER = (
    "# Same-side duplicate genotypes (lane-loading-record split tube: one embryo-half\n"
    "# sequenced in two libraries).  One row per candidate genotype with >= 100\n"
    "# nuclei in both libraries of a design duplicate pair; is_duplicate = the call.\n"
    "# Evidence: cell-type composition chi-square across the two lanes\n"
    "# (p_composition, tv_distance; no PITX2 input) and the one-tube nuclei-ratio\n"
    "# test (p_tube_ratio).  Read by d2l / d2n; recomputed from\n"
    "# figures/data/e3_composition.tsv by pipeline/validate_same_side_duplicates.py,\n"
    "# which raises on any disagreement.  Regenerate only with --write-reference.\n"
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--composition", type=Path, default=COMPOSITION_PATH)
    ap.add_argument("--reference", type=Path, default=DUPLICATE_REFERENCE)
    ap.add_argument("--design", type=Path, default=DESIGN_PATH)
    ap.add_argument("--write-reference", action="store_true",
                    help="write the reference from e3 instead of validating it")
    args = ap.parse_args()

    design = load_collapsed_design(args.design)
    if args.write_reference:
        ev = identify_duplicate_genotypes(design, args.composition, verbose=True)
        with open(args.reference, "w") as fh:
            fh.write(HEADER)
            ev[EVIDENCE_COLS].to_csv(fh, sep="\t", index=False, float_format="%.10g")
        print(f"wrote {args.reference} ({len(ev)} rows) -- review the diff before committing")
        return

    fresh = validate_reference(design, args.reference, args.composition)
    dups = sorted(fresh.loc[fresh["is_duplicate"], "genotype"])
    print(f"same-side duplicate reference OK: {len(fresh)} candidate genotype(s) "
          f"recomputed from {args.composition} match {args.reference} "
          f"(duplicates: {', '.join(dups) if dups else 'none'})")


if __name__ == "__main__":
    main()
