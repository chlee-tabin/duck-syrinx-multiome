#!/usr/bin/env python3
"""Assert a planted-BRE sensitivity floor actually recovered its injected SMAD core.

A planted-BRE floor is the load-bearing positive control for ANY "0 SMAD" verdict in
the right-peak BMP-SMAD probe: it injects a clustered BMP BRE (GGCGCC GC-SBE + CAGAC
SBE, PLANT_COPIES x) into the positives, retrains the gkm-SVM (h9) and reads it out
(h3b).  The companion "no SMAD" arms (foc_peak / dyn_peak / dynlost_peak) are
interpretable as a TRUE biological negative ONLY IF the matching floor proves the
readout can detect a SMAD footprint when one is present.  If a floor trains but does
NOT recover the planted core, the readout is BLIND there and the "0 SMAD" arms are
method-blindness, not absence -> this exits nonzero (fail loud, no
silent fallbacks) so a broken floor can never read as a pass.

This is the single source of truth for the floor gate, called by
run_rightpeak_smad_probe.sh, run_rightpeak_dynamic_lost.sh, and
run_h9_sensitivity_floor.sh.  It was extracted from the inline gate that previously
lived ONLY in run_h9_sensitivity_floor.sh and is NEWLY applied to the two rightpeak
runners, which until now trained planted-BRE floor arms without ever asserting recovery.

Usage:
    python3 pipeline/_floor_gate.py LABEL OUT_DIR [LABEL2 OUT_DIR2 ...]

For each (LABEL, OUT_DIR) floor arm it reads, from that arm's own output dir:
    OUT_DIR/h9_summary_*.tsv                  -> cv_auc                 (need >= MIN_CV_AUC)
    OUT_DIR/h3b_summary_strict_genomicbg.tsv  -> motifs_with_SMAD15_BMP (need >= 1)
The h3b path is EXACT (not a glob): a glob could sort-pick a stale
h3b_summary_strict_gc.tsv ('gc' sorts before 'genomicbg') and validate the wrong run.

Exits 1 if ANY floor arm fails the gate; 0 only if ALL pass.  MIN_CV_AUC defaults to
0.90 and is overridable via FLOOR_MIN_CV_AUC (the floors are saturated positive
controls with cv_auc ~ 1.0, so 0.90 is a generous lower bound).
"""
import csv
import glob
import os
import sys

MIN_CV_AUC = float(os.environ.get("FLOOR_MIN_CV_AUC", "0.90"))


def _cell(path, name):
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return None
    return rows[0].get(name) if rows else None


def gate(out_dir):
    """Return a list of failure strings for one floor arm (empty list == PASS)."""
    h3b = os.path.join(out_dir, "h3b_summary_strict_genomicbg.tsv")
    # h9 summary via glob is safe: each floor arm has its own OUT_DIR holding exactly one
    # h9_summary_<TAG>.tsv (same idiom the runners' summary printers use).  h3b is read by
    # EXACT path instead — a glob could sort-pick a stale h3b_summary_strict_gc.tsv.
    h9_matches = sorted(glob.glob(os.path.join(out_dir, "h9_summary_*.tsv")))
    h9 = h9_matches[0] if h9_matches else None
    cv = _cell(h9, "cv_auc")
    smad = _cell(h3b, "motifs_with_SMAD15_BMP")
    fails = []
    try:
        if cv is None or float(cv) < MIN_CV_AUC:
            fails.append(f"cv_auc={cv} (need >= {MIN_CV_AUC})")
    except ValueError:
        fails.append(f"cv_auc unparseable: {cv!r}")
    try:
        if smad is None or int(smad) < 1:
            fails.append(f"motifs_with_SMAD15_BMP={smad} (need >= 1)")
    except ValueError:
        fails.append(f"SMAD core count unparseable: {smad!r}")
    return fails


def main(argv):
    if len(argv) < 2 or len(argv) % 2 != 0:
        sys.stderr.write("usage: _floor_gate.py LABEL OUT_DIR [LABEL2 OUT_DIR2 ...]\n")
        return 2
    pairs = list(zip(argv[0::2], argv[1::2]))
    any_fail = False
    for label, out_dir in pairs:
        fails = gate(out_dir)
        if fails:
            any_fail = True
            sys.stderr.write(f"FLOOR FAIL [{label}]: " + "; ".join(fails) + "\n")
        else:
            print(f"FLOOR PASS [{label}]")
    if any_fail:
        sys.stderr.write(
            "A planted-BRE floor did NOT recover the injected SMAD core -> the "
            "gkm-SVM/h3b readout is BLIND for that arm, so the companion '0 SMAD' "
            "arms are NOT interpretable. Do NOT cite the no-SMAD result.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
