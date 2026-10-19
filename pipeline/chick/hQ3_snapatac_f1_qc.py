#!/usr/bin/env python3
"""hQ3 — SnapATAC2 per-library f1 QC for the chick snATAC QC panel.

Chick analog of the duck f1 step, matching `f1_coupled_atac.py`'s METHOD exactly
(not just its gate thresholds).  For each staged chick ARC library:

  1. import the cellranger-arc fragments via SnapATAC2 (galGal6a) at
     ``min_num_fragments=100`` — the 100-999-fragment cells are kept so
     ``tsse()`` computes per-cell TSS enrichment on the SAME broad denominator
     as the duck (a ≥1000-only import would shift the tsse() background and make
     the chick TSSe values not numerically comparable to the duck's);
  2. run ``snap.metrics.tsse()``;
  3. apply the snATAC QC gate

         f1-passed  ==  TSS enrichment >= 2  AND  n_fragment >= 1000

This makes the chick per-cell TSSe values numerically comparable to the duck's
(for cross-species QC comparison), with the identical gate.  Writes:

  - figdata/chick_atac_f1_passed_barcodes.tsv  (arc_library, barcode)
        consumed by hQ2 to build the "f1" qc_set reservoir (E2-E6)
  - figdata/chick_atac_f1_qc_summary.tsv       (arc_library, n_imported,
        n_qc_pass, median_tsse, median_n_fragment)
        consumed by figS04_supp E1 (light bars = n_qc_pass)

Self-contained (no repo import); libraries discovered by globbing the staged
tree.  Runs in the ``snapatac2`` conda env (SnapATAC2 is not in vireo).
SnapATAC2 2.9 obs is Polars-backed — see ``snap_obs_to_pandas``
(feedback_snapatac2_polars_quirks).

  sbatch pipeline/chick/hQ3_snapatac_f1_qc.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import glob
import os
import sys
import time

import numpy as np
import pandas as pd

QC_BASE = os.environ.get("CHICK_QC_BASE",
                         f"{SCRATCH_USER}/chick_demux/atac_qc")
REF_DIR = os.environ.get("CHICK_REF_DIR",
                         f"{SCRATCH_USER}/chick_demux/ref")
FIGDIR = os.path.join(QC_BASE, "figdata")
os.makedirs(FIGDIR, exist_ok=True)

GENOME_FAI = os.path.join(REF_DIR, "genome.fa.fai")
CHROMS_TXT = os.path.join(REF_DIR, "chroms.txt")
GTF = os.path.join(REF_DIR, "genes.gtf.gz")

# f1 QC, matching the duck f1_coupled_atac.py method exactly:
#   import at min_num_fragments=100 (keeps 100-999-frag cells so tsse() computes
#   on the same broad denominator as the duck), run tsse(), THEN apply the gate
#   (TSS enrichment >= 2 AND n_fragment >= 1000).  This makes the chick per-cell
#   TSSe values numerically comparable to the duck's (for cross-species QC
#   comparison), not just the gate thresholds.
IMPORT_MIN_FRAGMENTS = 100   # duck-parity import floor (tsse() denominator)
F1_MIN_FRAGMENTS = 1000      # gate: kept cells need >= 1000 fragments
F1_MIN_TSSE = 2.0            # gate: kept cells need TSS enrichment >= 2


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def discover_libs():
    return [os.path.basename(os.path.dirname(o))
            for o in sorted(glob.glob(os.path.join(QC_BASE, "*", "outs")))]


def snap_obs_to_pandas(obs, obs_names=None, columns=("tsse", "n_fragment")):
    """Extract columns from a SnapATAC2 2.9 Polars-backed obs into pandas.
    pd.DataFrame(obs) raises and there is no .to_pandas(); columns ARE
    iterable, so iterate per name.  Missing columns are skipped."""
    data = {}
    for col in columns:
        try:
            data[col] = list(obs[col])
        except Exception as e:  # noqa: BLE001
            # A genuinely absent optional column is rare; the two required
            # columns (tsse, n_fragment) fail loud at the caller's
            # obs["..."].astype(float), so this only swallows truly-optional
            # ones — but log it so a real extraction failure is visible.
            ts(f"  NOTE: obs column {col!r} unavailable, skipped ({e})")
            continue
    idx = list(obs_names) if obs_names is not None else None
    return pd.DataFrame(data, index=idx)


def load_chrom_sizes():
    if not os.path.isfile(GENOME_FAI):
        sys.exit(f"ERROR: galGal6a .fai not found at {GENOME_FAI}")
    chrom_sizes = {}
    with open(GENOME_FAI) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    # Optional restriction to the canonical chrom set (drops unplaced scaffolds
    # so TSS/import stays on chr1-28/30-33 + Z/W).
    if os.path.isfile(CHROMS_TXT):
        with open(CHROMS_TXT) as fh:
            names = set(x for x in fh.read().replace("\n", ",").split(",") if x.strip())
        kept = {c: s for c, s in chrom_sizes.items() if c in names}
        if kept:
            ts(f"  restricted to {len(kept)} canonical chroms (chroms.txt)")
            chrom_sizes = kept
    ts(f"  {len(chrom_sizes)} chroms; e.g. {list(chrom_sizes.items())[:3]}")
    return chrom_sizes


def qc_library(name, chrom_sizes, snap):
    """Import one library's fragments, compute TSS + n_fragment, apply f1 gate.
    Returns (passed_barcodes list, summary dict)."""
    frag = os.path.join(QC_BASE, name, "outs", "atac_fragments.tsv.gz")
    if not os.path.isfile(frag):
        ts(f"  [{name}] fragments missing ({frag}) — skipping")
        return [], {"library": name, "n_imported": 0, "n_qc_pass": 0,
                    "median_tsse": float("nan"), "median_n_fragment": float("nan"),
                    "status": "missing_fragments"}

    h5ad = os.path.join(QC_BASE, name, "outs", "_hQ3_snap.h5ad")
    import_fn = (getattr(snap.pp, "import_fragments", None)
                 or getattr(snap.pp, "import_data", None))
    if import_fn is None:
        sys.exit("ERROR: SnapATAC2 exposes neither pp.import_fragments nor pp.import_data.")
    ts(f"  [{name}] importing {frag} (min_num_fragments={IMPORT_MIN_FRAGMENTS}, "
       f"duck-parity tsse() denominator)")
    try:
        adata = import_fn(frag, chrom_sizes=chrom_sizes, file=h5ad,
                          min_num_fragments=IMPORT_MIN_FRAGMENTS, sorted_by_barcode=False)
    except TypeError:
        adata = import_fn(frag, chrom_size=chrom_sizes, file=h5ad,
                          min_num_fragments=IMPORT_MIN_FRAGMENTS, sorted_by_barcode=False)
    n_imported = adata.n_obs
    ts(f"  [{name}] {n_imported:,} barcodes with >= {IMPORT_MIN_FRAGMENTS} fragments")

    try:
        try:
            snap.metrics.tsse(adata, gene_anno=GTF)
        except TypeError:
            snap.metrics.tsse(adata, GTF)
    except Exception as e:  # noqa: BLE001
        adata.close()
        ts(f"  [{name}] TSS enrichment FAILED ({e}) — cannot apply f1 gate")
        return [], {"library": name, "n_imported": n_imported, "n_qc_pass": 0,
                    "median_tsse": float("nan"), "median_n_fragment": float("nan"),
                    "status": "tsse_failed"}

    obs = snap_obs_to_pandas(adata.obs, obs_names=adata.obs_names)
    adata.close()
    tsse = obs["tsse"].astype(float)
    nfrag = obs["n_fragment"].astype(float)
    # Apply the gate AFTER tsse() (duck-parity): keep cells with TSS >= 2 AND
    # >= 1000 fragments, from the >= 100-fragment import population.
    keep = (tsse >= F1_MIN_TSSE) & (nfrag >= F1_MIN_FRAGMENTS)
    passed = obs.index[keep]
    # Medians reported over the PASSED cells (the kept population), matching the
    # duck summary (whose obs is post-mask) — NOT over the broad >= 100 import.
    tsse_pass = tsse[keep]
    nfrag_pass = nfrag[keep]
    summary = {
        "library":           name,
        "n_imported":        n_imported,   # >= 100-fragment raw import (duck n_raw analog)
        "n_qc_pass":         int(len(passed)),
        "median_tsse":       round(float(np.median(tsse_pass)), 3) if len(tsse_pass) else float("nan"),
        "median_n_fragment": round(float(np.median(nfrag_pass)), 1) if len(nfrag_pass) else float("nan"),
        "status":            "ok",
    }
    ts(f"  [{name}] f1-pass {summary['n_qc_pass']:,}/{n_imported:,} "
       f"(median TSSe={summary['median_tsse']}, median frag={summary['median_n_fragment']:,.0f})")
    # Clean the scratch h5ad (large, reproducible).
    try:
        os.remove(h5ad)
    except OSError:
        pass
    return list(passed), summary


def main():
    ts("=== hQ3: chick SnapATAC2 per-library f1 QC ===")
    import snapatac2 as snap
    ts(f"  SnapATAC2 version: {snap.__version__}")

    chrom_sizes = load_chrom_sizes()
    if not os.path.isfile(GTF):
        sys.exit(f"ERROR: gene annotation for TSS not found at {GTF}")
    ts(f"  gene annotation: {GTF}")

    libs = discover_libs()
    ts(f"  discovered {len(libs)} staged libraries: {libs}")

    bc_rows = []
    summaries = []
    for name in libs:
        passed, summary = qc_library(name, chrom_sizes, snap)
        summaries.append(summary)
        for bc in passed:
            bc_rows.append((name, bc))

    bc_out = os.path.join(FIGDIR, "chick_atac_f1_passed_barcodes.tsv")
    with open(bc_out, "w") as f:
        f.write("arc_library\tbarcode\n")
        for lib, bc in bc_rows:
            f.write(f"{lib}\t{bc}\n")
    ts(f"  Wrote {len(bc_rows):,} f1-passed barcodes -> {bc_out}")

    sum_out = os.path.join(FIGDIR, "chick_atac_f1_qc_summary.tsv")
    header = ["arc_library", "n_imported", "n_qc_pass",
              "median_tsse", "median_n_fragment", "status"]
    keys = ["library", "n_imported", "n_qc_pass",
            "median_tsse", "median_n_fragment", "status"]
    with open(sum_out, "w") as f:
        f.write("\t".join(header) + "\n")
        for s in summaries:
            f.write("\t".join(str(s[k]) for k in keys) + "\n")
    ts(f"  Wrote summary ({len(summaries)} libs) -> {sum_out}")
    ts("=== hQ3 complete ===")


if __name__ == "__main__":
    main()
