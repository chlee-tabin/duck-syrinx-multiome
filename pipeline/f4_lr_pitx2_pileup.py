#!/usr/bin/env python3
"""Phase F Step 4: L/R pseudobulk peak calling + PITX2 locus pileup.

User reframe (2026-04-24): instead of calling peaks at (d2g_subcluster x
stage) and relying on "reproducible across all genotype replicates"
(SnapATAC2 2.9.0's overly strict default, documented in
docs/snapatac2_upstream_findings.md Finding 2, not distributed), pool cells within
``POOLS[stage]["inference"]`` (laterality-showing LPM populations:
ventral splanchnic + syringeal + etc.) by the snRNA-seq-inferred
L/R side and call peaks per ``(stage, side)`` bucket.

Each stage has only 2 replicates (L and R), so the reproducibility
trap doesn't bite.  The PITX2c promoter peak, hypothesized to be
L-biased, will appear in the L pseudobulk even if it was filtered out
of f2's per-genotype "all 13 replicates" intersection.

Design rationale
----------------
(1) HH13 / HH21 / HH25: d2l has firm (``direct`` or ``placement``)
    per-genotype L/R inference based on 4-metric PITX2 congruence.
    Cells in ``POOLS[stage]["inference"]`` with a firm call contribute
    to either the L or R pseudobulk for that stage.

(2) HH17/16: d2l inference is weak pre-TES (D/V not resolved, PITX2
    pan-expressed).  Handled downstream via paired-library test at
    the PITX2c coordinate discovered in (1), not here.

Output
------
  phase2/atac/lr_pileup/
    fragments/<group>.bed.gz           per-group fragment BEDs
    peaks/<group>_peaks.bed            MACS3 narrowPeak per group
    step_f4_pitx2_coverage.tsv         per-base depth at PITX2 +/- 50 kb
                                        (format mirrors d2i: pos, <group>,
                                        <group>_percell columns).  Under
                                        --compartment-split --per-genotype-anchors
                                        <group> = f"{anchor}_{compartment}_
                                        {genotype}_{side}" (the fig03 panel C
                                        per-embryo-half pileup; --out-subdir
                                        redirects the output dir).
    step_f4_cell_counts.tsv            (group, n_cells, n_peaks_pitx2_locus)

Run on O2
---------
  sbatch pipeline/run_f4_lr_pileup.sh
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import gzip
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import snapatac2 as snap

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"

LR_PILEUP_DIR = f"{ATAC_DIR}/lr_pileup"
FRAG_DIR = f"{LR_PILEUP_DIR}/fragments"
PEAKS_DIR = f"{LR_PILEUP_DIR}/peaks"
os.makedirs(FRAG_DIR, exist_ok=True)
os.makedirs(PEAKS_DIR, exist_ok=True)

D2L_TSV = f"{SCRATCH}/figures/data/d2l_lr_inference.tsv"
# fig03-canonical L/R: snRNA at HH13/21/25, chromatin-canonical at HH17/16
# (docs/chromatin_canonical_hh1716.md).  Same key schema as d2l
# (genotype, library_m, inferred_side, firmness).  Select via --lr-tsv.
MULTIOME_LR_TSV = f"{SCRATCH}/figures/data/lr_inference_multiome.tsv"
PITX2_BED = f"{SCRATCH}/references/pitx2_promoter_regions.bed"

# Import POOLS from the shared subcluster helper
sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.lr_firmness import resolve_sides, side_map  # noqa: E402

# ── Parameters ────────────────────────────────────────────────────────────
STAGES_FIRM = ["HH13", "HH21", "HH25"]
# Inferred-side firmness classes trusted (extend with --firmness-also).  `weak`
# is excluded by rule (logged); `same_side_dup` (lane-loading-record split tube:
# one embryo-half in two libraries) is always kept on its design side -- every
# group label below is library-free, so its two libraries pool into ONE half
# (figures/_shared/lr_firmness.py).  Unknown values raise.
FIRMNESS_ACCEPT = {"direct", "placement"}
MIN_CELLS = 100  # minimum cells per pseudobulk to be meaningful (overridable via --min-cells)

# d2g_stage → stage_anchor (matches d2l's stage→stage_anchor mapping).
# HH17 and HH16 are pre-TES so they collapse into the HH17/16 anchor;
# POOLS["HH17/16"] is the only HH17/16-relevant POOLS key.  This mapping
# applies when --include-stages adds any of {HH16, HH17, HH17/16}.
STAGE_TO_ANCHOR = {
    "HH13": "HH13",
    "HH16": "HH17/16",
    "HH17": "HH17/16",
    "HH17/16": "HH17/16",
    "HH21": "HH21",
    "HH25": "HH25",
}

# f1 builds the AnnDataSet with add_key="arc_library" (ARC naming), and
# does NOT explicitly propagate gex_library to dataset.obs.  d2l's
# library_m column uses a further-collapsed GEX naming (GEX03-1A and
# GEX03-1B both map to "GEX03-1"; PILOTARC4A/B to "PILOTGEX4").  This
# map takes an arc_library value straight to the d2l library_m value.
ARC_TO_D2L_LIBRARY = {
    "ARC03-1A":    "GEX03-1",
    "ARC03-1B":    "GEX03-1",
    "ARC03-2":     "GEX03-2",
    "ARC04":       "GEX04",
    "ARC05":       "GEX05",
    "ARC06":       "GEX06",
    "ARC07":       "GEX07",
    "ARC08":       "GEX08",
    "ARC0S-duck":  "GEX0S-duck",
    "PILOTARC3":   "PILOTGEX3",
    "PILOTARC4A":  "PILOTGEX4",
    "PILOTARC4B":  "PILOTGEX4",
}

# MACS3 parameters (match f2 for cross-comparability with plan §3.2)
PEAK_QVALUE = 0.05   # relaxed vs f2's 0.01 because we have only 2 groups per stage
PEAK_SHIFT = -100
PEAK_EXTSIZE = 200

# PITX2 locus — 50 kb flanks around gene body (NC_051775.1:3,961,386-3,976,254)
PITX2_CHROM = "NC_051775.1"
PITX2_LOCUS_START = 3_911_386   # locus_50kb|PITX2 from pitx2_promoter_regions.bed
PITX2_LOCUS_END = 4_026_254


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sanitize(name: str) -> str:
    return name.replace("/", "_").replace("|", "__").replace(" ", "_")


# ──────────────────────────────────────────────────────────────────────────
# d2l lookup
# ──────────────────────────────────────────────────────────────────────────
def load_lr_map(tsv_path: str) -> dict:
    """(genotype, library_m) -> (side, sample library) for the usable sidecar rows.

    The sidecar has one row per (genotype, library) pair: an embryo's L and
    R halves typically land in different libraries so the side depends
    on WHICH library the cell came from, not just the genotype.  Rows in
    FIRMNESS_ACCEPT are used; same-side duplicates are kept on their design
    side (shared rule); unknown firmness / accepted rows without an L/R side
    raise (figures/_shared/lr_firmness.py).
    """
    df = pd.read_csv(tsv_path, sep="\t")
    return side_map(resolve_sides(df, FIRMNESS_ACCEPT, source=os.path.basename(tsv_path),
                                  log=ts))


# ──────────────────────────────────────────────────────────────────────────
# stage_side assignment
# ──────────────────────────────────────────────────────────────────────────
def build_stage_side(
    dataset, lr_map, include_stages: list[str],
    per_genotype_anchors: set[str],
    compartment_split: bool = False,
) -> list[str | None]:
    """Assign group label to each cell; None if excluded.

    Include a cell iff:
      - d2g_stage in include_stages
      - unified_genotype_id + arc_library both non-null (coupled with
        vireo genotype assignment)
      - d2g_subcluster in POOLS[anchor]["inference"]  (or, when
        compartment_split, also POOLS[anchor]["bilateral"])
      - L/R sidecar lookup (--lr-tsv; multiome by default) yields a firm L or R call

    Group label format depends on per_genotype_anchors / compartment_split:
      - default: f"{anchor}_{side}"  (e.g., "HH13_L", "HH17/16_R")
      - per-genotype anchors: f"{anchor}_{genotype}_{side}"  (one
        track per (genotype, side) at that anchor)
      - compartment_split: f"{anchor}_{compartment}_{side}" where
        compartment = ventral (POOLS inference) | dorsal (POOLS bilateral).
      - compartment_split AND per-genotype anchor: f"{anchor}_{compartment}_
        {genotype}_{side}" — the per-embryo-half granularity the fig03 panel C
        pileup rectification needs (one track per (compartment, genotype, side)
        so each embryo-half can be weighted equally).
      The side is a per-(genotype, library) physical L/R-half call, applied to
      whichever cells (dorsal or ventral) came from that half — so dorsal cells
      inherit the same side as ventral cells from the same library half.

    The d2g_stage value is mapped to its anchor via STAGE_TO_ANCHOR
    before the POOLS lookup and the group label, so HH17 and HH16 cells
    fold into the HH17/16 anchor (POOLS only has the HH17/16 key).
    """
    # Polars-backed obs: list() to extract as Python values.
    # f1's AnnDataSet uses add_key="arc_library" and does NOT propagate
    # gex_library to the merged obs — translate arc_library to d2l's
    # library_m convention via ARC_TO_D2L_LIBRARY.
    genotype = list(dataset.obs["unified_genotype_id"])
    arc_lib = list(dataset.obs["arc_library"])
    stage = list(dataset.obs["d2g_stage"])
    subcluster = list(dataset.obs["d2g_subcluster"])
    n = len(stage)
    assert (
        len(genotype) == len(arc_lib) == len(subcluster) == n
    ), "obs column length mismatch"

    out: list[str | None] = []
    drops: Counter = Counter()
    unknown_arc: Counter = Counter()
    for g, alib, s, sc in zip(genotype, arc_lib, stage, subcluster):
        if s not in include_stages:
            drops["stage_not_included"] += 1
            out.append(None)
            continue
        if g is None or alib is None:
            drops["uncoupled_or_unassigned"] += 1
            out.append(None)
            continue
        d2l_lib = ARC_TO_D2L_LIBRARY.get(str(alib))
        if d2l_lib is None:
            unknown_arc[str(alib)] += 1
            drops["unknown_arc_library"] += 1
            out.append(None)
            continue
        anchor = STAGE_TO_ANCHOR[s]
        if compartment_split:
            inf_pool = POOLS[anchor]["inference"]
            bil_pool = POOLS[anchor]["bilateral"]
            if sc in inf_pool:
                compartment = "ventral"
            elif sc in bil_pool:
                compartment = "dorsal"
            else:
                drops["not_in_inference_or_bilateral_pool"] += 1
                out.append(None)
                continue
        else:
            pool_defs = POOLS[anchor]["inference"]
            if sc not in pool_defs:
                drops["not_in_inference_pool"] += 1
                out.append(None)
                continue
        lr = lr_map.get((str(g), d2l_lib))
        if lr is None:
            # no usable side: excluded by the logged firmness rule (load_lr_map)
            # or no sidecar row
            drops["no_usable_side"] += 1
            out.append(None)
            continue
        side, _sample_lib = lr   # group labels are library-free: a duplicate pools
        # snapatac2's export_fragments uses the group label verbatim as
        # the output filename — "/" in HH17/16 trips "invalid filename".
        # Sanitize to HH17_16 here; the anchor-aware render scripts
        # detect it via prefix match against the known anchors.
        anchor_safe = anchor.replace("/", "_")
        if compartment_split and anchor in per_genotype_anchors:
            # Rectification (2026-06-17): per-(compartment x embryo-half) tracks
            # so the fig03 pileup can weight each embryo-HALF equally (mean over
            # halves), not each CELL (sum over cells -> high-cell embryos dominate
            # the shape). Mirrors the g4f per-embryo-half stats rectification.
            out.append(f"{anchor_safe}_{compartment}_{g}_{side}")
        elif compartment_split:
            out.append(f"{anchor_safe}_{compartment}_{side}")
        elif anchor in per_genotype_anchors:
            out.append(f"{anchor_safe}_{g}_{side}")
        else:
            out.append(f"{anchor_safe}_{side}")
    ts(f"  Drop reasons: {dict(drops)}")
    if unknown_arc:
        ts(f"  WARN unknown arc_library values (not in ARC_TO_D2L_LIBRARY): "
           f"{dict(unknown_arc)}")
    return out


# ──────────────────────────────────────────────────────────────────────────
# Per-base PITX2 locus coverage from a fragment BED
# ──────────────────────────────────────────────────────────────────────────
def pitx2_coverage_from_bed(bed_gz: Path, locus_size: int) -> np.ndarray:
    """Integer coverage vector over PITX2 locus (length = locus_size)."""
    cov = np.zeros(locus_size, dtype=np.int64)
    opener = gzip.open if str(bed_gz).endswith(".gz") else open
    with opener(bed_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] != PITX2_CHROM:
                continue
            frag_start = int(parts[1])
            frag_end = int(parts[2])
            if frag_end <= PITX2_LOCUS_START or frag_start >= PITX2_LOCUS_END:
                continue
            s = max(frag_start, PITX2_LOCUS_START) - PITX2_LOCUS_START
            e = min(frag_end, PITX2_LOCUS_END) - PITX2_LOCUS_START
            cov[s:e] += 1
    return cov


# ──────────────────────────────────────────────────────────────────────────
# External MACS3 (same pattern as f2b_gkmqc_rescue.py::run_macs3)
# ──────────────────────────────────────────────────────────────────────────
def run_macs3(bed_gz: Path, out_prefix: Path) -> Path:
    """Run MACS3 callpeak on a fragment BED.gz; return the narrowPeak path."""
    # gsize: use total size of main chromosomes from the .fai
    # Canonical /n/groups reference (cdw) — never purged, reachable from
    # compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror
    # (do not reintroduce the scratch reference path).
    fai = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa.fai"
    gsize = 0
    with open(fai) as fh:
        for line in fh:
            gsize += int(line.split("\t")[1])
    args = [
        "macs3", "callpeak",
        "-t", str(bed_gz),
        "-f", "BED",
        "-g", str(gsize),
        "--nomodel",
        "--shift", str(PEAK_SHIFT),
        "--extsize", str(PEAK_EXTSIZE),
        "-q", str(PEAK_QVALUE),
        "--call-summits",
        "-n", out_prefix.name,
        "--outdir", str(out_prefix.parent),
    ]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        ts(f"    MACS3 stderr: {result.stderr[-500:]}")
        raise RuntimeError(f"MACS3 failed for {bed_gz}")
    np_path = out_prefix.parent / f"{out_prefix.name}_peaks.narrowPeak"
    return np_path


def count_peaks_in_pitx2(narrowpeak: Path) -> int:
    if not narrowpeak.exists():
        return 0
    n = 0
    with open(narrowpeak) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != PITX2_CHROM:
                continue
            start = int(parts[1])
            end = int(parts[2])
            if end <= PITX2_LOCUS_START or start >= PITX2_LOCUS_END:
                continue
            n += 1
    return n


# ──────────────────────────────────────────────────────────────────────────
def main() -> None:
    global MIN_CELLS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-peaks", action="store_true",
                        help="Skip MACS3 peak calling (coverage-only).")
    parser.add_argument("--firmness-also", default="",
                        help="Comma-separated firmness values to accept in "
                             "addition to direct/placement (e.g. singleton_rank)")
    parser.add_argument("--include-stages", default="HH13,HH21,HH25",
                        help="Comma-separated d2g_stage values to include. "
                             "Default keeps backward compat with HH13/HH21/HH25; "
                             "for HH17/16 work pass HH13,HH21,HH25,HH17,HH16,HH17/16 "
                             "(HH17 + HH16 fold into the HH17/16 anchor).")
    parser.add_argument("--per-genotype-anchors", default="",
                        help="Comma-separated stage anchors (HH13, HH17/16, "
                             "HH21, HH25) where the groupby should be "
                             "(anchor, genotype, side) instead of (anchor, side). "
                             "When set, output goes to lr_pileup_per_genotype/ "
                             "subdirectory so the main lr_pileup/ stays intact.")
    parser.add_argument("--out-subdir", default="",
                        help="Override output subdirectory under phase2/atac/.  "
                             "Default = lr_pileup (or lr_pileup_per_genotype "
                             "when --per-genotype-anchors is set).")
    parser.add_argument("--compartment-split", action="store_true",
                        help="Bucket by (anchor, compartment, side) where "
                             "compartment = ventral (POOLS inference) | dorsal "
                             "(POOLS bilateral).  For the fig03 panel C "
                             "ventral-vs-dorsal L/R mirror pileup "
                             "(run_f4_perembryo.sh -> panel_C_pileup.R).  Pass "
                             "--min-cells to match the g4f floor (20).")
    parser.add_argument("--lr-tsv", default="multiome", choices=["d2l", "multiome"],
                        help="L/R side source: 'multiome' (default; lr_inference_multiome.tsv, "
                             "the fig03/d2n-canonical merger = snRNA at HH13/21/25 + "
                             "chromatin-canonical at HH17/16 — matches g1's default axis) or "
                             "'d2l' (legacy snRNA-only d2l_lr_inference.tsv; differs at HH17/16). "
                             "Default is multiome so a forgotten flag can't silently pick the "
                             "wrong axis — fig03 panel C always matches fig03's L/R sides.")
    parser.add_argument("--min-cells", type=int, default=None,
                        help=f"Minimum cells per pseudobulk to be retained.  "
                             f"Default depends on mode: {MIN_CELLS} for the "
                             f"L/R aggregate pileup (HH13/HH21/HH25), 30 for "
                             f"per-genotype HH17/16 anchors (pre-TES embryos "
                             f"have low yield; 30 recovers ~3x more paired "
                             f"tracks vs 100 with acceptable per-base noise).")
    args = parser.parse_args()

    # Parse anchor/stage sets first — the MIN_CELLS smart default branches on
    # per_genotype_anchors, so it must be bound before that block (was defined
    # after it → NameError whenever --min-cells was omitted).
    include_stages = [s.strip() for s in args.include_stages.split(",") if s.strip()]
    per_genotype_anchors = {s.strip() for s in args.per_genotype_anchors.split(",")
                            if s.strip()}
    ts(f"include_stages = {include_stages}")
    ts(f"per_genotype_anchors = {sorted(per_genotype_anchors) or '(none)'}")

    if args.min_cells is None:
        # Smart default: per-genotype HH17/16 anchor mode uses MIN_CELLS=30
        # (pre-TES embryos have low yield, validated 2026-04-27 against the
        # g52 case study) -- canonical mode keeps MIN_CELLS=100.
        MIN_CELLS = 30 if per_genotype_anchors else 100
    else:
        MIN_CELLS = args.min_cells
    ts(f"MIN_CELLS = {MIN_CELLS}")

    global FIRMNESS_ACCEPT
    if args.firmness_also:
        extra = {s.strip() for s in args.firmness_also.split(",") if s.strip()}
        FIRMNESS_ACCEPT = FIRMNESS_ACCEPT | extra   # validated in resolve_sides (raises)
    ts(f"FIRMNESS_ACCEPT = {sorted(FIRMNESS_ACCEPT)}")

    # Output subdirectory: per-genotype runs go to a sibling dir so the
    # canonical L/R-only lr_pileup/ outputs aren't clobbered.
    if args.out_subdir:
        subdir_name = args.out_subdir
    elif per_genotype_anchors:
        subdir_name = "lr_pileup_per_genotype"
    else:
        subdir_name = "lr_pileup"
    global LR_PILEUP_DIR, FRAG_DIR, PEAKS_DIR
    LR_PILEUP_DIR = f"{ATAC_DIR}/{subdir_name}"
    FRAG_DIR = f"{LR_PILEUP_DIR}/fragments"
    PEAKS_DIR = f"{LR_PILEUP_DIR}/peaks"
    os.makedirs(FRAG_DIR, exist_ok=True)
    os.makedirs(PEAKS_DIR, exist_ok=True)
    ts(f"output dir = {LR_PILEUP_DIR}")

    ts("=== Phase F Step 4: L/R pseudobulk peak calling + PITX2 pileup ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")

    lr_tsv = MULTIOME_LR_TSV if args.lr_tsv == "multiome" else D2L_TSV
    if not os.path.isfile(F1_DATASET):
        sys.exit(f"ERROR: f1 AnnDataSet missing at {F1_DATASET}")
    if not os.path.isfile(lr_tsv):
        sys.exit(
            f"ERROR: L/R TSV missing at {lr_tsv}.  Copy the repository "
            f"figures/data/ into the scratch working tree."
        )

    ts("Loading f1 dataset...")
    dataset = snap.read_dataset(F1_DATASET)
    ts(f"  {dataset.n_obs:,} cells x {dataset.n_vars:,} features")

    ts(f"Loading L/R map ({args.lr_tsv}: {os.path.basename(lr_tsv)})...")
    lr_map = load_lr_map(lr_tsv)
    ts(f"  {len(lr_map):,} (genotype, library) entries")

    ts("Assigning stage_side labels...")
    stage_side = build_stage_side(
        dataset, lr_map,
        include_stages=include_stages,
        per_genotype_anchors=per_genotype_anchors,
        compartment_split=args.compartment_split,
    )
    counts = Counter(s for s in stage_side if s is not None)
    ts(f"  Assigned: {sum(counts.values()):,} / {len(stage_side):,} cells")
    for grp in sorted(counts.keys()):
        marker = "" if counts[grp] >= MIN_CELLS else " (BELOW MIN)"
        ts(f"    {grp}: {counts[grp]:,} cells{marker}")

    valid_groups = sorted(g for g in counts if counts[g] >= MIN_CELLS)
    if not valid_groups:
        sys.exit(f"ERROR: no group meets MIN_CELLS={MIN_CELLS}")

    # §8 Polars quirk + SnapATAC2 export_fragments can't accept None
    # in group_by (Rust side: 'NoneType' object cannot be cast as 'str').
    # Replace None with a sentinel string that isn't in valid_groups; the
    # subsequent selections=valid_groups filter excludes it in Rust, so
    # no __unused__ file is written.  Same pattern as f2's
    # replicate_for_peaks fix.
    SENTINEL = "__skip__"
    stage_side_safe = [SENTINEL if v is None else v for v in stage_side]
    ts(f"Attaching stage_side column to dataset.obs (sentinel for "
       f"unassigned: {sum(1 for v in stage_side_safe if v == SENTINEL):,})")
    dataset.obs["stage_side"] = pl.Series(
        "stage_side", stage_side_safe, dtype=pl.Utf8
    )

    # ── Fragment export ──────────────────────────────────────────────────
    ts("")
    ts(f"--- Exporting fragment BEDs for {len(valid_groups)} groups ---")
    # Stale BED cleanup first (idempotency across retries)
    for grp in valid_groups:
        stale = Path(FRAG_DIR) / f"{sanitize(grp)}.bed.gz"
        if stale.exists():
            stale.unlink()
    snap.ex.export_fragments(
        dataset,
        groupby="stage_side",
        selections=valid_groups,
        out_dir=FRAG_DIR,
        suffix=".bed.gz",
        compression="gzip",
    )
    ts("  fragment export done")

    # SnapATAC2 writes files named "<prefix><group><suffix>"; prefix defaults
    # to "" so filenames are "<group>.bed.gz".  It also sanitizes group names
    # via snake_case — we need to find the actual output file per group.
    frag_files: dict[str, Path] = {}
    for grp in valid_groups:
        # Try the literal name and the sanitized name
        for cand in [f"{grp}.bed.gz", f"{sanitize(grp)}.bed.gz"]:
            p = Path(FRAG_DIR) / cand
            if p.exists():
                frag_files[grp] = p
                break
        if grp not in frag_files:
            # Fall back: find any .bed.gz that contains the group's name (safer
            # against SnapATAC2's internal filename sanitization)
            matches = list(Path(FRAG_DIR).glob("*.bed.gz"))
            ts(f"  WARN: no exact match for {grp!r} in {FRAG_DIR}; "
               f"available: {[p.name for p in matches]}")
            sys.exit(f"ERROR: fragment BED for {grp!r} not found")
        ts(f"    {grp}: {frag_files[grp].name} "
           f"({frag_files[grp].stat().st_size/1e6:.1f} MB)")

    # ── Per-group PITX2 locus coverage ───────────────────────────────────
    ts("")
    ts(f"--- Computing per-base PITX2 +/- 50 kb coverage ({len(valid_groups)} groups) ---")
    locus_size = PITX2_LOCUS_END - PITX2_LOCUS_START
    cov_by_group: dict[str, np.ndarray] = {}
    for grp in valid_groups:
        cov = pitx2_coverage_from_bed(frag_files[grp], locus_size)
        cov_by_group[grp] = cov
        ts(f"    {grp}: sum={cov.sum():,} max={cov.max():,} "
           f"density={cov.sum() / locus_size:.2f}/bp")

    # Build d2i-style TSV: chrom, pos, <group>, <group>_percell columns
    positions = np.arange(PITX2_LOCUS_START, PITX2_LOCUS_END, dtype=np.int64)
    out_df = pd.DataFrame({"chrom": PITX2_CHROM, "pos": positions})
    for grp in valid_groups:
        out_df[grp] = cov_by_group[grp]
        out_df[f"{grp}_percell"] = cov_by_group[grp] / counts[grp]
    out_tsv = Path(LR_PILEUP_DIR) / "step_f4_pitx2_coverage.tsv.gz"
    out_df.to_csv(out_tsv, sep="\t", index=False, compression="gzip")
    ts(f"  Saved: {out_tsv}  ({out_tsv.stat().st_size/1e6:.1f} MB)")

    # ── Peak calling ─────────────────────────────────────────────────────
    pitx2_peaks_per_group: dict[str, int] = {}
    total_peaks_per_group: dict[str, int] = {}
    if args.skip_peaks:
        ts("")
        ts("--- Skipping peak calling (--skip-peaks) ---")
    else:
        ts("")
        ts(f"--- Calling peaks via MACS3 subprocess ({len(valid_groups)} groups) ---")
        for grp in valid_groups:
            bed = frag_files[grp]
            out_prefix = Path(PEAKS_DIR) / sanitize(grp)
            np_path = run_macs3(bed, out_prefix)   # a MACS3 failure raises (no -1 sentinel)
            total_peaks_per_group[grp] = sum(
                1 for _ in open(np_path) if _.strip()
            )
            pitx2_peaks_per_group[grp] = count_peaks_in_pitx2(np_path)
            ts(f"    {grp}: {total_peaks_per_group[grp]:,} peaks total, "
               f"{pitx2_peaks_per_group[grp]} in PITX2 locus → "
               f"{np_path.name}")

    # ── Provenance ───────────────────────────────────────────────────────
    prov = pd.DataFrame({
        "group": valid_groups,
        "n_cells": [counts[g] for g in valid_groups],
        "frag_bed": [str(frag_files[g]) for g in valid_groups],
        "n_peaks_total": [total_peaks_per_group.get(g, np.nan) for g in valid_groups],
        "n_peaks_pitx2_locus": [
            pitx2_peaks_per_group.get(g, np.nan) for g in valid_groups
        ],
        "coverage_sum_pitx2_locus": [
            int(cov_by_group[g].sum()) for g in valid_groups
        ],
        "coverage_sum_per_cell": [
            cov_by_group[g].sum() / counts[g] for g in valid_groups
        ],
    })
    prov_tsv = Path(LR_PILEUP_DIR) / "step_f4_cell_counts.tsv"
    prov.to_csv(prov_tsv, sep="\t", index=False)
    ts(f"  Provenance: {prov_tsv}")

    ts("")
    ts("=== Summary ===")
    ts(f"  Groups: {len(valid_groups)}  (include_stages={include_stages}, "
       f"per_genotype_anchors={sorted(per_genotype_anchors) or '(none)'})")
    ts(f"  Total cells in pileup: {sum(counts[g] for g in valid_groups):,}")
    ts(f"  Per-base PITX2 coverage TSV: {out_tsv.name}")
    if not args.skip_peaks:
        ts("  Peaks in PITX2 locus (per group):")
        for g in valid_groups:
            ts(f"    {g}: {pitx2_peaks_per_group.get(g, 'n/a')} / "
               f"{total_peaks_per_group.get(g, 'n/a'):,} total")
    ts("=== Done ===")


if __name__ == "__main__":
    main()
