#!/usr/bin/env python3
"""Phase F Step 2: Peak calling per (d2g_subcluster × stage) via MACS3.

Plan §3.2 SnapATAC2 track + §3.3 granularity.  Calls peaks per
``peak_pool`` — the composite label set in f1 as
``f'{d2g_bucket}|{d2g_subcluster}'`` for LPM-lineage cells and
``f'{d2g_bucket}|{negative_control_cell_type}'`` for paraxial and
epithelial cells used as negative controls.  SnapATAC2's native
``snap.tl.macs3`` handles fragment export and delegates to MACS3;
when ``replicate=unified_genotype_id`` is supplied, peaks are filtered
to those reproducible across replicate pseudobulks within a group
(plan §3.2 reproducibility requirement).

Pool size floors (plan §3.3 and §3.4):
  * ``MIN_CELLS_PEAKS = 500``.  Pools at or above the floor are called
    here.  Pools below get emitted to ``step_f2_pools_for_gkmqc.tsv``
    for the f2b gkmQC rescue step (§3.4; not implemented by this file).
  * Pharyngeal endoderm at HH21/HH25 is explicitly skipped (§3.3, gut
    L/R confound).

Peak width: fixed 501 bp (``half_width=250`` in ``snap.tl.merge_peaks``)
per plan §3.2, inherited from Corces et al. 2018 Science TCGA pan-
cancer ATAC peak convention.

Input
-----
  phase2/atac/step_f1_combined.h5ads           — f1 AnnDataSet
  /n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/fasta/genome.fa.fai
  references/pitx2_promoter_regions.bed        — PITX2 spot-check (optional)

Output
------
  phase2/atac/peaks/<pool>_peaks.bed           — per-pool MACS3 calls
  phase2/atac/peaks/union_peaks.bed            — merged 501 bp union
  phase2/atac/step_f2_peak_matrix.h5ad         — cells × peaks matrix
  phase2/atac/step_f2_pools_called.tsv         — (peak_pool, n_cells, n_peaks)
  phase2/atac/step_f2_pools_for_gkmqc.tsv      — (peak_pool, n_cells) sub-500

Run on O2
---------
  sbatch --mem=48G -t 0-6:00 -c 4 -p short pipeline/run_phase_fg.sh --from f2 --to f2
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import snapatac2 as snap

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
PEAK_DIR = f"{ATAC_DIR}/peaks"
os.makedirs(PEAK_DIR, exist_ok=True)

F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
# Canonical /n/groups reference (cdw) — never purged, reachable from
# compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror.
GENOME_FAI = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa.fai"
PITX2_BED = f"{SCRATCH}/references/pitx2_promoter_regions.bed"

# ── Plan §3.2 / §3.3 parameters ──────────────────────────────────────────
MIN_CELLS_PEAKS = 500          # pool floor (below → gkmQC rescue in §3.4)
PEAK_QVALUE = 0.01             # MACS3 peak q-value
PEAK_SHIFT = -100              # MACS3 shift (Tn5 adjustment)
PEAK_EXTSIZE = 200             # MACS3 extension size
PEAK_HALF_WIDTH = 250          # merge_peaks half-width → 501 bp fixed width
REPLICATE_COL = "unified_genotype_id"  # biological replicate dimension
GROUP_COL = "peak_pool"

# Pools to skip even if >=MIN_CELLS_PEAKS (plan §3.3 gut L/R confound)
SKIP_POOLS = {
    "HH21|Pharyngeal endoderm",
    "HH25|Pharyngeal endoderm",
}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_chrom_sizes(fai_path: str) -> dict[str, int]:
    """Parse genome.fa.fai into {chrom: size} for SnapATAC2."""
    sizes: dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def sanitize_pool_name(pool: str) -> str:
    """Make a peak_pool label safe as a filename component."""
    out = pool.replace("/", "_").replace("|", "__").replace(" ", "_")
    return "".join(c if c.isalnum() or c in "._-+" else "_" for c in out)


def read_pitx2_region() -> tuple[str, int, int] | None:
    """Return (chrom, start, end) of the PITX2 ±50 kb locus from the BED.

    Looks for the ``locus_50kb|PITX2`` entry committed in
    references/pitx2_promoter_regions.bed.  Returns None if the file
    or entry is missing (spot-check becomes a warning, not an error).
    """
    if not os.path.isfile(PITX2_BED):
        ts(f"  WARN: {PITX2_BED} not found; PITX2 spot-check skipped")
        return None
    with open(PITX2_BED) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            if fields[3].startswith("locus_50kb|"):
                return fields[0], int(fields[1]), int(fields[2])
    ts("  WARN: no locus_50kb|PITX2 entry in BED; PITX2 spot-check skipped")
    return None


def read_anndataset_obs(dataset, cols: list[str]) -> pd.DataFrame:
    """Extract requested columns from a SnapATAC2 AnnDataSet.obs.

    §8.1 quirk: pd.DataFrame(dataset.obs) raises "DataFrame constructor
    not properly called!" on Polars-backed obs.  Pull column-by-column
    using the iterable protocol like f1's snap_obs_to_pandas helper.
    Unavailable columns are silently skipped (caller checks).
    """
    data = {}
    for col in cols:
        try:
            data[col] = list(dataset.obs[col])
        except Exception:  # noqa: BLE001
            continue
    return pd.DataFrame(data)


def _macs3_with_selections_fix(
    adata,
    *,
    groupby,
    qvalue=0.05,
    replicate=None,
    replicate_qvalue=None,
    max_frag_size=None,
    selections=None,
    shift=-100,
    extsize=200,
    min_len=None,
    blacklist=None,
    tempdir=None,
    n_jobs=8,
):
    """Fork of ``snap.tl.macs3`` (SnapATAC2 2.9.0) with a ``selections`` bug fix.

    Upstream bug (``snapatac2/tools/_call_peaks.py`` lines 161-170 at
    v2.9.0): ``groupby`` is translated to string indices via
    ``group_idx = {g: str(i) for i, g in enumerate(group_names)}`` before
    being passed to ``_snapatac2.export_tags``, but ``selections`` is
    passed through untranslated.  The Rust-side filter in
    ``src/call_peaks.rs:474`` then executes
    ``unique_keys.retain(|x| select.contains(x.0))`` where ``x.0`` is a
    string index like ``"0"`` / ``"1"`` and ``select`` still holds the
    original names like ``"HH21|Splanchnic"`` — the intersection is
    empty, zero fragment files are written, ``export_tags`` returns an
    empty HashMap, and ``snap.tl.macs3`` returns an empty peaks dict
    after wasting ~22 min iterating every fragment unfiltered.

    Symptom observed in job 37823135: ``0it [00:00, ?it/s]`` iteration
    immediately after a 22-minute "Exporting fragments..." phase.

    This fork is a line-for-line copy of ``snap.tl.macs3`` restricted to
    the call path we use (``inplace=False``, narrow peaks, no
    broad-peak or bulk path) plus a single corrective line that
    translates ``selections`` through ``group_idx`` before passing to
    ``_snapatac2.export_tags``.  Upstream PR planned once the fix is
    validated end-to-end.
    """
    import logging
    import tempfile
    from math import log

    import snapatac2._snapatac2 as _snapatac2
    from MACS3.Signal.PeakDetect import PeakDetect
    from snapatac2.tools._call_peaks import _par_map

    if isinstance(groupby, str):
        groupby = list(adata.obs[groupby])
    if replicate is not None and isinstance(replicate, str):
        replicate = list(adata.obs[replicate])

    options = type("MACS3_OPT", (), {})()
    options.info = lambda _: None
    options.debug = lambda _: None
    options.warn = logging.warn
    options.name = "MACS3"
    options.bdg_treat = "t"
    options.bdg_control = "c"
    options.cutoff_analysis = False
    options.cutoff_analysis_file = "a"
    options.store_bdg = False
    options.do_SPMR = False
    options.trackline = False
    options.log_pvalue = None
    options.log_qvalue = log(qvalue, 10) * -1
    options.PE_MODE = False
    options.gsize = adata.uns["reference_sequences"][
        "reference_seq_length"
    ].sum()
    options.maxgap = 30
    options.minlen = extsize if min_len is None else min_len
    options.shift = shift
    options.nolambda = False
    options.smalllocal = 1000
    options.largelocal = 10000
    options.call_summits = True
    options.broad = False
    options.fecutoff = 1.0
    options.d = extsize
    options.scanwindow = 2 * options.d

    with tempfile.TemporaryDirectory(dir=tempdir) as tmpdirname:
        logging.info("Exporting fragments...")
        group_names = list(set(groupby))
        group_idx = {g: str(i) for i, g in enumerate(group_names)}

        # BUGFIX vs snap.tl.macs3 v2.9.0: translate selections from group
        # names to string indices so the Rust-side filter actually matches.
        if selections is not None:
            missing = [s for s in selections if s not in group_idx]
            if missing:
                raise ValueError(
                    f"selections referenced unknown group(s): {missing}"
                )
            selections_translated = {group_idx[s] for s in selections}
        else:
            selections_translated = None

        fragments = _snapatac2.export_tags(
            adata,
            tmpdirname,
            [group_idx[x] for x in groupby],
            replicate,
            max_frag_size,
            selections_translated,
        )

        def _call_peaks(tags):
            import tempfile as _tempfile

            _tempfile.tempdir = tmpdirname
            merged, reps = _snapatac2.create_fwtrack_obj(tags)
            options.log_qvalue = log(qvalue, 10) * -1
            logging.getLogger().setLevel(logging.CRITICAL + 1)
            peakdetect = PeakDetect(treat=merged, opt=options)
            peakdetect.call_peaks()
            peakdetect.peaks.filter_fc(fc_low=options.fecutoff)
            merged = peakdetect.peaks

            others = []
            if replicate_qvalue is not None:
                options.log_qvalue = log(replicate_qvalue, 10) * -1
            for x in reps:
                peakdetect = PeakDetect(treat=x, opt=options)
                peakdetect.call_peaks()
                peakdetect.peaks.filter_fc(fc_low=options.fecutoff)
                others.append(peakdetect.peaks)

            logging.getLogger().setLevel(logging.INFO)
            return _snapatac2.find_reproducible_peaks(merged, others, blacklist)

        logging.info("Calling peaks...")
        if n_jobs == 1:
            peaks = [_call_peaks(x) for x in fragments.values()]
        else:
            peaks = _par_map(
                _call_peaks, [(x,) for x in fragments.values()], n_jobs
            )
        return {
            group_names[int(k)]: v for k, v in zip(fragments.keys(), peaks)
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pools",
        default=None,
        help=(
            "Comma-separated list of peak_pool names to restrict to.  "
            "Pilot mode: outputs gain a _pilot suffix so they don't "
            "clobber a later full-cohort run."
        ),
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help="Override output suffix (default: _pilot when --pools given, empty otherwise).",
    )
    args = parser.parse_args()

    if args.suffix is not None:
        suffix = args.suffix
    elif args.pools is not None:
        suffix = "_pilot"
    else:
        suffix = ""
    pilot_pools = (
        set(p.strip() for p in args.pools.split(",") if p.strip())
        if args.pools else None
    )

    ts("=== Phase F Step 2: MACS3 peak calling per (d2g_subcluster × stage) ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")
    if pilot_pools is not None:
        ts(f"  PILOT MODE: restricting to {len(pilot_pools)} pool(s): "
           f"{sorted(pilot_pools)}")
        ts(f"  Output suffix: {suffix!r}")
    ts(f"  MIN_CELLS_PEAKS={MIN_CELLS_PEAKS} (pools below → gkmQC rescue, §3.4)")
    ts(f"  peak qvalue={PEAK_QVALUE}; fixed-width peaks="
       f"{2 * PEAK_HALF_WIDTH + 1} bp")

    if not os.path.isfile(F1_DATASET):
        sys.exit(f"ERROR: f1 AnnDataSet missing at {F1_DATASET}")

    ts("Loading AnnDataSet...")
    dataset = snap.read_dataset(F1_DATASET)
    ts(f"  {dataset.n_obs:,} cells × {dataset.n_vars:,} features")

    # §8.1 quirk: AnnDataSet.obs is a PyDataFrameElem that pandas cannot
    # consume directly; pull named columns into a plain DataFrame.
    obs_df = read_anndataset_obs(
        dataset,
        [GROUP_COL, REPLICATE_COL, "arc_library", "gex_library",
         "cell_type", "sex", "stage",
         "d2g_stage", "d2g_cluster", "d2g_subcluster"],
    )

    if GROUP_COL not in obs_df.columns:
        sys.exit(
            f"ERROR: {GROUP_COL!r} not in AnnDataSet.obs.  "
            f"f1 must be re-run with the peak_pool transfer to dataset.obs."
        )
    if REPLICATE_COL not in obs_df.columns:
        sys.exit(
            f"ERROR: {REPLICATE_COL!r} not in AnnDataSet.obs."
        )

    # ── Pool size distribution ───────────────────────────────────────────
    pool_counts = obs_df[GROUP_COL].value_counts(dropna=True)
    n_orphan = int(obs_df[GROUP_COL].isna().sum())
    ts(f"  {len(pool_counts)} distinct peak_pool values; "
       f"{n_orphan:,} cells without a peak_pool (excluded from calling)")

    callable_pools = [p for p in pool_counts.index
                      if pool_counts[p] >= MIN_CELLS_PEAKS]
    callable_pools = [p for p in callable_pools if p not in SKIP_POOLS]
    if pilot_pools is not None:
        unknown = pilot_pools - set(pool_counts.index)
        if unknown:
            sys.exit(
                f"ERROR: --pools argument referenced unknown pool(s): {sorted(unknown)}.  "
                f"Available: {sorted(pool_counts.index)[:10]}..."
            )
        callable_pools = [p for p in callable_pools if p in pilot_pools]
        ts(f"  Pilot filter: {len(callable_pools)} of requested {len(pilot_pools)} "
           f"pools pass the MIN_CELLS_PEAKS floor")
    rescue_pools = pool_counts[
        (pool_counts > 0) & (pool_counts < MIN_CELLS_PEAKS)
    ].sort_values(ascending=False)
    skipped_hits = sorted(SKIP_POOLS & set(pool_counts.index))

    ts(f"  Callable (≥{MIN_CELLS_PEAKS} cells, "
       f"{len(callable_pools)} pools):")
    for p in sorted(callable_pools):
        # For callable pools, also check replicate availability.
        sub = obs_df[obs_df[GROUP_COL] == p]
        n_rep = sub[REPLICATE_COL].nunique()
        flag = " !=2rep" if n_rep < 2 else ""
        ts(f"    {p}: {pool_counts[p]:,} cells, {n_rep} replicates{flag}")

    if len(rescue_pools) > 0:
        ts(f"  Sub-{MIN_CELLS_PEAKS} pools deferred to gkmQC rescue "
           f"({len(rescue_pools)} pools):")
        for p, n in rescue_pools.items():
            ts(f"    {p}: {n:,} cells")

    if skipped_hits:
        ts(f"  Pools explicitly skipped (plan §3.3): {skipped_hits}")

    if not callable_pools:
        sys.exit("ERROR: no pools passed the MIN_CELLS_PEAKS threshold")

    # ── Call peaks ───────────────────────────────────────────────────────
    # snap.tl.macs3 scans the entire replicate column even when
    # selections= is passed; uncoupled cells (None for unified_genotype_id)
    # would trigger
    #   TypeError: argument 'replicates': 'NoneType' object cannot be cast as 'str'
    # inside _snapatac2.export_tags.  Populate a sentinel-safe
    # replicate column on dataset.obs (None → "uncoupled") and groupby
    # restricts the actual reproducibility filter to coupled cells via
    # selections=.
    ts("")
    ts("--- Adding sentinel-safe replicate column for macs3 ---")
    raw_rep = list(dataset.obs[REPLICATE_COL])
    rep_safe = [v if v is not None else "uncoupled" for v in raw_rep]
    dataset.obs["replicate_for_peaks"] = rep_safe
    ts(f"  replicate_for_peaks: "
       f"{sum(1 for v in rep_safe if v != 'uncoupled'):,} coupled, "
       f"{sum(1 for v in rep_safe if v == 'uncoupled'):,} uncoupled→sentinel")

    ts("")
    ts(f"--- Calling peaks: groupby={GROUP_COL!r} "
       f"replicate='replicate_for_peaks' ---")
    ts("  (using _macs3_with_selections_fix — see docstring for the "
       "SnapATAC2 v2.9.0 upstream bug this works around)")
    peaks = _macs3_with_selections_fix(
        dataset,
        groupby=GROUP_COL,
        replicate="replicate_for_peaks",
        selections=set(callable_pools),
        qvalue=PEAK_QVALUE,
        shift=PEAK_SHIFT,
        extsize=PEAK_EXTSIZE,
        n_jobs=4,
    )

    if not peaks:
        sys.exit("ERROR: _macs3_with_selections_fix returned no peaks")

    # ── Per-pool BED exports ─────────────────────────────────────────────
    # Pilot outputs go under peaks_pilot/ so they don't clobber the
    # production peaks/ directory (and vice versa on a real run).
    pool_bed_dir = PEAK_DIR if not suffix else f"{PEAK_DIR}{suffix}"
    os.makedirs(pool_bed_dir, exist_ok=True)
    for pool_name, peak_df in peaks.items():
        safe = sanitize_pool_name(pool_name)
        peak_df.write_csv(
            os.path.join(pool_bed_dir, f"{safe}_peaks.bed"),
            separator="\t",
            include_header=False,
        )
        ts(f"  {pool_name}: {len(peak_df):,} peaks → {safe}_peaks.bed")

    # ── Union peak set at fixed 501 bp ──────────────────────────────────
    ts("")
    ts("--- Merging peaks into union set (fixed-width 501 bp) ---")
    chrom_sizes = parse_chrom_sizes(GENOME_FAI)
    merged_peaks = snap.tl.merge_peaks(
        peaks,
        chrom_sizes=chrom_sizes,
        half_width=PEAK_HALF_WIDTH,
    )
    n_union = len(merged_peaks)
    ts(f"  Union peaks: {n_union:,}")

    # SnapATAC2 2.9.0 merge_peaks returns a DataFrame with columns
    # ["Peaks", <group>, ...] where "Peaks" is a "chrom:start-end"
    # string (see src/call_peaks.rs::py_merge_peaks).  Parse it into
    # Chromosome/Start/End for BED output and downstream spot-checks.
    mp = merged_peaks.to_pandas() if hasattr(merged_peaks, "to_pandas") else merged_peaks
    parts = mp["Peaks"].str.extract(r"^([^:]+):(\d+)-(\d+)$")
    mp["Chromosome"] = parts[0]
    mp["Start"] = parts[1].astype(int)
    mp["End"] = parts[2].astype(int)

    union_bed = os.path.join(pool_bed_dir, "union_peaks.bed")
    mp[["Chromosome", "Start", "End"]].to_csv(
        union_bed, sep="\t", index=False, header=False,
    )
    ts(f"  Saved: {union_bed}")

    # ── Build cell × peak matrix ─────────────────────────────────────────
    ts("")
    ts("--- Building peak matrix ---")
    peak_list = mp["Peaks"].tolist()
    peak_matrix_path = f"{ATAC_DIR}/step_f2_peak_matrix{suffix}.h5ad"
    try:
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, use_rep=peak_list, file=peak_matrix_path,
        )
    except TypeError:
        # Older API variant
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, peak_file=union_bed, file=peak_matrix_path,
        )
    ts(f"  Peak matrix: {peak_matrix.shape[0]:,} cells × "
       f"{peak_matrix.shape[1]:,} peaks → {peak_matrix_path}")

    # Transfer metadata from source AnnDataSet to the peak matrix.
    # §8 Polars quirk: peak_matrix.obs is Polars-backed.  Writing a raw
    # numpy object array (from pandas .values) is interpreted as
    # ``Binary`` dtype and rejected with "Unsupported series data type:
    # Binary".  Wrap values in ``pl.Series`` with an explicit dtype,
    # inferred from the first non-null value (Int/Float vs. Utf8).
    import polars as pl
    for col in ["arc_library", "gex_library", "unified_genotype_id",
                "cell_type", "sex", "stage",
                "d2g_stage", "d2g_cluster", "d2g_subcluster", "peak_pool"]:
        if col not in obs_df.columns:
            continue
        vals = obs_df[col].tolist()
        non_null = next((v for v in vals if v is not None), None)
        if isinstance(non_null, (int, np.integer)) and not isinstance(non_null, bool):
            pl_dtype = pl.Int64
        elif isinstance(non_null, (float, np.floating)):
            pl_dtype = pl.Float64
        else:
            pl_dtype = pl.Utf8
        peak_matrix.obs[col] = pl.Series(col, vals, dtype=pl_dtype)

    # ── Provenance TSVs ──────────────────────────────────────────────────
    pools_called_df = pd.DataFrame({
        "peak_pool": list(peaks.keys()),
        "n_cells": [int(pool_counts.get(p, 0)) for p in peaks.keys()],
        "n_peaks_called": [len(peaks[p]) for p in peaks.keys()],
    }).sort_values("n_cells", ascending=False)
    pools_called_df.to_csv(
        f"{ATAC_DIR}/step_f2_pools_called{suffix}.tsv", sep="\t", index=False,
    )

    rescue_df = rescue_pools.rename_axis("peak_pool").reset_index(name="n_cells")
    rescue_df.to_csv(
        f"{ATAC_DIR}/step_f2_pools_for_gkmqc{suffix}.tsv", sep="\t", index=False,
    )
    ts(f"  Provenance TSVs: step_f2_pools_called{suffix}.tsv, "
       f"step_f2_pools_for_gkmqc{suffix}.tsv")

    # ── PITX2 locus spot-check (positive control) ────────────────────────
    ts("")
    ts("--- PITX2 locus spot-check (positive control) ---")
    region = read_pitx2_region()
    if region is not None:
        chrom, start, end = region
        hits = mp[
            (mp["Chromosome"] == chrom)
            & (mp["Start"] < end)
            & (mp["End"] > start)
        ]
        ts(f"  Union peaks overlapping PITX2 ±50 kb "
           f"({chrom}:{start:,}-{end:,}): {len(hits)}")
        for _, r in hits.head(5).iterrows():
            ts(f"    {r['Chromosome']}:{r['Start']}-{r['End']}")

    # ── Summary ──────────────────────────────────────────────────────────
    ts("")
    ts("=== Summary ===")
    ts(f"  Pools called: {len(peaks)}")
    ts(f"  Pools deferred to gkmQC rescue: {len(rescue_pools)}")
    ts(f"  Pools explicitly skipped: {len(skipped_hits)}")
    ts(f"  Union peaks (501 bp fixed): {n_union:,}")
    ts(f"  Peak matrix: {peak_matrix.shape}")
    ts("=== Done ===")


if __name__ == "__main__":
    main()
