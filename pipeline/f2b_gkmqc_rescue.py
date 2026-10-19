#!/usr/bin/env python3
"""Phase F Step 2b: gkmQC rescue for sub-500 cell peak pools.

Plan §3.4 implementation.  Reads ``step_f2_pools_for_gkmqc.tsv`` from
f2 (pools that fell below the 500-cell MIN_CELLS_PEAKS floor) and
runs the gkmQC pipeline per pool:

    fragments → MACS3 at q=0.01 (orig)    gkmqc evaluate
             → MACS3 at p=0.005 (relaxed) gkmqc evaluate → optimize → <pool>.e300.optz.bed

Rescued peaks are then folded into the existing union set at 501 bp
fixed width and the peak matrix is rebuilt so downstream g1 → g4
pick them up transparently.

STATUS: PRE-STAGED (2026-04-23).  This driver has NOT been executed
end-to-end.  Three verification gates to resolve on first real run:

    (G1) gkmQC install on o2 — conda env ``gkmqc`` via the upstream
         README (https://github.com/Dongwon-Lee/gkmQC), C source build.
         Falls back to a skip-with-warning if the import fails.
    (G2) Duck ZJU1.0 null-seq index — one-shot build via
         ``pipeline/build_duck_gkmqc_index.sh``; writes to
         ``$SCRATCH/references/gkmqc_idx_zju1/``.
    (G3) SnapATAC2 ``export_fragments`` signature — 2.9 accepts
         ``groupby=`` and ``out_dir=`` but emits one BED.gz per group;
         verify on first real call and adjust the per-pool file naming.

References
----------
  Han et al. 2022 PNAS, PMID 36508674
  https://github.com/Dongwon-Lee/gkmQC

Input
-----
  phase2/atac/step_f1_combined.h5ads           f1 AnnDataSet (fragments)
  phase2/atac/step_f2_pools_for_gkmqc.tsv      pool list from f2
  phase2/atac/peaks/union_peaks.bed            f2 union set (pre-rescue)
  /n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/fasta/genome.fa.fai  chrom sizes
  references/gkmqc_idx_zju1/                   null-seq index (G2)

Output
------
  phase2/atac/peaks/rescue/<pool>/<pool>.orig.narrowPeak
  phase2/atac/peaks/rescue/<pool>/<pool>.relaxed.narrowPeak
  phase2/atac/peaks/rescue/<pool>/<pool>.e300.optz.bed
  phase2/atac/peaks/rescue/<pool>/<pool>_eval.tsv   per-peak gkmQC score
  phase2/atac/peaks/union_peaks_with_rescue.bed    post-rescue union
  phase2/atac/step_f2b_rescue_summary.tsv          per-pool tally
  phase2/atac/step_f2_peak_matrix.h5ad             REBUILT (in place)

Run on O2
---------
  sbatch --mem=32G -t 0-8:00 -c 10 -p short pipeline/run_phase_fg.sh --from f2b --to f2b
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import snapatac2 as snap

SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
PEAK_DIR = f"{ATAC_DIR}/peaks"
RESCUE_DIR = f"{PEAK_DIR}/rescue"
FRAG_DIR = f"{RESCUE_DIR}/_fragments"
os.makedirs(RESCUE_DIR, exist_ok=True)
os.makedirs(FRAG_DIR, exist_ok=True)

F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
UNION_BED = f"{PEAK_DIR}/union_peaks.bed"
POOLS_TSV = f"{ATAC_DIR}/step_f2_pools_for_gkmqc.tsv"
# Canonical /n/groups reference (cdw) — never purged, reachable from
# compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror.
GENOME_FAI = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa.fai"
GENOME_FA = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa"
# gkmQC's own `data/<genome>/` convention is where `gkmqc.py evaluate
# -g zju1` resolves the null-seq index; buildidx wrote it there.
GKMQC_IDX_DIR = f"{SCRATCH}/references/_gkmqc_src/data/zju1"
GKMQC_GENOME = "zju1"

PEAK_HALF_WIDTH = 250
MACS3_ORIG_QVALUE = 0.01
MACS3_RELAXED_PVALUE = 0.005
MACS3_SHIFT = -100
MACS3_EXTSIZE = 200
MACS3_GSIZE = "1.1e9"
N_THREADS = 10


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sanitize_pool_name(pool: str) -> str:
    out = pool.replace("/", "_").replace("|", "__").replace(" ", "_")
    return "".join(c if c.isalnum() or c in "._-+" else "_" for c in out)


def parse_chrom_sizes(fai_path: str) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def check_gkmqc_available() -> bool:
    """(G1) Probe for gkmqc.py on PATH and the null-seq index on disk."""
    gkmqc = shutil.which("gkmqc.py")
    if gkmqc is None:
        ts("  WARN: gkmqc.py not on PATH — gkmQC rescue skipped")
        ts("        install via https://github.com/Dongwon-Lee/gkmQC")
        return False
    if not os.path.isdir(GKMQC_IDX_DIR):
        ts(f"  WARN: null-seq index missing at {GKMQC_IDX_DIR}")
        ts("        build via pipeline/build_duck_gkmqc_index.sh")
        return False
    return True


def export_pool_fragments(dataset, pool: str, out_bed_gz: str) -> int:
    """Export fragments for a single peak_pool to a BED.gz file.

    (G3) SnapATAC2 2.9 ``snap.ex.export_fragments`` emits one BED.gz
    per selected group into ``out_dir`` — mirroring ``snap.tl.macs3``'s
    ``selections=`` convention (used in f2).

    snapatac2 v2.9 quirk: passing ``groupby="peak_pool"`` (column name)
    fails on cells where the column is ``None`` (Rust binding cannot
    cast None → str).  Build the per-cell group list explicitly with a
    sentinel replacement for None.
    """
    safe = sanitize_pool_name(pool)
    pool_dir = f"{FRAG_DIR}/{safe}"
    os.makedirs(pool_dir, exist_ok=True)
    # snapatac2 v2.9 uses the group label as the output filename, so pool
    # names with "/" or "|" (e.g., "HH17/16|Splanchnic") break.  Sanitize
    # both the groupby list and the selections key so they match.
    peak_pools_raw = list(dataset.obs["peak_pool"])
    NONE_SENTINEL = "_NO_PEAK_POOL_"
    peak_pools = [sanitize_pool_name(p) if p is not None else NONE_SENTINEL
                  for p in peak_pools_raw]
    # snapatac2 v2.9: the Rust binding ignores the Python `compression`
    # kwarg and always writes Zstandard-compressed .bed even when the
    # filename ends in .gz.  We let it write whatever it writes, then
    # detect the magic bytes and re-encode to actual gzip below so
    # MACS3 (which expects gzip or plain BED) can read the file.
    snap.ex.export_fragments(
        dataset,
        groupby=peak_pools,
        selections={safe},
        out_dir=pool_dir,
        suffix=".bed.zst",
    )

    produced = [p for p in Path(pool_dir).iterdir()
                if p.is_file() and ".bed" in p.name]
    if not produced:
        return 0
    produced.sort(key=lambda p: p.stat().st_size, reverse=True)
    src = produced[0]

    # Detect actual compression via magic bytes and re-encode to gzip
    # for MACS3.  v2.9 hands us a .bed.zst regardless of args.  No
    # Python `zstandard` in the snapatac2 env, so use the `zstd` CLI
    # (available at /usr/bin/zstd) via subprocess streaming.
    import gzip as _gz
    with open(src, "rb") as fh:
        magic = fh.read(4)
    if magic[:4] == b"\x28\xb5\x2f\xfd":  # zstd magic
        ts(f"    decompressing zstd → gzip: {src.name} → {out_bed_gz}")
        with _gz.open(out_bed_gz, "wb", compresslevel=4) as fout:
            proc = subprocess.Popen(["zstd", "-d", "-c", str(src)],
                                    stdout=subprocess.PIPE)
            shutil.copyfileobj(proc.stdout, fout)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"zstd decompression failed: rc={proc.returncode}")
        try:
            os.remove(src)
        except OSError:
            pass
    elif magic[:2] == b"\x1f\x8b":  # gzip magic
        shutil.move(str(src), out_bed_gz)
    else:  # plain BED → gzip it
        ts(f"    gzipping plain BED: {src.name} → {out_bed_gz}")
        with open(src, "rb") as fin, _gz.open(out_bed_gz, "wb",
                                              compresslevel=4) as fout:
            shutil.copyfileobj(fin, fout)
        try:
            os.remove(src)
        except OSError:
            pass

    pool_counts = pd.DataFrame({"peak_pool": list(dataset.obs["peak_pool"])})
    return int((pool_counts["peak_pool"] == pool).sum())


def run_macs3(frag_bed_gz: str, out_prefix: str, *, relaxed: bool) -> str:
    """Call MACS3 at either the plan's orig or relaxed threshold.

    Returns path to the .narrowPeak file.  Fragment input is BED format;
    snATAC adjustments (``--nomodel --shift -100 --extsize 200``) follow
    SnapATAC2's internal defaults so orig-threshold output here matches
    f2's peaks when it did call them.
    """
    args = [
        "macs3", "callpeak",
        "-t", frag_bed_gz,
        "-f", "BED",
        "-g", MACS3_GSIZE,
        "--nomodel",
        "--shift", str(MACS3_SHIFT),
        "--extsize", str(MACS3_EXTSIZE),
        "--keep-dup", "all",
        "--outdir", str(Path(out_prefix).parent),
        "-n", Path(out_prefix).name,
    ]
    if relaxed:
        args += ["-p", str(MACS3_RELAXED_PVALUE)]
    else:
        args += ["-q", str(MACS3_ORIG_QVALUE)]
    ts(f"    macs3 {'relaxed' if relaxed else 'orig'}: {' '.join(args)}")
    subprocess.run(args, check=True)
    return f"{out_prefix}_peaks.narrowPeak"


def run_gkmqc_evaluate(narrow_peak: str, out_prefix: str) -> str:
    args = [
        "gkmqc.py", "evaluate",
        "-i", narrow_peak,
        "-g", GKMQC_GENOME,
        "-n", out_prefix,
        "-@", str(N_THREADS),
    ]
    ts(f"    gkmqc evaluate: {narrow_peak}")
    subprocess.run(args, check=True, env={**os.environ,
                                          "GKMQC_IDX_DIR": GKMQC_IDX_DIR})
    return f"{out_prefix}.gkmqc/{Path(out_prefix).name}.gkmqc.eval.out"


def run_gkmqc_optimize(orig_prefix: str, relaxed_prefix: str) -> str:
    args = [
        "gkmqc.py", "optimize",
        "-p1", orig_prefix,
        "-p2", relaxed_prefix,
    ]
    ts(f"    gkmqc optimize: {orig_prefix} vs {relaxed_prefix}")
    subprocess.run(args, check=True, env={**os.environ,
                                          "GKMQC_IDX_DIR": GKMQC_IDX_DIR})
    return f"{relaxed_prefix}.e300.optz.bed"


def rescue_one_pool(dataset, pool: str, n_cells_expected: int) -> dict | None:
    safe = sanitize_pool_name(pool)
    pool_dir = f"{RESCUE_DIR}/{safe}"
    os.makedirs(pool_dir, exist_ok=True)
    frag_bed_gz = f"{FRAG_DIR}/{safe}.bed.gz"

    ts(f"  [{pool}] exporting fragments ({n_cells_expected:,} cells expected)")
    n_cells = export_pool_fragments(dataset, pool, frag_bed_gz)
    if n_cells == 0:
        ts(f"  [{pool}] no fragments exported — skipped")
        return None

    orig_prefix = f"{pool_dir}/{safe}.orig"
    relaxed_prefix = f"{pool_dir}/{safe}.relaxed"
    orig_peaks = run_macs3(frag_bed_gz, orig_prefix, relaxed=False)
    relaxed_peaks = run_macs3(frag_bed_gz, relaxed_prefix, relaxed=True)

    n_orig = sum(1 for _ in open(orig_peaks)) if os.path.isfile(orig_peaks) else 0
    n_relaxed = sum(1 for _ in open(relaxed_peaks)) if os.path.isfile(relaxed_peaks) else 0
    ts(f"  [{pool}] MACS3 orig={n_orig}  relaxed={n_relaxed}")

    if n_orig < 10:
        ts(f"  [{pool}] orig peak count too low for gkm-SVM training — skipping rescue")
        return {"peak_pool": pool, "n_cells": n_cells,
                "n_orig": n_orig, "n_relaxed": n_relaxed,
                "n_rescued": 0, "rescued_bed": None, "reason": "orig_too_few"}

    run_gkmqc_evaluate(orig_peaks, orig_prefix)
    run_gkmqc_evaluate(relaxed_peaks, relaxed_prefix)
    optz_bed = run_gkmqc_optimize(orig_prefix, relaxed_prefix)

    n_rescued = sum(1 for _ in open(optz_bed)) if os.path.isfile(optz_bed) else 0
    ts(f"  [{pool}] gkmQC optimized → {n_rescued} peaks")
    return {"peak_pool": pool, "n_cells": n_cells,
            "n_orig": n_orig, "n_relaxed": n_relaxed,
            "n_rescued": n_rescued, "rescued_bed": optz_bed,
            "reason": "ok"}


def merge_rescue_into_union(rescue_records: list[dict], chrom_sizes: dict[str, int]) -> str:
    """Fold rescued peaks into the f2 union at 501 bp fixed width."""
    if not os.path.isfile(UNION_BED):
        sys.exit(f"ERROR: union_peaks.bed missing at {UNION_BED}")

    combined = f"{PEAK_DIR}/union_peaks_with_rescue.bed"
    pre = pd.read_csv(UNION_BED, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chrom", "start", "end"])
    extras = []
    for rec in rescue_records:
        if not rec.get("rescued_bed") or not os.path.isfile(rec["rescued_bed"]):
            continue
        df = pd.read_csv(rec["rescued_bed"], sep="\t", header=None,
                         usecols=[0, 1, 2], names=["chrom", "start", "end"])
        extras.append(df)
    if extras:
        pre = pd.concat([pre, *extras], ignore_index=True)

    pre["center"] = ((pre["start"] + pre["end"]) // 2).astype(int)
    pre["start"] = (pre["center"] - PEAK_HALF_WIDTH).clip(lower=0)
    pre["end"] = pre["center"] + PEAK_HALF_WIDTH + 1
    pre["chrom_size"] = pre["chrom"].map(chrom_sizes).fillna(0).astype(int)
    pre = pre[pre["chrom_size"] > 0].copy()
    pre["end"] = pre[["end", "chrom_size"]].min(axis=1)
    pre = pre.sort_values(["chrom", "start", "end"]).drop_duplicates(["chrom", "start", "end"])

    merged_rows = []
    for chrom, grp in pre.groupby("chrom", sort=False):
        grp = grp.sort_values(["start", "end"])
        cur_s, cur_e = None, None
        for _, r in grp.iterrows():
            s, e = int(r["start"]), int(r["end"])
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s < cur_e:
                cur_e = max(cur_e, e)
            else:
                merged_rows.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            merged_rows.append((chrom, cur_s, cur_e))

    merged = pd.DataFrame(merged_rows, columns=["chrom", "start", "end"])
    merged.to_csv(combined, sep="\t", header=False, index=False)
    ts(f"  Post-rescue union: {len(merged):,} peaks → {combined}")
    return combined


def rebuild_peak_matrix(combined_bed: str) -> None:
    """Rebuild step_f2_peak_matrix.h5ad using the post-rescue union."""
    dataset = snap.read_dataset(F1_DATASET)
    peak_matrix_path = f"{ATAC_DIR}/step_f2_peak_matrix.h5ad"
    backup = f"{peak_matrix_path}.pre_rescue.bak"
    if os.path.isfile(peak_matrix_path) and not os.path.isfile(backup):
        shutil.copy2(peak_matrix_path, backup)
        ts(f"  Backed up pre-rescue peak matrix → {backup}")

    try:
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, peak_file=combined_bed, file=peak_matrix_path,
        )
    except TypeError:
        bed_df = pd.read_csv(combined_bed, sep="\t", header=None,
                             usecols=[0, 1, 2], names=["chrom", "start", "end"])
        peak_list = [f"{r.chrom}:{r.start}-{r.end}" for r in bed_df.itertuples()]
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, use_rep=peak_list, file=peak_matrix_path,
        )
    ts(f"  Rebuilt peak matrix: {peak_matrix.shape[0]:,} × {peak_matrix.shape[1]:,}")

    # Transfer metadata from the source AnnDataSet (same approach as
    # f2_peak_calling: read_anndataset_obs + the Polars metadata block).
    # §8.1 quirk: pd.DataFrame(dataset.obs) fails on Polars-backed obs, so
    # pull column by column.  §8 quirk: writing a raw numpy object array into
    # the Polars-backed peak_matrix.obs is read as Binary and rejected, so
    # wrap each column in a pl.Series with an explicit dtype.
    import polars as pl
    cols = ["arc_library", "gex_library", "unified_genotype_id",
            "cell_type", "sex", "stage",
            "d2g_stage", "d2g_cluster", "d2g_subcluster", "peak_pool"]
    data = {}
    for col in cols:
        try:
            data[col] = list(dataset.obs[col])
        except Exception:  # noqa: BLE001 — absent optional column
            continue
    if "peak_pool" not in data:
        raise RuntimeError(f"{F1_DATASET}: obs has no 'peak_pool' column — "
                           "cannot label the rebuilt peak matrix")
    for col, vals in data.items():
        non_null = next((v for v in vals if v is not None), None)
        if isinstance(non_null, (int, np.integer)) and not isinstance(non_null, bool):
            pl_dtype = pl.Int64
        elif isinstance(non_null, (float, np.floating)):
            pl_dtype = pl.Float64
        else:
            pl_dtype = pl.Utf8
        peak_matrix.obs[col] = pl.Series(col, vals, dtype=pl_dtype)
    peak_matrix.close()
    dataset.close()


def main() -> None:
    ts("=== Phase F Step 2b: gkmQC rescue for sub-500 peak pools ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")

    if not os.path.isfile(POOLS_TSV):
        sys.exit(f"ERROR: pool TSV missing at {POOLS_TSV} — run f2 first")
    pools = pd.read_csv(POOLS_TSV, sep="\t")
    ts(f"  {len(pools)} sub-500 pools from f2")
    for _, r in pools.iterrows():
        ts(f"    {r['peak_pool']}: {int(r['n_cells']):,} cells")

    if not check_gkmqc_available():
        summary = pools.copy()
        summary["n_orig"] = 0
        summary["n_relaxed"] = 0
        summary["n_rescued"] = 0
        summary["reason"] = "gkmqc_not_available"
        summary.to_csv(f"{ATAC_DIR}/step_f2b_rescue_summary.tsv",
                       sep="\t", index=False)
        ts("  f2b skipped — summary written with reason=gkmqc_not_available")
        return

    ts("Loading AnnDataSet...")
    dataset = snap.read_dataset(F1_DATASET)

    records: list[dict] = []
    for _, row in pools.iterrows():
        rec = rescue_one_pool(dataset, row["peak_pool"], int(row["n_cells"]))
        if rec is not None:
            records.append(rec)

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(f"{ATAC_DIR}/step_f2b_rescue_summary.tsv",
                      sep="\t", index=False)
    ts(f"  Summary → step_f2b_rescue_summary.tsv")

    chrom_sizes = parse_chrom_sizes(GENOME_FAI)
    combined_bed = merge_rescue_into_union(records, chrom_sizes)
    rebuild_peak_matrix(combined_bed)

    ts("")
    ts("=== Summary ===")
    ts(f"  Pools attempted: {len(pools)}")
    ts(f"  Pools with rescue: {sum(1 for r in records if r['n_rescued'] > 0)}")
    ts(f"  Total rescued peaks: {sum(r['n_rescued'] for r in records):,}")
    ts("=== Done ===")


if __name__ == "__main__":
    main()
