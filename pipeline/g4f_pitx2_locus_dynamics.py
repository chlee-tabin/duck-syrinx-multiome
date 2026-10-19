#!/usr/bin/env python3
"""g4f_pitx2_locus_dynamics.py — depth-confound-controlled PITX2-locus peak dynamics
per (compartment × stage), over the 359-peak Pitx2-TAD focal set.

REPLICATE UNIT = EMBRYO-HALF (rectified 2026-06-16).  An embryo-half is one
(compartment × stage × genotype × side) pseudobulk built by SUMMING the genome-
wide counts over ALL inference-pool subclusters that the half clears the floor
in.  This fixes a pseudoreplication: previously the laterality/stage tests fit
DESeq2 over per-(genotype × d2g_subcluster × side) pseudobulks, so one embryo-half
entered as 2–3 *correlated* "L" (or "R") replicates — inflating significance AND
distorting the apparent L/R balance (e.g. HH25 looked 12L/21R in pseudobulks but
is 8L/8R in embryo-halves).  Dorsal was already 1× (single subcluster per stage);
HH13 ventral was 1× (single subcluster) — both unchanged.  Pooling to the embryo-
half is the fig02 panel C / figS02D convention (one pseudobulk per embryo-half).

Two orthogonal dynamics, each reported as EFFECT SIZES (log2FC distributions), never
raw sig-peak counts alone — because sig counts scale with depth / n_cells /
n_halves (power), not just biology:

  LATERALITY  : L vs R within a compartment, per stage.  Computed TWO ways:
      unpaired  ~ inferred_side                         (ALL embryo-halves; max power,
                                                          but side may be confounded
                                                          with embryo identity when the
                                                          L- and R-halves come from
                                                          different embryos)
      paired    ~ unified_genotype_id + inferred_side   (both-sided genotypes only:
                                                          within-embryo L-vs-R, removes
                                                          between-embryo variance; loses
                                                          one-sided embryos — at HH21/HH25
                                                          this loses NOTHING (fully paired),
                                                          at HH17/16 it keeps 4 of 11)
    Both are emitted (lr_long carries a `design` column).  The panel leads with
    unpaired (primary, uniform, all data); paired confirms.
  STAGE       : gained / lost within a compartment, consecutive stage pairs, over
                embryo-halves (~ [inferred_side +] commit_stage, contrast late/early).

Compartment pools per stage (g2c SSOT, identical to the earlier focal-locus L/R-by-
compartment analysis — validated by prior use):
  ventral (= inference / laterality pool):
    HH13   : Early LPM, Cardiopharyngeal LPM
    HH17/16: Early LPM, Cardiopharyngeal LPM, BARX1+
    HH21   : Ventral spl., Syringeal, BARX1+
    HH25   : Ventral spl., Syringeal, Chondrogenic LPM
  dorsal:
    HH17/16: Dorsal-biased LPM        HH21/HH25: Dorsal sub-epi.   (HH13: none)

DESeq2 size factors are fit GENOME-WIDE (full count>=MIN matrix of the POOLED
embryo-halves) so depth normalization is not driven by the 359 locus peaks;
results are then restricted to the locus + locus-FDR.  A global poscounts size-
factor vector + the 359-peak count slice (over embryo-halves) are exported so the
depth-matched subsampling / fragment-downsampling re-runs can be done locally
without re-touching the cluster.  The per-subcluster slice + meta are ALSO exported
(`*_bysubcluster.*`) for any subcluster-resolved supporting analysis (which must
represent N as CELLS/EMBRYOS, never pseudobulk count).

Outputs (all to $OUT_DIR, default $DA_DIR/pitx2_locus_dynamics):
  g4f_census.tsv              per (stage,compartment): n_pb (per-subcluster, for
                              transparency) / n_half / n_L_half / n_R_half / n_paired
                              / Sum n_cells / median libsize   [the power panel]
  g4f_cpm_trajectory_long.tsv per (peak, compartment, stage, side[all/L/R]): mean CPM
                              over EMBRYO-HALVES, n_half  [depth-normalized dynamics]
  g4f_lr_long.tsv             per (peak, stage, compartment, design): full L/R log2FC
                              + p + padj + padj_locus (ALL tested locus peaks)
  g4f_lr_grid.tsv             per (stage,compartment,design): n_L/R_half, n_paired,
                              n_locus_tested, n_L_sig/n_R_sig, median_log2FC, frac_L
  g4f_stage_long.tsv          per (peak, compartment, stage_pair): late/early log2FC...
  g4f_stage_grid.tsv          per (compartment, stage_pair): n_early/late_half,
                              gw_gained/gw_lost/gw_median_log2FC, locus_*,  design
  g4f_locus_count_slice.tsv.gz  359 peaks x EMBRYO-HALVES raw counts (peaks on rows)
  g4f_sample_meta.tsv         per EMBRYO-HALF: libsize, size_factor_global, n_cells,
                              total_fragments, n_subclusters, stage, side, genotype, compartment
  g4f_locus_count_slice_bysubcluster.tsv.gz / g4f_sample_meta_bysubcluster.tsv
                              the pre-pooling per-(genotype×subcluster×side) substrate

Run (snapatac2 env; pydeseq2):  sbatch pipeline/run_pitx2_locus_dynamics.sh
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys
import time

import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
COUNTS_TSV = os.environ.get("COUNTS_TSV", f"{DA_DIR}/pseudobulk_counts.tsv.gz")
META_TSV = os.environ.get("META_TSV", f"{DA_DIR}/pseudobulk_metadata.tsv")
LOCUS_TSV = os.environ.get("LOCUS_TSV", f"{SCRATCH}/figures/data/pitx2_locus_tad_peaks.tsv")
OUT_DIR = os.environ.get("OUT_DIR", f"{DA_DIR}/pitx2_locus_dynamics")
MIN_PEAK = int(os.environ.get("MIN_PEAK_COUNT", "10"))
PADJ = 0.1

# Pool population lists come from the POOLS single-source-of-truth
# (figures/_shared/subcluster_map.py) — the same source f4/h1/h8 import, NOT a
# hardcoded copy.  A POOLS rename would otherwise silently desync g4f/g4g from the
# fig02/figS02 grouping (project convention: pool membership is POOLS-owned).  Fail loud if
# POOLS can't be imported — never fall back to a stale hardcoded copy.
# Use a REPO-RELATIVE figures/ path (not $SCRATCH/figures like f4/h8): g4f is also
# imported LOCALLY by g4g + the _test_/_check_ diagnostics, where the scratch path
# does not exist.  pipeline/../figures resolves both locally and on the cluster
# (pipeline/ + figures/ are synced side-by-side to $SCRATCH).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from _shared.subcluster_map import POOLS  # noqa: E402

STAGES = ["HH13", "HH17/16", "HH21", "HH25"]
VENTRAL = {st: POOLS[st]["inference"] for st in STAGES if POOLS.get(st, {}).get("inference")}
DORSAL = {st: POOLS[st]["bilateral"] for st in STAGES if POOLS.get(st, {}).get("bilateral")}
COMPARTMENTS = (("ventral", VENTRAL), ("dorsal", DORSAL))


def ts(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _det(df, keys):
    """Sort a per-peak output by its key columns before writing, for a deterministic
    row ORDER: rows are appended in set(locus_ids) iteration order (hash-randomized
    per process), so without this the order churns every run on identical inputs —
    bloating diffs and breaking bit-reproducibility (values unchanged; consumers read
    by merge-key, so this is pure hygiene). Shared by g4f + g4g (g4g imports it) so the
    two engines stay byte-identical; the committed ``*_long`` base files reorder once on
    their next regen. Gates on row-count (len), NOT column membership: the empty-result
    path returns a column-less ``pd.DataFrame()`` so ``sort_values`` would KeyError —
    len==0 is the real trigger. Do NOT "improve" this to a
    ``set(keys) <= df.columns`` check: that would silently no-op a non-empty frame that
    legitimately lost a key column.
    """
    return df.sort_values(keys) if len(df) else df


def load_counts_meta():
    ts("loading metadata + counts...")
    meta = pd.read_csv(META_TSV, sep="\t")
    t0 = time.time()
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
    ids = set(meta["sample_id"])
    # orient so rows = samples, cols = peaks
    if len(set(counts.index) & ids) < len(set(counts.columns) & ids):
        counts = counts.T
    ts(f"  counts {counts.shape[0]} samples x {counts.shape[1]:,} peaks in {time.time()-t0:.1f}s")
    common = [s for s in meta["sample_id"] if s in set(counts.index)]
    meta = meta[meta["sample_id"].isin(common)].copy()
    counts = counts.loc[meta["sample_id"]].copy()
    return counts, meta


def assign_compartment(meta: pd.DataFrame) -> pd.DataFrame:
    """Tag each pseudobulk ventral/dorsal/none by (stage, subcluster)."""
    comp = pd.Series("none", index=meta.index)
    for label, table in COMPARTMENTS:
        for stage, pops in table.items():
            sel = (meta["d2g_stage"] == stage) & (meta["d2g_subcluster"].isin(pops))
            comp[sel] = label
    meta = meta.copy()
    meta["compartment"] = comp.values
    return meta


def build_embryo_half(counts: pd.DataFrame, meta: pd.DataFrame):
    """Pool (SUM genome-wide counts) over all inference-pool subclusters within each
    (compartment × stage × genotype × side) -> ONE pseudobulk per EMBRYO-HALF.

    This is the laterality/stage replicate unit (fixes the per-subcluster
    pseudoreplication).  'none'-compartment pseudobulks (Somatic / Mesothelium /
    Venous-pole spl. / paraxial / ambiguous) are dropped — they are neither the
    inference pool nor the bilateral control.

    Returns (pooled_counts [halves x peaks], half_meta [one row per half]).
    """
    keep = meta[meta["compartment"].isin(["ventral", "dorsal"])].copy()
    keep["half_id"] = (
        keep["compartment"] + "|" + keep["d2g_stage"].astype(str) + "|"
        + keep["unified_genotype_id"].astype(str) + "|" + keep["inferred_side"].astype(str)
    )
    s2h = keep.set_index("sample_id")["half_id"]
    c = counts.loc[keep["sample_id"]].copy()
    c.index = s2h.loc[c.index].to_numpy()
    pooled = c.groupby(level=0).sum()
    half_meta = (
        keep.groupby("half_id")
        .agg(
            compartment=("compartment", "first"),
            d2g_stage=("d2g_stage", "first"),
            unified_genotype_id=("unified_genotype_id", "first"),
            inferred_side=("inferred_side", "first"),
            n_cells=("n_cells", "sum"),
            total_fragments=("total_fragments", "sum"),
            n_subclusters=("d2g_subcluster", "nunique"),
            subclusters=("d2g_subcluster", lambda s: ",".join(sorted(set(s)))),
        )
        .reset_index()
        .rename(columns={"half_id": "sample_id"})
    )
    pooled = pooled.loc[half_meta["sample_id"]]
    half_meta["libsize"] = pooled.sum(axis=1).to_numpy()
    ts(f"  pooled {len(meta)} per-subcluster pseudobulks (ventral+dorsal+none) -> "
       f"{len(half_meta)} embryo-halves "
       f"({int((half_meta.compartment=='ventral').sum())} ventral, "
       f"{int((half_meta.compartment=='dorsal').sum())} dorsal)")
    return pooled, half_meta


def global_size_factors(counts: pd.DataFrame) -> pd.Series:
    """poscounts median-of-ratios size factors over the genome-wide count>=MIN matrix."""
    mat = counts.loc[:, counts.sum(axis=0) >= MIN_PEAK].to_numpy(dtype=float)  # samples x peaks
    with np.errstate(divide="ignore"):
        logc = np.where(mat > 0, np.log(mat), np.nan)
    loggeo = np.nanmean(logc, axis=0)            # per-peak log geomean over positive samples
    ok = np.isfinite(loggeo)
    sf = np.empty(mat.shape[0])
    for i in range(mat.shape[0]):
        r = logc[i, ok] - loggeo[ok]
        sf[i] = np.exp(np.nanmedian(r))
    sf = sf / np.exp(np.nanmean(np.log(sf)))     # center (geomean 1)
    return pd.Series(sf, index=counts.index, name="size_factor_global")


def _deseq(counts_sub: pd.DataFrame, md: pd.DataFrame, design, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    # Parallelize dispersion/LFC fitting ACROSS PEAKS (pydeseq2's joblib
    # inference) over the allocated cores — the genome-wide fit (283K peaks ×
    # n_contrasts) is the slow tail, and historically NO pipeline step set
    # n_cpus, so this parallelism was unused (canonical g4f ~100 min at -c 2).
    # Built defensively: tries (poscounts + inference) first, then drops
    # whichever kwarg an older pydeseq2 rejects.  It can ONLY speed things up —
    # never break (one combo always works) and never change results
    # (parallelism is numerically identical).  Set N_CPUS (run script) to enable.
    base = dict(counts=counts_sub.astype(int), metadata=md,
                design_factors=design, refit_cooks=True)
    extra = {}
    n_cpus = int(os.environ.get("N_CPUS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    if n_cpus > 1:
        from pydeseq2.default_inference import DefaultInference
        extra["inference"] = DefaultInference(n_cpus=n_cpus)
    # pydeseq2==0.5.4 (envs/o2_snapatac2.yml) takes these kwargs directly; one
    # constructor call, no TypeError cascade that could drop poscounts silently.
    dds = DeseqDataSet(**base, size_factors_fit_type="poscounts", **extra)
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast)
    st.summary()
    return st.results_df.copy()


def _half_side_coverage(h: pd.DataFrame) -> dict:
    """genotype -> set of sides present among that genotype's halves in `h`."""
    return (h.groupby("unified_genotype_id")["inferred_side"]
            .apply(lambda s: set(s)).to_dict())


def census(meta: pd.DataFrame, half_meta: pd.DataFrame) -> pd.DataFrame:
    """Per (stage, compartment): pseudobulk counts (transparency) AND the honest
    embryo-half counts that drive the DESeq2 replication + the panel axis."""
    rows = []
    for label, table in COMPARTMENTS:
        for stage in STAGES:
            pops = table.get(stage)
            if not pops:
                continue
            pb = meta[(meta["d2g_stage"] == stage) & (meta["d2g_subcluster"].isin(pops))]
            h = half_meta[(half_meta["compartment"] == label)
                          & (half_meta["d2g_stage"] == stage)]
            nLh = int((h["inferred_side"] == "L").sum())
            nRh = int((h["inferred_side"] == "R").sum())
            cov = _half_side_coverage(h)
            n_paired = sum(1 for v in cov.values() if {"L", "R"} <= v)
            rows.append({
                "stage": stage, "compartment": label,
                # per-subcluster pseudobulk counts (kept ONLY for transparency)
                "n_pb": len(pb),
                "n_L_pb": int((pb["inferred_side"] == "L").sum()),
                "n_R_pb": int((pb["inferred_side"] == "R").sum()),
                # embryo-half counts (THE sample size)
                "n_half": len(h), "n_L_half": nLh, "n_R_half": nRh,
                "n_paired": n_paired,
                "sum_n_cells": int(h["n_cells"].sum()) if len(h) else 0,
                "sum_total_fragments": int(h["total_fragments"].sum()) if len(h) else 0,
                "median_n_cells_half": float(h["n_cells"].median()) if len(h) else np.nan,
                "median_libsize_half": float(h["libsize"].median()) if len(h) else np.nan,
                # column gates for the panel (derived, never hardcoded)
                "lr_testable": int(nLh >= 2 and nRh >= 2),          # unpaired column
                "paired_testable": int(n_paired >= 2),               # paired column
            })
    return pd.DataFrame(rows)


def cpm_trajectory(cpm: pd.DataFrame, half_meta: pd.DataFrame, locus_ids) -> pd.DataFrame:
    """mean CPM per (peak, compartment, stage, side) over EMBRYO-HALVES."""
    loc_cols = [p for p in locus_ids if p in cpm.columns]
    cpm_loc = cpm[loc_cols]
    out = []
    for label, table in COMPARTMENTS:
        for stage in STAGES:
            if stage not in table:
                continue
            base = half_meta[(half_meta["compartment"] == label)
                             & (half_meta["d2g_stage"] == stage)]
            for side, sel in (("all", base),
                              ("L", base[base["inferred_side"] == "L"]),
                              ("R", base[base["inferred_side"] == "R"])):
                if len(sel) == 0:
                    continue
                m = cpm_loc.loc[sel["sample_id"]].mean(axis=0)
                for pid, v in m.items():
                    out.append({"peak_id": pid, "compartment": label, "stage": stage,
                                "side": side, "mean_cpm": round(float(v), 5),
                                "n_half": len(sel)})
    return pd.DataFrame(out)


def _lr_one(counts_h, h, design_cols, dname, stage, label, locus_ids):
    """Run one laterality DESeq2 fit (unpaired or paired) over embryo-halves;
    return (locus_long_df, grid_row); a failed fit raises."""
    c = counts_h.loc[h["sample_id"]]
    c = c.loc[:, c.sum(axis=0) >= MIN_PEAK]
    md = h.set_index("sample_id")[design_cols].copy()
    md["inferred_side"] = pd.Categorical(md["inferred_side"], categories=["R", "L"])
    if "unified_genotype_id" in design_cols:
        md["unified_genotype_id"] = pd.Categorical(md["unified_genotype_id"])
    nL = int((h["inferred_side"] == "L").sum())
    nR = int((h["inferred_side"] == "R").sum())
    n_paired = sum(1 for v in _half_side_coverage(h).values() if {"L", "R"} <= v)
    res = _deseq(c, md, design_cols, ["inferred_side", "L", "R"])   # a failed fit raises
    res.index.name = "peak_id"
    res = res.reset_index()
    loc = res[res["peak_id"].isin(locus_ids) & res["pvalue"].notna()].copy()
    from statsmodels.stats.multitest import multipletests
    loc["padj_locus"] = (multipletests(loc["pvalue"], method="fdr_bh")[1]
                         if len(loc) else [])
    loc["stage"], loc["compartment"], loc["design"] = stage, label, dname
    nLs = int(((loc["log2FoldChange"] > 0) & (loc["padj_locus"] < PADJ)).sum())
    nRs = int(((loc["log2FoldChange"] < 0) & (loc["padj_locus"] < PADJ)).sum())
    grid = {"stage": stage, "compartment": label, "design": dname,
            "n_L_half": nL, "n_R_half": nR, "n_paired": n_paired,
            "status": "ok", "n_locus_tested": len(loc), "n_L_sig": nLs, "n_R_sig": nRs,
            "median_log2FC": round(float(loc["log2FoldChange"].median()), 4) if len(loc) else np.nan,
            "frac_L_of_tested": round(float((loc["log2FoldChange"] > 0).mean()), 4) if len(loc) else np.nan}
    ts(f"  LR {stage} {label} [{dname}]: L={nL}/R={nR} paired={n_paired}, "
       f"tested={len(loc)}, L_sig={nLs} R_sig={nRs}, "
       f"medL2FC={loc['log2FoldChange'].median():+.3f}" if len(loc) else
       f"  LR {stage} {label} [{dname}]: L={nL}/R={nR} paired={n_paired}, tested=0")
    return loc, grid


def lr_dynamics(counts_h: pd.DataFrame, half_meta: pd.DataFrame, locus_ids):
    """Laterality over embryo-halves, BOTH designs:
       unpaired  ~ inferred_side                       (all halves; >=2L & >=2R)
       paired    ~ unified_genotype_id + inferred_side (both-sided genotypes only; >=2 pairs)
    """
    long_rows, grid_rows = [], []
    for label, table in COMPARTMENTS:
        for stage in STAGES:
            if stage not in table:
                continue
            h = half_meta[(half_meta["compartment"] == label)
                          & (half_meta["d2g_stage"] == stage)].copy()
            nL = int((h["inferred_side"] == "L").sum())
            nR = int((h["inferred_side"] == "R").sum())
            cov = _half_side_coverage(h)
            paired_g = [g for g, v in cov.items() if {"L", "R"} <= v]

            # ── unpaired (primary): all embryo-halves ──
            if nL >= 2 and nR >= 2:
                loc, grid = _lr_one(counts_h, h, ["inferred_side"], "unpaired",
                                    stage, label, locus_ids)
                if loc is not None:
                    long_rows.append(loc)
                    grid_rows.append(grid)
            else:
                ts(f"  LR {stage} {label} [unpaired]: SKIP (L={nL}, R={nR})")
                grid_rows.append({"stage": stage, "compartment": label, "design": "unpaired",
                                  "n_L_half": nL, "n_R_half": nR, "n_paired": len(paired_g),
                                  "status": "underpowered", "n_locus_tested": 0, "n_L_sig": 0,
                                  "n_R_sig": 0, "median_log2FC": np.nan, "frac_L_of_tested": np.nan})

            # ── paired (confirmatory): both-sided genotypes only ──
            if len(paired_g) >= 2:
                hp = h[h["unified_genotype_id"].isin(paired_g)].copy()
                loc, grid = _lr_one(counts_h, hp, ["unified_genotype_id", "inferred_side"],
                                    "paired", stage, label, locus_ids)
                if loc is not None:
                    long_rows.append(loc)
                    grid_rows.append(grid)
            else:
                ts(f"  LR {stage} {label} [paired]: SKIP (n_paired={len(paired_g)})")
                grid_rows.append({"stage": stage, "compartment": label, "design": "paired",
                                  "n_L_half": len(paired_g), "n_R_half": len(paired_g),
                                  "n_paired": len(paired_g), "status": "underpowered",
                                  "n_locus_tested": 0, "n_L_sig": 0, "n_R_sig": 0,
                                  "median_log2FC": np.nan, "frac_L_of_tested": np.nan})
    long = (pd.concat(long_rows, ignore_index=True)[
        ["peak_id", "stage", "compartment", "design", "log2FoldChange",
         "pvalue", "padj", "padj_locus"]]
        if long_rows else pd.DataFrame())
    return long, pd.DataFrame(grid_rows)


def stage_dynamics(counts_h: pd.DataFrame, half_meta: pd.DataFrame, locus_ids):
    """Gained/lost within a compartment over consecutive stages, EMBRYO-HALVES.
    design ~ [inferred_side +] commit_stage, contrast late/early."""
    long_rows, grid_rows = [], []
    for label, table in COMPARTMENTS:
        sts = [s for s in STAGES if s in table]
        for early, late in zip(sts[:-1], sts[1:]):
            e = half_meta[(half_meta["compartment"] == label)
                          & (half_meta["d2g_stage"] == early)].copy()
            l = half_meta[(half_meta["compartment"] == label)
                          & (half_meta["d2g_stage"] == late)].copy()
            if len(e) < 2 or len(l) < 2:
                ts(f"  STAGE {label} {early}->{late}: SKIP (e={len(e)}, l={len(l)})")
                continue
            e["commit_stage"], l["commit_stage"] = "early", "late"
            sub = pd.concat([e, l], ignore_index=True)
            side_ok = all(((sub["commit_stage"] == g) & (sub["inferred_side"] == s)).sum() >= 1
                          for g in ("early", "late") for s in ("L", "R"))
            design = (["inferred_side", "commit_stage"] if side_ok else ["commit_stage"])
            c = counts_h.loc[sub["sample_id"]]
            c = c.loc[:, c.sum(axis=0) >= MIN_PEAK]
            md = sub.set_index("sample_id").copy()
            md["commit_stage"] = pd.Categorical(md["commit_stage"], categories=["early", "late"])
            if "inferred_side" in design:
                md["inferred_side"] = pd.Categorical(md["inferred_side"])
            res = _deseq(c, md[design], design, ["commit_stage", "late", "early"])  # raises on failure
            res.index.name = "peak_id"
            res = res.reset_index()
            ok = res["pvalue"].notna()
            gw_gain = int(((res["padj"] < PADJ) & (res["log2FoldChange"] > 0)).sum())
            gw_lost = int(((res["padj"] < PADJ) & (res["log2FoldChange"] < 0)).sum())
            loc = res[res["peak_id"].isin(locus_ids) & res["pvalue"].notna()].copy()
            from statsmodels.stats.multitest import multipletests
            loc["padj_locus"] = (multipletests(loc["pvalue"], method="fdr_bh")[1]
                                 if len(loc) else [])
            loc["compartment"], loc["stage_pair"] = label, f"{early}->{late}"
            long_rows.append(loc)
            grid_rows.append({
                "compartment": label, "stage_pair": f"{early}->{late}", "design": "+".join(design),
                "n_early_half": len(e), "n_late_half": len(l),
                "gw_tested": int(ok.sum()), "gw_gained": gw_gain, "gw_lost": gw_lost,
                "gw_median_log2FC": round(float(res.loc[ok, "log2FoldChange"].median()), 4),
                "locus_tested": len(loc),
                "locus_gained": int(((loc["padj_locus"] < PADJ) & (loc["log2FoldChange"] > 0)).sum()),
                "locus_lost": int(((loc["padj_locus"] < PADJ) & (loc["log2FoldChange"] < 0)).sum()),
                "locus_median_log2FC": round(float(loc["log2FoldChange"].median()), 4) if len(loc) else np.nan})
            ts(f"  STAGE {label} {early}->{late} [{'+'.join(design)}]: "
               f"e={len(e)} l={len(l)} | gw gain={gw_gain} lost={gw_lost} "
               f"med={res.loc[ok,'log2FoldChange'].median():+.3f} | locus med={loc['log2FoldChange'].median():+.3f}")
    long = (pd.concat(long_rows, ignore_index=True)[
        ["peak_id", "compartment", "stage_pair", "log2FoldChange", "pvalue", "padj", "padj_locus"]]
        if long_rows else pd.DataFrame())
    return long, pd.DataFrame(grid_rows)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    ts("=== g4f_pitx2_locus_dynamics (replicate unit = EMBRYO-HALF) ===")
    for p in (COUNTS_TSV, META_TSV, LOCUS_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing input {p}")
    locus = pd.read_csv(LOCUS_TSV, sep="\t")
    locus_ids = set(locus["peak_id"].astype(str))
    ts(f"locus peaks: {len(locus_ids)}")

    counts, meta = load_counts_meta()
    meta["libsize"] = counts.sum(axis=1).reindex(meta["sample_id"]).values
    meta = assign_compartment(meta)

    # ── pool to EMBRYO-HALF (the laterality/stage replicate unit) ──
    counts_h, half_meta = build_embryo_half(counts, meta)

    # pool sanity: every (stage,compartment) cell must be non-empty
    cen = census(meta, half_meta)
    empties = cen[cen["n_half"] == 0]
    if len(empties):
        raise SystemExit(f"ERROR: {len(empties)} empty pool cells (pool-name drift?):\n{empties}")
    cen.to_csv(f"{OUT_DIR}/g4f_census.tsv", sep="\t", index=False)
    ts(f"  -> g4f_census.tsv\n{cen.to_string(index=False)}")

    # global size factors (genome-wide, over the POOLED embryo-half matrix) + CPM
    sf = global_size_factors(counts_h)
    cpm = counts_h.div(half_meta.set_index("sample_id")["libsize"], axis=0) * 1e6

    # _det (module scope, shared with g4g) sorts each per-peak output by its key
    # columns for a deterministic row order — see its docstring for the len-gate rationale.
    traj = cpm_trajectory(cpm, half_meta, locus_ids)
    _det(traj, ["peak_id", "compartment", "stage", "side"]).to_csv(
        f"{OUT_DIR}/g4f_cpm_trajectory_long.tsv", sep="\t", index=False)
    ts(f"  -> g4f_cpm_trajectory_long.tsv ({len(traj)} rows)")

    lr_long, lr_grid = lr_dynamics(counts_h, half_meta, locus_ids)
    _det(lr_long, ["peak_id", "stage", "compartment", "design"]).to_csv(
        f"{OUT_DIR}/g4f_lr_long.tsv", sep="\t", index=False, float_format="%.5g")
    lr_grid.to_csv(f"{OUT_DIR}/g4f_lr_grid.tsv", sep="\t", index=False)
    ts(f"  -> g4f_lr_long.tsv / g4f_lr_grid.tsv\n{lr_grid.to_string(index=False)}")

    st_long, st_grid = stage_dynamics(counts_h, half_meta, locus_ids)
    _det(st_long, ["peak_id", "compartment", "stage_pair"]).to_csv(
        f"{OUT_DIR}/g4f_stage_long.tsv", sep="\t", index=False, float_format="%.5g")
    st_grid.to_csv(f"{OUT_DIR}/g4f_stage_grid.tsv", sep="\t", index=False)
    ts(f"  -> g4f_stage_long.tsv / g4f_stage_grid.tsv\n{st_grid.to_string(index=False)}")

    # locus count slice (peaks x EMBRYO-HALVES) + per-half meta, for LOCAL depth-matched work
    loc_cols = [p for p in locus_ids if p in counts_h.columns]
    slice_df = counts_h[loc_cols].T  # peaks x halves
    slice_df.index.name = "peak_id"
    slice_df.to_csv(f"{OUT_DIR}/g4f_locus_count_slice.tsv.gz", sep="\t")
    smeta = half_meta[["sample_id", "libsize", "n_cells", "total_fragments", "d2g_stage",
                       "n_subclusters", "subclusters", "inferred_side",
                       "unified_genotype_id", "compartment"]].copy()
    smeta = smeta.merge(sf.rename("size_factor_global"), left_on="sample_id", right_index=True)
    smeta.to_csv(f"{OUT_DIR}/g4f_sample_meta.tsv", sep="\t", index=False, float_format="%.6g")
    ts(f"  -> g4f_locus_count_slice.tsv.gz ({slice_df.shape[0]} peaks x {slice_df.shape[1]} halves)"
       f" + g4f_sample_meta.tsv ({len(smeta)} halves)")

    # per-subcluster substrate (pre-pooling) for subcluster-resolved supporting analyses
    sf_pb = global_size_factors(counts)
    loc_cols_pb = [p for p in locus_ids if p in counts.columns]
    slice_pb = counts[loc_cols_pb].T
    slice_pb.index.name = "peak_id"
    slice_pb.to_csv(f"{OUT_DIR}/g4f_locus_count_slice_bysubcluster.tsv.gz", sep="\t")
    smeta_pb = meta[["sample_id", "libsize", "n_cells", "total_fragments", "d2g_stage",
                     "d2g_subcluster", "inferred_side", "unified_genotype_id", "compartment"]].copy()
    smeta_pb = smeta_pb.merge(sf_pb.rename("size_factor_global"), left_on="sample_id", right_index=True)
    smeta_pb.to_csv(f"{OUT_DIR}/g4f_sample_meta_bysubcluster.tsv", sep="\t", index=False, float_format="%.6g")
    ts(f"  -> g4f_locus_count_slice_bysubcluster.tsv.gz ({slice_pb.shape[0]} peaks x {slice_pb.shape[1]} pseudobulks)"
       f" + g4f_sample_meta_bysubcluster.tsv ({len(smeta_pb)} pseudobulks)")
    ts("=== g4f DONE ===")


if __name__ == "__main__":
    main()
