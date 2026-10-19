#!/usr/bin/env python3
"""Vireo seed QC: rank seeds, validate against baseline, compute GT distance.

For each group's 20-seed rerun:
  1. Parse logLik from _log.txt → rank seeds
  2. Compare best seed's donor_ids.tsv against backed-up baseline → concordance %
  3. Parse GT_donors.vireo.vcf.gz → compute pairwise GT distance within group
  4. Flag anomalies (wrong N_donors, logLik regression, low baseline concordance)

Input:  vireo_seeded/{g6_libpair,g9_dp500,g2_old,g3_denovo}/seed_*/
        vireo_baseline_20260330/  (backed-up baseline on standby, rsynced to scratch)
Output: vireo_seeded/seed_qc_report.tsv
        vireo_seeded/gt_distance/  (pairwise correlation matrices per group)

Run:  sbatch FC_08656/run_vireo_seed_qc.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import os
import sys
import glob
import gzip
import time
import numpy as np
import pandas as pd
from collections import Counter

SCRATCH = SCRATCH_ROOT
SEEDED_DIR = f"{SCRATCH}/vireo_seeded"
# Baselines are the existing vireo outputs on scratch (no separate backup needed)
BASELINE_DIR = SCRATCH
OUTDIR = f"{SEEDED_DIR}/seed_qc"
os.makedirs(OUTDIR, exist_ok=True)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ════════════════════════════════════════════════════════════════════════
# 1. Parse vireo outputs
# ════════════════════════════════════════════════════════════════════════
def parse_loglik(log_path):
    """Extract logLik from vireo _log.txt."""
    with open(log_path) as f:
        for line in f:
            if line.startswith("logLik:"):
                return float(line.split(":")[1].strip())
    return None


def parse_summary(summary_path):
    """Parse summary.tsv → dict of donor_id: n_cells."""
    df = pd.read_csv(summary_path, sep="\t")
    return dict(zip(df["Var1"], df["Freq"]))


def count_donors(summary):
    """Count singlet donors (excluding doublet/unassigned)."""
    return sum(1 for k in summary if k not in ("doublet", "unassigned"))


# ════════════════════════════════════════════════════════════════════════
# 2. Seed ranking
# ════════════════════════════════════════════════════════════════════════
def rank_seeds(group_dir, group_name, expected_donors=None):
    """Rank seeds by logLik for a single group."""
    ts(f"\n=== Ranking seeds: {group_name} ===")

    seed_dirs = sorted(glob.glob(f"{group_dir}/seed_*"))
    if not seed_dirs:
        ts(f"  No seed directories found in {group_dir}")
        return pd.DataFrame()

    rows = []
    for sd in seed_dirs:
        seed = os.path.basename(sd).replace("seed_", "")
        log_path = f"{sd}/_log.txt"
        summary_path = f"{sd}/summary.tsv"
        gt_path = f"{sd}/GT_donors.vireo.vcf.gz"

        if not os.path.isfile(log_path):
            ts(f"  SKIP {seed}: no _log.txt (incomplete run)")
            continue

        loglik = parse_loglik(log_path)
        summary = parse_summary(summary_path) if os.path.isfile(summary_path) else {}
        n_donors = count_donors(summary)
        n_singlets = sum(v for k, v in summary.items()
                         if k not in ("doublet", "unassigned"))
        has_gt = os.path.isfile(gt_path)

        rows.append({
            "group": group_name,
            "seed": seed,
            "logLik": loglik,
            "n_donors": n_donors,
            "n_singlets": n_singlets,
            "has_GT_vcf": has_gt,
            "dir": sd,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df = df.sort_values("logLik", ascending=False)  # higher logLik = better
    df["rank"] = range(1, len(df) + 1)

    ts(f"  {len(df)} seeds completed")
    ts(f"  logLik range: {df['logLik'].max():.6e} to {df['logLik'].min():.6e}")
    ts(f"  N_donors: {sorted(df['n_donors'].unique())}")

    # Flag seeds with wrong donor count — use configured target, not top seed
    if expected_donors is not None:
        wrong_n = df[df["n_donors"] != expected_donors]
        if len(wrong_n) > 0:
            ts(f"  WARNING: {len(wrong_n)}/{len(df)} seeds found "
               f"{sorted(wrong_n['n_donors'].unique())} donors instead of "
               f"expected {expected_donors}")
        if df.iloc[0]["n_donors"] != expected_donors:
            ts(f"  *** CRITICAL: best seed has {df.iloc[0]['n_donors']} donors, "
               f"expected {expected_donors} ***")

    return df


# ════════════════════════════════════════════════════════════════════════
# 3. Baseline concordance
# ════════════════════════════════════════════════════════════════════════
def baseline_concordance(new_donor_ids_path, baseline_donor_ids_path):
    """Compare cell-level assignments between new and baseline runs.

    Returns (concordance_rate, n_both_singlet, concordance_matrix).
    Donor labels may differ between runs; finds best 1:1 mapping.
    """
    new = pd.read_csv(new_donor_ids_path, sep="\t")
    baseline = pd.read_csv(baseline_donor_ids_path, sep="\t")

    # Align by cell barcode
    merged = new[["cell", "donor_id"]].merge(
        baseline[["cell", "donor_id"]],
        on="cell", suffixes=("_new", "_base"), how="inner"
    )

    # Filter to both-singlet
    singlet_mask = (
        ~merged["donor_id_new"].isin(["doublet", "unassigned"]) &
        ~merged["donor_id_base"].isin(["doublet", "unassigned"])
    )
    both_sing = merged[singlet_mask]

    if len(both_sing) == 0:
        return 0.0, 0, pd.DataFrame()

    # Build concordance matrix
    conc = pd.crosstab(both_sing["donor_id_new"], both_sing["donor_id_base"])

    # Find best 1:1 mapping (greedy: assign each new donor to its best baseline match)
    mapping = {}
    used_base = set()
    for new_d in conc.index:
        row = conc.loc[new_d].drop(labels=list(used_base), errors="ignore")
        if len(row) == 0:
            continue
        best_base = row.idxmax()
        mapping[new_d] = best_base
        used_base.add(best_base)

    # Count concordant cells
    concordant = 0
    for _, row in both_sing.iterrows():
        if mapping.get(row["donor_id_new"]) == row["donor_id_base"]:
            concordant += 1

    rate = concordant / len(both_sing)
    return rate, len(both_sing), conc


# ════════════════════════════════════════════════════════════════════════
# 4. GT distance from VCF
# ════════════════════════════════════════════════════════════════════════
def parse_gt_vcf(vcf_path):
    """Parse GT_donors.vireo.vcf.gz → donor × variant GT matrix (0/1/2 dosage).

    Returns (donors, positions, gt_matrix) where gt_matrix is n_donors × n_variants.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open

    donors = []
    gt_rows = []

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                donors = fields[9:]  # sample names after FORMAT column
                continue

            fields = line.strip().split("\t")
            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            dosages = []
            for sample_field in fields[9:]:
                gt_str = sample_field.split(":")[gt_idx]
                alleles = gt_str.replace("|", "/").split("/")
                try:
                    dosage = sum(int(a) for a in alleles)
                except ValueError:
                    dosage = np.nan
                dosages.append(dosage)
            gt_rows.append(dosages)

    if not gt_rows:
        return donors, np.array([]), np.array([])

    gt_matrix = np.array(gt_rows, dtype=float).T  # donors × variants
    return donors, gt_matrix


def pairwise_gt_correlation(gt_matrix, donors):
    """Compute pairwise Pearson correlation on GT dosage vectors.

    Returns DataFrame (donors × donors) of correlation values.
    """
    n = len(donors)
    corr = np.full((n, n), np.nan)

    for i in range(n):
        for j in range(i, n):
            vi, vj = gt_matrix[i], gt_matrix[j]
            # Use only positions where both donors have non-NaN GT
            mask = ~(np.isnan(vi) | np.isnan(vj))
            if mask.sum() < 100:
                continue
            r = np.corrcoef(vi[mask], vj[mask])[0, 1]
            corr[i, j] = r
            corr[j, i] = r

    return pd.DataFrame(corr, index=donors, columns=donors)


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════
def main():
    ts("=== Vireo Seed QC ===")

    groups = {
        "g6_libpair": {
            "seeded_dir": f"{SEEDED_DIR}/g6_libpair",
            "baseline_dir": f"{BASELINE_DIR}/libpair_g6/reference_vireo",
            "expected_donors": 37,
        },
        "g9_dp500": {
            "seeded_dir": f"{SEEDED_DIR}/g9_dp500",
            "baseline_dir": f"{BASELINE_DIR}/universal_genotyping/group9/vireo_filt_dp500",
            "expected_donors": 19,
        },
        "g3_denovo": {
            "seeded_dir": f"{SEEDED_DIR}/g3_denovo",
            "baseline_dir": f"{BASELINE_DIR}/universal_genotyping/group3/vireo_N11_denovo_best",
            "expected_donors": 11,
        },
    }

    all_rankings = []

    for group_name, cfg in groups.items():
        # 1. Rank seeds (validated against configured expected_donors)
        ranking = rank_seeds(cfg["seeded_dir"], group_name,
                             expected_donors=cfg["expected_donors"])
        if ranking.empty:
            ts(f"  Skipping {group_name} (no completed seeds)")
            continue

        # 2. Baseline concordance for best seed
        best_seed_dir = ranking.iloc[0]["dir"]
        baseline_donor_ids = f"{cfg['baseline_dir']}/donor_ids.tsv"
        new_donor_ids = f"{best_seed_dir}/donor_ids.tsv"

        if os.path.isfile(baseline_donor_ids) and os.path.isfile(new_donor_ids):
            conc_rate, n_cells, conc_matrix = baseline_concordance(
                new_donor_ids, baseline_donor_ids
            )
            ranking.loc[ranking.index[0], "baseline_concordance"] = conc_rate
            ts(f"  Best seed baseline concordance: {conc_rate:.4f} "
               f"({n_cells} both-singlet cells)")

            if conc_rate < 0.95:
                ts(f"  *** WARNING: concordance {conc_rate:.4f} < 0.95 ***")

            # Save concordance matrix
            conc_out = f"{OUTDIR}/{group_name}_concordance_matrix.tsv"
            conc_matrix.to_csv(conc_out, sep="\t")
        else:
            ts(f"  Baseline comparison skipped (file not found)")

        # 3. GT distance from best seed's VCF
        best_gt_vcf = f"{best_seed_dir}/GT_donors.vireo.vcf.gz"
        if not os.path.isfile(best_gt_vcf):
            best_gt_vcf = f"{best_seed_dir}/GT_donors.vireo.vcf"

        if os.path.isfile(best_gt_vcf):
            ts(f"  Computing GT distance from {os.path.basename(best_gt_vcf)}...")
            donors, gt_matrix = parse_gt_vcf(best_gt_vcf)[:2]
            if gt_matrix.size > 0:
                gt_corr = pairwise_gt_correlation(gt_matrix, donors)
                gt_out = f"{OUTDIR}/{group_name}_gt_pairwise_corr.tsv"
                gt_corr.to_csv(gt_out, sep="\t")
                ts(f"  Saved: {gt_out} ({len(donors)} donors × "
                   f"{gt_matrix.shape[1]} variants)")
        else:
            ts(f"  No GT VCF found for best seed")

        all_rankings.append(ranking)

    # G2: chunked — handle separately
    g2_base = f"{SEEDED_DIR}/g2_old"
    g2_chunks = sorted(glob.glob(f"{g2_base}/chunk_*"))
    for chunk_dir in g2_chunks:
        chunk_name = os.path.basename(chunk_dir)
        ranking = rank_seeds(chunk_dir, f"g2_{chunk_name}", expected_donors=16)
        if not ranking.empty:
            all_rankings.append(ranking)

            # GT distance for best seed
            best_dir = ranking.iloc[0]["dir"]
            for ext in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf"):
                gt_path = f"{best_dir}/{ext}"
                if os.path.isfile(gt_path):
                    donors, gt_matrix = parse_gt_vcf(gt_path)[:2]
                    if gt_matrix.size > 0:
                        gt_corr = pairwise_gt_correlation(gt_matrix, donors)
                        gt_out = f"{OUTDIR}/g2_{chunk_name}_gt_pairwise_corr.tsv"
                        gt_corr.to_csv(gt_out, sep="\t")
                        ts(f"  Saved: {gt_out}")
                    break

    # Combine and save
    if all_rankings:
        combined = pd.concat(all_rankings, ignore_index=True)
        report_path = f"{OUTDIR}/seed_qc_report.tsv"
        combined.to_csv(report_path, sep="\t", index=False)
        ts(f"\nSaved combined report: {report_path}")
        ts(f"Total: {len(combined)} seed runs across {combined['group'].nunique()} groups")

        # Summary
        ts("\n=== Summary ===")
        for group in combined["group"].unique():
            g = combined[combined["group"] == group]
            best = g.iloc[0]
            ts(f"  {group}: best seed={best['seed']}, logLik={best['logLik']:.6e}, "
               f"donors={best['n_donors']}, singlets={best['n_singlets']}")
            bc = best.get("baseline_concordance")
            if pd.notna(bc):
                ts(f"    baseline concordance: {bc:.4f}")
    else:
        ts("No completed seed runs found!")

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
