#!/usr/bin/env python3
"""Phase H Step 3: per-nucleus detection-fraction export for NODAL,
CFC1/Cryptic and PITX2 (positive control).

Committed producer behind the manuscript statement that NODAL transcripts sit
at the multiome's detection floor ("~0.1% of nuclei"); this TSV makes that
number reproducible, with explicit denominators.

Scopes
------
  all_nuclei      every post-doublet singlet nucleus in the DEG-substrate
                  h5ad (step02_preCB_clean — the h1/h2 universe; raw
                  CellBender-corrected integer counts, all genes)
  inference_pool  nuclei whose d2g population is in the stage's laterality
                  inference pool (POOLS[stage]["inference"]) — the population
                  the letter's HH13 paired DE runs on

Stage anchors come from the FIXED e1b sidecar per genotype (same discipline
as h1); "all" rows aggregate every scope nucleus regardless of stage.

Memory note (O2 RC over-request notice, 2026-08-25): doublets are MASKED, not
subset — no `.copy()` of the full matrix — so peak RAM ≈ one h5ad load. Run
via pipeline/run_h3.sh (right-sized #SBATCH), not run_pipeline.sh's 64G
default header.

Output
------
  figures/data/h3_gene_detection.tsv
    scope, stage, gene_label, gene_id, n_nuclei, n_detected,
    frac_detected, total_umi

Run on O2
---------
  sbatch pipeline/run_h3.sh
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"),
)
from _shared.subcluster_map import POOLS  # noqa: E402

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
H5AD = f"{PHASE2}/step02_preCB_clean.h5ad"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
ANNOT_TSV = f"{SCRATCH}/references/d2g_subcluster_annotations.tsv"
D2G_TSVS = {
    "HH13":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}
LEIDEN_COL_PER_STAGE = {
    "HH13":    "leiden_0.2",
    "HH17/16": "leiden_0.3",
    "HH21":    "leiden_0.3",
    "HH25":    "leiden_0.3",
}
STAGE_TO_ANCHOR = {
    "HH13": "HH13", "HH16": "HH17/16", "HH17": "HH17/16",
    "HH17/16": "HH17/16", "HH21": "HH21", "HH25": "HH25",
}
OUT_TSV = f"{SCRATCH}/figures/data/h3_gene_detection.tsv"

# (label, var id in the h5ad). Missing id = hard error, never skip:
# a renamed/absent gene id means the wrong assembly or a stale h5ad.
GENES = [
    ("NODAL", "LOC101798108"),
    ("CFC1_Cryptic", "LOC101799976"),
    ("PITX2", "PITX2"),
]

STAGE_ANCHORS = ["HH13", "HH17/16", "HH21", "HH25"]


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_e1b_stage_map(tsv_path: str) -> dict[str, str]:
    df = pd.read_csv(tsv_path, sep="\t")
    return dict(zip(df["genotype_id"].astype(str),
                    df["inferred_stage"].astype(str)))


def load_inference_pool_cells() -> dict[tuple[str, str, str], str]:
    """(genotype, library, barcode) -> stage anchor, for inference-pool cells."""
    annot = pd.read_csv(ANNOT_TSV, sep="\t")
    out: dict[tuple[str, str, str], str] = {}
    for anchor, tsv_path in D2G_TSVS.items():
        if not os.path.isfile(tsv_path):
            sys.exit(f"ERROR: missing d2g sidecar {tsv_path}")
        leiden_col = LEIDEN_COL_PER_STAGE[anchor]
        df = pd.read_csv(tsv_path, sep="\t", compression="gzip",
                         usecols=["genotype", "library", "barcode", leiden_col])
        annot_sub = annot[(annot["stage"] == anchor) &
                          (annot["leiden_col"] == leiden_col)]
        leiden_to_pop = dict(zip(annot_sub["cluster"].astype(int),
                                 annot_sub["population"]))
        df["population"] = df[leiden_col].astype(int).map(leiden_to_pop)
        df = df.dropna(subset=["population"])
        inference = set(POOLS[anchor]["inference"])
        sub = df[df["population"].isin(inference)]
        ts(f"  d2g[{anchor}]: {len(sub):,} inference-pool cells")
        for r in sub.itertuples():
            out[(str(r.genotype), str(r.library), str(r.barcode))] = anchor
    return out


def main() -> None:
    ts("=== Phase H Step 3: gene detection-fraction export ===")
    if not os.path.isfile(H5AD):
        sys.exit(f"ERROR: {H5AD} missing")
    ts("Loading DEG-substrate h5ad ...")
    adata = sc.read_h5ad(H5AD)
    ts(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    for col in ("unified_genotype_id", "library", "scDblFinder_class"):
        if col not in adata.obs.columns:
            sys.exit(f"ERROR: obs missing column '{col}'. "
                     f"Available: {sorted(adata.obs.columns.tolist())}")

    # Mask doublets — no subset/copy (keeps peak RAM at one h5ad load).
    obs = adata.obs
    singlet = (obs["scDblFinder_class"].astype(str) != "doublet").to_numpy()
    ts(f"  {int((~singlet).sum()):,} doublets masked; "
       f"{int(singlet.sum()):,} singlets in scope")

    missing = [gid for _, gid in GENES if gid not in adata.var_names]
    if missing:
        sys.exit(f"ERROR: gene id(s) not in var_names: {missing}")

    bc = (obs["barcode"].astype(str) if "barcode" in obs.columns
          else obs.index.astype(str))
    stage_map = load_e1b_stage_map(E1B_TSV)
    stage_anchor = (obs["unified_genotype_id"].astype(str).map(stage_map)
                    .map(STAGE_TO_ANCHOR))
    n_no_stage = int((stage_anchor.isna() & singlet).sum())
    ts(f"  {n_no_stage:,} singlet nuclei without an e1b stage anchor "
       f"(kept in 'all' rows, absent from per-stage rows)")

    pool_cells = load_inference_pool_cells()
    keys = list(zip(obs["unified_genotype_id"].astype(str),
                    obs["library"].astype(str), bc))
    pool_anchor = pd.Series([pool_cells.get(k) for k in keys],
                            index=obs.index, dtype=object)
    n_pool = int((pool_anchor.notna().to_numpy() & singlet).sum())
    ts(f"  {n_pool:,} singlet nuclei in the inference pool")
    if n_pool == 0:
        sys.exit("ERROR: zero inference-pool nuclei — key join failed")

    X = adata[:, [gid for _, gid in GENES]].X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X)

    scope_masks = {
        "all_nuclei": {"all": singlet,
                       **{a: singlet & (stage_anchor == a).to_numpy()
                          for a in STAGE_ANCHORS}},
        "inference_pool": {"all": singlet & pool_anchor.notna().to_numpy(),
                           **{a: singlet & (pool_anchor == a).to_numpy()
                              for a in STAGE_ANCHORS}},
    }

    rows = []
    for scope, stages in scope_masks.items():
        for stage_label, mask in stages.items():
            n = int(mask.sum())
            if n == 0:
                sys.exit(f"ERROR: empty scope×stage bucket: "
                         f"{scope} × {stage_label}")
            for (label, gid), col in zip(GENES, X[mask].T):
                rows.append({
                    "scope": scope,
                    "stage": stage_label,
                    "gene_label": label,
                    "gene_id": gid,
                    "n_nuclei": n,
                    "n_detected": int((col > 0).sum()),
                    "frac_detected": float((col > 0).mean()),
                    "total_umi": float(col.sum()),
                })

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    ts(f"Wrote {len(out)} rows -> {OUT_TSV}")
    nod = out[(out.gene_label == "NODAL") & (out.scope == "inference_pool")
              & (out.stage == "all")].iloc[0]
    ts(f"  NODAL inference-pool overall: {nod.n_detected}/{nod.n_nuclei} "
       f"nuclei = {100 * nod.frac_detected:.3f}%")


if __name__ == "__main__":
    main()
