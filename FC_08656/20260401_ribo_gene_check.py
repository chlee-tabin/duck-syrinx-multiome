#!/usr/bin/env python3
"""Check ribosomal gene coverage in duck h5ad var_names.

GO/NO-GO gate for pct_ribo panel: if <80 ribosomal protein genes found,
defer the panel.

Expected: ~130 ribosomal protein genes in vertebrates
  - ~80 RPL* (large subunit: RPL3-RPL41, some with A/B suffixes)
  - ~33 RPS* (small subunit: RPS2-RPS30, some with A suffixes)
  - Possibly LOC* gene IDs for unresolved orthologs

Input:  phase2/step04_preCB_clustered.h5ad (only reads var_names)
Output: phase2/figures/step04b_ribo_gene_audit.tsv (gene list + classification)

Run via sbatch:
    sbatch --mem=4G -t 0-0:05 -p short pipeline/run_fig02_phase1.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import scanpy as sc
import pandas as pd
import re
import os
import time


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Known vertebrate ribosomal protein genes (canonical names)
EXPECTED_RPL = [
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9",
    "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14",
    "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22",
    "RPL22L1", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27",
    "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34",
    "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL37", "RPL37A", "RPL38",
    "RPL39", "RPL40", "RPL41", "RPLP0", "RPLP1", "RPLP2",
]
EXPECTED_RPS = [
    "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8",
    "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15",
    "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21",
    "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28",
    "RPS29", "RPS30", "RPSA",
]
EXPECTED_ALL = set(EXPECTED_RPL + EXPECTED_RPS)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Ribosomal gene orthology check ===")

    ts("Loading h5ad (var_names only)...")
    adata = sc.read_h5ad(f"{PHASE2}/step04_preCB_clustered.h5ad", backed="r")
    var_names = list(adata.raw.var_names) if adata.raw is not None else list(adata.var_names)
    ts(f"  {len(var_names):,} genes in var_names")

    var_set = set(var_names)

    # 1. Direct RPS*/RPL* matches
    rps_found = sorted([g for g in var_names if g.startswith("RPS")])
    rpl_found = sorted([g for g in var_names if g.startswith("RPL")])
    ribo_direct = rps_found + rpl_found

    ts(f"  RPS* genes found: {len(rps_found)}")
    ts(f"  RPL* genes found: {len(rpl_found)}")
    ts(f"  Total direct matches: {len(ribo_direct)}")

    # 2. Check against expected canonical set
    found_canonical = [g for g in EXPECTED_ALL if g in var_set]
    missing_canonical = [g for g in EXPECTED_ALL if g not in var_set]
    extra = [g for g in ribo_direct if g not in EXPECTED_ALL]

    ts(f"  Canonical hits: {len(found_canonical)}/{len(EXPECTED_ALL)}")
    if missing_canonical:
        ts(f"  Missing canonical: {', '.join(sorted(missing_canonical)[:20])}")
    if extra:
        ts(f"  Extra (non-canonical): {', '.join(extra[:20])}")

    # 3. Check for LOC* genes that might be ribosomal
    loc_genes = [g for g in var_names if g.startswith("LOC")]
    ts(f"  LOC* genes total: {len(loc_genes)}")

    # Build audit table
    rows = []
    for g in ribo_direct:
        rows.append({
            "gene": g,
            "subunit": "large" if g.startswith("RPL") else "small",
            "canonical": g in EXPECTED_ALL,
            "source": "direct_match",
        })
    for g in missing_canonical:
        rows.append({
            "gene": g,
            "subunit": "large" if g.startswith("RPL") else "small",
            "canonical": True,
            "source": "missing",
        })

    audit = pd.DataFrame(rows)
    out = f"{FIGDIR}/step04b_ribo_gene_audit.tsv"
    audit.to_csv(out, sep="\t", index=False)
    ts(f"  Saved audit: {out}")

    # GO/NO-GO decision
    n_ribo = len(ribo_direct)
    threshold = 80
    if n_ribo >= threshold:
        ts(f"\n  GO: {n_ribo} ribosomal genes found (>= {threshold} threshold)")
        ts(f"  pct_ribo panel is viable")
    else:
        ts(f"\n  NO-GO: only {n_ribo} ribosomal genes found (< {threshold} threshold)")
        ts(f"  Defer pct_ribo panel — check LOC* genes or duck genome annotation")

    ts("=== Ribo gene check complete ===")


if __name__ == "__main__":
    main()
