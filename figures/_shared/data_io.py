"""Data loading utilities and path constants."""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import os
import pandas as pd

# ── Scratch paths ────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
PHASE2_DIR = f"{SCRATCH}/phase2"

# Figures data: the committed figures/data/ (single path; no scratch fallback)
_LOCAL_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
if not os.path.isdir(_LOCAL_DATA_DIR):
    raise FileNotFoundError(f"figures/data/ not found at {os.path.abspath(_LOCAL_DATA_DIR)}")
FIGURES_DATA_DIR = os.path.abspath(_LOCAL_DATA_DIR)
W_GENES_FILE = f"{PHASE2_DIR}/20220910_Wgenes.tsv"
Z_GENES_FILE = f"{PHASE2_DIR}/20220910_Zgenes.tsv"

# Demux result paths
LIBPAIR_DONOR_IDS = f"{SCRATCH}/libpair_g6/reference_vireo/result_unified/donor_ids.tsv"
G9_DONOR_IDS = f"{SCRATCH}/dp500_vireo_results/result_A_unfilt/donor_ids.tsv"
# Old G2 vireo: multi-chunk result, best-seed from QC.
# For the unified cell-level result, use cb_genotype_map_v2.tsv.gz instead.
# This path points to chunk_0 best seed for single-chunk analysis scripts.
G2_DONOR_IDS = (
    f"{SCRATCH}/GEX03-1A_GEX/outs/pseudobulk_group2_merged/"
    "vireo_output/n_donor_16_40000/chunk_0_overlap50/seed_600914/donor_ids.tsv"
)
CROSS_ANCHOR_DIR = f"{SCRATCH}/cross_anchor_output"
GENOTYPE_DISTANCE_DIR = f"{SCRATCH}/genotype_distance_output"


def ts(msg):
    """Timestamped log message."""
    import time
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_tsv(filename, subdir=None, **kwargs):
    """Load a TSV file from the figures data directory."""
    if subdir:
        path = os.path.join(subdir, filename)
    else:
        path = os.path.join(FIGURES_DATA_DIR, filename)
    return pd.read_csv(path, sep="\t", **kwargs)


def load_tsv_gz(filename, subdir=None, **kwargs):
    """Load a gzipped TSV file."""
    return load_tsv(filename, subdir=subdir, compression="gzip", **kwargs)


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)
    return path


