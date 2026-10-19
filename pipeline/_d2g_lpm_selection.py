"""Shared LPM cell-selection logic for d2g-family pipeline steps.

Three pipeline scripts (d2g_stage_restricted_subclustering.py,
d2g_signaling_export.py, d2e_silhouette_export.py) share the same three
primitives: resolve LPM cluster IDs from the authoritative annotation TSV,
fold HH16+HH17 into an HH17/16 bucket, and build the LPM × stage mask.

Keeping them in one place ensures renames in `cluster_annotations.tsv`
propagate uniformly, and eliminates a latent NameError in the
`f"{SCRATCH}/references/..."` fallback branch that used to live in each
caller (SCRATCH wasn't always defined).
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
from pathlib import Path

import pandas as pd

# Annotation names that compose the "LPM" subset at res=0.2.  Both
# populations co-cluster at res=0.2 and are what d2g / d2e drill into.
LPM_ANNOTATIONS = frozenset({"LPM mesenchyme", "Splanchnic mesenchyme"})

# Default search order for references/cluster_annotations.tsv.
# Callers inside the repo hit the first path; cluster jobs whose repo
# layout is mirrored at /n/scratch fall through to the second.
DEFAULT_ANN_CANDIDATES: tuple[str, ...] = (
    # Repo-relative (relative to pipeline/)
    str(Path(__file__).resolve().parent.parent
        / "references" / "cluster_annotations.tsv"),
    # Scratch mirror used by cluster jobs when the repo path isn't present
    f"{SCRATCH_ROOT}/references/cluster_annotations.tsv",
)


def load_cluster_annotations(
    candidates: tuple[str, ...] = DEFAULT_ANN_CANDIDATES,
) -> dict[str, str]:
    """Read {cluster_id_str: annotation} from cluster_annotations.tsv.

    Checks `candidates` in order, loading the first that exists.  Raises
    FileNotFoundError (not the previous silent-fallback-to-a-missing-path
    NameError) if none are readable.
    """
    for path in candidates:
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            return dict(zip(df["cluster"].astype(str), df["annotation"]))
    raise FileNotFoundError(
        f"cluster_annotations.tsv not found at any of: {list(candidates)}. "
        f"Run pipeline/d2b_validate_clusters.py first, or rsync "
        f"references/ to scratch before launching this step."
    )


def lpm_cluster_ids(
    annotations: dict[str, str] | None = None,
) -> list[int]:
    """Return sorted list of integer cluster IDs whose annotation is LPM.

    `annotations` is the dict returned by `load_cluster_annotations`;
    pass None to load it implicitly from the default search path.
    """
    if annotations is None:
        annotations = load_cluster_annotations()
    ids = [int(cl) for cl, name in annotations.items()
           if name in LPM_ANNOTATIONS]
    if not ids:
        raise ValueError(
            f"No clusters in cluster_annotations.tsv match "
            f"{sorted(LPM_ANNOTATIONS)}.  Current alphabet: "
            f"{sorted(set(annotations.values()))!r}.  Update "
            f"LPM_ANNOTATIONS in pipeline/_d2g_lpm_selection.py if "
            f"d2b renamed these populations."
        )
    return sorted(ids)


def fold_hh17_hh16(adata, stage_col: str = "stage") -> None:
    """Fold HH16 + HH17 into HH17/16, in place.

    Both stages are from the same pre-TES developmental window and are
    co-clustered in d2g.  The design file preserves the HH16/HH17
    distinction (fig01A is wetlab ground truth), but downstream
    subclustering treats them as one bucket.
    """
    for s in ("HH16", "HH17"):
        adata.obs.loc[adata.obs[stage_col] == s, stage_col] = "HH17/16"
