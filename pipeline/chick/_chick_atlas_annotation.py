"""Shared cell-type annotation for the chick whole-embryo atlases.

Single source of truth for the marker panel + annotation method used by BOTH
`cAE_full_atlas.py` (pooled Group A+B, figS04 featured) and the Group-B-only
robustness atlas (not distributed) — so the two atlases label identically.

Method (mirrors duck d2b scoring): per-cluster score = mean fraction-expressing across a
cell type's markers; ARGMAX per cluster (NOT the duck greedy one-per-type, which mislabels
the chick LPM blob "Splanchnic"); <0.01 floor -> Unassigned; a "Low-quality" override for
ribo/mito-dominated clusters; and a Sclerotome margin-promotion (surfaces the sclerotome that
the broad "Mesenchyme (early)" type otherwise out-scores). Detail: docs/chick_control/ANNOTATION_DIFF.md (not distributed).
"""
import numpy as np
import pandas as pd

# ── Cell-type marker panel (duck d2b ortholog set + chick-specific "Mesenchyme (early)"). ──
# ARGMAX assignment (per cluster, independent) — dict order does NOT matter. The chick res=0.2
# mesenchyme splits ~99% by STAGE, not sublineage: an early/broad pool ("Mesenchyme (early)";
# EBF1/2/SOX5/BNC2/CDH11 + DCN, low TBX5/TBX4) and an HH25 TBX5/TBX4+ specified "LPM mesenchyme".
# "Splanchnic mesenchyme" is kept but the chick LPM does not separate splanchnic, so it should not
# win argmax on the LPM blob.
CELLTYPE_MARKERS = {
    "RBC": ["HBZ", "HBA1", "LOC101798111", "FTH1"],
    "Notochord": ["COL8A1", "SHH", "KIAA1217"],
    "Cardiomyocyte": ["RYR2", "RBM20", "GATA5"],
    "Otic placode": ["AGR3", "EYA1", "PAX2"],
    "Endothelial": ["FLI1", "EGFL7", "PTPRM", "PECAM1"],
    "Pharyngeal endoderm": ["ADAMTS17", "TFAP2E", "FOXI2"],
    "Ventral neural prog.": ["SPON1", "RFX4", "ARX"],
    "Dorsal neural prog.": ["LMX1B", "SLIT2", "DMD"],
    "Diff. neuron": ["BRSK2", "MYT1L", "KALRN", "ELAVL4"],
    "Dermomyotome": ["CDH4", "TMTC2", "PAX3", "SPOCK3"],
    "Epithelial": ["EPCAM", "KIAA1217", "PTPRF"],
    "Neural crest": ["COL25A1", "TFAP2B"],
    "Neural progenitor": ["PRTG", "NPAS3", "PTPRN2", "SOX2"],
    "Mesenchyme (early)": ["EBF1", "EBF2", "CDH11", "SOX5", "BNC2"],
    "Sclerotome": ["PAX1", "PAX9", "FOXC2", "MEOX1"],
    "Splanchnic mesenchyme": ["FOXF1", "GATA4", "ZFPM2", "ISL1"],
    "LPM mesenchyme": ["KIF26B", "CACNA1C", "TBX5", "TBX4", "NFIA"],
}

# Paraxial / presomitic mesoderm probe (score_genes) — orthogonal corroboration for the
# Sclerotome promotion; NOT in the argmax.
PARAXIAL_PROBE = ["TBX6", "MSGN1", "DLL1", "DLL3", "HES7", "LFNG", "RIPPLY2",
                  "TCF15", "MEOX1", "MEOX2", "FOXC1", "FOXC2", "ALDH1A2", "NKX1-2", "FGF8"]

# Chick mitochondrial genes (galGal6a; exact set avoids matching nuclear NDUF*/COX-assembly).
MITO_EXACT = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
              "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB",
              "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
              "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB"}

LOWQ_TOP10_MIN = 6        # >= this many of a cluster's top-10 DE genes are ribo/mito -> Low-quality
SCLERO_PROMOTE_MIN = 0.30  # mean frac-expressing of PAX1/PAX9/FOXC2/MEOX1 to promote -> Sclerotome
PARAXIAL_PROMOTE_MIN = 0.05  # mean paraxial/PSM score_genes probe (orthogonal corroboration)


def _is_ribo(g):
    return str(g).upper().startswith(("RPL", "RPS", "MRPL", "MRPS"))


def annotate_clusters(adata, leiden_key, lowq=frozenset()):
    """ARGMAX per-cluster annotation. Each cluster gets the cell type whose markers it expresses
    most (mean fraction-expressing across that type's markers). <0.01 floor -> Unassigned_<cl>.
    Clusters in `lowq` -> "Low-quality". Returns (assignments {cl->celltype}, score_mat)."""
    raw_var = adata.raw.var_names if adata.raw is not None else adata.var_names
    clusters = sorted(adata.obs[leiden_key].unique(), key=lambda x: int(x))
    scores = {}
    for ct, markers in CELLTYPE_MARKERS.items():
        present = [m for m in markers if m in raw_var]
        if not present:
            print(f"  WARNING: no markers found for {ct}", flush=True)
            continue
        midx = [np.where(raw_var == m)[0][0] for m in present]
        for cl in clusters:
            idx = np.where((adata.obs[leiden_key] == cl).values)[0]
            x = adata.raw.X[idx][:, midx] if adata.raw is not None else adata.X[idx][:, midx]
            if hasattr(x, "toarray"):
                x = x.toarray()
            scores[(str(cl), ct)] = float((x > 0).mean(axis=0).mean())
    score_mat = pd.DataFrame(
        {ct: {str(cl): scores.get((str(cl), ct), 0.0) for cl in clusters}
         for ct in CELLTYPE_MARKERS}
    )
    assignments = {}
    for cl in clusters:
        cl = str(cl)
        if cl in lowq:
            assignments[cl] = "Low-quality"
            continue
        row = score_mat.loc[cl]
        best_ct, best_score = row.idxmax(), float(row.max())
        assignments[cl] = best_ct if best_score >= 0.01 else f"Unassigned_{cl}"
    return assignments, score_mat


def lowq_from_de(rgg, top_n=10, min_hk=LOWQ_TOP10_MIN):
    """Return the set of cluster ids whose top-`top_n` rank_genes_groups DE genes are
    dominated (>= min_hk) by ribosomal/mito genes (no real cell-type signal)."""
    lowq = set()
    if rgg is None:
        return lowq
    for cl in rgg["names"].dtype.names:
        top = [rgg["names"][cl][r] for r in range(min(top_n, len(rgg["names"][cl])))]
        n_hk = sum(1 for g in top if _is_ribo(g) or str(g).upper() in MITO_EXACT)
        if n_hk >= min_hk:
            lowq.add(str(cl))
    return lowq


def promote_sclerotome(assignments, score_mat, paraxial_by_cl, log=print):
    """In place: relabel a non-paraxial cluster to Sclerotome when its sclerotome-marker score is
    genuinely high (>= SCLERO_PROMOTE_MIN) AND its paraxial-mesoderm probe is positive
    (> PARAXIAL_PROMOTE_MIN) — argmax + the broad "Mesenchyme (early)" type otherwise suppress it.
    Resolution-dependent by construction (fires where sclerotome concentrates). Returns n promoted."""
    n = 0
    for cl in list(assignments):
        if assignments[cl] in ("Dermomyotome", "Sclerotome", "Low-quality"):
            continue
        if (score_mat.loc[cl, "Sclerotome"] >= SCLERO_PROMOTE_MIN
                and float(paraxial_by_cl.get(cl, 0.0)) > PARAXIAL_PROMOTE_MIN):
            log(f"  [promote] cl{cl} {assignments[cl]!r} -> Sclerotome "
                f"(score={score_mat.loc[cl, 'Sclerotome']:.3f}, paraxial={paraxial_by_cl[cl]:.3f})")
            assignments[cl] = "Sclerotome"
            n += 1
    return n
