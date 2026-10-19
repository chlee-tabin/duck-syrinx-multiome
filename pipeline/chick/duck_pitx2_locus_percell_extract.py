#!/usr/bin/env python3
"""(i) Duck PITX2-locus per-cell ArchR-peak counts → chick schema, for the cross-species
orientation-free half-vs-half (duck + chick in the SAME metric).

Reads the already-built duck ArchR-peak cell×peak matrix
(`phase2/atac/step_f2_peak_matrix_archr.h5ad` — the 283K ArchR union peaks projected onto
f1's per-cell AnnData, SAME peak-set family as the chick build, NOT SnapATAC2-sparse), so no
ArchR rebuild is needed.

Why this is a FAIR positive control (non-circular for magnitude):
  g1b_pitx2_side_inference assigns `inferred_side` per (genotype, library): top_library→L,
  bot_library→R. The two libraries are the embryo's two PHYSICAL halves — the split is
  physical, independent of PITX2. Only the L/R *label* is PITX2-derived. The orientation-free
  half-vs-half discards the label (compares the two halves' magnitudes), so the divergence it
  measures is non-circular. We tag each cell `cell = "<library_m>#<barcode>"` so the chick
  `orientation_free_halfvshalf` (lane = the half) runs UNCHANGED.

Side join: firmness∈{direct, placement} + pair_constraint_ok, the production g1/g3 firm set,
via the shared rule in figures/_shared/lr_firmness.py.  Same-side duplicates (lane-loading-record tube
group 17-B(2) = embryos 17-5/17-6: one LEFT half sequenced in GEX04 AND GEX0S-duck) are kept once on their design
side with both libraries POOLED into one lane `GEX04+GEX0S-duck` -- so the orientation-free
half-vs-half (lane = half) can never pair the two aliquots -- and are not gated by
pair_constraint_ok (a property of INFERRED sides), the same rule g1 applies.

Output (written directly to $SCRATCH/figures/data/duck_locus/ — the
git-synced tree, shipped by derived_data_manifest.tsv; chick schema):
  locus_peak_coords.tsv         peak_id, chrom, start, end, in_promoter_c, in_gene_tss
  locus_percell_counts.tsv.gz   cell, stage, population, side, group, genotype, nFrags, <peaks>

Run: sbatch pipeline/chick/duck_pitx2_locus_percell_extract.sh
"""
from __future__ import annotations
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")
SCRATCH_ROOT = (_os.environ.get("DUCK_SCRATCH_ROOT") or f"{SCRATCH_USER}/FC_08656")


import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import snapatac2 as snap

# Resolve the figures/ helper dir (pipeline/chick/ -> repo/figures) so the canonical
# d2g subcluster SSOT loader imports cleanly on scratch and in dev.  The sync script
# rsyncs figures/_shared/ to scratch, so it is reachable in either context.
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures"))
from _shared.subcluster_map import get_domain_map  # noqa: E402
from _shared.lr_firmness import resolve_sides  # noqa: E402

SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
# F2_H5AD / OUT_DIR default to the production peak matrix + duck_locus sidecar; override both to
# re-extract from an alternative reprojected matrix (e.g. the genotype-reproducibility build)
# into an isolated staging dir without touching the committed sidecar.
F2_H5AD = os.environ.get("F2_H5AD", f"{ATAC_DIR}/step_f2_peak_matrix_archr.h5ad")
# f1 per-cell metadata sidecar carries the leiden id (d2g_cluster) needed to re-map the
# STALE baked d2g_subcluster NAME in the f2 h5ad through the current SSOT (see main()).
F1_META_TSV = os.environ.get("F1_META_TSV", f"{ATAC_DIR}/step_f1_coupled_atac_metadata.tsv")
# d2l RAW snRNA L/R sidecar.  This extract is PLUMBING, not a laterality conclusion — it
# emits per-cell counts plus whatever `side` the sidecar carries — so the laterality-sidecar routing rule's "may
# read either sidecar" applies.  ⚠ The emitted `side` column is therefore RAW-d2L and is
# STALE AT HH17/16, where PITX2 RNA sits at the detection floor and 4 genotypes flip
# against the chromatin (d2n) call.  Any consumer drawing a DIRECTIONAL conclusion MUST
# re-map HH17/16 first via per_stage_laterality.apply_d2n_hh1716_override — the ONE
# shared mechanism (gated d2n map, "d2n direct" rather than the d2o merge;
# docs/laterality_sidecar_routing.md §3, not distributed).  (The cross-species
# half-vs-half is safe without that: it discards the L/R label entirely and reports
# orientation-free magnitude.)
# Read the COMMITTED figures/data/ copy — d2l's actual output path
# (pipeline/d2l_lr_inference.py) + g1's convention — not the purge-prone, now-absent
# phase2/figures/ scratch path.
D2L_TSV = os.environ.get("D2L_TSV", f"{SCRATCH}/figures/data/d2l_lr_inference.tsv")
# 2026-08-10: write to figures/data/duck_locus (the git-synced tree the manifest
# ships from), not the purge-prone phase2/figures — the INPUT paths moved off
# phase2/figures earlier but the OUTPUT did not, so every regeneration had to be
# hand-ferried and the manifest had no duck_locus row at all.
OUT_DIR = os.environ.get("LOCUS_OUT_DIR", f"{SCRATCH}/figures/data/duck_locus")
os.makedirs(OUT_DIR, exist_ok=True)

# duck PITX2 locus (references/pitx2_promoter_regions.bed; ZJU1.0 NC_051775.1, + strand)
CHROM = "NC_051775.1"
GENE_START, GENE_END = 3_961_386, 3_976_254
LOCUS_START, LOCUS_END = GENE_START - 500_000, GENE_END + 500_000   # ±500 kb
PROMC_START, PROMC_END = 3_971_084, 3_971_685     # promoter|PITX2c (TSS 3,971,385 ±300)
GENE_TSS = 3_961_387                              # PITX2a TSS = gene start (chick in_gene_tss convention)

ACCEPTED_FIRMNESS = {"direct", "placement"}
LIBRARY_COLLAPSE = {"GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1",
                    "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4"}
ARC_TO_GEX = {"PILOTARC3": "PILOTGEX3", "PILOTARC4A": "PILOTGEX4A", "PILOTARC4B": "PILOTGEX4B",
              "ARC05": "GEX05", "ARC06": "GEX06", "ARC07": "GEX07", "ARC08": "GEX08",
              "ARC03-1A": "GEX03-1A", "ARC03-1B": "GEX03-1B", "ARC03-2": "GEX03-2",
              "ARC04": "GEX04", "ARC0S-duck": "GEX0S-duck"}


def ts(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def collapse(r):
    return r if pd.isna(r) else LIBRARY_COLLAPSE.get(r, r)


def parse_peak(p):
    """'<chrom>:<start>-<end>' -> (chrom, start, end).  Every var name of the ArchR peak
    matrix has this form; one that does not is a wrong input, so raise and name it
    (this used to return None and skip the peak without a word)."""
    try:
        c, se = str(p).rsplit(":", 1)
        s, e = se.split("-")
        return c, int(s), int(e)
    except ValueError as exc:
        raise ValueError(f"peak name {p!r} is not '<chrom>:<start>-<end>'") from exc


def snap_obs_to_pandas(obs, required, optional=()):
    """SnapATAC2 2.9 obs is Polars-backed; pull columns via the iterable protocol.

    A REQUIRED column that cannot be read raises, naming it.  An OPTIONAL one is left out
    and logged with the error type (main() handles each optional column's absence
    explicitly).  This used to swallow every exception for every column and drop it
    without a word.  (The Rust-side PyDataFrameElem's missing-column exception type is
    not documented, hence the broad catch; it is re-raised or logged, never silent.)"""
    d = {}
    for c in list(required) + list(optional):
        try:
            d[c] = list(obs[c])
        except Exception as exc:  # noqa: BLE001 -- re-raised (required) or logged (optional)
            if c not in optional:
                raise KeyError(f"peak-matrix obs: cannot read required column {c!r} "
                               f"({type(exc).__name__}: {exc})") from exc
            ts(f"  obs has no optional column {c!r} ({type(exc).__name__}); handled below")
    return pd.DataFrame(d)


def main():
    if not os.path.isfile(F2_H5AD):
        sys.exit(f"ERROR: duck ArchR peak matrix missing at {F2_H5AD}")
    ts(f"reading {F2_H5AD}")
    adata = snap.read(F2_H5AD)
    ts(f"  {adata.n_obs:,} cells × {adata.n_vars:,} peaks")
    var = list(adata.var_names)

    locus_idx, coords = [], []
    for i, p in enumerate(var):
        pp = parse_peak(p)
        if pp is None:
            continue
        c, s, e = pp
        if c == CHROM and s < LOCUS_END and e > LOCUS_START:
            locus_idx.append(i)
            coords.append({"peak_id": p, "chrom": c, "start": s, "end": e,
                           "in_promoter_c": (s < PROMC_END and e > PROMC_START),
                           "in_gene_tss": (s <= GENE_TSS < e)})
    ts(f"  locus peaks (chr {CHROM} ±500 kb): {len(locus_idx)} "
       f"(promoter_c={sum(c['in_promoter_c'] for c in coords)}, "
       f"gene_tss={sum(c['in_gene_tss'] for c in coords)})")
    if not locus_idx:
        sys.exit("no peaks in PITX2 ±500 kb window")
    coords_df = pd.DataFrame(coords)

    # gex_library is optional (derived from arc_library below when absent); "stage" is not
    # used downstream and is no longer read.
    obs = snap_obs_to_pandas(
        adata.obs,
        required=["arc_library", "unified_genotype_id", "d2g_subcluster", "d2g_stage"],
        optional=["gex_library"])
    obs["cell_idx"] = np.arange(len(obs))
    # obs_names IS the cell barcode and is load-bearing (it builds the per-cell key below);
    # fail loud if it can't be read rather than fabricating synthetic 'cellN' barcodes.
    obs["barcode"] = list(adata.obs_names)
    if "gex_library" not in obs.columns or obs["gex_library"].isna().all():
        if "arc_library" in obs.columns:
            obs["gex_library"] = obs["arc_library"].map(ARC_TO_GEX)
    obs["library_m"] = obs["gex_library"].map(collapse)

    # ── Re-derive d2g_subcluster from the CANONICAL SSOT via the f1 sidecar ──────
    # The d2g_subcluster NAME baked into the f2 ArchR peak-matrix obs is written at
    # f1-time from references/d2g_subcluster_annotations.tsv and goes STALE when that
    # TSV is re-labelled without a fragment-level f1 rebuild (the 2026-07 HH25
    # Chondrogenic LPM→Syringeal / cl5 Syringeal→Ventral spl. vet, on top of the older
    # Somatic→Skeletogenic LPM and Splanchnic→Venous-pole spl. renames).  The leiden id
    # needed to re-map is NOT in the peak matrix, but IS in f1's per-cell sidecar
    # (step_f1_coupled_atac_metadata.tsv: arc_library, cb16, d2g_stage, d2g_cluster).
    # Join that leiden by (arc_library, cb16) and re-map leiden -> fresh label via
    # get_domain_map (f1's own rule) so an annotation change flows into this artifact
    # WITHOUT a rebuild.  Mirrors pipeline/g1_pseudobulk_atac.py exactly.  Fail loud on
    # a coverage gap — never silently drop.
    if not os.path.isfile(F1_META_TSV):
        sys.exit(f"ERROR: f1 metadata sidecar missing at {F1_META_TSV}; "
                 f"restore from standby (atac_backup/atac/) to re-derive labels.")
    if "arc_library" not in obs.columns:
        sys.exit("ERROR: obs missing 'arc_library'; cannot join the leiden "
                 "sidecar by (arc_library, cb16).")

    def _leiden_str(v):
        if pd.isna(v):
            return None
        try:
            return str(int(v))       # int / numpy-int / float-integer / "8"
        except (ValueError, TypeError):
            return str(v)

    # cb16 = the 16 bp barcode (obs_name minus the "-1" GEM suffix).
    # rsplit, matching f1's cb16 derivation — first-hyphen split diverges on a
    # double-suffixed barcode (obs_names_make_unique can emit "...-1-1")
    obs["cb16"] = [str(b).rsplit("-", 1)[0] if "-" in str(b) else str(b)
                   for b in obs["barcode"]]

    _meta = pd.read_csv(
        F1_META_TSV, sep="\t",
        usecols=lambda c: c in {"arc_library", "cb16", "d2g_stage", "d2g_cluster"},
        dtype={"arc_library": str, "cb16": str, "d2g_stage": str},
    )
    # `usecols=lambda` silently tolerates a missing column (a list would raise naming it),
    # so check explicitly — a restored/partial sidecar otherwise dies as a bare KeyError.
    _need = {"arc_library", "cb16", "d2g_stage", "d2g_cluster"}
    if not _need <= set(_meta.columns):
        sys.exit(f"ERROR: {F1_META_TSV} missing column(s) {sorted(_need - set(_meta.columns))}; "
                 f"restore the f1 sidecar from standby (atac_backup/atac/) or re-run f1.")
    _meta = _meta[_meta["d2g_cluster"].notna() & _meta["d2g_stage"].notna()].copy()
    if _meta.empty:
        sys.exit(f"ERROR: {F1_META_TSV} has no rows with both d2g_stage and d2g_cluster set; "
                 f"the sidecar carries no LPM assignments. Re-run f1.")
    # (arc_library, cb16) is the join key; a duplicate would let an arbitrary last-written
    # label win silently. f1 itself enforces this (f1_coupled_atac.py drop_duplicates on
    # ["library","cb16"]), so a duplicate here means a concatenated/partial sidecar.
    _dup = int(_meta.duplicated(["arc_library", "cb16"]).sum())
    if _dup:
        sys.exit(f"ERROR: {F1_META_TSV} has {_dup:,} duplicate (arc_library, cb16) rows — the "
                 f"join key is not unique, so labels would be assigned arbitrarily. The sidecar "
                 f"is likely a concatenation of multiple f1 runs; regenerate it.")
    _fresh = []
    for _st, _g in _meta.groupby("d2g_stage", dropna=True):
        _dm = get_domain_map(_st)              # {str(cluster_id): population}
        _fresh.append(_g["d2g_cluster"].map(_leiden_str).map(_dm))
    _meta["fresh_label"] = pd.concat(_fresh).reindex(_meta.index)
    # An unmapped (stage, leiden) means the sidecar's d2g run and the current
    # annotations TSV disagree — a whole subcluster would vanish from the artifact.
    # Raise naming the pairs (project convention: never silently substitute/drop on a failed
    # lookup); same discipline as archr_metadata_sidecar._cal_pop.
    _unmapped = _meta.loc[_meta["fresh_label"].isna(), ["d2g_stage", "d2g_cluster"]]
    if len(_unmapped):
        _pairs = sorted({(s, _leiden_str(c)) for s, c in
                         zip(_unmapped["d2g_stage"], _unmapped["d2g_cluster"])})
        sys.exit(
            f"ERROR: {len(_unmapped):,} f1-sidecar cells carry a (stage, leiden) with no "
            f"entry in references/d2g_subcluster_annotations.tsv: {_pairs}. Re-run "
            f"pipeline/d2g_validate_subclusters.py, or re-run f1 against the current d2g."
        )
    _lut = {
        (lib, cb): lab
        for lib, cb, lab in zip(_meta["arc_library"], _meta["cb16"],
                                _meta["fresh_label"])
    }
    ts(f"  f1 sidecar: {len(_lut):,} LPM cells carry a fresh SSOT label")

    _old_label = obs["d2g_subcluster"].astype("object")
    obs["d2g_subcluster"] = [
        _lut.get((lib, cb))
        for lib, cb in zip(obs["arc_library"], obs["cb16"])
    ]
    _n_old = int(_old_label.notna().sum())
    _n_new = int(obs["d2g_subcluster"].notna().sum())
    _n_changed = int(
        (_old_label.fillna("∅").astype(str)
         != obs["d2g_subcluster"].fillna("∅").astype(str)).sum()
    )
    ts(f"  SSOT relabel via sidecar: baked LPM={_n_old:,} fresh LPM={_n_new:,} "
       f"changed={_n_changed:,}")
    # Coverage check on the CELLS THEMSELVES, not on a count ratio.  A ratio compares
    # two different cell sets, so a cell that gains a label (baked NaN → fresh) offsets
    # a cell that loses one 1:1 and the gap passes silently; a global 2% budget also
    # hides the total loss of one ARC library (= one physical embryo-half, which this
    # embryo-paired artifact cannot tolerate); and an all-NaN baked column made the old
    # `_n_old > 0` form skip the check entirely and write an empty artifact.  Every cell
    # in the peak matrix came from f1, so "had a baked label, has no sidecar entry" is
    # exactly zero when the two builds agree.
    _lost = _old_label.notna() & obs["d2g_subcluster"].isna()
    if _lost.any():
        _eg = list(zip(obs.loc[_lost, "arc_library"], obs.loc[_lost, "cb16"]))[:10]
        _by_lib = obs.loc[_lost, "arc_library"].value_counts().to_dict()
        sys.exit(
            f"ERROR: SSOT relabel coverage gap — {int(_lost.sum()):,} cells have a baked "
            f"d2g label but no f1-sidecar entry (by arc_library: {_by_lib}; e.g. {_eg}). "
            f"The sidecar and the peak matrix come from different f1 builds. Re-run f1 to "
            f"regenerate the sidecar."
        )

    d2l = resolve_sides(pd.read_csv(D2L_TSV, sep="\t"), ACCEPTED_FIRMNESS,
                        source=os.path.basename(D2L_TSV), require_pair_constraint=True,
                        log=ts)
    d2l = (d2l[["genotype", "library_m", "inferred_side", "library_key"]]
           .rename(columns={"genotype": "unified_genotype_id"}))
    ts(f"  d2l firm-oriented (genotype,library) rows: {len(d2l)}")

    m = obs.merge(d2l, how="left", on=["unified_genotype_id", "library_m"])
    m = m[m["inferred_side"].notna()
          & m["d2g_subcluster"].notna()
          & (m["d2g_subcluster"] != "ambiguous")].copy()
    ts(f"  sided LPM-lineage cells: {len(m):,} "
       f"({m['unified_genotype_id'].nunique()} genotypes)")

    ts("  materialising peak matrix (X[:]) and subsetting locus peaks...")
    X = adata.X[:]
    rows = m["cell_idx"].to_numpy()
    Xrows = X[rows]
    # Depth denominator: the h5ad obs carries NO total-fragment column, so use total in-peak
    # fragments per cell (row-sum over ALL 283,818 ArchR peaks). This is a valid per-cell depth
    # weight; because the orientation-free metric is a WITHIN-embryo half ratio compared to a
    # WITHIN-embryo shuffle null using the SAME denominator, the exact denominator choice
    # (this vs ArchR nFrags) is second-order — it rescales acc but not the off-vs-on-diagonal
    # divergence. (Chick used ArchR nFrags; FRiP is stable enough that the two track per cell.)
    depth = np.asarray(Xrows.sum(axis=1)).ravel().astype(float)
    sub = Xrows[:, locus_idx]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    ts(f"  count submatrix: {sub.shape}; median total-in-peak depth/cell: {np.median(depth):.0f}")

    # lane (before '#') = the half: library_m, or the pooled label of a same-side
    # duplicate.  A pooled lane keeps the source library in the barcode part so the
    # two aliquots' barcodes cannot collide.
    lane = m["library_key"].astype(str)
    bc = m["barcode"].astype(str)
    pooled = lane != m["library_m"].astype(str)
    cell_name = (lane + "#" + bc.where(~pooled, m["library_m"].astype(str) + ":" + bc)).to_numpy()
    if pd.Index(cell_name).duplicated().any():
        raise RuntimeError("duplicate duck cell names in the locus extract: "
                           f"{pd.Index(cell_name)[pd.Index(cell_name).duplicated()][:5].tolist()}")
    out = pd.DataFrame({"cell": cell_name,
                        "stage": m["d2g_stage"].to_numpy(),
                        "population": m["d2g_subcluster"].to_numpy(),
                        "side": m["inferred_side"].to_numpy(),
                        "group": "duck",
                        "genotype": m["unified_genotype_id"].to_numpy(),
                        "nFrags": depth})
    peak_ids = [var[i] for i in locus_idx]
    cnt = pd.DataFrame(sub, columns=peak_ids, index=out.index)
    out = pd.concat([out, cnt], axis=1)

    coords_df.to_csv(f"{OUT_DIR}/locus_peak_coords.tsv", sep="\t", index=False)
    out.to_csv(f"{OUT_DIR}/locus_percell_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    ts(f"  wrote {OUT_DIR}/locus_percell_counts.tsv.gz ({len(out):,} cells × {len(peak_ids)} peaks)")

    pv = (out.groupby(["stage", "side"])
          .agg(n_cells=("cell", "size"), n_geno=("genotype", "nunique"))
          .reset_index())
    ts("design preview (firm-oriented L/R):\n" + pv.to_string(index=False))
    ts("done")


if __name__ == "__main__":
    main()
