#!/usr/bin/env python3
"""Compare the reproduced Fig 5A / Fig 5C inputs with the committed tables.

Run inside the repro root (snapatac2 env, SLURM).  Every check writes one row to
<out>/summary.tsv (check, status PASS/DIFF/MISSING, detail); detail lines go to stdout.
Exits 1 unless every check is PASS, so the job's sacct STATE is the verdict.
Cells are keyed on (arc_library, barcode): a barcode alone recurs across libraries.

Fig 5A (figures/data/step_f5_pitx2_region_counts.tsv):
  A1  f5 --from-cached-matrix (root_cache) vs committed           byte md5 + keyed row diff
  A2  fresh f1 metadata vs the 2026-08-10 metadata (private copy)  cell set, coupled set, labels
  A3  fresh f5 4-region matrix vs the 2026-08-10 cache             barcode-aligned counts
  A4  fresh f5 default path (root) vs committed                    byte md5 + keyed row diff
Fig 5C (figures/data/f4_compartment_perembryo_pitx2_coverage.tsv.gz, split-tube aliquots summed):
  C2  shipped f4 vs the committed export                           every column, exact
  C2n shipped f4 cell counts vs f4_..._cell_counts.tsv             group, n_cells, sums
  C3  shipped f4 HH13 vs f4_hh13_perembryo_pitx2_coverage.tsv.gz   every column, exact
"""
import argparse
import atexit
import hashlib
import os
import sys

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--repro", required=True, help="repro root ($DUCK_REPRO_ROOT)")
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
X = a.repro
R = f"{X}/root"
D = f"{R}/figures/data"
PRIV_META = f"{X}/root_cache/phase2/atac/step_f1_coupled_atac_metadata.tsv"
PRIV_CACHE = f"{X}/root_cache/phase2/atac/step_f5_pitx2_region_matrix_tmp.h5ad"
MANIFEST_MD5 = {PRIV_META: "aec783bd637e96727433b0d26561b8a7",
                PRIV_CACHE: "ab6077658ab52d3c716e78a8aa8e0c76"}
rows = []


@atexit.register
def write_summary():
    """Runs on normal exit AND after an exception, so a late crash keeps the rows so far."""
    pd.DataFrame(rows, columns=["check", "status", "detail"]).to_csv(f"{a.out}/summary.tsv", sep="\t", index=False)
    print(f"wrote {a.out}/summary.tsv ({len(rows)} checks)", flush=True)


def rec(check, status, detail):
    rows.append((check, status, detail))
    print(f"[{status}] {check}: {detail}", flush=True)


def cells(meta):
    """Index a metadata frame by (arc_library, barcode): a barcode alone recurs across libraries."""
    return meta.set_index([meta["arc_library"].astype(str), meta.index.astype(str)])


def cache_keys(meta_path, barcodes):
    """(arc_library, barcode) per cache row, from the metadata TSV the cache is positionally
    aligned to (what _locus_region_counts.check_cached_matrix asserts); None if not aligned."""
    meta = pd.read_csv(meta_path, sep="\t", index_col=0, low_memory=False)
    if [str(b) for b in meta.index] != barcodes:
        return None
    return pd.MultiIndex.from_arrays([meta["arc_library"].astype(str), barcodes])


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


for p, m in MANIFEST_MD5.items():
    got = md5(p)
    rec(f"input_md5 {os.path.basename(p)}", "PASS" if got == m else "DIFF", f"{got} (manifest {m})")

# ---------------- Fig 5A ----------------
F5_COMMITTED = f"{D}/step_f5_pitx2_region_counts.tsv"
KEY = ["genotype", "library_d2l", "side", "stage", "stage_raw", "region"]


def cmp_f5(check, path):
    if not os.path.isfile(path):
        rec(check, "MISSING", path)
        return
    m_new, m_old = md5(path), md5(F5_COMMITTED)
    new = pd.read_csv(path, sep="\t")
    old = pd.read_csv(F5_COMMITTED, sep="\t")
    j = old.merge(new, on=KEY, how="outer", suffixes=("_old", "_new"), indicator=True)
    only_old = int((j["_merge"] == "left_only").sum())
    only_new = int((j["_merge"] == "right_only").sum())
    both = j[j["_merge"] == "both"]
    ndiff = {c: int((both[f"{c}_old"] != both[f"{c}_new"]).sum())
             for c in ("n_cells", "n_frags")}
    fpc = float(np.abs(both["frags_per_cell_old"] - both["frags_per_cell_new"]).max()) if len(both) else float("nan")
    ok = m_new == m_old
    rec(check, "PASS" if ok else "DIFF",
        f"md5 {m_new} vs committed {m_old}; rows {len(new)} vs {len(old)}; only_committed {only_old}, "
        f"only_new {only_new}; differing n_cells {ndiff['n_cells']}, n_frags {ndiff['n_frags']}; "
        f"max |d frags_per_cell| {fpc:.3g}")
    if not ok:
        diff = j[(j["_merge"] != "both") | (j["n_cells_old"] != j["n_cells_new"]) |
                 (j["n_frags_old"] != j["n_frags_new"])]
        diff.to_csv(f"{a.out}/{check}_rowdiff.tsv", sep="\t", index=False)


cmp_f5("A1_f5_cached_matrix_vs_committed", f"{X}/out/f5_cache/step_f5_pitx2_region_counts.tsv")

FRESH_META = f"{R}/phase2/atac/step_f1_coupled_atac_metadata.tsv"
if os.path.isfile(FRESH_META):
    old = cells(pd.read_csv(PRIV_META, sep="\t", index_col=0, low_memory=False))
    new = cells(pd.read_csv(FRESH_META, sep="\t", index_col=0, low_memory=False))
    unique = old.index.is_unique and new.index.is_unique
    same_order = list(old.index) == list(new.index)
    so, sn = set(old.index), set(new.index)
    common = old.index.intersection(new.index)
    det = [f"md5 {md5(FRESH_META)} vs {MANIFEST_MD5[PRIV_META]}", f"(arc_library, barcode) unique {unique}", f"cells {len(new):,} vs {len(old):,}",
           f"same order {same_order}", f"only_old {len(so - sn):,}", f"only_new {len(sn - so):,}"]
    co, cn = old["unified_genotype_id"].notna(), new["unified_genotype_id"].notna()
    det.append(f"coupled {int(cn.sum()):,} vs {int(co.sum()):,}")
    coldiff = {}
    for c in () if not unique else ("unified_genotype_id", "stage", "d2g_stage", "d2g_subcluster", "peak_pool", "cell_type"):
        if c not in old.columns or c not in new.columns:
            coldiff[c] = "missing"  # an absent label column is a difference, not a skip
            continue
        o = old.loc[common, c].astype(str)
        n = new.loc[common, c].astype(str)
        coldiff[c] = int((o != n).sum())
    det.append("differing labels on shared cells " + str(coldiff))
    status = "PASS" if (unique and same_order and so == sn and not any(coldiff.values())) else "DIFF"
    rec("A2_f1_metadata_vs_20260810", status, "; ".join(det))
else:
    rec("A2_f1_metadata_vs_20260810", "MISSING", FRESH_META)

FRESH_CACHE = f"{R}/phase2/atac/step_f5_pitx2_region_matrix_tmp.h5ad"
if os.path.isfile(FRESH_CACHE):
    import snapatac2 as snap
    po = snap.read(PRIV_CACHE, backed="r")
    pn = snap.read(FRESH_CACHE, backed="r")
    bo, bn = [str(b) for b in po.obs_names], [str(b) for b in pn.obs_names]
    vo, vn = [str(v) for v in po.var_names], [str(v) for v in pn.var_names]
    Xo, Xn = po.X[:], pn.X[:]
    Xo = Xo.toarray() if hasattr(Xo, "toarray") else np.asarray(Xo)
    Xn = Xn.toarray() if hasattr(Xn, "toarray") else np.asarray(Xn)
    po.close(); pn.close()
    det = [f"shape {Xn.shape} vs {Xo.shape}", f"regions equal {vo == vn}",
           f"totals {Xn.sum(0).astype(int).tolist()} vs {Xo.sum(0).astype(int).tolist()}"]
    if bo == bn and vo == vn:
        ncell = int((Xo != Xn).any(1).sum())
        det.append(f"barcodes equal True; cells with any differing count {ncell:,}")
        status = "PASS" if ncell == 0 else "DIFF"
    else:
        ko, kn = cache_keys(PRIV_META, bo), cache_keys(FRESH_META, bn) if os.path.isfile(FRESH_META) else None
        if ko is None or kn is None or not (ko.is_unique and kn.is_unique):
            det.append("barcodes equal False; caches not positionally aligned to unique (arc_library, barcode) "
                       "metadata, so no per-cell alignment")
        else:
            common = kn.intersection(ko)
            det.append(f"barcodes equal False; common cells {len(common):,}, only_old {len(ko.difference(kn)):,}, "
                       f"only_new {len(kn.difference(ko)):,}")
            if vo == vn and len(common):
                ro, rn = ko.get_indexer(common), kn.get_indexer(common)
                det.append(f"common cells with differing counts {int((Xo[ro] != Xn[rn]).any(1).sum()):,}")
        status = "DIFF"
    rec("A3_f5_region_matrix_vs_20260810_cache", status, "; ".join(det))
else:
    rec("A3_f5_region_matrix_vs_20260810_cache", "MISSING", FRESH_CACHE)

cmp_f5("A4_f5_fresh_default_vs_committed", f"{X}/out/f5_fresh/step_f5_pitx2_region_counts.tsv")

# ---------------- Fig 5C ----------------
F4 = f"{D}/f4_compartment_perembryo_pitx2_coverage.tsv.gz"
HH13 = f"{D}/f4_hh13_perembryo_pitx2_coverage.tsv.gz"
SHIP = f"{R}/phase2/atac/lr_pileup_compartment_perembryo"
SHIP13 = f"{R}/phase2/atac/lr_pileup_hh13_perembryo"


def load_cov(p):
    return pd.read_csv(p, sep="\t")


def groups(df):
    return sorted(c for c in df.columns if c not in ("chrom", "pos") and not c.endswith("_percell"))


def cmp_cov(check, new_p, ref_p, skip=()):
    if not os.path.isfile(new_p):
        rec(check, "MISSING", new_p)
        return None
    new, ref = load_cov(new_p), load_cov(ref_p)
    if not (new[["chrom", "pos"]].equals(ref[["chrom", "pos"]])):
        rec(check, "DIFF", "chrom/pos axis differs")
        return new
    gn, gr = set(groups(new)) - set(skip), set(groups(ref)) - set(skip)
    common = sorted(gn & gr)
    bad = []
    for g in common:
        if not np.array_equal(new[g].to_numpy(), ref[g].to_numpy()):
            bad.append(f"{g} raw")
        elif not np.allclose(new[f"{g}_percell"], ref[f"{g}_percell"], rtol=1e-9, atol=0):
            bad.append(f"{g} percell")
    ok = not bad and gn == gr
    rec(check, "PASS" if ok else "DIFF",
        f"groups {len(gn)} vs {len(gr)} (excluded {len(skip)}); exact {len(common) - len(bad)}/{len(common)}; "
        f"only_new {sorted(gn - gr)}; only_ref {sorted(gr - gn)}; differing {bad[:20]}")
    return new


cmp_cov("C2_f4_shipped_vs_committed", f"{SHIP}/step_f4_pitx2_coverage.tsv.gz", F4)

F4N = f"{D}/f4_compartment_perembryo_pitx2_cell_counts.tsv"  # panel C's nuclei-per-column sidecar
new_n = f"{SHIP}/step_f4_cell_counts.tsv"
if os.path.isfile(new_n):
    ref = pd.read_csv(F4N, sep="\t").sort_values("group").reset_index(drop=True)
    new = pd.read_csv(new_n, sep="\t")
    missing = [c for c in ref.columns if c not in new.columns]
    if missing:
        rec("C2n_f4_cell_counts_vs_committed", "DIFF", f"producer lacks {missing}")
    else:
        new = new[list(ref.columns)].sort_values("group").reset_index(drop=True)  # drops producer-only path columns
        same = (new["group"].equals(ref["group"]) and new["n_cells"].equals(ref["n_cells"])
                and new["coverage_sum_pitx2_locus"].equals(ref["coverage_sum_pitx2_locus"])
                and np.allclose(new["coverage_sum_per_cell"], ref["coverage_sum_per_cell"], rtol=1e-9, atol=0))
        rec("C2n_f4_cell_counts_vs_committed", "PASS" if same else "DIFF",
            f"groups {len(new)} vs {len(ref)}; group/n_cells/coverage_sum identical {same}")
else:
    rec("C2n_f4_cell_counts_vs_committed", "MISSING", new_n)
cmp_cov("C3_f4_shipped_hh13_vs_committed", f"{SHIP13}/step_f4_pitx2_coverage.tsv.gz", HH13)

sys.exit(1 if any(status != "PASS" for _, status, _ in rows) else 0)  # sacct STATE carries the verdict
