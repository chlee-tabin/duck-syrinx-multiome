#!/usr/bin/env python3
"""Group B cS-split doublet/unassigned cells: build the committed exclusion table + drop helpers.

The resolved-12 Group B overlay (``_groupB_resolved12.py``) splits two merged vireo donors with
an N=2 re-vireo (``cS_targeted_split.py``): donor11 -> donor11a/b and donor2 -> donor2a/b.  In
each split's ``N2/donor_ids.tsv`` a cell is either cleanly ``donor0``/``donor1`` or flagged
``doublet``/``unassigned``.  Until 2026-10 the overlay defaulted the flagged cells to
"{base}a" (``.fillna``), so 163 cells (donor11: 37 doublet + 8 unassigned = 45; donor2:
106 doublet + 12 unassigned = 118) were silently counted as embryo donor11a / donor2a.  These
cells are now EXCLUDED — consistent with doublet removal using the vireo calls only.

This module owns that rule:

* ``build`` (CLI): read the two split ``donor_ids.tsv`` files and write
  ``figures/data/chick_control/groupB_split_excluded_cells.tsv``
  (barcode, base_donor, split_status), with the source paths + md5s in ``#`` header lines.
  The header records the path actually read (``--split``) and, via ``--origin``, where that copy
  came from; the md5 ties the two.  The committed table was built from local copies of:
    STANDBY_FC/chick_backup/chick_demux_B/vireo_seeded/seed_qc/split_donor11/N2/donor_ids.tsv
        md5 3c1298d32574a98292bf2698d0b15de1  (1449 donor1 / 288 donor0 / 37 doublet / 8 unassigned)
    STANDBY_FC/chick_backup/chick_demux_B/vireo_seeded/seed_qc/split_donor2/N2/donor_ids.tsv
        md5 eb6a0dd3ee62fb95f18eb91f300339be  (2551 donor1 / 2704 donor0 / 106 doublet / 12 unassigned)
  (scratch originals: SCRATCH_USER/chick_demux_B/vireo_seeded/seed_qc/split_donor{11,2}/N2/).
  Each header line also records the CLEAN assignment signature (``clean_signature``: the
  donor0/donor1 counts + a sha256 over the sorted ``cell<TAB>donor0|donor1`` pairs).  It is
  content-level (independent of float formatting / extra columns), so ``load_split_map`` can
  require that a split file read at run time assigns exactly the same cells to the same clean
  sub-donor as the files the table was built from: a rerun that swaps donor0/donor1 (which would
  exchange the fixed {base}a/{base}b stage identities) or reassigns any clean cell raises.

* ``fold_split(dmap, base, path, group)``: the split fold for the cluster producers that read
  the split files directly (cAD / cZ / cP2 / cB / cAB): donor0 -> {base}b, donor1 -> {base}a,
  doublet/unassigned keep that status (the caller's doublet filter drops them); missing file,
  duplicated/absent cell, unexpected status, or an excluded cell no longer labelled ``base``
  raises; Group B donor11/donor2 asserted against the table.  ``load_split_map`` +
  ``check_excluded_labels`` are the shared validation, also used by ``resolve_donor``.

* ``load_split_excluded()``: the committed table as a DataFrame (the single reader).

* ``drop_split_excluded(df, table)``: drop those cells from a per-cell table that predates the
  fix, scoped to Group B.  The hit count must equal ``EXPECTED[table]`` (the count in the
  committed, pre-fix table) or 0 (a table regenerated with the fixed overlay); anything else
  raises, so a drifted table can never be filtered silently; EVERY Group B barcode must be in
  the exclusion ids' format (per barcode), so no excluded cell can escape under another format.  ``excluded_archr_cbs()`` gives
  the same cells as ArchR ids for ATAC-keyed tables.

The table path is ``$CHICK_SPLIT_EXCLUDED_TSV`` if set, else the repo copy (resolved lazily).

Usage:
  conda run -n duckfig python3 pipeline/chick/_groupB_split_excluded.py \\
      --split donor11=<path>/split_donor11/N2/donor_ids.tsv \\
      --split donor2=<path>/split_donor2/N2/donor_ids.tsv \\
      [--origin donor11=<where the copy came from> --origin donor2=...]
"""
import argparse
import hashlib
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


TABLE_NAME = "groupB_split_excluded_cells.tsv"
TABLE_ENV = "CHICK_SPLIT_EXCLUDED_TSV"   # override (cluster scratch copies have no repo tree)


def _repo_root():
    """Walk up from this file to the checkout root (dir containing pipeline/_paths.py), or None."""
    d = HERE
    while True:
        if os.path.exists(os.path.join(d, "pipeline", "_paths.py")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def table_path():
    """The committed exclusion table: $CHICK_SPLIT_EXCLUDED_TSV, else the repo copy. Resolved
    lazily (not at import) so a scratch copy of pipeline/chick/ without the repo tree can still
    import the modules that depend on this one; reading the table then fails loud."""
    env = os.environ.get(TABLE_ENV)
    if env:
        return env
    root = _repo_root()
    if root is None:
        raise RuntimeError(f"{TABLE_NAME}: no repo checkout above {HERE} — set {TABLE_ENV} to the "
                           f"committed figures/data/chick_control/{TABLE_NAME}")
    return os.path.join(root, "figures", "data", "chick_control", TABLE_NAME)


CLEAN = ("donor0", "donor1")          # the two clean sub-donors of an N=2 split
EXCLUDED_STATUSES = ("doublet", "unassigned")
BASES = ("donor11", "donor2")
GROUPS = ("A", "B")                   # the only chick control group labels in per-cell tables
SPLIT_EXCLUDED = "split_excluded"     # sentinel donor label written by resolve_donor

# Hits of the 163 excluded cells in the committed PRE-FIX per-cell tables (a table regenerated
# with the fixed overlay, or pruned with --prune-resolved12, has 0).  All 163 are in
# groupB_resolved12_percell (29,524 Group B cells, written by an earlier cAK run than the committed
# cAK/cAE atlas tables, whose Group B has 29,080 cells), hence 148 in cAK/cAE; 92 survive into the
# LPM-lineage cAM subset, 72 of those are coupled ATAC nuclei in the chick PITX2-locus per-cell
# table, and the ArchR sidecar (archr_metadata_sidecar.py, regenerated) carries 0.
EXPECTED = {
    "groupB_resolved12_percell": 163,
    "cAK_atlas_wz_percell": 148,
    "cAE_atlas_umap": 148,
    "cAM_percell_populations": 92,
    "chick_archr_coupled_metadata": 92,
    "chick_locus_percell_counts": 72,
}
# GEX barcode lane suffix -> ARC (ATAC fragment) library, Group B multiome lanes
_GEX_SUFFIX_TO_ARC = {"9": "ARC09", "10": "ARC10"}


def split_status_counts(path):
    d = pd.read_csv(path, sep="\t")
    bad = sorted(set(d["donor_id"].astype(str)) - set(CLEAN) - set(EXCLUDED_STATUSES))
    if bad:
        raise ValueError(f"{path}: unexpected vireo donor_id values {bad}")
    dup = d.loc[d["cell"].duplicated(keep=False), "cell"].astype(str)
    if len(dup):
        raise ValueError(f"{path}: {dup.nunique()} duplicated cell barcodes "
                         f"(e.g. {sorted(set(dup))[:3]})")
    return d


def excluded_from_split(base, path):
    """Rows of the split file that are not cleanly donor0/donor1 -> (barcode, base_donor, status)."""
    d = split_status_counts(path)
    ex = d[~d["donor_id"].astype(str).isin(CLEAN)]
    return pd.DataFrame({"barcode": ex["cell"].astype(str).values,
                         "base_donor": base,
                         "split_status": ex["donor_id"].astype(str).values})


def clean_signature(d):
    """(counts, sha256) of a split file's CLEAN assignment: {'donor0': n, 'donor1': n} and a sha256
    over the sorted 'cell<TAB>donor0|donor1' lines.  Catches a donor0/donor1 swap or any reassigned
    clean cell, independent of the file's float formatting or extra columns."""
    c = d[d["donor_id"].astype(str).isin(CLEAN)]
    pairs = sorted(zip(c["cell"].astype(str), c["donor_id"].astype(str)))
    h = hashlib.sha256("\n".join(f"{a}\t{b}" for a, b in pairs).encode()).hexdigest()
    counts = {k: int((c["donor_id"].astype(str) == k).sum()) for k in CLEAN}
    return counts, h


_SIG_RE = re.compile(r"^#\s+(donor\d+): .*clean_donor0=(\d+) clean_donor1=(\d+) "
                     r"clean_sha256=([0-9a-f]{64})")


def recorded_clean_signatures():
    """{base: (counts, sha256)} from the committed table's '#' header (written by ``main``).
    Raises if a base's signature is absent, so the clean-mapping check can never be skipped."""
    path = table_path()
    sig = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            m = _SIG_RE.match(line)
            if m:
                if m.group(1) in sig:
                    raise RuntimeError(f"{path}: two clean_sha256 header lines for {m.group(1)}")
                sig[m.group(1)] = ({"donor0": int(m.group(2)), "donor1": int(m.group(3))},
                                   m.group(4))
    missing = [b for b in BASES if b not in sig]
    if missing:
        raise RuntimeError(f"{path}: no clean_sha256 header line for {missing} — rebuild the table "
                           f"with pipeline/chick/_groupB_split_excluded.py")
    return sig


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_split_excluded():
    path = table_path()
    t = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    if list(t.columns) != ["barcode", "base_donor", "split_status"]:
        raise ValueError(f"{path}: unexpected columns {list(t.columns)}")
    if t["barcode"].duplicated().any():
        raise ValueError(f"{path}: duplicated barcodes")
    if not set(t["base_donor"]) <= set(BASES) or not set(t["split_status"]) <= set(EXCLUDED_STATUSES):
        raise ValueError(f"{path}: unexpected base_donor/split_status values")
    return t


def excluded_archr_cbs():
    """The excluded cells as ArchR cell ids ('ARC09#<raw>-1' / 'ARC10#<raw>-1').

    GEX barcode '<raw>-1-<9|10>' -> ATAC id, the same translation as archr_metadata_sidecar.py."""
    out = set()
    for bc in load_split_excluded()["barcode"]:
        raw, suf = bc.rsplit("-", 1)
        if suf not in _GEX_SUFFIX_TO_ARC:
            raise ValueError(f"excluded barcode {bc!r}: lane suffix {suf!r} is not a Group B lane")
        out.add(f"{_GEX_SUFFIX_TO_ARC[suf]}#{raw}")
    return out


_GEX_ID = re.compile(r"^[ACGTN]+-1-(9|10)$")       # Group B GEX cell id '<raw>-1-<9|10>'
_ARC_ID = re.compile(r"^ARC(09|10)#[ACGTN]+-1$")     # Group B ArchR cell id 'ARC09#<raw>-1'


def drop_split_excluded(df, table, barcode_col="barcode", group_col="group", ids=None,
                        log=print):
    """Return ``df`` without the Group B split-excluded cells.

    ``table`` names the per-cell table (a key of ``EXPECTED``); the hit count must equal its
    pre-fix count or 0 (a table regenerated with the fixed overlay) — anything else raises.
    ``ids`` defaults to the GEX barcodes; pass ``excluded_archr_cbs()`` for ATAC-keyed tables.
    If ``group_col`` is None the caller guarantees ``df`` is Group-B-only (or ids are
    Group-B-specific, as the ARC09/ARC10 ids are).

    Every Group B ``barcode_col`` value must be in the same id format as ``ids`` (checked per
    barcode) — a wrong column or a different suffix convention, for all rows or only some,
    would otherwise let excluded cells through silently.  With ``group_col``, its values must be
    exactly ``A``/``B`` and each row's label must agree with its barcode's lane/library identity
    (Group B <=> the exclusion-id format), so a drifted label (``"B "``, ``"b"``, a B cell
    relabelled ``A``) raises instead of skipping the drop.
    """
    expected = EXPECTED[table]
    ex = set(load_split_excluded()["barcode"]) if ids is None else set(ids)
    pat = _ARC_ID if all("#" in i for i in ex) else _GEX_ID
    if not all(pat.match(i) for i in ex):
        raise ValueError(f"[split-excluded] {table}: exclusion ids are not all {pat.pattern}")
    col = df[barcode_col].astype(str)
    if group_col is not None:
        grp = df[group_col].astype(str)
        bad_grp = sorted(set(grp) - set(GROUPS))
        if bad_grp:
            raise RuntimeError(f"[split-excluded] {table}: unexpected '{group_col}' values {bad_grp} "
                               f"(expected {GROUPS}) — a drifted group label would hide Group B rows")
        in_b = grp.eq("B")
        # the group label must agree with the barcode's own lane/library identity (Group B =
        # GEX lanes 9/10 or ArchR ARC09/ARC10), so a relabelled row can neither escape the
        # Group B scope nor smuggle a Group A cell into it
        by_id = col.str.match(pat.pattern)
        mism = col[by_id != in_b]
        if len(mism):
            raise RuntimeError(f"[split-excluded] {table}: {len(mism)} rows whose '{group_col}' label "
                               f"disagrees with their barcode lane/library (Group B ids match "
                               f"{pat.pattern}; e.g. {list(zip(mism.iloc[:3], grp[mism.index[:3]]))})")
    else:
        in_b = pd.Series(True, index=df.index)
    hit = col.isin(ex) & in_b
    n = int(hit.sum())
    if n not in (expected, 0):
        raise RuntimeError(f"[split-excluded] {table}: {n} Group B rows match {TABLE_NAME}, "
                           f"expected {expected} (pre-fix table) or 0 (regenerated) — the table "
                           f"and the exclusion list have drifted")
    bad_fmt = col[in_b & ~col.str.match(pat.pattern)]
    if len(bad_fmt):
        raise RuntimeError(f"[split-excluded] {table}: {len(bad_fmt)} Group B '{barcode_col}' values "
                           f"are not in the exclusion id format {pat.pattern} (e.g. "
                           f"{list(bad_fmt.iloc[:3])}) — wrong column or barcode convention; an "
                           f"excluded cell in another format would not be matched")
    log(f"[split-excluded] {table}: dropped {n} Group B split doublet/unassigned cells "
        f"(expected {expected} or 0)")
    return df[~hit].copy()


def load_split_map(base, path, group):
    """The shared split-file validation (fold_split and _groupB_resolved12.resolve_donor).

    Reads one cS N=2 ``donor_ids.tsv`` -> {cell: donor0/donor1/doublet/unassigned}.  A missing
    file, a duplicated cell barcode or an unexpected vireo status raises; for the Group B
    donor11/donor2 splits the clean assignment must match the recorded ``clean_signature``
    (counts + sha256 in the table header) and the non-clean set must equal the committed table."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"[split] {group}/{base}: split file missing: {path}")
    d = split_status_counts(path)
    m = dict(zip(d["cell"].astype(str), d["donor_id"].astype(str)))
    if group == "B" and base in BASES:
        counts, sha = clean_signature(d)
        exp_counts, exp_sha = recorded_clean_signatures()[base]
        if sha != exp_sha:
            raise RuntimeError(f"[split] {group}/{base}: clean donor0/donor1 assignment of {path} "
                               f"{counts} sha256={sha[:12]}... differs from the one {TABLE_NAME} "
                               f"was built from {exp_counts} sha256={exp_sha[:12]}... — the split "
                               f"was rerun (a donor0/donor1 swap would exchange the fixed "
                               f"{base}a/{base}b identities); rebuild the table and re-check the overlay")
        committed = load_split_excluded()
        table_ex = set(committed.loc[committed["base_donor"] == base, "barcode"])
        split_ex = {c for c, s in m.items() if s not in CLEAN}
        if split_ex != table_ex:
            raise RuntimeError(f"[split] {group}/{base}: split file excludes {len(split_ex)} cells, the "
                               f"committed {TABLE_NAME} lists {len(table_ex)} "
                               f"(sym. diff {len(split_ex ^ table_ex)})")
    return m


def check_excluded_labels(labels, base, group):
    """Every committed excluded cell of ``base`` present in ``labels`` (barcode -> current donor
    label, a dict or Series) must still carry the label ``base``; a rerun that relabelled one
    (e.g. a donor11 doublet now called donor4) would otherwise keep it silently.  Group B
    donor11/donor2 only; per barcode, raises naming the offending cells."""
    if not (group == "B" and base in BASES):
        return
    committed = load_split_excluded()
    ex = committed.loc[committed["base_donor"] == base, "barcode"]
    drift = {bc: str(labels[bc]) for bc in ex if bc in labels and str(labels[bc]) != base}
    if drift:
        shown = dict(list(drift.items())[:5])
        raise RuntimeError(f"[split] {group}/{base}: {len(drift)} committed split-excluded cells are "
                           f"no longer labelled {base!r} in the donor map (e.g. {shown}) — donor "
                           f"labels and {TABLE_NAME} have drifted")


def fold_split(dmap, base, path, group, log=print):
    """Fold one cS N=2 split into a cell -> donor dict (the overlay used by cAD/cZ/cP2/cB/cAB).

    Cells labelled ``base`` become ``{base}b`` (split vireo donor0) / ``{base}a`` (donor1); a
    cell the split calls ``doublet``/``unassigned`` takes that status as its label, so the
    caller's existing doublet/unassigned filter drops it (no default to ``{base}a``).  A
    missing split file, a ``base`` cell absent from it, or an unexpected status raises.  For
    the Group B donor11/donor2 splits the non-clean set must equal the committed exclusion
    table.  Returns a new dict.
    """
    m = load_split_map(base, path, group)
    check_excluded_labels(dmap, base, group)
    cells = [c for c, v in dmap.items() if v == base]
    absent = [c for c in cells if c not in m]
    if absent:
        raise RuntimeError(f"[split] {group}/{base}: {len(absent)} of {len(cells)} cells absent "
                           f"from {path} — donor map and split have drifted")
    lab = {"donor0": f"{base}b", "donor1": f"{base}a"}
    out = dict(dmap)
    for c in cells:
        out[c] = lab.get(m[c], m[c])
    vc = pd.Series([out[c] for c in cells]).value_counts().to_dict()
    log(f"[split] folded {group}/{base} ({len(cells)} cells) from {path}: {vc} "
        f"(doublet/unassigned are dropped by the caller)")
    return out


def prune_resolved12(path):
    """Drop the excluded cells from the committed resolved-12 per-cell record, in place.

    That table is written by cAK on the cluster (no local producer); pruning it here gives the
    rows a fixed cAK would write (it drops SPLIT_EXCLUDED cells before writing). Deterministic
    gzip (mtime=0) so a re-run is byte-identical."""
    import gzip
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    out = drop_split_excluded(df, "groupB_resolved12_percell", group_col=None)
    left = out[out["base_donor"].isin(BASES)]
    bad = left[~left["donor"].isin([f"{b}{s}" for b in BASES for s in "ab"])]
    if len(bad):
        raise RuntimeError(f"{path}: {len(bad)} split-base rows still not {{base}}a/b after pruning")
    with gzip.GzipFile(path, "wb", mtime=0) as g:
        g.write(out.to_csv(sep="\t", index=False).encode())
    print(f"pruned {path}: {len(df)} -> {len(out)} rows")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--split", action="append", metavar="BASE=PATH",
                    help="base donor and its N2 donor_ids.tsv (give once per base: donor11, donor2)")
    ap.add_argument("--origin", action="append", default=[], metavar="BASE=WHERE",
                    help="optional: where the --split file was copied from (recorded in the "
                         "header next to the path actually read; the md5 ties the two)")
    ap.add_argument("--out", default=None, help="default: the committed repo table")
    ap.add_argument("--prune-resolved12", metavar="TSV_GZ", default=None,
                    help="instead of building: drop the excluded cells from this "
                         "groupB_resolved12_percell.tsv.gz in place (no local reader re-derives it)")
    a = ap.parse_args()
    if a.prune_resolved12:
        prune_resolved12(a.prune_resolved12)
        return
    if not a.split:
        raise SystemExit("--split is required when building the table")
    a.out = a.out or table_path()
    splits = dict(s.split("=", 1) for s in a.split)
    origins = dict(s.split("=", 1) for s in a.origin)
    if set(splits) != set(BASES):
        raise SystemExit(f"need --split for exactly {BASES}, got {sorted(splits)}")
    parts, header = [], [
        "# Group B cS-split cells NOT cleanly assigned to donor0/donor1 by the N=2 re-vireo;",
        "# EXCLUDED from all resolved-12 analyses (formerly defaulted to {base}a).",
        "# built by pipeline/chick/_groupB_split_excluded.py from:",
    ]
    for base in BASES:
        p = splits[base]
        ex = excluded_from_split(base, p)
        parts.append(ex)
        counts = ex["split_status"].value_counts().to_dict()
        cc, sha = clean_signature(split_status_counts(p))
        line = (f"#   {base}: read {p} md5={md5(p)} excluded={len(ex)} {counts} "
                f"clean_donor0={cc['donor0']} clean_donor1={cc['donor1']} clean_sha256={sha}")
        if base in origins:
            line += f" (copy of {origins[base]})"
        header.append(line)
    out = pd.concat(parts, ignore_index=True)
    with open(a.out, "w") as f:
        f.write("\n".join(header) + "\n")
        out.to_csv(f, sep="\t", index=False)
    print(f"wrote {a.out}: {len(out)} cells {out.groupby(['base_donor', 'split_status']).size().to_dict()}")


if __name__ == "__main__":
    main()
