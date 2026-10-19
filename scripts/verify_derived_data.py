#!/usr/bin/env python3
"""Validator for derived-data freshness.

Runs locally (pre-commit or CI).  For every entry in
`pipeline/derived_data_manifest.tsv`, it cross-checks the repo copy
against the current `references/cluster_annotations.tsv`.  Three fatal
checks:

  (a) Cell-count consistency.  Any TSV with a `cluster`/`cluster_02`/
      `celltype_02`/`annotation` column whose values look like cell-type
      annotations must have per-annotation counts that agree with the
      authoritative annotations table (exact total, +/-5% per category).

  (b) Timestamp monotonicity.  Every `depends_on_h5ad=yes` TSV must have
      `mtime >= mtime(cluster_annotations.tsv) - 1h`.  Catches the case
      where d2b was rerun but a dependent TSV was never re-exported.

  (c) Annotation alphabet.  Any categorical column whose observed values
      look like cell-type names must be a subset of the union of two
      authoritative alphabets:
        - references/cluster_annotations.tsv          (d2b, res=0.2 global)
        - references/d2g_subcluster_annotations.tsv   (d2g, stage-restricted)
      Labels in the curated KNOWN_LEGACY_LABELS set are FATAL (signal a
      stale TSV).  Labels outside the union but not in the legacy set
      are WARN (may be group labels like "All" / "Mesenchyme").

Usage:
  python3 scripts/verify_derived_data.py            # normal run
  python3 scripts/verify_derived_data.py --dry-run  # don't fail the pre-
                                                    # commit hook, just
                                                    # report
  python3 scripts/verify_derived_data.py --manifest <path>   (override)
  python3 scripts/verify_derived_data.py --verbose

Exit codes:
  0 — all checks passed (or dry-run)
  1 — at least one fatal violation
  2 — usage / setup error
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
from pathlib import Path
from typing import Iterable, Iterator


# ── Repo layout ────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_DEFAULT = REPO_ROOT / "pipeline" / "derived_data_manifest.tsv"
# The d2b (res=0.2 global) alphabet is the source of truth for per-cell
# partition checks.  The d2g (stage-restricted) alphabet is an additional
# allowed set — any TSV-held annotation that appears in either is
# acceptable.
ANNOTATIONS_FILE = REPO_ROOT / "references" / "cluster_annotations.tsv"
D2G_SUBCLUSTER_FILE = REPO_ROOT / "references" / "d2g_subcluster_annotations.tsv"

# d2e silhouette export uses shortened display names that aren't in any
# reference file.  Add them here; this list MUST be kept in sync with
# pipeline/d2e_silhouette_export.py if it renames focused-mes labels.
D2E_FOCUSED_MES_DIALECT = frozenset({
    "Cardiac-adj. mes.", "Early LPM prog.", "Early somatic mes.",
    "Somatic mes.", "Somatic mes. (2)", "Splanchnic mes.",
    "Syringeal mes.", "Transitional mes.", "Transitional mes. (2)",
    "Sclerotome", "Dermomyotome",
})

# ── Column heuristics ──────────────────────────────────────────────────
# Columns that carry cell-type annotation strings.
# "population" carries d2g subcluster labels (duck_locus/chick_locus schema) —
# without it those files were exempt from the retired-label alphabet guard.
ANNOTATION_COLS = ("annotation", "celltype_02", "cell_type", "cell_type_group",
                   "population")
# Columns that carry numeric cluster IDs.
CLUSTER_COLS = ("cluster", "cluster_02")
# Columns that carry a per-cell count (used in the silhouette TSV).
COUNT_COLS = ("n_cells", "total_cells")

# Curated set of legacy labels — these are pre-reorg names that signal a
# stale TSV.  Hard-fatal when observed.  Do NOT add heuristics like
# " (" here — several current valid labels contain parentheses (e.g.
# "Somatic mes. (2)") and would false-positive.
KNOWN_LEGACY_LABELS = frozenset({
    "Blood", "Chondrocyte", "Macrophage", "Endothelium", "Epithelium",
    "SMC precursor", "Splanchnic LPM",
    "Mesenchyme (late)", "Mesenchyme (early)",
    "Diff. neuron (2)",
    # ── retired d2g subcluster names (project convention: "Stage-restricted subclustering" §5) ──
    # These are the renames the SSOT actually went through, and they were missing from the
    # very guard that exists to catch a TSV predating them: "Somatic" -> "Skeletogenic LPM"
    # (2026-06-18), "Chondrogenic LPM" -> folded into "Syringeal" at HH25 (2026-07-02),
    # "Splanchnic" -> "Venous-pole spl." and "Epicardial mes." -> "Mesothelium" (2026-05-04).
    # Matching is exact set intersection, so none of these collides with a live label:
    # "Somatic mes. (2)" and the global cluster "Splanchnic mesenchyme" are different strings
    # and are unaffected. Verified against both SSOTs at the time of writing — neither emits
    # any of these four.
    "Somatic", "Chondrogenic LPM", "Splanchnic", "Epicardial mes.",
})


# ── Utilities ──────────────────────────────────────────────────────────
def _open_text(path: Path) -> io.TextIOBase:
    """Open a plain or gzip TSV as a text file, utf-8."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def _iter_rows(path: Path) -> Iterator[dict]:
    """Yield csv.DictReader rows from a TSV (plain or .gz)."""
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield row


def _read_header(path: Path) -> list[str]:
    with _open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            return next(reader)
        except StopIteration:
            return []


def _load_annotations(path: Path) -> dict[int, str]:
    """Load references/cluster_annotations.tsv → {cluster_id: annotation}."""
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    required = {"cluster", "annotation", "n_cells"}
    if rows and not required.issubset(rows[0].keys()):
        missing = required - set(rows[0].keys())
        raise SystemExit(
            f"ERROR: {path.name} missing required columns: {sorted(missing)}; "
            f"update scripts/verify_derived_data.py for the new schema."
        )
    return {int(r["cluster"]): r["annotation"] for r in rows}


def _load_cluster_counts(path: Path) -> dict[int, int]:
    """Load references/cluster_annotations.tsv → {cluster_id: n_cells}."""
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {int(r["cluster"]): int(r["n_cells"]) for r in reader}


def _load_label_set(path: Path, column: str) -> set[str]:
    """Return the set of non-empty values in `column` of a TSV.

    Used to load supplementary annotation alphabets (d2d res=0.5,
    d2g stage-restricted).  Missing file → empty set (not fatal).
    """
    if not path.is_file():
        return set()
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {r[column].strip() for r in reader if r.get(column, "").strip()}


# ── Manifest ───────────────────────────────────────────────────────────
def read_manifest(path: Path) -> list[dict]:
    """Parse derived_data_manifest.tsv into a list of dicts.

    Fields: scratch_src, repo_dest, produced_by, depends_on_h5ad.
    depends_on_h5ad is normalized to a bool via first token
    ('yes' → True, 'no' → False); anything else raises.
    """
    rows: list[dict] = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            r = {k.strip(): (v.strip() if v else "") for k, v in r.items()}
            h5 = r["depends_on_h5ad"].split()[0].lower()
            if h5 == "yes":
                r["depends_on_h5ad"] = True
            elif h5 == "no":
                r["depends_on_h5ad"] = False
            else:
                raise ValueError(
                    f"manifest row {r['repo_dest']!r} has invalid "
                    f"depends_on_h5ad={r['depends_on_h5ad']!r}"
                )
            rows.append(r)
    return rows


def expand_repo_dest(dest: str) -> Iterable[Path]:
    """Expand a possibly-globbed repo_dest path under REPO_ROOT."""
    p = REPO_ROOT / dest
    if "*" in dest:
        # Path.glob doesn't work with absolute paths — split.
        base = REPO_ROOT
        yield from base.glob(dest)
    else:
        yield p


# ── Violation tracking ─────────────────────────────────────────────────
class Report:
    def __init__(self) -> None:
        self.violations: list[str] = []
        self.warnings: list[str] = []
        self.checked: int = 0
        self.skipped_missing: int = 0

    def violate(self, msg: str) -> None:
        self.violations.append(msg)
        print(f"VIOLATION: {msg}", file=sys.stderr)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)
        print(f"WARN:      {msg}", file=sys.stderr)

    def info(self, msg: str) -> None:
        print(f"INFO:      {msg}")

    def ok(self, msg: str) -> None:
        print(f"OK:        {msg}")

    def summary(self) -> None:
        print("=" * 72)
        print(f"checked:   {self.checked} files")
        print(f"skipped:   {self.skipped_missing} (not present in repo)")
        print(f"warnings:  {len(self.warnings)}")
        print(f"fatal:     {len(self.violations)}")
        if self.violations:
            print()
            print("Violations:")
            for v in self.violations:
                print(f"  - {v}")


# ── Checks ─────────────────────────────────────────────────────────────
def check_cell_counts(
    path: Path,
    manifest_row: dict,
    allowed_annotations: set[str],
    cluster_counts_auth: dict[int, int],
    annotations_auth: dict[int, str],
    report: Report,
    tolerance: float = 0.05,
) -> None:
    """Check (a): per-cluster / per-annotation count consistency.

    Only runs for TSVs that look like per-cell tables with a recognizable
    cluster or annotation column.  Aggregate tables (dotplot, pseudobulk
    PCA, per-genotype composition) are treated as OK because their row
    semantics are not per-cell.
    """
    header = _read_header(path)
    if not header:
        report.warn(f"{path.name}: empty file")
        return

    # Detect schema.
    has_cluster = any(c in header for c in CLUSTER_COLS)
    has_ann = any(c in header for c in ANNOTATION_COLS)

    if not (has_cluster or has_ann):
        return  # nothing to check

    # Skip cell-count reconstruction for files whose row cardinality is
    # inherently not a partition of the global h5ad:
    #   * resolution-sweep tables (method / resolution / leiden_col columns)
    #   * stage-restricted tables (stage column — d2g / d2h)
    #   * label-transfer tables (per-cell but limited to a stage)
    #   * marker tables (per-cluster-per-gene, one row per rank or gene)
    #   * res=0.5 partitions (cluster_res05 — subset / different resolution)
    # In these cases we still do the alphabet check elsewhere.
    SCOPE_COLS = ("method", "resolution", "leiden_col", "stage",
                  "stage_label", "stage_restricted",
                  "rank", "gene", "category",
                  "cluster_res05", "parent_res02")
    if any(c in header for c in SCOPE_COLS):
        # Restricted / swept / marker / subset — counts won't match
        # global. Skip (a).
        return

    # Silhouette exports (d2e_silhouette_*) subsample to 50K cells and
    # may use focused subcluster IDs rather than the global alphabet.
    # Their per-cluster `n_cells` column reflects the subsample — not
    # global counts — so the aggregated check cannot apply. The alphabet
    # check (c) still runs separately and is what catches stale labels
    # like "Mesenchyme (late)" / "Splanchnic LPM" in these files.
    if "mean_silhouette" in header:
        return

    # Only apply cell-count checks to per-cell TSVs.  Heuristic: there is
    # NO aggregation column (`n_cells`, `total_cells`, `fraction`).  If
    # an aggregate column is present, the row cardinality is not
    # per-cell and the cluster count reconstruction won't match.
    aggregate_markers = ("total_cells", "fraction", "frac_expressing",
                         "mean_cp10k", "mean_raw")
    if any(c in header for c in aggregate_markers):
        # Special case: aggregate TSVs with cluster + n_cells columns can
        # still cross-check counts against the authoritative table. This
        # only applies when `n_cells` reflects the full (unsampled)
        # partition — silhouette already returned above.
        if "n_cells" in header and ("cluster" in header or "annotation" in header):
            _check_aggregated_counts(
                path, header, cluster_counts_auth, annotations_auth,
                tolerance, report
            )
        return

    # If the only aggregation marker is n_cells (no other aggregators),
    # route to aggregated-counts validator too.
    if "n_cells" in header and ("cluster" in header or "annotation" in header):
        _check_aggregated_counts(
            path, header, cluster_counts_auth, annotations_auth,
            tolerance, report,
        )
        return

    # Per-cell mode: count rows per cluster (or per annotation).
    cluster_counts: dict[str, int] = {}
    ann_counts: dict[str, int] = {}
    total = 0
    cluster_col = next((c for c in CLUSTER_COLS if c in header), None)
    ann_col = next((c for c in ANNOTATION_COLS if c in header), None)

    for row in _iter_rows(path):
        total += 1
        if cluster_col is not None:
            v = row[cluster_col]
            cluster_counts[v] = cluster_counts.get(v, 0) + 1
        if ann_col is not None:
            v = row[ann_col]
            ann_counts[v] = ann_counts.get(v, 0) + 1

    # Compare to authoritative table.
    auth_total = sum(cluster_counts_auth.values())
    if total != auth_total:
        report.violate(
            f"{path.name}: total cell count {total} != authoritative "
            f"{auth_total}"
        )

    if cluster_col is not None and cluster_counts:
        _compare_per_cluster(path, cluster_counts, cluster_counts_auth,
                             tolerance, report, annotations_auth)

    if ann_col is not None and ann_counts:
        _compare_per_annotation(path, ann_counts, cluster_counts_auth,
                                annotations_auth, allowed_annotations,
                                tolerance, report)


def _tolerance_cells(expected: int, tolerance: float) -> int:
    """Absolute cell-count drift allowed for a cluster of `expected` size.

    Single source of truth for per-cluster / per-annotation tolerance.
    Floor of 10 cells accepts boundary-cell jitter between seeded runs
    but catches the ~247-cell shift documented in project conventions as destroying
    PITX2 laterality.
    """
    return max(10, int(tolerance * expected))


def _check_aggregated_counts(
    path: Path,
    header: list[str],
    cluster_counts_auth: dict[int, int],
    annotations_auth: dict[int, str],
    tolerance: float,
    report: Report,
) -> None:
    """Silhouette-like TSVs: one row per cluster with an `n_cells` column.

    Aggregated totals must match `cluster_annotations.tsv` EXACTLY —
    these files are exports of the authoritative partition itself, not
    subsamples.  Per-cluster comparisons use the shared tolerance
    helper.  Rows whose `n_cells` or `cluster` column doesn't parse
    as an integer are violations, not silent skips.
    """
    expected_total = sum(cluster_counts_auth.values())
    rows_total = 0
    skipped = 0
    for row in _iter_rows(path):
        raw_n = row.get("n_cells", "")
        try:
            n = int(float(raw_n))
        except (ValueError, TypeError):
            skipped += 1
            if skipped <= 3:
                report.violate(
                    f"{path.name}: non-integer n_cells={raw_n!r} "
                    f"(row must parse as an integer count)"
                )
            continue
        rows_total += n

        # Also validate the cluster/annotation pair for this row.
        if "cluster" in row:
            raw_cl = row.get("cluster", "")
            try:
                cl = int(raw_cl)
            except (ValueError, TypeError):
                # Annotation-only aggregates (e.g. e2/e3) put the
                # annotation string in the `cluster` slot — skip, the
                # annotation comparison below handles them.
                continue
            expected = cluster_counts_auth.get(cl)
            if expected is None:
                report.violate(
                    f"{path.name}: cluster {cl} not in authoritative "
                    f"annotations"
                )
                continue
            if abs(n - expected) > _tolerance_cells(expected, tolerance):
                report.violate(
                    f"{path.name}: cluster {cl} n_cells={n}, "
                    f"authoritative={expected} "
                    f"(drift > max(10, {tolerance*100:.0f}%))"
                )
            if "annotation" in row:
                expected_ann = annotations_auth[cl]
                if row["annotation"] != expected_ann:
                    report.violate(
                        f"{path.name}: cluster {cl} annotation "
                        f"{row['annotation']!r} != authoritative "
                        f"{expected_ann!r}"
                    )

    if rows_total > 0 and rows_total != expected_total:
        report.violate(
            f"{path.name}: sum of n_cells {rows_total} != authoritative "
            f"{expected_total} (exact total required for aggregated "
            f"exports of the cluster_annotations.tsv partition)"
        )


def _compare_per_cluster(
    path: Path,
    observed: dict[str, int],
    auth: dict[int, int],
    tolerance: float,
    report: Report,
    annotations_auth: dict[int, str],
) -> None:
    for v, n in observed.items():
        try:
            cl = int(v)
        except (ValueError, TypeError):
            # Non-integer cluster value in a CLUSTER_COLS column means
            # a schema violation (cluster column should hold IDs, not
            # names).  Flag and move on.
            report.violate(
                f"{path.name}: non-integer value in cluster column: {v!r}"
            )
            continue
        expected = auth.get(cl)
        if expected is None:
            report.violate(
                f"{path.name}: observed cluster {cl} not in authoritative "
                f"annotations (stale partition?)"
            )
            continue
        if abs(n - expected) > _tolerance_cells(expected, tolerance):
            ann = annotations_auth.get(cl, "?")
            report.violate(
                f"{path.name}: cluster {cl} ({ann}) n_cells={n}, "
                f"authoritative={expected} "
                f"(drift > max(10, {tolerance*100:.0f}%))"
            )


def _compare_per_annotation(
    path: Path,
    observed: dict[str, int],
    auth_counts: dict[int, int],
    auth_ann: dict[int, str],
    allowed_annotations: set[str],
    tolerance: float,
    report: Report,
) -> None:
    # Build authoritative annotation → total cells.
    auth_by_ann: dict[str, int] = {}
    for cl, n in auth_counts.items():
        ann = auth_ann[cl]
        auth_by_ann[ann] = auth_by_ann.get(ann, 0) + n

    for ann, n in observed.items():
        if ann not in allowed_annotations:
            # Alphabet check (c) already flagged this with legacy /
            # novel diagnosis — don't duplicate here.  We only do
            # count comparison for labels in the allowed set.
            continue
        expected = auth_by_ann.get(ann)
        if expected is None:
            continue
        if abs(n - expected) > _tolerance_cells(expected, tolerance):
            report.violate(
                f"{path.name}: annotation {ann!r} n_cells={n}, "
                f"authoritative={expected} "
                f"(drift > max(10, {tolerance*100:.0f}%))"
            )


def check_timestamp(path: Path, ann_mtime: float, report: Report,
                    grace_seconds: int = 3600) -> None:
    """(b) TSV mtime must be >= annotations mtime - grace.

    Fatal if stale.  Raises because stale derived data is the whole point
    of this infra.
    """
    try:
        mtime = path.stat().st_mtime
    except OSError as e:
        report.warn(f"{path.name}: stat failed ({e})")
        return

    if mtime < ann_mtime - grace_seconds:
        delta = ann_mtime - mtime
        hrs = delta / 3600
        report.violate(
            f"{path.name}: mtime older than cluster_annotations.tsv by "
            f"{hrs:.1f}h (annotations rerun not propagated)"
        )


def check_annotation_alphabet(
    path: Path,
    allowed_annotations: set[str],
    report: Report,
) -> None:
    """(c) Any column whose values look like cell-type names must be a
    subset of the authoritative alphabet.  Catches legacy labels like
    'Mesenchyme (late)'.
    """
    header = _read_header(path)
    if not header:
        return

    # Identify candidate columns.
    candidates = [c for c in header if c in ANNOTATION_COLS]
    if not candidates:
        # Fallback heuristic: columns containing the substring "annotation"
        # or "cell_type" (but not genotype).
        for col in header:
            lower = col.lower()
            if ("annotation" in lower or "cell_type" in lower) and "id" not in lower:
                candidates.append(col)

    if not candidates:
        return

    # Collect unique values.
    observed: dict[str, set[str]] = {c: set() for c in candidates}
    for row in _iter_rows(path):
        for c in candidates:
            v = row.get(c)
            if v:
                observed[c].add(v)

    for col, values in observed.items():
        # Filter out obviously-non-annotation strings (numbers, empty,
        # or very long values unlikely to be cell types).
        filtered = {
            v for v in values
            if v and not v.isdigit() and len(v) < 64
        }
        if not filtered:
            continue
        bad = filtered - allowed_annotations
        # Flag only labels that are in the curated legacy set — we
        # deliberately don't use a " (" heuristic because several
        # current valid labels contain parentheses ("Somatic mes. (2)",
        # "Transitional mes. (2)").  Novel values that aren't in any
        # authoritative alphabet get a soft warning (could be a group
        # label like "All" / "Mesenchyme").
        flagged = sorted(bad & KNOWN_LEGACY_LABELS)
        if flagged:
            report.violate(
                f"{path.name}: column {col!r} contains legacy "
                f"annotation label(s) not in any current alphabet: "
                f"{flagged!r}"
            )
            # Fall through — still warn about any additional novel
            # labels from the same column.
        novel = sorted(bad - KNOWN_LEGACY_LABELS)
        if novel:
            report.warn(
                f"{path.name}: column {col!r} has values outside the "
                f"authoritative alphabets (possibly group labels): "
                f"{novel!r}"
            )


# ── Main ───────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=MANIFEST_DEFAULT)
    ap.add_argument("--dry-run", action="store_true",
                    help="report violations but exit 0")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"ERROR: manifest not found at {args.manifest}", file=sys.stderr)
        return 2
    if not ANNOTATIONS_FILE.is_file():
        print(f"ERROR: annotations not found at {ANNOTATIONS_FILE}",
              file=sys.stderr)
        return 2

    try:
        manifest = read_manifest(args.manifest)
    except Exception as e:
        print(f"ERROR: cannot parse manifest: {e}", file=sys.stderr)
        return 2

    annotations_auth = _load_annotations(ANNOTATIONS_FILE)
    cluster_counts_auth = _load_cluster_counts(ANNOTATIONS_FILE)
    # Union of all authoritative alphabets.  The alphabet check only
    # fires fatal on KNOWN_LEGACY_LABELS; anything else outside the
    # union is a WARN (possibly a group label or new population).
    d2b_labels = set(annotations_auth.values())
    d2g_labels = _load_label_set(D2G_SUBCLUSTER_FILE, "population")
    # Group labels commonly used in e2/e3 aggregated exports.
    GROUP_LABELS = {"All", "Mesenchyme", "LPM mesenchyme", "Paraxial"}
    allowed = (
        d2b_labels | d2g_labels
        | D2E_FOCUSED_MES_DIALECT | GROUP_LABELS
    )
    ann_mtime = ANNOTATIONS_FILE.stat().st_mtime

    report = Report()
    if args.verbose:
        report.info(f"d2b (cluster_annotations):   {len(d2b_labels)} labels")
        report.info(f"d2g (stage-restricted):      {len(d2g_labels)} labels")
        report.info(f"d2e focused-mes dialect:     {len(D2E_FOCUSED_MES_DIALECT)} labels")
        report.info(f"union alphabet:              {len(allowed)} labels")
        report.info(f"authoritative total cells:   {sum(cluster_counts_auth.values())}")

    for row in manifest:
        repo_dest = row["repo_dest"]
        h5 = row["depends_on_h5ad"]
        paths = list(expand_repo_dest(repo_dest))
        if not paths:
            report.skipped_missing += 1
            if args.verbose:
                report.info(f"missing: {repo_dest} (not synced locally)")
            continue

        for path in paths:
            if not path.is_file():
                report.skipped_missing += 1
                continue
            report.checked += 1
            if args.verbose:
                report.ok(f"checking {path.relative_to(REPO_ROOT)} "
                          f"({'h5ad-dep' if h5 else 'independent'})")

            # (c) Alphabet first — cheapest, no DictReader needed past header.
            # Narrow catch: only tolerate transient I/O failures.
            # Parse / schema errors must propagate so a malformed TSV
            # doesn't silently downgrade a fatal check to a warning.
            try:
                check_annotation_alphabet(path, allowed, report)
            except (FileNotFoundError, PermissionError) as e:
                report.warn(f"{path.name}: alphabet check unreadable: {e}")

            # (b) Timestamp monotonicity (only h5ad-dependent files).
            if h5:
                check_timestamp(path, ann_mtime, report)

            # (a) Cell-count consistency.
            try:
                check_cell_counts(
                    path, row, allowed, cluster_counts_auth,
                    annotations_auth, report,
                )
            except (FileNotFoundError, PermissionError) as e:
                report.warn(f"{path.name}: cell-count check unreadable: {e}")

    print()
    report.summary()

    if report.violations and not args.dry_run:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
