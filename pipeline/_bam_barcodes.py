"""Shared BAM-locator + barcode-suffix helpers for d2i / d2j scripts.

Three pipeline scripts (d2i_pitx2_rna_pileup, d2j_pitx2_isoform_per_genotype,
d2j_pitx2_isoform_percell) previously carried byte-similar copies of
`find_bam`, `detect_bam_suffix`, and `adjust_barcodes_for_bam`.  Keeping
them in one place prevents divergence — e.g., the PILOTARC cellranger-arc
GEX paths used to live only in d2j_percell and were missing from the
other two.

The unified `find_bam` checks both cellranger-count and cellranger-arc
GEX candidates by default; d2i and d2j_per_genotype only need the count
candidates in practice, but running the arc check costs nothing (fast
os.path.isfile calls) and removes the divergence.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import subprocess
import tempfile
from collections import Counter

SCRATCH = os.environ.get(
    "FC08656_SCRATCH", SCRATCH_ROOT
)


def _base_barcode(bc: str) -> str:
    """Strip every trailing `-<suffix>` segment down to the 16-bp base.

    h5ad barcodes are sometimes double-suffixed (`ATCG...-1-1`) after a
    round-trip rebuild; `bc.split("-", 1)[0]` reduces both single- and
    double-suffix forms to the 16-bp base sequence.  Shared between
    `adjust_barcodes_list` and `write_adjusted_barcode_file` so the two
    helpers cannot drift (previously one used `split("-")[0]` and the
    other `rsplit("-", 1)[0]`, which disagreed on double-suffix cases
    and silently produced non-matching BAM CB queries).
    """
    return bc.split("-", 1)[0]


def find_bam(library: str, scratch: str = SCRATCH) -> str | None:
    """Find BAM for a library across cellranger-count + cellranger-arc paths.

    Checks, in order:
      1. cellranger-count standard BAM
      2. cellranger-count modsuffix2 variant (re-suffixed barcodes)
      3. cellranger-count modsuffix2 in bam_staging/ (PILOT-era BAMs)
      4. staged cellranger-arc v2.1.0 GEX BAM (post-purge recovery, all libs)
      5. cellranger-arc GEX BAM in bare arc dirs (stale v2.0.0 PILOT copies —
         kept LAST so a staged v2.1.0 BAM always wins if both exist)
      6. cellranger-arc GEX modsuffix2 variant (same bare dirs)

    A companion .bai index is required alongside any candidate.  Returns
    None (not raises) so callers can decide whether a missing BAM is
    fatal or skippable — d2j currently collects the skips and raises
    at end-of-run once every library has been attempted.
    """
    candidates = [
        f"{scratch}/{library}_GEX/outs/possorted_genome_bam.bam",
        f"{scratch}/{library}_GEX/outs/possorted_genome_bam_modsuffix2.bam",
        f"{scratch}/{library}_GEX/outs/bam_staging/possorted_genome_bam_modsuffix2.bam",
    ]
    # cellranger-arc v2.1.0 GEX BAMs staged from the standby _v210 mirror after
    # the 2026-08 scratch purge removed every cellranger-count BAM.  Same GEX
    # reads (the multiome GEX modality of the same library, aligned by
    # cellranger-arc 2.1.0); CB suffixes are handled by detect_bam_suffix.
    # Staging: each library's cellranger-arc GEX BAM (GEO: <ARC>_gex_possorted_bam.bam)
    #   -> <scratch>/staged_v210_gex_bams/<ARC>.gex_possorted_bam.bam{,.bai}
    # where <ARC> = library.replace("GEX", "ARC") (ARC0S-duck covers GEX0S-duck).
    arc_v210 = library.replace("GEX", "ARC")
    if arc_v210 != library:
        candidates.append(
            f"{scratch}/staged_v210_gex_bams/{arc_v210}.gex_possorted_bam.bam"
        )
    arc_name = library.replace("PILOTGEX", "PILOTARC")
    if arc_name != library:
        candidates += [
            f"{scratch}/{arc_name}/outs/gex_possorted_bam.bam",
            f"{scratch}/{arc_name}/outs/gex_possorted_bam_modsuffix2.bam",
        ]
    for p in candidates:
        if os.path.isfile(p) and os.path.isfile(p + ".bai"):
            return p
    return None


def detect_bam_suffix(
    bam_path: str,
    n_reads: int = 50_000,
    region: str | None = None,
) -> str:
    """Detect the barcode suffix used in a BAM by majority vote over CB tags.

    Barcodes in cellranger BAMs are tagged as `CB:Z:XXXXXXXXXXXXXXXX-N`
    where `N` is a library-specific suffix (often `-1`, but modsuffix2
    BAMs use values like `-2`, `-3`, `-6`, etc.).  Streams the first
    `n_reads` records (optionally within `region`) and returns the most
    common suffix.

    `n_reads=50_000` (up from 200) prevents a minority leading-contig
    contamination — e.g., after a PILOTARC cellranger-arc merge — from
    winning the vote: older cellranger BAMs often have small blocks of
    unbarcoded-or-differently-suffixed reads at the file head.  Passing
    `region` (e.g., `"NC_051775.1:3950000-3980000"` for PITX2) further
    hardens the vote by sampling from a cell-dense locus.

    Raises `RuntimeError` (not a silent `-1` default) when no CB tags
    are found — previously a sampled-empty-region BAM would silently
    fall through to `-1` and produce a zero-overlap barcode list against
    the actual BAM suffix.  The aggregate zero-count guards in d2i / d2j
    would catch all-libraries-broken, but a single misrouted library
    could slip through.
    """
    region_spec = f" {region}" if region else ""
    cmd = (
        f"samtools view {bam_path}{region_spec} | "
        f"head -{n_reads} | "
        f"grep -oP 'CB:Z:\\S+'"
    )
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=120,
    )
    suffixes: Counter[str] = Counter()
    for line in result.stdout.strip().split("\n"):
        if not line.startswith("CB:Z:"):
            continue
        bc = line[5:]
        parts = bc.rsplit("-", 1)
        if len(parts) == 2:
            suffixes[f"-{parts[1]}"] += 1
    if not suffixes:
        raise RuntimeError(
            f"detect_bam_suffix: no CB:Z tags found in the first "
            f"{n_reads} reads of {bam_path}"
            f"{' ' + region if region else ''}.  BAM may be truncated, "
            f"un-barcoded, or pointing at a region with no CB-tagged "
            f"reads.  Inspect with `samtools view {bam_path} | head`."
        )
    return suffixes.most_common(1)[0][0]


def adjust_barcodes_list(barcodes: list[str], bam_suffix: str) -> list[str]:
    """Rewrite a list of h5ad barcodes (-1 suffix) to match a BAM's suffix.

    h5ad barcodes typically carry `-1` (and sometimes double suffixes
    like `-1-1` from round-trip rebuilds); `_base_barcode` strips all
    trailing suffixes down to the 16-bp base sequence and then
    `bam_suffix` is re-appended.
    """
    return [f"{_base_barcode(bc)}{bam_suffix}" for bc in barcodes]


def write_adjusted_barcode_file(
    barcode_file: str, bam_suffix: str, tmp_dir: str = "/tmp",
) -> str:
    """Materialize an adjusted-suffix barcode list to a temp file.

    For callers (d2i) that need a file path to hand to
    `samtools view -D CB:<file>`.  Returns the original `barcode_file`
    unchanged when `bam_suffix == "-1"` (no rewrite needed).  Otherwise
    returns a new temp file path.  Caller is responsible for unlinking
    the temp file when done (d2i already tracks these in `tmp_files`).

    Uses `_base_barcode` (same normalization as `adjust_barcodes_list`)
    so both helpers produce identical base sequences for double-suffix
    inputs; divergence between them previously caused silently-empty
    BAM CB queries on double-suffix h5ad rebuilds.
    """
    if bam_suffix == "-1":
        return barcode_file
    adjusted = tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", dir=tmp_dir, delete=False, prefix="bc_"
    )
    with open(barcode_file) as f:
        for line in f:
            bc = line.strip()
            if bc:
                adjusted.write(f"{_base_barcode(bc)}{bam_suffix}\n")
    adjusted.close()
    return adjusted.name
