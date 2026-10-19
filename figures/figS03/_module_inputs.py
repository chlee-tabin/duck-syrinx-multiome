"""_module_inputs.py — the HH21 module positive control's input fingerprint, computed in ONE place.

`module_inputs_sha` identifies the module-test run that produced a set of artifacts. It covers
BOTH files the NB-GLMM consumes — `hh21_hh1716_response_counts.tsv` (the counts) and
`hh21_hh1716_halfmeta.tsv` (every offset and the paired flags) — because fingerprinting the
counts alone let a stale same-tier halfmeta pass every other check while supplying another run's
offsets.

Consumers must do TWO things, and the second is the easy one to miss:

  1. check that the artifacts they read agree with each other on the stamp, AND
  2. RECOMPUTE the stamp from the input files on disk and check the artifacts against that.

Without (2) the guard compares two derived outputs to each other: a wholly stale set — stale
scores, stale GLMM, and a label-swap table derived from those stale scores — all carry the same
old digest and satisfy every cross-check, even though the real inputs have since changed. That
is the "two identically-wrong values" failure mode.

`hh21_hh1716_nb_glmm.R` recomputes the same digest independently in R (digest on the
concatenated raw bytes); the two implementations are verified to agree.
"""
import hashlib
import os

# A bare 12-hex-char digest is sometimes NUMERIC-LOOKING (`13038807e193` parses as 1.3e+200),
# and both pandas and R's read.delim then coerce the whole column to float — so a perfectly
# fresh, perfectly consistent artifact set is REFUSED, with an error no rerun can clear. That
# happens for roughly one digest in 280. Tagging the value at the source keeps every column
# textual in both languages, instead of needing dtype= pinned at every read site.
SHA_PREFIX = "sha:"


def _tag(hexdigest):
    return SHA_PREFIX + hexdigest


COUNTS_NAME = "hh21_hh1716_response_counts.tsv"
HALFMETA_NAME = "hh21_hh1716_halfmeta.tsv"
SCAN_NAME = "grn_memory_node_promoter_scan.tsv"


def scan_sha(data_dir, digits=12):
    """sha256 of the committed promoter scan — the SOURCE the module's arms are projected from.

    The module artifacts are consistent with each other and with their immediate inputs; that
    says nothing about whether the scan those inputs were projected from is still current. A
    scan re-run updates panel H while the module panel keeps a coherent OLD projection, and every
    module_inputs_sha check passes because its inputs are consistently old.
    """
    path = os.path.join(data_dir, SCAN_NAME)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{path} is missing — cannot confirm the module artifacts were "
                                "projected from the current scan.")
    with open(path, "rb") as fh:
        return _tag(hashlib.sha256(fh.read()).hexdigest()[:digits])


def module_inputs_sha(data_dir, digits=12):
    """sha256 over the two GLMM input files, byte-concatenated (counts first)."""
    parts = []
    for name in (COUNTS_NAME, HALFMETA_NAME):
        path = os.path.join(data_dir, name)
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"{path} is missing — cannot verify that the committed module artifacts match "
                "their inputs. Regenerate with pipeline/run_hh21_module_v2.sh.")
        with open(path, "rb") as fh:
            parts.append(fh.read())
    return _tag(hashlib.sha256(b"".join(parts)).hexdigest()[:digits])


# Which committed artifacts must carry which stamps. Keeping this HERE rather than in each
# caller is the point: the previous signature took `stamps` and an OPTIONAL `scan_stamps`, so a
# caller could pass one and silently halve the guard — which one of them did.
_REQUIRE = {
    "hh21_hh1716_module_scores.tsv": ("module_inputs_sha", "scan_sha"),
    "hh21_hh1716_module_halfscores.tsv": ("module_inputs_sha", "scan_sha"),
    "hh21_hh1716_module_labelswap_stats.tsv": ("module_inputs_sha", "scan_sha"),
    "hh21_power_calibration.tsv": ("scan_sha",),
    "hh21_hh1716_nb_glmm.tsv": ("module_inputs_sha",),   # the R model never reads the scan
}


def assert_fresh(data_dir, artifacts, label):
    """Raise unless every artifact's stamps match digests RECOMPUTED from their sources.

    `artifacts` maps a display name -> a DataFrame already read from one of the files in
    `_REQUIRE` (keyed by basename via the `_source` attribute set here, or passed as a
    (frame, basename) pair). Comparing artifacts only to each other passes a wholly stale set;
    the digests below come from the files those artifacts describe.
    """
    import pandas as pd  # local: keeps this module importable without pandas for the digests

    current = {"module_inputs_sha": module_inputs_sha(data_dir), "scan_sha": scan_sha(data_dir)}
    problems = []
    for name, item in artifacts.items():
        frame, basename = item if isinstance(item, tuple) else (item, name)
        required = _REQUIRE.get(basename)
        if required is None:
            raise KeyError(f"{basename!r} is not a stamped module artifact — add it to "
                           "_module_inputs._REQUIRE or do not pass it here")
        for col in required:
            if col not in frame.columns:
                raise KeyError(f"{name} has no `{col}` column — regenerate it "
                               "(pipeline/run_hh21_module_v2.sh, then the local figS03 steps)")
            got = set(pd.Series(frame[col]).astype(str))
            if got != {current[col]}:
                problems.append(f"{name}.{col}={sorted(got)} (current {current[col]})")
    if problems:
        raise ValueError(
            f"{label}: artifact(s) do not match their sources — " + "; ".join(problems) +
            f". `module_inputs_sha` is sha256({COUNTS_NAME} + {HALFMETA_NAME}); `scan_sha` is "
            f"sha256({SCAN_NAME}). Re-run the module test on the current scan, then "
            "Rscript figures/figS03/hh21_hh1716_nb_glmm.R and the figS03 panels."
        )


if __name__ == "__main__":  # ponytail: one runnable check for the cross-language contract
    import subprocess

    here = os.path.dirname(os.path.abspath(__file__))
    data = os.path.join(here, "..", "data")
    py = module_inputs_sha(data)
    r_expr = (
        '.rb <- function(p) readBin(p, "raw", file.info(p)$size); '
        'cat(paste0("sha:", substr(digest::digest(c(.rb(commandArgs(TRUE)[1]), '
        '.rb(commandArgs(TRUE)[2])), algo = "sha256", serialize = FALSE), 1, 12)))'
    )
    r = subprocess.run(["Rscript", "-e", r_expr,
                        os.path.join(data, COUNTS_NAME), os.path.join(data, HALFMETA_NAME)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"R side failed (needs Rscript + digest): {r.stderr.strip()}")
    # serialize=FALSE is load-bearing: serialize=TRUE hashes an R object header too and the two
    # languages silently disagree, which would lock every consumer out with a false positive.
    assert r.stdout.strip() == py, f"R {r.stdout.strip()} != Python {py}"
    print(f"OK module_inputs_sha: Python and R agree ({py}); scan_sha {scan_sha(data)}")
