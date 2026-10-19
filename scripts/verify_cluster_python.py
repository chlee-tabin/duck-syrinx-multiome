#!/usr/bin/env python3
"""verify_cluster_python.py — compile cluster-executed Python against an OLDER interpreter.

WHY THIS EXISTS. On 2026-08-02 a cluster job died in seconds on

    f"NB2, alpha={float(nb.params["alpha"]):.4g}"

Nested same-quote f-strings are PEP 701, new in Python 3.12. The workstation runs 3.12 and the
cluster's `snapatac2` env runs 3.11.15, so `python3 -m py_compile` passed locally and the code was
synced, submitted, queued and executed before anything complained. A local syntax check that uses
the local interpreter cannot, even in principle, catch a version skew in that direction.

WHAT IT DOES. Byte-compiles every Python file in the trees that the repo-to-scratch sync
transfers to scratch — i.e. exactly the files the cluster runs — using the OLDEST interpreter it
can find that is still >= the cluster's. Anything that interpreter accepts, the cluster accepts,
so the check is conservative: it can produce a false alarm on genuinely-new syntax the cluster
would in fact support, never a false pass. A false alarm is a one-line fix; a false pass costs a
queue slot and a round trip.

`ast.parse(..., feature_version=(3, 11))` does NOT work for this — f-strings are handled by the
tokenizer before feature_version is consulted, so it accepts the very construct above. Only a real
older interpreter catches it.

Configure with CLUSTER_PY (path to an interpreter) if the autodetected one is wrong. If no
suitable older interpreter exists the check SKIPS with a warning rather than failing the commit —
it is a safety net, not a dependency.

Run:  make check-cluster-python      (also runs in .githooks/pre-commit on staged .py)
Exit: 0 = all clear or skipped; 1 = a file the cluster would reject.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The trees the repo-to-scratch sync copies to scratch. Keep in step with its SYNC_TREES.
# Recursive trees. FC_08656 is deliberately NOT here: the rsync takes only its TOP-LEVEL
# `FC_08656/*.py` glob, so recursing would gate archive_genotyping/ and other subdirectories that
# never reach the cluster. (Doing so surfaced a genuinely broken archived script — a real
# IndentationError on every interpreter, not a version skew — which is worth fixing on its own
# terms but is not this gate's business.)
CLUSTER_TREES = ["pipeline", "figures/figS03", "figures/_shared",
                 "figures/fig01", "figures/fig01_supp", "figures/fig02"]
CLUSTER_GLOBS = ["FC_08656/*.py"]        # non-recursive, matching the rsync

# The cluster's snapatac2 env, measured 2026-08-02. Anything <= this is a valid conservative gate.
CLUSTER_VERSION = (3, 11)


def find_interpreter() -> tuple[str, tuple[int, int]] | tuple[None, None]:
    """Oldest interpreter that is <= the cluster's, preferring the closest match."""
    if env := os.environ.get("CLUSTER_PY"):
        return env, _version_of(env)
    cands: list[tuple[tuple[int, int], str]] = []
    for pat in ("/usr/local/bin/python3.*", "/opt/homebrew/bin/python3.*", "/usr/bin/python3.*"):
        for path in glob.glob(pat):
            if path.endswith("-config") or not os.access(path, os.X_OK):
                continue
            v = _version_of(path)
            if v and v <= CLUSTER_VERSION:
                cands.append((v, path))
    if not cands:
        return None, None
    v, path = max(cands)          # closest to the cluster from below
    return path, v


def _version_of(path: str) -> tuple[int, int] | None:
    try:
        out = subprocess.run([path, "-c", "import sys;print(sys.version_info[0],sys.version_info[1])"],
                             capture_output=True, text=True, timeout=15)
        a, b = out.stdout.split()
        return int(a), int(b)
    except Exception:
        return None


def main(argv: list[str]) -> int:
    files = [f for f in argv[1:] if f.endswith(".py")]
    if not files:
        for tree in CLUSTER_TREES:
            files += glob.glob(os.path.join(ROOT, tree, "**", "*.py"), recursive=True)
        for g in CLUSTER_GLOBS:
            files += glob.glob(os.path.join(ROOT, g))
    files = sorted({os.path.abspath(f) for f in files if os.path.exists(f)})
    # only gate files that actually live in a cluster-executed tree
    roots = tuple(os.path.join(ROOT, t) + os.sep for t in CLUSTER_TREES)
    flat = {os.path.abspath(f) for g in CLUSTER_GLOBS for f in glob.glob(os.path.join(ROOT, g))}
    files = [f for f in files if f.startswith(roots) or f in flat]
    if not files:
        print("verify_cluster_python: no cluster-executed Python in scope")
        return 0

    interp, ver = find_interpreter()
    if not interp:
        print(f"verify_cluster_python: SKIPPED — no interpreter <= {CLUSTER_VERSION} found. "
              f"The cluster runs {CLUSTER_VERSION[0]}.{CLUSTER_VERSION[1]}; a newer local "
              f"interpreter cannot catch syntax it accepts but the cluster rejects. "
              f"Set CLUSTER_PY=/path/to/python3.11 to enable this check.", file=sys.stderr)
        return 0

    bad = []
    for f in files:
        r = subprocess.run([interp, "-m", "py_compile", f], capture_output=True, text=True)
        if r.returncode != 0:
            bad.append((f, (r.stderr or r.stdout).strip().splitlines()[-1] if (r.stderr or r.stdout) else "?"))

    label = f"{os.path.basename(interp)} ({ver[0]}.{ver[1]}, cluster runs "\
            f"{CLUSTER_VERSION[0]}.{CLUSTER_VERSION[1]})"
    if bad:
        print(f"verify_cluster_python: FAILED against {label}")
        for f, msg in bad:
            print(f"  {os.path.relpath(f, ROOT)}: {msg}")
        print("\nThese files run ON THE CLUSTER. A newer local interpreter accepts them; the "
              "cluster will not, and the job dies after it is queued and scheduled.")
        return 1
    print(f"verify_cluster_python: OK ({len(files)} file(s) checked against {label})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
