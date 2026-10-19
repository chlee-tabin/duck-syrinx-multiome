#!/usr/bin/env python3
"""Centralized, portable path roots for duck-syrinx-multiome (Python).

Why: ~127 pipeline scripts hardcode ``SCRATCH = "/n/scratch/users/c/cl266/FC_08656"``
(the username ``cl266`` is baked into 1,407 path literals), plus the standby and
canonical-reference roots.  That makes the repo non-portable: no other user — and not
even the same user on relocated/renamed scratch — can run it without editing source.
Full audit + rollout plan: ``docs/audit_hardcoded_paths.md``.

This module centralizes the roots.  **Defaults equal the historically hardcoded
paths, so existing cluster runs are byte-for-byte unchanged**; set the env var to run
as a different user or on relocated storage.

    REPO_ROOT          this repo checkout    derived from __file__, never hardcoded (Path)
    SCRATCH_USER       user scratch parent   $DUCK_SCRATCH_USER                      (str)
    SCRATCH_ROOT       cluster working tree  $DUCK_SCRATCH_ROOT                      (str)
    STANDBY_DATA_ROOT  standby raw-data tree $DUCK_STANDBY_DATA_ROOT                 (str)
    STANDBY_REPO       standby git repo      $DUCK_STANDBY_REPO                      (str)
    REFERENCE_ROOT     canonical genome refs $DUCK_REFERENCE_ROOT (the stable 'cdw') (str)
    HOME_ROOT          cluster user home     $DUCK_HOME                              (str)

Usage (drop-in for ``SCRATCH = "..."``)::

    from _paths import SCRATCH_ROOT
    SCRATCH = SCRATCH_ROOT
    df = pd.read_csv(f"{SCRATCH}/figures/data/foo.tsv")

``_paths`` lives in ``pipeline/`` and resolves because a script run as
``python3 pipeline/<step>.py`` has its own directory on ``sys.path[0]``.

From a **subdirectory that the scratch rsync carries** (``FC_08656/``,
``figures/_shared/``) ``pipeline/`` is not on ``sys.path``; use the depth-agnostic
bootstrap, which walks up to whichever tree the file is executing from (repo checkout
*or* the scratch mirror, where ``pipeline/`` is rsynced alongside).  **NOT
``pipeline/chick/``** — see below::

    import sys
    from pathlib import Path as _P
    _r = _P(__file__).resolve()
    while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
        _r = _r.parent
    sys.path.insert(0, str(_r / "pipeline"))
    from _paths import SCRATCH_ROOT

**When NOT to import this module.** Some Python is ferried to standalone scratch trees
that do not carry ``pipeline/``, so the walk-up cannot find ``_paths.py`` and the import
would be a hard, new failure.  Two such trees exist:

  * ``old_duck_cr9/analysis/`` — ferried per file;
  * **all of ``pipeline/chick/``** — the chick steps run from
    ``$DUCK_SCRATCH_USER/chick_demux{,_B}/pipeline[/chick]/``, a *sibling* of
    ``FC_08656/`` that the scratch sync never reaches.  Despite living
    under ``pipeline/`` in the repo, ``pipeline/chick/**`` is tier B;
    tier A raises ``ModuleNotFoundError`` there.

Those files use the *same* env-with-default idiom inline instead — identical semantics,
no file dependency::

    SCRATCH = (os.environ.get("DUCK_SCRATCH_ROOT")
               or "/n/scratch/users/c/cl266/FC_08656")

Write it in exactly that shape.  ``os.environ.get(name, default)`` is NOT equivalent and
is rejected by the guardrail: it returns ``""`` for a set-but-EMPTY variable.  The
parentheses are load-bearing — ``or`` binds looser than ``+``, so the unparenthesised
form would silently drop the tail of ``get(...) + "/sub"`` whenever the variable is set.
The R equivalent is ``.duck_env("DUCK_X", "<literal>")`` with the one-line ``nzchar``
body from ``scripts/_paths.R``.

See ``docs/audit_hardcoded_paths.md`` ("Two adoption tiers") for the rule.
"""
from __future__ import annotations

import os
from pathlib import Path

# Repo root is ALWAYS derived from this file's location (pipeline/_paths.py -> repo/).
# Never an env var: the running file's location is the authoritative checkout path.
REPO_ROOT: Path = Path(__file__).resolve().parent.parent

# Environment-specific roots: env override, behavior-preserving default.
#
# `or default` (not `os.environ.get(name, default)`) so that a set-but-EMPTY variable
# falls back to the default instead of yielding "" — which would silently turn
# f"{SCRATCH_ROOT}/phase2" into "/phase2", a wrong path rather than an error.  This
# matches `: "${DUCK_X:=default}"` in _paths.sh (`:=` also treats null as unset) and
# nzchar() in _paths.R.  Since the tier-B inline sites use the SAME
# normalizing forms — `(os.environ.get(...) or …)` and `.duck_env(...)` — and the
# guardrail rejects the two-argument variants *whose default is a path literal*, which
# is every occurrence in this tree (re-measured: 0 non-normalizing inline sites
# in 653).  A two-argument form whose default is a VARIABLE carries no literal, so the
# guardrail's line scan cannot see it — the contract "a DUCK_* root is either unset or a
# non-empty absolute path" is enforced for the literal forms and conventional otherwise.
# SCRATCH_USER is the *parent* of SCRATCH_ROOT; the sibling chick / old-duck scratch
# trees hang off it, so it is a root in its own right (audit doc, mechanism item 3).
SCRATCH_USER: str = os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266"
SCRATCH_ROOT: str = os.environ.get("DUCK_SCRATCH_ROOT") or f"{SCRATCH_USER}/FC_08656"
# The standby *raw-data* tree is one level above the standby git checkout; backups and
# the cellranger-arc _v210 mirrors live there.  Off-scratch by definition: parameterize,
# never relativize.
STANDBY_DATA_ROOT: str = (
    os.environ.get("DUCK_STANDBY_DATA_ROOT")
    or "/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data"
)
STANDBY_REPO: str = (
    os.environ.get("DUCK_STANDBY_REPO")
    or f"{STANDBY_DATA_ROOT}/R4.4.0/duck-syrinx-multiome"
)
# Category C — the stable `cdw` reference home.  Parameterized for portability but it
# stays ABSOLUTE: it is deliberately NOT scratch-relative (docs/audit_hardcoded_paths.md).
REFERENCE_ROOT: str = (
    os.environ.get("DUCK_REFERENCE_ROOT")
    or "/n/groups/tabin/chlee/singlecell/reference"
)
# Cluster user home (conda installs, R user lib).  NOT `os.path.expanduser("~")`: many
# call sites build a *remote* path from a local machine, where `~` is the wrong home.
HOME_ROOT: str = os.environ.get("DUCK_HOME") or "/home/cl266"

__all__ = [
    "REPO_ROOT",
    "SCRATCH_USER",
    "SCRATCH_ROOT",
    "STANDBY_DATA_ROOT",
    "STANDBY_REPO",
    "REFERENCE_ROOT",
    "HOME_ROOT",
]


if __name__ == "__main__":  # quick self-report / sanity check
    print(f"REPO_ROOT         = {REPO_ROOT}")
    print(f"SCRATCH_USER      = {SCRATCH_USER}")
    print(f"SCRATCH_ROOT      = {SCRATCH_ROOT}")
    print(f"STANDBY_DATA_ROOT = {STANDBY_DATA_ROOT}")
    print(f"STANDBY_REPO      = {STANDBY_REPO}")
    print(f"REFERENCE_ROOT    = {REFERENCE_ROOT}")
    print(f"HOME_ROOT         = {HOME_ROOT}")
