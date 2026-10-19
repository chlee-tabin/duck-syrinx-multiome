#!/usr/bin/env python3
"""Guardrail: fail on a NEW single-user hardcoded root literal.

The repo was re-rooted onto env-overridable roots (``pipeline/_paths.py``,
``scripts/_paths.sh``, ``scripts/_paths.R``; audit + rationale in
``docs/audit_hardcoded_paths.md``).  Without a gate the single-user lock-in creeps
straight back the next time someone writes ``SCRATCH="/n/scratch/users/c/cl266/..."``,
so this check is wired into ``.githooks/pre-commit`` (staged files) and ``make
check-portable-paths`` (whole repo).

It does NOT ban the literal outright — several contexts legitimately keep it, and
banning those would either break behavior or be unenforceable:

  * the ``_paths.*`` helpers            — they ARE the defaults
  * a tier-B inline default             — ``: "${DUCK_X:=<literal>}"`` /
    ``(os.environ.get("DUCK_X") or "<literal>")`` / ``.duck_env("DUCK_X", "<literal>")``
    is the sanctioned portable idiom for files that cannot reach a helper.  Only these
    EMPTY-NORMALIZING forms are sanctioned — see ``SANCTIONED_RE``.
  * ``#SBATCH`` directives              — SLURM parses them before any shell
    expansion, so a variable there simply does not resolve.  (`the cluster submission wrapper`
    passes ``sbatch --output/--error`` built from the resolved roots, which override
    the in-file directive; a direct ``sbatch`` still needs ``-o``/``-e``.)
  * comments and python docstrings      — documentation, not a resolved path
  * ``archive*/`` and ``legacy_2024_notebooks/`` — frozen historical trees
  * ``NARROW_EXEMPTIONS`` — a hand-audited, file+line-content allowlist.  Today it holds
    exactly the three ``sed`` **patterns** in ``FC_08656/fix_hardcoded_paths.sh`` (not distributed): the
    literals there are the search side of a path-rewriting tool, so replacing them
    breaks what the tool does.

Anything else is a regression: a raw root literal in executable code.

**Known limit.** This is a LINE regex, so a literal split across concatenation
(``"/n/scratch/users/c/" + "cl266/FC_08656"``) is not detected.  Catching that needs a
per-language constant-folding pass.

**Deliberately NOT exempt** (both were blanket passes once, and both protected nothing
while leaving a real hole open):

  * ``<<'EOF'`` quoted heredocs — the shell does not expand there, but the *embedded*
    interpreter (python/R) still reads the sanctioned ``os.environ.get("DUCK_…", …)`` /
    ``Sys.getenv("DUCK_…", …)`` idiom perfectly well.  Exempting the whole heredoc let a
    literal split-brain sit inside embedded code unnoticed.
  * ``'single quotes'`` in shell — see ``NARROW_EXEMPTIONS`` above; the blanket version
    waved through any literal that merely happened to be single-quoted.

Usage:
    scripts/verify_portable_paths.py            # staged files (pre-commit; reads the INDEX)
    scripts/verify_portable_paths.py --all      # whole repo (working tree)
    scripts/verify_portable_paths.py FILE...    # explicit files (working tree)
"""
from __future__ import annotations

import ast
import io
import os
import re
import subprocess
import sys
import tokenize

# Any user's roots, not just the historical one — the point is to stop ANY
# single-user path from being baked in again.  The last alternative catches the
# *local-Mac* roots (a developer's checkout / Dropbox tree): 9 of those pointed at a
# git worktree that no longer exists, i.e. they were already broken.
ROOT_RE = re.compile(
    r"/n/scratch/users/[a-z0-9]/[a-z0-9_]+"
    r"|/n/standby/hms/[a-z]+/[a-z]+/compute/[a-z0-9_]+"
    # the whole per-user /n/groups dir, not only its singlecell/reference tail — the
    # tail is what REFERENCE_ROOT covers, but anything else under it is equally
    # single-user.  The LAB component is generic too: pinning it to `tabin` let
    # /n/groups/otherlab/<user>/… through.  `shared_databases` is excluded — it is the
    # institution-wide read-only namespace (bowtie2 indexes etc.), shared by everyone
    # and therefore already portable.
    r"|/n/groups/(?!shared_databases/)[a-z0-9_]+/[a-z0-9_]+"
    # raw BCL/FASTQ home (transfer-node only): /n/files/Genetics/Tabin_Lab-<user>
    r"|/n/files/[A-Za-z_]+/[A-Za-z_]+-[A-Za-z0-9_]+"
    r"|/home/[A-Za-z][A-Za-z0-9._-]{1,30}"
    r"|/Users/[A-Za-z][A-Za-z0-9._-]{0,30}"
)

# `.sbatch` is included: those ARE job scripts (16 of them carry root literals), they
# were simply invisible to an extension filter built around .sh.
MISSING = object()      # tracked, but not present in the worktree
EXTS = (".py", ".sh", ".sbatch", ".R", ".r")
HELPERS = {"_paths.py", "_paths.sh", "_paths.R"}
FROZEN_RE = re.compile(r"(^|/)(archive[A-Za-z0-9_]*|legacy_2024_notebooks)(/|$)")
SHELL_EXTS = (".sh", ".sbatch")

# The sanctioned portable idioms, per language.  Each one must NORMALIZE a set-but-EMPTY
# variable back to the default (`:=`, `or`, `nzchar`), so `DUCK_SCRATCH_ROOT=` can never
# silently turn f"{SCRATCH_ROOT}/phase2" into "/phase2".  The two-argument forms
# `os.environ.get("DUCK_X", <lit>)` / `Sys.getenv("DUCK_X", <lit>)` are therefore NOT
# sanctioned: they return "" for a set-but-empty variable.
SANCTIONED_RE = re.compile(
    r':\s+"\$\{DUCK_[A-Z_]+:=[^"]*\}"'              # bash   : "${DUCK_X:=...}"
    # python  (os.environ.get("DUCK_X") or <default>)  — the `)` before `or` is what
    # distinguishes the normalizing form from the positional-default one.
    # The trailing lookahead (not a consuming `\S`) matters: strip_sanctioned() swallows
    # the default that begins exactly at the end of this match, so the match must stop
    # BEFORE the default's opening quote rather than eating it.
    r'|_?os\.environ\.get\(\s*["\']DUCK_[A-Z_]+["\']\s*\)\s*or\s+(?=\S)'
    r'|\.duck_env\(\s*"DUCK_[A-Z_]+"\s*,'           # R      .duck_env("DUCK_X", ...)
)

# Hand-audited exemptions, matched on (path, substring-of-the-line) so they cannot
# silently widen: a literal that is the *search pattern* of a path-rewriting tool.
# ONE entry today — the file key maps to a substring that covers the three `sed -i` lines
# in it.  A SECOND path key is a blanket exemption for another whole file; add one
# only with explicit justification.
# Keyed on the SPECIFIC literals, not just the `sed -i 's|` shape: keying on the shape
# alone exempted any NEW line added to that file, whatever path it carried.
NARROW_EXEMPTIONS = {
    # This file: the table below necessarily QUOTES the literals it exempts.  Keyed on
    # the `sed -i 's|` shape plus the three leading-quote refbuild shapes below, each of
    # which appears verbatim on its own defining line (self-covering) and on the
    # refbuild entry's lines — so none can widen to anything else here.
    "scripts/verify_portable_paths.py": (
        "\"sed -i 's|",
        '"cd /n/groups/tabin/chlee/singlecell/reference/',
        '"--fasta=/n/groups/tabin/chlee/singlecell/reference/',
        '"--gene=/n/groups/tabin/chlee/singlecell/reference/',
    ),
    "FC_08656/fix_hardcoded_paths.sh": (  # (not distributed)
        "sed -i 's|cd /n/scratch/users/c/cl266/FC_08656|",
        "sed -i 's|cd /n/scratch/users/c/cl266/FC_08858|",
        "sed -i 's|source \"/home/cl266/anaconda3/etc/profile.d/conda.sh\"|",
    ),
    # Archival provenance: the VERBATIM 2023 mkref record (SLURM job 19590176) for the
    # chick Group-B GEX reference cr7.reference.gg6a.ncbi.EGFP.  A historical record,
    # not runnable pipeline code — rewriting its paths would falsify provenance
    # (see pipeline/chick/refbuild/README.md).  Keyed on the three executable literals;
    # any NEW hardcoded line added to the file is still checked.
    "pipeline/chick/refbuild/20231014_01_build_gg6a_ncbi.sh": (
        "cd /n/groups/tabin/chlee/singlecell/reference/",
        "--fasta=/n/groups/tabin/chlee/singlecell/reference/",
        "--gene=/n/groups/tabin/chlee/singlecell/reference/",
    ),
}


def py_docstring_lines(src: str) -> set:
    """1-based lines covered by a genuine DOCSTRING (module / class / def head).

    Deliberately NOT "every triple-quoted string": a triple-quoted string can be an
    executable payload (an embedded shell or R program handed to `subprocess.run`), and
    exempting those is the same mistake as the blanket quoted-heredoc exemption that was
    removed in — it would let a literal split-brain sit inside embedded code.
    Falls back to the token scan only when the file does not parse.
    """
    covered: set = set()
    try:
        tree = ast.parse(src)
    except (SyntaxError, ValueError):
        try:
            for tok in tokenize.generate_tokens(io.StringIO(src).readline):
                if tok.type == tokenize.STRING and (
                    tok.string.count('"""') or tok.string.count("'''")
                ):
                    covered.update(range(tok.start[0], tok.end[0] + 1))
        except (tokenize.TokenError, IndentationError, SyntaxError):
            pass
        return covered
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef)):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        first = body[0]
        if (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)):
            covered.update(range(first.lineno, (first.end_lineno or first.lineno) + 1))
    return covered


def strip_sanctioned(line: str) -> str:
    """Remove each sanctioned idiom AND its default literal from `line`.

    The check then re-scans what is left, so a SECOND raw root literal sitting beside a
    sanctioned one on the same line is still caught.  A blanket `continue` on the whole
    line waved those through — and the tier-B preamble line is exactly where a hurried
    edit appends another path.
    """
    out, last = [], 0
    for m in SANCTIONED_RE.finditer(line):
        if m.start() >= last:
            out.append(line[last:m.start()])
            last = m.end()
        # Swallow the idiom's OWN default literal — and only that.  It must start
        # IMMEDIATELY after the idiom (modulo whitespace); anything else on the line is
        # somebody else's path and must survive to be re-scanned.  An unanchored
        # "skip to the next quoted string" swallowed a second, unrelated literal and
        # made the gate green on exactly the shape it was written to catch.
        # The BASH form's default is already inside the match (`:="…"}"`), so a swallow
        # after it can only consume somebody else's path — skip it there.
        if not m.group(0).endswith('}"'):
            k = re.match(r'\s*(["\'][^"\']*["\'])', line[last:])
            if k:
                last += k.end()
    out.append(line[last:])
    return "".join(out)


def scan(rel: str, src: str):
    """Yield (lineno, text) for every offending line in one file."""
    if os.path.basename(rel) in HELPERS or FROZEN_RE.search(rel):
        return
    exempt_subs = NARROW_EXEMPTIONS.get(rel, ())
    doc = py_docstring_lines(src) if rel.endswith(".py") else set()
    for i, line in enumerate(src.splitlines(), 1):
        if not ROOT_RE.search(line):
            continue
        s = line.strip()
        if s.startswith("#SBATCH") or s.startswith("#"):
            continue
        if i in doc:
            continue
        if any(sub in line for sub in exempt_subs):
            continue
        if not ROOT_RE.search(strip_sanctioned(line)):
            continue
        yield i, s


def scan_sbatch_order(src: str):
    """Flag a ``#SBATCH`` directive that sbatch will IGNORE because an executable
    line precedes it — sbatch stops parsing directives at the first non-comment,
    non-blank line.  Guards the tier-B portable-root idiom's PLACEMENT (the
    2026-08-26 incident: ``: "${DUCK_SCRATCH_ROOT:=…}"`` inserted at line 2
    silently killed every directive in 4 runners).  Scanning stops at the first
    heredoc opener (``<<``) so job scripts embedded in a writer script are not
    false-flagged; a conservative early stop on a literal ``<<`` is acceptable —
    directive blocks belong at the top of a file, above any heredoc.
    Returns (sbatch_lineno, first_exec_lineno) or None."""
    first_exec = None
    for i, line in enumerate(src.splitlines(), 1):
        if "<<" in line:
            return None
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if s.startswith("#SBATCH") and first_exec is not None:
                return (i, first_exec)
            continue
        if first_exec is None:
            first_exec = i
    return None


def target_files(argv):
    """Return (mode, files).  mode is 'index' when reading staged content."""
    args = [a for a in argv if not a.startswith("--")]
    if args:
        return "worktree", args
    root = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True, check=True
    ).stdout.strip()
    os.chdir(root)
    if "--all" in argv:
        out = subprocess.run(["git", "ls-files"], capture_output=True, text=True,
                             check=True).stdout
        mode = "worktree"
    else:
        out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR",
                              "-z"], capture_output=True, text=True, check=True).stdout
        out = out.replace("\0", "\n")
        mode = "index"
    return mode, [f for f in out.splitlines() if f.endswith(EXTS)]


def read_source(rel: str, mode: str):
    """Content of `rel`, from the INDEX in staged mode and from disk otherwise.

    Staged mode must read the index: `git diff --cached` selects what is *about to be
    committed*, so checking the working tree would both miss a violation that is staged
    but since edited away on disk, and reject a violation that only exists unstaged.
    """
    if mode == "index":
        p = subprocess.run(["git", "show", f":{rel}"], capture_output=True)
        if p.returncode != 0:
            return None
        try:
            return p.stdout.decode("utf-8")
        except UnicodeDecodeError:
            return None
    if not os.path.isfile(rel):
        # Tracked but absent from the worktree — a normal mid-edit state (deleted, not
        # yet `git rm`'d).  NOT the same as undecodable, and reporting it as
        # "unreadable … hardcoded single-user root(s)" failed `make check` with a false
        # banner and a false reason.  Sentinel, so main() can tell the two apart.
        return MISSING
    try:
        with open(rel, encoding="utf-8") as fh:
            return fh.read()
    except (UnicodeDecodeError, OSError):
        return None


def main(argv) -> int:
    mode, files = target_files(argv)
    files = [f for f in files if f.endswith(EXTS)]
    bad = []
    checked = 0
    for rel in files:
        src = read_source(rel, mode)
        if src is MISSING:
            continue
        if src is None:
            # NOT a silent skip.  A code file that cannot be decoded cannot be checked,
            # and one stray non-UTF-8 byte would otherwise drop the whole file out of
            # the gate while the run still printed OK and exited 0.  (The migration
            # harness already fails on this; the STANDING gate must too.)
            bad.append((rel, 0, "unreadable or not UTF-8 — cannot be checked"))
            continue
        checked += 1
        for i, text in scan(rel, src):
            bad.append((rel, i, text))
        if rel.endswith(SHELL_EXTS) and not FROZEN_RE.search(rel):
            hit = scan_sbatch_order(src)
            if hit:
                bad.append((rel, hit[0],
                            f"#SBATCH directive AFTER executable line {hit[1]} — sbatch "
                            "ignores it; move the executable line below the directive block"))
    if not bad:
        print(f"verify_portable_paths: OK ({checked} file(s) checked, source={mode})")
        return 0
    print("verify_portable_paths: FAILED — hardcoded single-user root(s) "
          "and/or dead #SBATCH directives:\n")
    for rel, i, text in bad:
        print(f"  {rel}:{i}: {text[:160]}")
    print(
        "\nUse a portable root instead (docs/audit_hardcoded_paths.md):\n"
        "  python (under pipeline/):  from _paths import SCRATCH_ROOT\n"
        "  python (elsewhere):        SCRATCH_ROOT = (os.environ.get(\"DUCK_SCRATCH_ROOT\")\n"
        "                                              or \"<the old literal>\")\n"
        "  bash:                      : \"${DUCK_SCRATCH_ROOT:=<the old literal>}\"\n"
        "  R:                         .duck_env <- function(n, d) {\n"
        "                                 v <- Sys.getenv(n); if (nzchar(v)) v else d }\n"
        "                             SCRATCH_ROOT <- .duck_env(\n"
        "                                 \"DUCK_SCRATCH_ROOT\", \"<the old literal>\")\n"
        "  local Mac / repo root:     derive it — no literal, no env var:\n"
        "                             _r = Path(__file__).resolve(); walk up to the\n"
        "                             checkout (see docs/audit_hardcoded_paths.md)\n"
        "\nBypass (only for a deliberate exception):  SKIP_PORTABLE_PATHS_CHECK=1 git commit ...\n"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
