## Centralized, portable path roots for duck-syrinx-multiome R scripts.
##
## Every R script here is launched by a Bash runner (`Rscript ...` from a sbatch
## wrapper).  That runner cannot reliably source scripts/_paths.sh (see the note in
## that file), so this helper does NOT depend on it: it reads the same DUCK_* env with
## the same absolute defaults.  If the runner *did* export DUCK_*, those values win.
##
## Why + full audit/rollout: docs/audit_hardcoded_paths.md.  Defaults equal the
## historically hardcoded paths, so existing runs are unchanged; override via env to
## run as another user / on relocated storage.
##
## Source it with the depth-agnostic walk-up (works from any subdirectory, in the repo
## checkout or the scratch mirror):
##
##     .dp <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
##     .dp <- if (length(.dp)) dirname(normalizePath(.dp[1])) else normalizePath(getwd())
##     while (!file.exists(file.path(.dp, "scripts", "_paths.R")) && dirname(.dp) != .dp)
##       .dp <- dirname(.dp)
##     source(file.path(.dp, "scripts", "_paths.R"))
##
## Scripts that are ferried to a standalone scratch tree (no scripts/ alongside) use the
## same env-with-default idiom INLINE instead — see "Two adoption tiers" in the audit doc.

.duck_env <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) v else default
}

## SCRATCH_USER is the *parent* of SCRATCH_ROOT; the sibling chick_demux / old_duck_*
## scratch trees hang off it, so it is a root in its own right.
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(SCRATCH_USER, "FC_08656"))
## The standby *raw-data* tree is one level above the standby git checkout.
STANDBY_DATA_ROOT <- .duck_env(
  "DUCK_STANDBY_DATA_ROOT",
  "/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data"
)
STANDBY_REPO <- .duck_env(
  "DUCK_STANDBY_REPO", file.path(STANDBY_DATA_ROOT, "R4.4.0", "duck-syrinx-multiome")
)
## Category C — the stable `cdw` reference home.  Parameterized for portability but it
## stays ABSOLUTE: deliberately NOT scratch-relative (docs/audit_hardcoded_paths.md).
REFERENCE_ROOT <- .duck_env(
  "DUCK_REFERENCE_ROOT", "/n/groups/tabin/chlee/singlecell/reference"
)
## Cluster user home (conda installs, R user lib).  NOT path.expand("~"): call sites may
## build a *remote* path from a local machine, where ~ is the wrong home.
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")
## NOTE: no REPO_ROOT here on purpose.  Defaulting it to "" would make file.path()
## silently yield "/figures/..." whenever DUCK_REPO_ROOT is unset — a silent wrong
## path, not an error.  The documented walk-up above already binds the checkout root
## as `.dp`; use that (or Sys.getenv("DUCK_REPO_ROOT") with your own explicit check).

if (!interactive() && identical(Sys.getenv("DUCK_PATHS_SELFTEST"), "1")) {
  cat("SCRATCH_USER      =", SCRATCH_USER, "\n")
  cat("SCRATCH_ROOT      =", SCRATCH_ROOT, "\n")
  cat("STANDBY_DATA_ROOT =", STANDBY_DATA_ROOT, "\n")
  cat("STANDBY_REPO      =", STANDBY_REPO, "\n")
  cat("REFERENCE_ROOT    =", REFERENCE_ROOT, "\n")
  cat("HOME_ROOT         =", HOME_ROOT, "\n")
}
