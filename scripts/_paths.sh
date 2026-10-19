# shellcheck shell=bash
# Centralized, portable path roots for duck-syrinx-multiome shell scripts.
#
# Source near the top of a script:
#     . "$(dirname "${BASH_SOURCE[0]}")/_paths.sh"   # when the script lives in scripts/
#     . "$DUCK_REPO_ROOT/scripts/_paths.sh"          # from elsewhere, once REPO root is known
#
# Why + full audit/rollout: docs/audit_hardcoded_paths.md.  Defaults equal the
# historically hardcoded paths, so existing runs are unchanged; override via env to run
# as another user / on relocated storage.  `: "${VAR:=default}"` means a caller's
# pre-set env value wins (e.g. the cluster submission wrapper's own overrides).
#
# IMPORTANT — sourcing this file is NOT possible from most cluster-side scripts:
#   * a SLURM batch script executes from the node spool copy, so `dirname
#     "${BASH_SOURCE[0]}"` is /var/spool/slurmd/..., not the repo; and
#   * scripts/ is not part of the scratch rsync (the repo-to-scratch sync ships
#     FC_08656/*, references/, pipeline/, figures/{_shared,data} only).
# Those scripts therefore inline the SAME `: "${DUCK_X:=default}"` lines instead.
# Because the idiom is identical, DUCK_* present in the JOB's environment makes every
# inline default a no-op.  Getting it there is NOT automatic, though:
#   * sbatch propagates its own environment (--export=ALL), but the cluster submission wrapper
#     invokes sbatch on o2 THROUGH THE RELAY — a fresh login shell that does not
#     inherit this machine's env.  It therefore embeds the roots in the remote
#     command (DUCK_SCRATCH_ROOT=... sbatch ...); any other submission path must too.
#   * `#SBATCH -o/-e` directives are parsed before any shell runs, so SLURM log paths
#     stay pinned to the literal and need `sbatch --output=... --error=...` instead.
# See docs/audit_hardcoded_paths.md, "Two adoption tiers".
#
# Exports: DUCK_REPO_ROOT DUCK_SCRATCH_USER DUCK_SCRATCH_ROOT DUCK_STANDBY_DATA_ROOT
#          DUCK_STANDBY_REPO DUCK_REFERENCE_ROOT DUCK_HOME DUCK_LOCAL_DATA_ROOT
#
# NAMESPACE NOTE: `DUCK_*` is ALSO used in this repo for per-script duck-vs-chick data
# variables (DUCK_FA, DUCK_SIZES, DUCK_GTF, DUCK_CHR, …).  None of them collide with the
# root names above, and the cluster submission wrapper propagates only the roots it names explicitly,
# so a script-local DUCK_* never leaks into a job.  Check for a collision before adding a
# root:  git grep -hoE 'DUCK_[A-Z_]+' -- '*.py' '*.sh' '*.sbatch' '*.R' | sort -u

# Repo root derived from this file's location (scripts/_paths.sh -> repo/).
DUCK_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# DUCK_SCRATCH_USER is the *parent* of DUCK_SCRATCH_ROOT; the sibling chick_demux /
# old_duck_* scratch trees hang off it, so it is a root in its own right.
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER}/FC_08656}"
# The standby *raw-data* tree is one level above the standby git checkout (backups,
# cellranger-arc _v210 mirrors).  Off-scratch by definition: parameterize, never relativize.
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_STANDBY_REPO:=${DUCK_STANDBY_DATA_ROOT}/R4.4.0/duck-syrinx-multiome}"
# Category C — the stable `cdw` reference home.  Parameterized for portability but it
# stays ABSOLUTE: deliberately NOT scratch-relative (docs/audit_hardcoded_paths.md).
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
# Cluster user home (conda installs, R user lib).  NOT "$HOME": many call sites build a
# *remote* path from a local machine, where $HOME is the wrong home.
: "${DUCK_HOME:=/home/cl266}"
# The LOCAL (this-machine) data tree the laptop<->cluster sync/pull helpers mirror
# cluster data into.  No default: it is personal to each machine, so those helpers
# require it to be set (`: "${DUCK_LOCAL_DATA_ROOT:?…}"`).  Nothing in the pipeline or
# the figure code reads it.  Distinct from DUCK_REPO_ROOT (the git checkout).
DUCK_LOCAL_DATA_ROOT="${DUCK_LOCAL_DATA_ROOT:-}"

export DUCK_REPO_ROOT DUCK_SCRATCH_USER DUCK_SCRATCH_ROOT DUCK_STANDBY_DATA_ROOT
export DUCK_STANDBY_REPO DUCK_REFERENCE_ROOT DUCK_HOME DUCK_LOCAL_DATA_ROOT
