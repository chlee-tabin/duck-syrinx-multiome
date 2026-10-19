# Portable, env-overridable cluster path roots

**How paths work in this repo.** Every cluster storage root is an environment variable
whose default is the path the pipeline originally ran under (the `cl266` tree on the HMS
O2 cluster). A run with nothing set behaves exactly as the original run did; one variable
re-roots a whole tree. The mechanism lives in `pipeline/_paths.py`, `scripts/_paths.sh`
and `scripts/_paths.R`; the standing gate is `scripts/verify_portable_paths.py`
(`make check-portable-paths`, also run by the pre-commit hook on staged content after
`make install-hooks`).

**If you are adding a script, read "Two adoption tiers" — that is the only rule you need.**

## Roots

| Variable | Default | What lives there |
|---|---|---|
| `DUCK_SCRATCH_USER` | `/n/scratch/users/c/cl266` | parent work tree; the chick (`chick_demux*`) and E6/E9/E11 (`old_duck_cr9`) work trees hang off it |
| `DUCK_SCRATCH_ROOT` | `$DUCK_SCRATCH_USER/FC_08656` | multiome pipeline outputs (h5ads, intermediates, `phase2/`) |
| `DUCK_REFERENCE_ROOT` | `/n/groups/tabin/chlee/singlecell/reference` | genome references (FASTA, `.fai`, GTF) |
| `DUCK_HOME` | `/home/cl266` | conda installs and the R user library |
| `DUCK_STANDBY_DATA_ROOT` | `/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data` | backup / raw-data tree |
| `DUCK_STANDBY_REPO` | `$DUCK_STANDBY_DATA_ROOT/R4.4.0/duck-syrinx-multiome` | a second checkout of this repository next to the data |

The repository checkout root is never a literal and never a variable: code walks up from
its own file to the directory that contains `pipeline/_paths.py` and raises if there is
none.

## Two adoption tiers

Importing or sourcing the helper is the better idiom only where the helper is reachable at
run time. Three run-time facts make a blanket "everyone sources the helper" rule unsafe:

1. A SLURM batch script runs from the node's spool copy, so `dirname "${BASH_SOURCE[0]}"`
   inside a `#SBATCH` script is not the repository.
2. The scratch working tree carries `pipeline/`, `references/`, `figures/{_shared,data}`
   and the top level of `FC_08656/`, but not `scripts/`, so `scripts/_paths.sh` and
   `scripts/_paths.R` are absent at run time.
3. The chick (`pipeline/chick/`) and E6/E9/E11 steps run from their own work trees under
   `$DUCK_SCRATCH_USER`, which carry no `pipeline/_paths.py`.

| Tier | Where | Idiom |
|---|---|---|
| **A — helper** | Python under `pipeline/*.py` (not `pipeline/chick/`), `FC_08656/*.py`, `figures/_shared/`, `scripts/*.py`; shell at the top level of `scripts/` without `#SBATCH` | `from _paths import SCRATCH_ROOT` (walk-up bootstrap in subdirectories) / `. "$(dirname "${BASH_SOURCE[0]}")/_paths.sh"` |
| **B — inline** | every `#SBATCH` script, `pipeline/*.sh`, all of `pipeline/chick/**`, `FC_08656/*.sh`, shell in `scripts/` subdirectories, every `.R` | the same env-with-default line the helper uses, in the empty-normalizing form below |

Tier-B forms (the gate accepts exactly these):

```bash
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
```
```python
SCRATCH = (os.environ.get("DUCK_SCRATCH_ROOT")
           or (os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266") + "/FC_08656")
```
```r
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_ROOT <- .duck_env("DUCK_SCRATCH_ROOT", file.path(.duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266"), "FC_08656"))
```

Both tiers use the same variable names and the same "environment wins, else default"
rule, so once `DUCK_*` is set in a job's environment every tier-B default is inert.

**Contract: a `DUCK_*` root is either unset or a non-empty absolute path.** All three
languages treat a set-but-empty value as unset. The two-argument
`os.environ.get(name, default)` / `Sys.getenv(name, default)` forms are rejected by the
gate because they return `""` for an empty variable, which turns
`f"{SCRATCH_ROOT}/phase2"` into `"/phase2"`. The Python form is parenthesised on purpose:
`or` binds looser than `+`.

## Running as another user

```bash
export DUCK_SCRATCH_USER=/path/to/work        # re-roots every scratch code path
export DUCK_REFERENCE_ROOT=/path/to/references
export DUCK_HOME=$HOME
sbatch --export=ALL \
    --output=$DUCK_SCRATCH_USER/FC_08656/logs/%x_%j.out \
    --error=$DUCK_SCRATCH_USER/FC_08656/logs/%x_%j.err \
    pipeline/run_pipeline.sh --from d2 --to d2
```

- `sbatch` passes the submitting shell's environment to the job (`--export=ALL`); a job
  submitted from a shell without `DUCK_*` falls back to the defaults.
- `#SBATCH -o/-e` lines are read by SLURM before any shell runs, so they stay literals.
  Pass `--output=` and `--error=` on the command line for every direct `sbatch`.
- The cellranger-arc `--libraries` CSVs (`pipeline/nlibraries_v210/*.csv`,
  `pipeline/chick/nlibraries0S-chick.csv`, `FC_08656/nlibraries*.csv`) are read by
  cellranger itself, so their FASTQ paths are not re-rooted by a variable; edit or
  regenerate them for a new location.

## What still carries a literal

- Tier-B inline defaults (one per root per file): this is the portable idiom.
- `#SBATCH -o/-e` directives (see above).
- `pipeline/_paths.py`, `scripts/_paths.sh`, `scripts/_paths.R`: the single source of the
  defaults.
- Docstrings and comments that show usage.
- Records of past runs (provenance columns in committed tables); no script reads a path
  back out of them.

## How to verify

```bash
python3 pipeline/_paths.py
bash -c '. scripts/_paths.sh && echo "$DUCK_SCRATCH_ROOT | $DUCK_REPO_ROOT"'
DUCK_PATHS_SELFTEST=1 Rscript scripts/_paths.R
DUCK_SCRATCH_USER=/tmp/x python3 pipeline/_paths.py     # -> SCRATCH_ROOT = /tmp/x/FC_08656
make check-portable-paths
```
