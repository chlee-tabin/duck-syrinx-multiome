#!/bin/bash
# Fig 5A / Fig 5C input reproduction: build the scratch roots and stage the
# cellranger-arc v2.1.0 fragment files from the raw-data tree.  Runs on a cluster transfer
# node (file / git ops only; no python, no modules).
#
# Layout (everything under $DUCK_REPRO_ROOT; nothing is written to the production tree):
#   root/         shallow clone of $REPRO_REF (default main) -- the FRESH chain
#                 f1 (from fragments) -> f5 default path -> f4 (shipped rule)
#     <LIB>       -> ../arc_v210/<LIB>         (staged v2.1.0 outs, f1's ARC_LIBRARIES paths)
#     phase2/step04_preCB_clustered.h5ad -> production RNA master (read-only input)
#   root_cache/   the in-repo --from-cached-matrix path: code trees -> root/, and
#                 phase2/atac/{step_f1_coupled_atac_metadata.tsv, step_f5_pitx2_region_matrix_tmp.h5ad}
#                 -> archived copies of the August 2026 run's inputs (not deposited in GEO, so
#                 the cache mode and checks A1-A3 need them; the fresh chain does not).
#                 Kept apart from root/ so the fresh f5 run cannot "reuse" the old cache.
#   arc_v210/<LIB>/outs/atac_fragments.tsv.gz(.tbi) + _staged_<LIB> markers + _staged_done
#   out/{f5_cache,f5_fresh,compare}/ , logs/
set -euo pipefail
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_REPRO_ROOT:=${DUCK_SCRATCH_ROOT}/repro_20261016}"
S="$DUCK_SCRATCH_ROOT"; B="$DUCK_STANDBY_DATA_ROOT"; X="$DUCK_REPRO_ROOT"
BRANCH="${REPRO_REF:-main}"
# clone source: this checkout's origin unless REPRO_REMOTE is set
REMOTE="${REPRO_REMOTE:-$(git -C "$(dirname "${BASH_SOURCE[0]}")" remote get-url origin)}"
PRIV="$S/_arc0sL_inputs/phase2/atac"
LIBS=(ARC0S-duck PILOTARC4A PILOTARC3 PILOTARC4B ARC05 ARC06 ARC07 ARC08 ARC03-1A ARC03-1B ARC03-2 ARC04)

echo "=== setup_and_stage $(date) on $(hostname) ==="
mkdir -p "$X"/arc_v210 "$X"/out/f5_cache "$X"/out/f5_fresh "$X"/out/compare "$X"/logs

# ---- root/: clone (or fast-forward an existing clone of the same branch) ----
if [[ -d "$X/root/.git" ]]; then
    git -C "$X/root" pull --ff-only origin "$BRANCH"
else
    [[ -e "$X/root" ]] && { echo "ERROR: $X/root exists and is not a clone" >&2; exit 1; }
    git clone --depth 1 -b "$BRANCH" "$REMOTE" "$X/root"
fi
echo "root HEAD $(git -C "$X/root" rev-parse HEAD)"
mkdir -p "$X/root/phase2/atac"
ln -sfn "$S/phase2/step04_preCB_clustered.h5ad" "$X/root/phase2/step04_preCB_clustered.h5ad"
for L in "${LIBS[@]}"; do ln -sfn "$X/arc_v210/$L" "$X/root/$L"; done
# The fresh root must NOT hold an f5 region-matrix cache (f5 would reuse it).
[[ ! -e "$X/root/phase2/atac/step_f5_pitx2_region_matrix_tmp.h5ad" ]] || {
    echo "ERROR: fresh root already holds an f5 cache" >&2; exit 1; }

# ---- root_cache/ ----
mkdir -p "$X/root_cache/phase2/atac"
for d in pipeline figures references design; do ln -sfn "$X/root/$d" "$X/root_cache/$d"; done
for f in step_f1_coupled_atac_metadata.tsv step_f5_pitx2_region_matrix_tmp.h5ad; do
    [[ -f "$PRIV/$f" ]] || { echo "ERROR: missing private input $PRIV/$f" >&2; exit 1; }
    ln -sfn "$PRIV/$f" "$X/root_cache/phase2/atac/$f"
done

echo "roots ready $(date)"

# ---- stage v2.1.0 fragments (additive; size-verified per file) ----
src_of() {
    case "$1" in
        PILOTARC3)             echo "$B/FC_08429/PILOTARC3_v210_PITX2/outs" ;;
        PILOTARC4A|PILOTARC4B) echo "$B/FC_08470/${1}_v210/outs" ;;
        *)                     echo "$B/FC_08656/${1}_v210/outs" ;;
    esac
}
: > "$X/arc_v210/_sizes.tsv"
for L in "${LIBS[@]}"; do
    src=$(src_of "$L"); dst="$X/arc_v210/$L/outs"
    mkdir -p "$dst"
    echo "[$(date +%T)] rsync $L from $src"
    rsync -a --partial "$src/atac_fragments.tsv.gz" "$src/atac_fragments.tsv.gz.tbi" "$dst/"
    for f in atac_fragments.tsv.gz atac_fragments.tsv.gz.tbi; do
        a=$(stat -c %s "$src/$f"); b=$(stat -c %s "$dst/$f")
        [[ "$a" == "$b" ]] || { echo "ERROR: size mismatch $L/$f $a vs $b" >&2; exit 1; }
        printf '%s\t%s\t%s\t%s\n' "$L" "$f" "$b" "$src/$f" >> "$X/arc_v210/_sizes.tsv"
    done
    touch "$X/arc_v210/_staged_$L"
    echo "[$(date +%T)] staged $L"
done
date > "$X/arc_v210/_staged_done"
echo "=== STAGE_OK $(date) ==="
