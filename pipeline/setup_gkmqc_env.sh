#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/setup_gkmqc_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/setup_gkmqc_%j.err
#SBATCH -J setup_gkmqc

# Create conda environment for gkmQC (Han et al. 2022 PNAS).  Upstream
# install instructions: https://github.com/Dongwon-Lee/gkmQC
# Pipeline consumer: pipeline/f2b_gkmqc_rescue.py (plan §3.4).
#
# Workflow: clone repo → create env from environment.yml → build C
# source → install `gkmqc.py` onto the env's PATH.  Idempotent: a
# second invocation reuses the existing env and skips the clone if the
# repo is already on scratch.
#
# Run: sbatch pipeline/setup_gkmqc_env.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date
echo "=== Setting up gkmqc conda environment ==="

SCRATCH="${DUCK_SCRATCH_ROOT}"
GKMQC_DIR="$SCRATCH/references/_gkmqc_src"

module purge
module load gcc/14.2.0 conda/miniforge3

if [ ! -d "$GKMQC_DIR/.git" ]; then
    echo "Cloning gkmQC..."
    git clone https://github.com/Dongwon-Lee/gkmQC.git "$GKMQC_DIR"
else
    echo "gkmQC source already present at $GKMQC_DIR — pulling latest"
    git -C "$GKMQC_DIR" pull --ff-only   # a stale checkout must not pass silently
fi

cd "$GKMQC_DIR"

if ! conda env list 2>/dev/null | grep -q '^gkmqc '; then
    echo "Creating gkmqc env from environment.yml..."
    conda env create -f environment.yml -n gkmqc
else
    echo "gkmqc env already exists, skipping create"
fi

conda activate gkmqc

echo ""
echo "=== Building C source ==="
cd "$GKMQC_DIR/src"
make clean || true
make
make install

echo ""
echo "=== Exposing gkmQC binaries on PATH ==="
# The upstream Makefile copies into $GKMQC_DIR/bin/, not the conda env.
# Symlink every executable into $CONDA_PREFIX/bin/ so `conda activate
# gkmqc` is sufficient to run gkmqc.py / gkmkern from anywhere.
for tool in "$GKMQC_DIR"/bin/*; do
    [ -f "$tool" ] || continue
    link="$CONDA_PREFIX/bin/$(basename "$tool")"
    ln -sf "$tool" "$link"
    echo "  $link -> $tool"
done

echo ""
echo "=== Patching pyfasta for Python 3 ==="
# environment.yml pins python=3.7 + pip-installs pyfasta, which is
# abandoned Python-2 code (uses implicit relative imports + cPickle +
# .iteritems).  Job 37810426 crashed at "No module named 'fasta'"
# then 'records' then 'cPickle'.  Run 2to3 over the installed package
# once — fixes all the legacy idioms in one pass.  Idempotent: 2to3
# on already-refactored code is a no-op.
PYFASTA_PKG="$CONDA_PREFIX/lib/python3.7/site-packages/pyfasta"
if [ -d "$PYFASTA_PKG" ]; then
    python -m lib2to3 -wn "$PYFASTA_PKG" 2>&1 | tail -10
else
    echo "WARN: pyfasta not found at $PYFASTA_PKG"
fi

echo ""
echo "=== Cloning + building lsgkm (Lee 2016) into gkmqc env ==="
# lsgkm provides gkmtrain + gkmpredict — the SVM training/scoring binaries
# behind gkm-SVM motif discovery (item c of overnight chromatin queue
# 2026-05-12).  gkmQC bundles only gkmkern; lsgkm is the upstream source.
# Idempotent: skip clone + make if binaries already present.
LSGKM_DIR="$SCRATCH/references/_lsgkm_src"
if [ ! -d "$LSGKM_DIR/.git" ]; then
    echo "Cloning lsgkm..."
    git clone https://github.com/Dongwon-Lee/lsgkm.git "$LSGKM_DIR"
else
    echo "lsgkm source already present — pulling latest"
    git -C "$LSGKM_DIR" pull --ff-only   # a stale checkout must not pass silently
fi
if [ ! -x "$LSGKM_DIR/bin/gkmtrain" ] || [ ! -x "$LSGKM_DIR/bin/gkmpredict" ]; then
    (cd "$LSGKM_DIR/src" && make && make install)
fi
for tool in "$LSGKM_DIR"/bin/*; do
    [ -f "$tool" ] || continue
    link="$CONDA_PREFIX/bin/$(basename "$tool")"
    ln -sf "$tool" "$link"
    echo "  $link -> $tool"
done

echo ""
echo "=== Patching preprocess.py chromosome regex for NCBI accessions ==="
# gkmQC was authored for UCSC-style names; gkmQC/scripts/preprocess.py:96
# filters peaks with `$1 ~ /^chr[0-9XY]+$/`, which rejects every duck
# chromosome (NC_051772.1, …).  Job 39834221 hit this: 11,013 MACS2
# peaks → empty .e300.bed → empty .qc.bed → IndexError in split_posset.
# gkmQC's own zju1 bit-data is already keyed by NC accessions, so
# widening the regex is consistent, not a hack.  Idempotent.
PREPROCESS_PY="$GKMQC_DIR/scripts/preprocess.py"
OLD_REGEX='/^chr[0-9XY]+$/'
NEW_REGEX='/^(chr[0-9XYM]+|NC_[0-9]+[.][0-9]+)$/'
if [ -f "$PREPROCESS_PY" ]; then
    if grep -qF "$OLD_REGEX" "$PREPROCESS_PY"; then
        python -c "
src = '$PREPROCESS_PY'
c = open(src).read()
open(src, 'w').write(c.replace('$OLD_REGEX', '$NEW_REGEX'))
"
        echo "  patched: $PREPROCESS_PY"
    elif grep -qF "$NEW_REGEX" "$PREPROCESS_PY"; then
        echo "  already patched: $PREPROCESS_PY"
    else
        echo "WARN: neither old nor new regex found in $PREPROCESS_PY"
    fi
else
    echo "WARN: preprocess.py not found at $PREPROCESS_PY"
fi

echo ""
echo "=== Registering gkmQC/scripts on PYTHONPATH via .pth ==="
# gkmqc.py imports helper modules from scripts/ (seqs_nullgen,
# gkmsvm, optimize, preprocess, visualize) — none of them are packaged
# as installable modules, so drop a .pth file into the env's
# site-packages to prepend scripts/ to sys.path whenever the env is
# active.  Job 37809463 crashed on ModuleNotFoundError: seqs_nullgen.
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
PTH="$SITE_PACKAGES/gkmqc_scripts.pth"
echo "$GKMQC_DIR/scripts" > "$PTH"
echo "  $PTH -> $GKMQC_DIR/scripts"

echo ""
echo "=== Verification ==="
command -v gkmqc.py >/dev/null && echo "gkmqc.py: $(which gkmqc.py)" || {
    echo "ERROR: gkmqc.py not on PATH after install"
    exit 1
}
gkmqc.py --help 2>&1 | head -5 || true

date
echo "=== Done ==="
echo ""
echo "Next: sbatch pipeline/build_duck_gkmqc_index.sh    # one-shot index"
