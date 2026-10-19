#!/bin/bash
# Stage-F conda env build — run on an INTERACTIVE-partition compute node via srun
# (o2 login is sbatch/srun-only; installs go to interactive compute, which has internet).
# Invoke (from o2, inside a detached tmux so it survives SSH drops):
#   tmux new-session -d -s stagef_envs \
#     'srun -p interactive -c 4 --mem 16G -t 0-04:00:00 bash -l \
#        /n/scratch/users/c/cl266/FC_08656/stagef/scripts/stagef_build_envs.sh \
#        > /n/scratch/users/c/cl266/FC_08656/stagef/logs/build_envs.log 2>&1'
# Two envs (independent solves):
#   stagef   = phast (phyloP/phastCons/phyloFit/msa_view/tree_doctor/maf_parse) + multiz/roast + lastz
#   phyloacc = phyloacc (PhyloAcc-ST/-GT + phyloacc.py)
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_ROOT DUCK_SCRATCH_USER

mkdir -p ${DUCK_SCRATCH_ROOT}/stagef/logs
module load gcc/14.2.0 conda/miniforge3
export CONDA_PKGS_DIRS=${DUCK_SCRATCH_USER}/conda_pkgs_stagef
mkdir -p "$CONDA_PKGS_DIRS"
echo "=== host=$(hostname) conda=$(conda --version 2>&1) START $(date) ==="
echo "=== net reachability (compute node) ==="
curl -fsS -m 15 -o /dev/null -w "bioconda repodata HTTP %{http_code}\n" https://conda.anaconda.org/bioconda/noarch/repodata.json 2>&1 | head -1 || echo "curl to bioconda FAILED"

echo "=== [1/2] create env: stagef (phast + multiz + lastz) $(date) ==="
conda create -y -n stagef -c bioconda -c conda-forge phast multiz lastz
echo "rc_stagef_create=$?"
echo "--- stagef tool check ---"
conda run -n stagef bash -lc 'for t in phyloP phastCons phyloFit msa_view tree_doctor maf_parse multiz roast maf_project lastz; do printf "%-12s %s\n" "$t" "$(command -v $t || echo MISSING)"; done'

echo "=== [2/2] create env: phyloacc $(date) ==="
conda create -y -n phyloacc -c bioconda -c conda-forge phyloacc
echo "rc_phyloacc_create=$?"
echo "--- phyloacc tool check ---"
conda run -n phyloacc bash -lc 'for t in phyloacc.py PhyloAcc-ST PhyloAcc-GT; do printf "%-14s %s\n" "$t" "$(command -v $t || echo MISSING)"; done; phyloacc.py --version 2>&1 | head -3'

echo "=== ALL ENV BUILDS DONE $(date) __STAGEF_ENVS_OK__ ==="
