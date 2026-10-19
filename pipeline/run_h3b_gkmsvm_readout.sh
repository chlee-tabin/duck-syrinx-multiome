#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3b_kmer_readout_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3b_kmer_readout_%j.err
#SBATCH -J h3b_kmer_readout

# NATIVE gkm-SVM k-mer readout (Lee deltaSVM-style): gkmpredict every
# non-redundant 10-mer against the trained model -> rank weights -> cluster
# top k-mers into motifs -> scan for TF cores.  This is the real "what did the
# SVM learn" readout, replacing the contradicted HOMER "SMAD5" de novo call.
#
# Runs BOTH models from the HH17/16 strict diagnostic:
#   gc      — GC-matched negatives (PRIMARY: top k-mers are GC-independent)
#   random  — legacy negatives (compare: which top k-mers were GC-driven)
# The k-mer FASTA is built once and reused (SKIP_PREDICT only reuses scores;
# the enumeration is cheap either way).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate gkmqc

echo "gkmpredict: $(which gkmpredict)"

for NM in gc random; do
  echo "=== h3b native readout: NEG_MATCH=${NM} ==="
  MODE=strict NEG_MATCH="${NM}" python3 pipeline/h3b_gkmsvm_kmer_readout.py
done

date
echo "=== DONE ==="
