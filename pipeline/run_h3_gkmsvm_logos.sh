#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3_gkmsvm_logos_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3_gkmsvm_logos_%j.err
#SBATCH -J h3_gkmsvm_logos

# De novo motif discovery on h2 strict-mode L-biased peaks via HOMER
# findMotifs.pl (Path A in the spec).  Closes the loop opened in h2's
# docstring.  Runs on 99 positives × 990 negatives — small input, ~10-30 min
# wall.
#
# Pre-defined motif scans (h3_jaspar_scan.R, g3b_*) already covered the
# "do PITX/FOX/SMAD2-3 enrich?" question.  This script asks the
# COMPLEMENTARY "what motifs emerge UNSUPERVISED — including motifs we did
# not ask about?".
#
# MEME suite is NOT available on O2; HOMER 5.1 is.  HOMER's compareMotifs
# runs internally vs the bundled vertebrate motif library (TomTom analog).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3 homer/5.1
conda activate gkmqc  # for python3 + pandas + numpy

echo "findMotifs.pl: $(which findMotifs.pl)"
echo "python3: $(which python3)"

# NEG_MATCH selects the h2 negative FASTA variant (random = legacy bilateral
# draw the committed h3 known-scan was produced from; gc = GC-matched negatives).
MODE=${MODE:-strict} NEG_MATCH=${NEG_MATCH:-random} python3 pipeline/h3_gkmsvm_logos.py

date
echo "=== DONE ==="
