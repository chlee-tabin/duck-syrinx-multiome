#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=192G
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_peaks_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_peaks_%j.err
#SBATCH -J archr_peaks
#
# ArchR track A peak calling (plan §3.2, §3.3).  Parallel to the
# SnapATAC2 f2 track; writes to phase2/atac_archr/ (separate from
# phase2/atac/peaks/ which belongs to f2).
#
# Usage:
#   sbatch pipeline/run_archr_peaks.sh                          # full 12 libs
#   sbatch pipeline/run_archr_peaks.sh --libraries ARC07         # pilot 1 lib
#   sbatch pipeline/run_archr_peaks.sh --libraries ARC07,ARC03-2 # pilot 2 libs
#
# Prerequisites:
#   * ~/R/4.4.2/ has ArchR + BSgenome.Anasplatyrhynchos.ZJU1.0 installed
#     (run pipeline/run_archr_env.sh first)
#   * /n/scratch/.../FC_08656/phase2/atac/step_f1_coupled_atac_metadata.tsv
#     exists (f1 must have run on all 12 libs, or at least the subset
#     requested via --libraries)
#   * MACS3 binary at /home/cl266/.conda/envs/snapatac2/bin/macs3 (from
#     the snapatac2 conda env)
#
# Resource model (post-pilot, full cohort):
#   * locking=TRUE eliminates intra-library HDF5 race (no subThreading
#     within a single library's chunks).
#   * threads=4 parallelises ACROSS libraries via mclapply — 4 libs
#     processed concurrently, each single-threaded internally.
#   * Per-library MaxRSS at threads=1 = ~43 GB (pilot 37895184).
#     4 × 43 GB ≈ 172 GB → 192 GB allocation (12% headroom).
#   * Wall: empirical from job 37913139 (cancelled at 5h):
#       - Parallel Arrow build (8 libs / 4 threads): ~4h ✓
#       - subsetArchRProject Arrow copy: ~30 min ✓
#       - addGroupCoverages with 21 callable pools = 50 (pool×rep)
#         groups × ~12 min/group ≈ 10 h (ALL HH13/HH17 took 4-13 min,
#         HH21/HH25 will be larger).  Blew the 12h short ceiling.
#     Switched to medium partition (5d max) with -t 1-00:00.  Most
#     jobs end up ~12-15h actual; medium gives us comfortable headroom.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
OUT_DIR="$SCRATCH/phase2/atac_archr"
mkdir -p "$OUT_DIR"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== ArchR peak-calling track ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Args: $*"

Rscript pipeline/archr_build_peaks.R --threads 4 "$@"

echo ""
echo "=== done ==="
echo "Finished: $(date)"
