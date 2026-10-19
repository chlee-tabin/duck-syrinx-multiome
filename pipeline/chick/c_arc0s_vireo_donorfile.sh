#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 4
#SBATCH --mem 24G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/c_arc0s_vireo_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/c_arc0s_vireo_%j.err
#SBATCH -J chk_arc0s_v
#
# c_arc0s_vireo — REFERENCE-VCF vireo: project the ARC0S-chick cells onto the
# FIXED Group A donor genotypes (N18/seed_117), mirroring the duck Phase-2
# strategy (pipeline/a3_sub_vireo_g6_libpair.sh). Each ARC0S cell is assigned to
# one of the 18 Group A donors or left unassigned/doublet.
#
# PRECONDITION: c_arc0s_cellsnp_genotype.sh has produced
#   $WORK/arc0s/cellsnp/  (cellData at the SAME 283,435 positions as the donor VCF).
#
# INTERPRETATION (the science):
#   Donors that receive SUBSTANTIAL ARC0S cells == the HH15 + HH25 embryos whose
#   LEFT half was sorted into the ARC0S lane in 2024. ARC0S holds chick-15-1L..4L
#   + chick-25-1L..4L (8 embryos). This gives an INDEPENDENT, PHYSICAL (non-PITX2)
#   LEFT-side anchor for those 8 Group A donors -> breaks the L/R-assignment
#   circularity for Group A HH15/HH25.
#   Duck-contaminant ARC0S cells (2024 analyst flagged heavy duck contamination)
#   do not match any chick donor genotype -> land as unassigned / low prob_max.
#
# Donor VCF carries all 18 N18 donors (donor0..17); ~6 are phantoms (2-11 cells)
# with flat PL (5,5,5) genotypes -> uninformative, won't spuriously attract cells.
# The 12 real embryos are the ~450-970-cell donors. We DELIBERATELY keep all 18
# (no pre-subset) so a cell matching a real donor isn't forced onto a phantom.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux
CELL="$WORK/arc0s/cellsnp"
DONOR_VCF="$WORK/vireo_seeded/N18/seed_117/GT_donors.vireo.vcf.gz"
SEED="${SEED:-117}"
OUTDIR="$WORK/arc0s/vireo_donorfile/seed_${SEED}"
mkdir -p "$OUTDIR" "$WORK/logs"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cat > "$OUTDIR/run_params.txt" <<EOF
GROUP=chick_ARC0S_projected_onto_GroupA_N18seed117
MODE=reference (--donorFile, fixed donor GT; mirrors a3_sub_vireo_g6_libpair.sh)
GENO_TAG=PL
SEED=$SEED
CELLDATA=$CELL
DONOR_VCF=$DONOR_VCF
SLURM_JOB_ID=$SLURM_JOB_ID
DATE=$(date -Iseconds)
EOF

echo "=== c_arc0s reference-VCF vireo (project ARC0S onto Group A donors) ==="
echo "cellData:  $CELL"
echo "donor VCF: $DONOR_VCF ($(zcat "$DONOR_VCF" | grep '^#CHROM' | tr '\t' '\n' | tail -n +10 | wc -l) donors)"
test -s "$CELL/cellSNP.samples.tsv" || { echo "ERROR: ARC0S cellData missing (run c_arc0s_cellsnp_genotype.sh first)"; exit 1; }

# Reference mode: vireo infers donor count from the VCF (no --nDonor), as in the
# duck G6 libpair reference run. PL genoTag (default) = soft donor genotypes.
vireo \
    --cellData "$CELL" \
    --donorFile "$DONOR_VCF" \
    --genoTag PL \
    --randSeed "$SEED" \
    -p 4 \
    --outDir "$OUTDIR"

echo "=== c_arc0s vireo done -> $OUTDIR ==="
echo "--- donor assignment counts (summary.tsv) ---"
cat "$OUTDIR/summary.tsv"
