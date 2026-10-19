#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:15
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/prep_g2_edge_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/prep_g2_edge_%j.err
#SBATCH -J prep_edge

# Build edge cellSNP directory for G2 by subsetting genotyping_all matrices.
# Edge = first 20K + last 20K cells (from 75,162 total with chunk_size=40K).
# This creates a contiguous-indexed cellSNP that vireo can consume directly.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH=${DUCK_SCRATCH_ROOT}
CELLSNP="$SCRATCH/GEX03-1A_GEX/outs/pseudobulk_group2_merged/genotyping_all"
EDGE_DIR="$SCRATCH/GEX03-1A_GEX/outs/pseudobulk_group2_merged/genotyping_edge_40000"
CELLS_PER_CHUNK=40000
EDGE_SIZE=$((CELLS_PER_CHUNK / 2))  # 20000

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

mkdir -p "$EDGE_DIR"

echo "=== Building G2 Edge CellSNP ==="
echo "Source: $CELLSNP"
echo "Output: $EDGE_DIR"

TOTAL_CELLS=$(wc -l < "$CELLSNP/cellSNP.samples.tsv")
echo "Total cells: $TOTAL_CELLS, Edge size: $EDGE_SIZE per side"

# Extract edge barcodes: first EDGE_SIZE + last EDGE_SIZE
head -n $EDGE_SIZE "$CELLSNP/cellSNP.samples.tsv" > "$EDGE_DIR/barcodes_edges.txt"
tail -n $EDGE_SIZE "$CELLSNP/cellSNP.samples.tsv" >> "$EDGE_DIR/barcodes_edges.txt"
EDGE_CELLS=$(wc -l < "$EDGE_DIR/barcodes_edges.txt")
echo "Edge barcodes: $EDGE_CELLS"

# Copy VCF (same variants)
cp "$CELLSNP/cellSNP.base.vcf.gz" "$EDGE_DIR/"

# Subset matrices via Python
python3 - "$CELLSNP" "$EDGE_DIR/barcodes_edges.txt" "$EDGE_DIR" << 'PYEOF'
import sys, scipy.io as sio, pandas as pd, numpy as np, os

genotype_dir, edge_bc_file, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]

print(f"Loading matrices from {genotype_dir}")
ad = sio.mmread(f"{genotype_dir}/cellSNP.tag.AD.mtx").tocsr()
dp = sio.mmread(f"{genotype_dir}/cellSNP.tag.DP.mtx").tocsr()
print(f"AD: {ad.shape}, DP: {dp.shape}")

oth_path = f"{genotype_dir}/cellSNP.tag.OTH.mtx"
oth = sio.mmread(oth_path).tocsr() if os.path.exists(oth_path) else None

full_bc = pd.read_csv(f"{genotype_dir}/cellSNP.samples.tsv", header=None)[0].tolist()
edge_bc = pd.read_csv(edge_bc_file, header=None)[0].tolist()

bc_to_idx = {bc: i for i, bc in enumerate(full_bc)}
indices = [bc_to_idx[bc] for bc in edge_bc if bc in bc_to_idx]
print(f"Matched {len(indices)} of {len(edge_bc)} edge barcodes")

ad_sub = ad[:, indices].tocoo()
dp_sub = dp[:, indices].tocoo()
sio.mmwrite(f"{output_dir}/cellSNP.tag.AD.mtx", ad_sub)
sio.mmwrite(f"{output_dir}/cellSNP.tag.DP.mtx", dp_sub)

if oth is not None:
    sio.mmwrite(f"{output_dir}/cellSNP.tag.OTH.mtx", oth[:, indices].tocoo())

pd.Series(edge_bc).to_csv(f"{output_dir}/cellSNP.samples.tsv", header=False, index=False)
print(f"Done: {ad_sub.shape[1]} cells, {ad_sub.shape[0]} variants")
PYEOF

# Verify
if [ ! -f "$EDGE_DIR/cellSNP.tag.AD.mtx" ] || [ ! -f "$EDGE_DIR/cellSNP.tag.DP.mtx" ]; then
    echo "ERROR: Edge matrix subsetting failed"
    exit 1
fi

echo "=== Edge CellSNP ready: $EDGE_DIR ==="
