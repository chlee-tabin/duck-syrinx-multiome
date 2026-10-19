#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4c_motif_agg_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4c_motif_agg_%j.err
#SBATCH -J g4c_motif_agg

# Run motif enrichment on the aggregated descriptive L vs R peak sets
# (g2c output) for both SnapATAC2 and ArchR peak sets.  Adapter:
# rename log2FC -> log2FoldChange and binom_padj -> padj so the
# existing g4_hh1716_motif_enrichment.R script reads it as a DESeq2
# table.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== g4c motif on aggregated ==="
date

# Adapter: convert g2c aggregated TSV to DESeq2-compatible schema
adapt_table() {
  local src="$1"
  local dst="$2"
  python3 - <<EOF
import pandas as pd
df = pd.read_csv("$src", sep="\t")
df = df.rename(columns={"log2FC": "log2FoldChange", "binom_padj": "padj"})
df["lfcSE"] = 0.0  # placeholder for downstream tools
df.to_csv("$dst", sep="\t", index=False, float_format="%.6f")
print(f"Wrote {len(df)} rows to {{'$dst'}}", flush=True)
EOF
}

# Variant 1: SnapATAC2 peak set
echo ""
echo "--- adapter: SnapATAC2 g2c -> DESeq2 schema ---"
mkdir -p phase2/atac/da_agg_snapatac
adapt_table \
  phase2/atac/da/g2c_aggregated_HH1716_INFERENCE_POOL.tsv \
  phase2/atac/da_agg_snapatac/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv

module unload conda/miniforge3
module load R/4.4.2

echo ""
echo "--- g4 motif enrichment (SnapATAC2 peaks, aggregated L vs R) ---"
DA_DIR="$SCRATCH/phase2/atac/da_agg_snapatac" \
  PEAKS_BED="$SCRATCH/phase2/atac/peaks/union_peaks.bed" \
  Rscript pipeline/g4_hh1716_motif_enrichment.R

# Variant 2: ArchR peak set
module unload R/4.4.2
module load conda/miniforge3
conda activate snapatac2

echo ""
echo "--- adapter: ArchR g2c -> DESeq2 schema ---"
mkdir -p phase2/atac/da_agg_archr
adapt_table \
  phase2/atac/da_archr/g2c_aggregated_HH1716_INFERENCE_POOL.tsv \
  phase2/atac/da_agg_archr/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv

module unload conda/miniforge3
module load R/4.4.2

echo ""
echo "--- g4 motif enrichment (ArchR peaks, aggregated L vs R) ---"
DA_DIR="$SCRATCH/phase2/atac/da_agg_archr" \
  PEAKS_BED="$SCRATCH/phase2/atac_archr/peaks/union_peaks.bed" \
  Rscript pipeline/g4_hh1716_motif_enrichment.R

echo ""
echo "=== DONE ==="
date
