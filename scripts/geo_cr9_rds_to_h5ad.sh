#!/bin/bash
# GEO submission helper — convert the cr9 integrated scRNA master (Seurat v5 .rds)
# to AnnData .h5ad for deposit alongside the .rds (open, language-agnostic format).
#
# Robust two-step (no sceasy/SeuratDisk/zellkonverter in the user R lib):
#   Phase 1 (R/4.4.2 + Seurat user lib): export RNA counts (MTX) + lognorm + obs + embeddings.
#   Phase 2 (conda `jupyter`: scanpy/anndata/scipy): assemble + write .h5ad.
#
# Submit on o2 (SLURM only):  sbatch scripts/geo_cr9_rds_to_h5ad.sh
# Reads/writes under the GEO staging work dir on scratch.
#SBATCH -J cr9_h5ad
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-2:00
#SBATCH -o /n/scratch/users/c/cl266/geo_upload/work/cr9_h5ad_%j.log
#SBATCH -e /n/scratch/users/c/cl266/geo_upload/work/cr9_h5ad_%j.log
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


W=${DUCK_SCRATCH_USER}/geo_upload/work
RDS="$W/cr9_03_integrated_CB.rds"
H5AD="$W/cr9_03_integrated_CB.h5ad"
echo "[$(date)] START cr9 rds->h5ad on $(hostname); RDS=$RDS"
test -s "$RDS" || { echo "MISSING_RDS $RDS"; exit 2; }

# ---- Phase 1: R exports portable components --------------------------------
module purge; module load gcc/14.2.0 R/4.4.2
Rscript --vanilla /dev/stdin "$RDS" "$W" <<'RS'
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
.libPaths(c(file.path(.duck_env("DUCK_HOME", "/home/cl266"), "R", "4.4.2"), .libPaths()))
suppressMessages({library(Seurat); library(Matrix)})
a <- commandArgs(trailingOnly=TRUE); rds <- a[1]; W <- a[2]
obj <- readRDS(rds)
cat("class:", class(obj)[1], "| assays:", paste(Assays(obj), collapse=","),
    "| default:", DefaultAssay(obj), "\n")
asy <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cat("using assay:", asy, "\n")
# Seurat v5 may keep per-sample split layers; join them so counts/data are single matrices.
obj[[asy]] <- tryCatch(JoinLayers(obj[[asy]]),
                       error=function(e){cat("JoinLayers skipped:", conditionMessage(e), "\n"); obj[[asy]]})
cts <- tryCatch(GetAssayData(obj, assay=asy, layer="counts"),
                error=function(e) GetAssayData(obj, assay=asy, slot="counts"))
cat("counts dim:", paste(dim(cts), collapse=" x "), "\n")
dat <- tryCatch(GetAssayData(obj, assay=asy, layer="data"),
                error=function(e) GetAssayData(obj, assay=asy, slot="data"))   # v4/v5 API shim; no NULL tier
Matrix::writeMM(cts, file.path(W, "cr9_counts.mtx"))
writeLines(rownames(cts), file.path(W, "cr9_features.tsv"))
writeLines(colnames(cts), file.path(W, "cr9_barcodes.tsv"))
if (!is.null(dat) && nrow(dat) == nrow(cts) && ncol(dat) == ncol(cts) && sum(dim(dat)) > 0) {
  Matrix::writeMM(as(dat, "CsparseMatrix"), file.path(W, "cr9_data.mtx"))
  cat("wrote lognorm data layer\n")
}
write.csv(obj@meta.data, file.path(W, "cr9_obs.csv"))
for (rd in Reductions(obj)) write.csv(Embeddings(obj, rd), file.path(W, paste0("cr9_obsm_", rd, ".csv")))
cat("R_EXPORT_OK ncells=", ncol(cts), " ngenes=", nrow(cts),
    " reductions=", paste(Reductions(obj), collapse=","), "\n", sep="")
RS

# ---- Phase 2: Python assembles AnnData ------------------------------------
module purge; module load gcc/14.2.0 conda/miniforge3
eval "$(conda shell.bash hook)"; conda activate jupyter
python3 - "$W" "$H5AD" <<'PY'
import sys, os, re, glob
import numpy as np, pandas as pd, anndata as ad
from scipy.io import mmread
W, H5AD = sys.argv[1], sys.argv[2]
cts = mmread(os.path.join(W, "cr9_counts.mtx")).tocsr()           # genes x cells
features = [l.strip() for l in open(os.path.join(W, "cr9_features.tsv"))]
barcodes = [l.strip() for l in open(os.path.join(W, "cr9_barcodes.tsv"))]
obs = pd.read_csv(os.path.join(W, "cr9_obs.csv"), index_col=0)
obs.index = obs.index.astype(str)
obs = obs.loc[barcodes]
X = cts.T.tocsr()                                                 # cells x genes
adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(features)))
adata.var_names_make_unique()
adata.layers["counts"] = X.copy()
dpath = os.path.join(W, "cr9_data.mtx")
if os.path.exists(dpath):
    adata.layers["lognorm"] = mmread(dpath).T.tocsr()
for f in sorted(glob.glob(os.path.join(W, "cr9_obsm_*.csv"))):
    name = re.sub(r".*cr9_obsm_(.+)\.csv$", r"\1", f)
    emb = pd.read_csv(f, index_col=0); emb.index = emb.index.astype(str)
    adata.obsm["X_" + name] = emb.loc[barcodes].values.astype("float32")
adata.write_h5ad(H5AD)
print(f"H5AD_OK shape={adata.shape} layers={list(adata.layers)} obsm={list(adata.obsm)} -> {H5AD}")
PY

ls -lh "$H5AD"
echo "[$(date)] DONE cr9 rds->h5ad"
echo "CR9_H5AD_PIPELINE_OK"
