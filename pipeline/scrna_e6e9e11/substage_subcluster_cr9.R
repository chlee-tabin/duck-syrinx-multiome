#!/usr/bin/env Rscript
# =====================================================================
# cr9 SUBSTAGE subclustering — the cr9-CB analogue of substage_subcluster_v2.R.
# IDENTICAL recipe (so cr9-vs-cr6 substage differences reflect the cr9-CB master,
# NOT a changed method): CP10K log-norm · 2000 HVG · regress pct.mito+S+G2M
# (+ sex_numeric, the §4a addition) · PCA 30 · Harmony(library) · Louvain res-sweep.
#
# DIFFS vs substage_subcluster_v2.R (the only changes):
#   1. RDS -> cr9 master  (old_duck_cr9/master/cr9_03_integrated_CB.rds)
#   2. OUT -> old_duck_cr9/substage  (SEPARATE dir; cr6 outputs sit beside it for comparison)
#   3. CC  -> + "sex_numeric"  (§4a: regress the sex axis like cell-cycle; CATEGORICAL
#             M=0/F=1/unknown=0.5 from the v2 genotype demux — multiome d1 + chick cAE
#             parity, REPLACING the old continuous sex.score)
#   4. demux pre-filter after load  (§4b; defensive — the cr9 master is already
#      singlet-filtered at integration, so this removes ~0, but documents intent
#      and protects if the master ever carries raw calls)
# Exploratory only.  R/4.4.2, o2 SLURM.
# =====================================================================
## portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
# empty-normalizing env read (same body as scripts/_paths.R): a set-but-empty DUCK_* falls back to the default
.duck_env <- function(n, d) { v <- Sys.getenv(n); if (nzchar(v)) v else d }
SCRATCH_USER <- .duck_env("DUCK_SCRATCH_USER", "/n/scratch/users/c/cl266")
HOME_ROOT <- .duck_env("DUCK_HOME", "/home/cl266")

suppressWarnings(suppressMessages({
  .libPaths(c(paste0(HOME_ROOT, "/R/4.4.2"), .libPaths()))
  library(Seurat); library(Matrix); library(jsonlite); library(harmony)
}))
set.seed(1)
options(future.globals.maxSize = 8 * 1024^3)
try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)

RDS <- paste0(SCRATCH_USER, "/old_duck_cr9/master/cr9_03_integrated_CB.rds")   # (1) cr9 master
OUT <- paste0(SCRATCH_USER, "/old_duck_cr9/substage")                          # (2) separate dir
# Repo destination for the figure-facing substage_E{09,11}_percell.tsv.gz (copied back from
# the cluster): figures/data/scrna_e6e9e11/.  Other OUT products stay scratch-only.
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
t0 <- Sys.time()

say("== loading %s", RDS)
obj <- readRDS(RDS); DefaultAssay(obj) <- "RNA"
# cr9 master (Seurat v5) carries per-library SPLIT RNA layers (counts.D1/.D2/.D3,
# data.*) from the merge — GetAssayData(layer=...) errors on multiple layers.
# Join into single counts/data (no-op if already joined; the cr6 master was single-layer).
obj <- JoinLayers(obj)
say("loaded+joined in %.1f min; RNA genes=%d cells=%d",
    as.numeric(difftime(Sys.time(), t0, "mins")), nrow(obj[["RNA"]]), ncol(obj))

# ---- (4) demux pre-filter (defensive — master already applies genotype carry-forward) ----
# Keep genotype-valid cells; drop true doublets only. HTO-Negative cells are KEPT
# (genotype carry-forward: stage from genotype, LR=NA) — they are valuable for
# clustering and only laterality-ambiguous.  No-op on the cr9 master (already filtered).
md0 <- obj@meta.data; keep <- rep(TRUE, ncol(obj))
if ("donor_id" %in% colnames(md0))
  keep <- keep & !(md0$donor_id %in% c("doublet", "unassigned"))
if ("scDblFinder_class" %in% colnames(md0))
  keep <- keep & (is.na(md0$scDblFinder_class) | md0$scDblFinder_class != "doublet")
if ("final.calls" %in% colnames(md0))   # NOTE: the v2 master no longer carries final.calls
  keep <- keep & (is.na(md0$final.calls) | md0$final.calls != "Doublet")   # (dead guard; doublet removal is via donor_id + scDblFinder above). Kept for a hypothetical raw-call master.
say("demux pre-filter: keep %d / %d cells (HTO-Negative kept via genotype)", sum(keep), ncol(obj))
obj <- subset(obj, cells = colnames(obj)[keep])

genes_all <- rownames(obj[["RNA"]]); NKX61 <- "LOC113843492"
md <- obj@meta.data
C_stage <- if ("stage.called" %in% colnames(md)) "stage.called" else "stage"
C_lib   <- if ("library" %in% colnames(md)) "library" else "orig.ident"
C_glob  <- if ("SCT_snn_res.0.4" %in% colnames(md)) "SCT_snn_res.0.4" else "seurat_clusters"
# (3) cell-cycle / mito / SEX regression covariates (must exist for ccreg)
# sex_numeric (CATEGORICAL M=0/F=1/unknown=0.5 from the v2 genotype demux) REPLACES the
# old continuous sex.score — multiome d1 / chick cAE parity.
CC <- intersect(c("percent.mito", "S.Score", "G2M.Score", "sex_numeric"), colnames(md))
say("regress covariates available: %s", paste(CC, collapse=","))
stopifnot(length(CC) >= 1)

# ---- mesenchyme filter (marker-based, identical to v2) --------------
obj <- NormalizeData(obj, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
logm <- GetAssayData(obj, assay="RNA", layer="data")
CONTAM <- list(immune=intersect(c("PTPRC","LCP1","CD74"),genes_all),
               endo=intersect(c("CDH5","PECAM1","EGFL7","KDR"),genes_all),
               epith=intersect(c("EPCAM","CDH1","CLDN1","WNT7B"),genes_all),
               smc=intersect(c("ACTA2","MYH11","MYL9","TAGLN","DES"),genes_all))
LPMP <- intersect(c("FOXF1","PRRX1","TWIST1","GATA4","TBX5","COL2A1","WT1","OSR1","ISL1"),genes_all)
glob <- as.character(md[[C_glob]]); cls <- sort(unique(glob))
fp <- function(gs, cc) if(length(gs)==0||length(cc)==0) NA_real_ else mean(Matrix::colSums(logm[gs,cc,drop=FALSE]>0)>0)
ctab <- data.frame(cluster=cls)
for(k in names(CONTAM)) ctab[[k]] <- sapply(cls, function(c) fp(CONTAM[[k]], colnames(obj)[glob==c]))
ctab$lpm <- sapply(cls, function(c) fp(LPMP, colnames(obj)[glob==c]))
drop <- with(ctab,(immune>0.5)|(endo>0.5)|(epith>0.5)|(smc>0.6 & lpm<0.25))
kept <- ctab$cluster[!drop]
say("mesenchyme filter: drop=%s keep=%s", paste(ctab$cluster[drop],collapse=","), paste(kept,collapse=","))
mes <- subset(obj, cells = colnames(obj)[glob %in% kept])
say("kept %d / %d cells", ncol(mes), ncol(obj))

# ---- POOLS-aligned module-score panels (kept for the table) ---------
PANELS <- list(
  Ventral_spl=intersect(c("FOXF1","GATA4","ZFPM2","HAND2","NKX2-3","TBX18"),genes_all),
  Dorsal_subepi=intersect(c(NKX61,"WNT2B","WNT2","MSX1","MSX2"),genes_all),
  Chondrogenic=intersect(c("SOX9","COL2A1","ACAN","RUNX2","WNT5A","DPP10","FNDC1","CNMD"),genes_all),
  Somatic=intersect(c("TBX5","TBX4","PITX1","PRRX1","TWIST1"),genes_all),
  Mesothelium=intersect(c("WT1","TCF21"),genes_all),
  VenousPole=intersect(c("TBX20","GATA4","ZFPM2"),genes_all),
  BARX1=intersect(c("BARX1"),genes_all),
  EarlyLPM=intersect(c("ISL1","OSR1","TBX5","FOXF1"),genes_all))

EXPORT_GENES <- intersect(unique(c(
  "PITX2","LOC106016802","LOC106016804","LOC106017192","FAM241A","LOC113843539",NKX61,"LOC101798108",
  "FOXF1","GATA4","GATA5","GATA6","ZFPM2","TBX18","TBX20","HAND2","NKX2-3","NKX2-8","ISL1","OSR1","FOXD1",
  "WNT2B","WNT2","MSX1","MSX2","SOX9","COL2A1","ACAN","RUNX2","WNT5A","WNT5B","DPP10","FNDC1","CNMD","SCX",
  "TBX5","TBX4","PITX1","PRRX1","TWIST1","PAX1","PAX9","PAX3","PAX7","FOXC2","WT1","TCF21","NR2F2","BARX1",
  "PTPRC","CDH5","PECAM1","EPCAM","ACTA2","MYH11","BMP4","BMP7","BMP5","SHISA2","SFRP2","GREM1","GREM2",
  "ESR1","ESR2","ESRRG","FOXL2","AMH","DMRT1","MKI67","TOP2A","PCNA","CCNB2")), genes_all)

RES_SWEEP <- c(0.1,0.2,0.3,0.5)
STAGES <- sort(unique(as.character(mes@meta.data[[C_stage]])))
say("stages: %s", paste(STAGES, collapse=", "))

process_stage <- function(s){
  say("\n#### STAGE %s (d2g-mirror: CP10K+ccreg+Harmony) ####", s)
  ss <- subset(mes, cells = colnames(mes)[as.character(mes@meta.data[[C_stage]])==s])
  DefaultAssay(ss) <- "RNA"
  ss <- NormalizeData(ss, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
  ss <- FindVariableFeatures(ss, nfeatures=2000, verbose=FALSE)
  say("[%s] cells=%d  regress=%s", s, ncol(ss), paste(CC,collapse="+"))
  ss <- ScaleData(ss, vars.to.regress=CC, verbose=FALSE)
  ss <- RunPCA(ss, npcs=30, verbose=FALSE, seed.use=1)
  nlib <- length(unique(as.character(ss@meta.data[[C_lib]])))
  rd <- "pca"
  if(nlib >= 2){
    ss <- RunHarmony(ss, group.by.vars=C_lib, dims.use=1:30, verbose=FALSE)
    rd <- "harmony"; say("[%s] Harmony on %d libraries", s, nlib)
  } else say("[%s] 1 library -> skip Harmony", s)
  ss <- FindNeighbors(ss, reduction=rd, dims=1:30, verbose=FALSE)
  algo <- 1L  # Louvain (robust, no python dep); d2g uses scanpy-leiden (2nd-order)
  for(r in RES_SWEEP) ss <- FindClusters(ss, resolution=r, algorithm=algo, random.seed=1, verbose=FALSE)
  ss <- RunUMAP(ss, reduction=rd, dims=1:30, seed.use=1, verbose=FALSE)
  ss <- AddModuleScore(ss, features=PANELS, assay="RNA", name="MS_", seed=1, ctrl=50)
  msc <- paste0("MS_", seq_along(PANELS))
  colnames(ss@meta.data)[match(msc, colnames(ss@meta.data))] <- paste0("ms_", names(PANELS))

  res_cols <- intersect(paste0("RNA_snn_res.", RES_SWEEP), colnames(ss@meta.data))
  logS <- as.matrix(GetAssayData(ss,assay="RNA",layer="data")[EXPORT_GENES,,drop=FALSE])
  cntS <- as.matrix(GetAssayData(ss,assay="RNA",layer="counts")[EXPORT_GENES,,drop=FALSE])
  um <- Embeddings(ss,"umap")
  mk <- intersect(c(C_stage,"donor_id","embryo_id","sex.called","sex_numeric","breed","LR.called",C_lib,"nCount_RNA","S.Score","G2M.Score","Phase"), colnames(ss@meta.data))
  pcdf <- data.frame(cell_id=colnames(ss),
    barcode=if("barcode"%in%colnames(ss@meta.data)) ss@meta.data$barcode else colnames(ss),
    ss@meta.data[,mk,drop=FALSE], ss@meta.data[,res_cols,drop=FALSE],
    ss@meta.data[,paste0("ms_",names(PANELS)),drop=FALSE],
    UMAP_1=um[,1], UMAP_2=um[,2], check.names=FALSE, stringsAsFactors=FALSE)
  logdf <- as.data.frame(t(logS)); colnames(logdf) <- paste0("log__",colnames(logdf))
  cntdf <- as.data.frame(t(cntS)); colnames(cntdf) <- paste0("cnt__",colnames(cntdf))
  pcdf <- cbind(pcdf, logdf, cntdf)
  gz <- gzfile(file.path(OUT,sprintf("substage_%s_percell.tsv.gz",s)),"w")
  write.table(pcdf, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)

  dotr <- list(); msr <- list()
  for(rc in res_cols){
    cl <- as.character(ss@meta.data[[rc]])
    for(g in EXPORT_GENES){ v <- logS[g,]
      for(cc in sort(unique(cl))){ i <- cl==cc
        dotr[[length(dotr)+1]] <- data.frame(resolution=rc,cluster=cc,gene=g,
          mean_log=mean(v[i]),frac_pos=mean(v[i]>0),n=sum(i),stringsAsFactors=FALSE)}}
    for(p in names(PANELS)){ msv <- ss@meta.data[[paste0("ms_",p)]]
      for(cc in sort(unique(cl))){ i <- cl==cc
        msr[[length(msr)+1]] <- data.frame(resolution=rc,cluster=cc,panel=p,
          mean_score=mean(msv[i]),n=sum(i),stringsAsFactors=FALSE)}}
  }
  write.table(do.call(rbind,dotr), file.path(OUT,sprintf("substage_%s_dotplot.tsv",s)),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,msr), file.path(OUT,sprintf("substage_%s_modscore.tsv",s)),sep="\t",quote=FALSE,row.names=FALSE)
  say("[%s] n clusters: %s", s, paste(sprintf("%s=%d",res_cols,sapply(res_cols,function(rc)length(unique(ss@meta.data[[rc]])))),collapse="  "))
  # Wilcoxon DE per subcluster (data-driven corroboration, like d2g) at the
  # candidate primary resolutions.
  for(r in c(0.2,0.3)){
    rc <- paste0("RNA_snn_res.", r); if(!rc %in% colnames(ss@meta.data)) next
    Idents(ss) <- factor(as.character(ss@meta.data[[rc]]))
    de <- tryCatch(FindAllMarkers(ss, assay="RNA", only.pos=TRUE, min.pct=0.25,
                                  logfc.threshold=0.25, verbose=FALSE),
                   error=function(e){say("[%s %s] DE err %s",s,rc,conditionMessage(e));NULL})
    if(!is.null(de) && nrow(de)>0){
      de <- de[order(de$cluster,-de$avg_log2FC),]
      de$rank <- ave(seq_len(nrow(de)), de$cluster, FUN=seq_along)
      de <- de[de$rank<=40,]
      write.table(de[,c("cluster","rank","gene","avg_log2FC","pct.1","pct.2","p_val_adj")],
        file.path(OUT,sprintf("substage_%s_de_res%s.tsv",s,r)),sep="\t",quote=FALSE,row.names=FALSE)
    }
  }
  saveRDS(ss, file.path(OUT, sprintf("substage_%s_obj.rds", s)))   # for reuse
  invisible(NULL)
}
for(s in STAGES){
  ok <- tryCatch({process_stage(s); TRUE}, error=function(e){say("[%s] ERROR %s",s,conditionMessage(e));FALSE})
  say("[%s] done=%s", s, ok)
}
writeLines(toJSON(list(method="cr9-CB substage d2g-mirror Seurat: CP10K+2000HVG+regress(pct.mito,S,G2M,sex_numeric)+PCA30+Harmony(library)+Louvain",
  rds=RDS, n_total=ncol(obj), mesenchyme_kept=ncol(mes), kept_clusters=kept,
  res_sweep=RES_SWEEP, regress=CC, stages=STAGES,
  seurat=as.character(packageVersion("Seurat")), date=as.character(Sys.time())),
  pretty=TRUE, auto_unbox=TRUE), file.path(OUT,"substage_provenance.json"))
say("\n== cr9 substage DONE in %.1f min ==", as.numeric(difftime(Sys.time(),t0,"mins")))
cat("SBATCH_DONE_OK\n")
