# Environment capture

Captured **2026-08-26** for the public code release.
Two sides: the **laptop** env that renders figures from the committed TSVs
(Tier 1), and the **O2 / SLURM** envs the pipeline activates (Tier 2).

## Files

| file | what | produced by |
|---|---|---|
| `duckfig.yml` | laptop figure env — python + R (ggplot2, DESeq2), 340 lines | `conda env export -n duckfig --no-builds` (local, macOS/osx-64), then curated: channels reordered to conda-forge, bioconda; osx-arm64 cross-compiler and Google API client packages (not used by any figure) removed |
| `duckfig_portable.yml` | platform-neutral version of `duckfig.yml` for Linux and other platforms: the explicitly requested conda packages (python 3.12, r-base 4.4, ggplot2/tidyverse, DESeq2 1.46.0, …) pinned to the versions in `duckfig.yml`, no build strings, no macOS-only packages, plus the pip packages the figure code imports and the R/Bioconductor packages of the figS03/figS17 mixed models and motif scans (lme4, glmmTMB, emmeans, motifmatchr, TFBSTools, universalmotif, JASPAR2020, BSgenome, BSgenomeForge). Dry-run solves for linux-64 and osx-64 (conda 26.1.1, `libmamba`, strict channel priority). The custom `BSgenome.Anasplatyrhynchos.ZJU1.0` is built separately with `pipeline/archr_build_bsgenome.R` | hand-curated from `duckfig.yml` + the cluster R library |
| `o2_vireo.yml` | main pipeline env: `cellsnp-lite=1.2.3`, `vireosnp=0.5.9`, `harmonypy==0.2.0`, `leidenalg==0.11.0`, python 3.9.23 — activated by `pipeline/run_pipeline.sh` and every `a*`/`b*` demux step. scanpy/anndata do **not** appear in this export; the snRNA analysis used scanpy 1.10.3 (paper Methods), which is the floor pinned in `pipeline/environment.yml` | `conda env export -n vireo --no-builds` (O2) |
| `o2_snapatac2.yml` | chromatin/ATAC env: `snapatac2==2.9.0`, `macs3==3.0.4` (MACS3 is the peak caller; `macs2==2.2.9.1` is also installed), `pydeseq2==0.5.4`, `scanpy==1.11.5`, python 3.11.15 — `f*`/`g*`/`h*` chromatin steps | idem, `-n snapatac2` |
| `o2_gkmqc.yml` | gkm-SVM host env (python 3.7.12 + numpy/pandas) — `h2`/`h3`/`h8`/`h9` gkm-SVM steps. **The `gkmQC` and `lsgkm` binaries (`gkmtrain`, `gkmpredict`, `gkmkern`) are built from source, not conda packages** — recipe in `pipeline/setup_gkmqc_env.sh` | idem, `-n gkmqc` |
| `o2_cellbender.yml` | CellBender env: `cellbender==0.3.2`, `pyro-ppl==1.9.0` (pinned — 1.9.1 crashes at checkpoint save), `torch==2.2.2+cu121`, python 3.10.20. The `+cu121` torch/torchvision wheels are not on PyPI; the pip block's first entry `--extra-index-url https://download.pytorch.org/whl/cu121` points pip at the PyTorch CUDA 12.1 index. The yml alone does not avoid a crash at checkpoint save under torch 2.2: build the env with `scripts/cellbender/cb_build_env.sh`, which installs the same pins and drops in the patched `checkpoint.py` vendored next to it (`cellbender_checkpoint_v0.3.2_patched.py`, BSD 3-Clause, see the repository README's License section) | `conda env export -p $DUCK_SCRATCH_USER/conda_envs/cellbender --no-builds` (prefix env on scratch, not a named env) |
| `o2_jupyter.yml` | general python env — `pipeline/run_g3_cfc1_deseq2.sh`, `scripts/geo_cr9_rds_to_h5ad.sh` | idem, `-n jupyter` |
| `o2_stagef.yml` | conservation toolchain: `phast=1.9.7`, `multiz=11.2`, `lastz=1.04.52` — `scripts/stagef/*.sbatch` | idem, `-n stagef` |
| `o2_phyloacc.yml` | `phyloacc=2.5.0`, python 3.13.14 — `scripts/stagef/stagef_0{6,8,9}*.sbatch` | idem, `-n phyloacc` |
| — | the legacy 2024 `jlab` env named by `FC_08656/2024*_{cellsnp-lite,vireo}*.sh` could not be exported (it lived under a retired conda install). Not a reproduction gap: the same cellSNP-lite + vireo steps are covered by the maintained `vireo` env (`o2_vireo.yml`), which `pipeline/run_pipeline.sh` and every current a*/b* demux step activate; the FC_08656 `jlab` scripts are the 2024 provenance record of the same tools | — |
| `o2_R442_sessionInfo.txt` | R 4.4.2 `sessionInfo()` + the `.libPaths()` actually used | `Rscript -e 'sessionInfo()'` |
| `o2_R442_installed_packages.tsv` | 335 packages, Package/Version (incl. ArchR 1.0.3, DESeq2 1.46.0, Seurat 5.4.0) | `Rscript -e 'installed.packages()[,c("Package","Version")]'` |

The machine-specific `prefix:` lines of every `*.yml` were stripped after
capture (and `o2_cellbender.yml`'s path-valued `name:` set to `cellbender`), so
`conda env create -f envs/<file>.yml` creates a named env on any machine.

## How the O2 side was produced

One SLURM job (`-p priority -c 1 --mem 4G -t 0:20`) — no login-node compute:

```bash
module purge && module load gcc/14.2.0 conda/miniforge3
for E in vireo snapatac2 gkmqc jupyter jlab stagef phyloacc; do
    conda env export -n "$E" --no-builds > o2_${E}.yml     # or an o2_${E}.MISSING marker
done
conda env export -p "$DUCK_SCRATCH_USER/conda_envs/cellbender" --no-builds > o2_cellbender.yml

module purge && module load gcc/14.2.0 R/4.4.2
Rscript -e '.libPaths(c(path.expand("~/R/4.4.2"), .libPaths()))' -e 'sessionInfo()'
Rscript -e '.libPaths(c(path.expand("~/R/4.4.2"), .libPaths()))' \
        -e 'write.table(installed.packages()[,c("Package","Version")], stdout(), sep="\t", quote=FALSE, row.names=FALSE)'
```

**R on O2 is `module load gcc/14.2.0 R/4.4.2`** plus a user library at `~/R/4.4.2`
prepended to `.libPaths()` — the system module lib alone carries only base R.
`R_LIBS_USER` is *not* the mechanism (a stale `~/.Renviron` pins it elsewhere on
this account); every cluster R script prepends the user lib explicitly. The
install recipe for that library is already recorded in
**`pipeline/install_o2_r_packages.R`**.

## Not captured here

**cellranger-arc 2.1.0** is a binary install put on `PATH` by
`pipeline/rerun_arc_v210.sh`, not a conda env, so it has no yml here. CellBender,
vireo and cellSNP-lite are pinned by the ymls above (`cellbender==0.3.2`,
`vireosnp=0.5.9`, `cellsnp-lite=1.2.3`).

Scope: only envs that a **shipped `pipeline/` or `scripts/` file activates**. Other envs
named in comments or in superseded scripts are not part of the released pipeline.
