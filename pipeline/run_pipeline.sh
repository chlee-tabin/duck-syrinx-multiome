#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-6:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/pipeline_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/pipeline_%j.err
#SBATCH -J duck_pipe

# Production pipeline: single source of truth for all phases.
# Phases: a0 (genotype map) → a1 (pre-CB filter) → b1 (consolidation) →
#         b1b (overmerge) → b2 (clean genotypes) → bv (validate) →
#         c1 (QC) → c2 (genotype assign) → c3 (sex) → c4 (embryo table) →
#         d0 (patch meta) → d1 (Harmony) → d2 (clustering) →
#         d2g (stage-restricted subclustering + compendium annotations) →
#         d2gs (signaling-gene sidecar re-export; standalone, covered by d2g) →
#         d2i (PITX2 RNA pileup) → d2j (PITX2 isoform per-genotype + per-cell) →
#         d2k (PITX2c laterality QC) →
#         e1 (pseudobulk) →
#         h1 (pseudobulk for L vs R DEG) → h2 (DESeq2 L vs R) → fig
#
# Input: CellBender h5 files + vireo outputs on scratch.
#
# Usage:
#   sbatch pipeline/run_pipeline.sh                    # full end-to-end
#   sbatch pipeline/run_pipeline.sh --from d1           # recluster only
#   sbatch pipeline/run_pipeline.sh --from b1           # genotype rebuild (includes c3, c4)
#   sbatch pipeline/run_pipeline.sh --from c3 --to c4   # re-export sex + embryo table
#   sbatch pipeline/run_pipeline.sh --from e1 --to fig  # e1 + figures
#   sbatch pipeline/run_pipeline.sh --from h1 --to h2   # rerun L/R DEG only
#
# Resource guide (override SBATCH defaults with sbatch CLI args):
#   Phase range     Memory   Time    Example
#   a0-a1 (prep)    32G      0:30    sbatch --mem=32G -t 0-0:30 pipeline/run_pipeline.sh --from a0 --to a1
#   b1-b2 (consol.) 32G      0:30    sbatch --mem=32G -t 0-0:30 pipeline/run_pipeline.sh --from b1 --to b2
#   c1-c4 (QC+sex)  32G      0:30    sbatch --mem=32G -t 0-0:30 pipeline/run_pipeline.sh --from c1 --to c4
#   d1 (Harmony)    32G      0:30    sbatch --mem=32G -t 0-0:30 pipeline/run_pipeline.sh --from d1 --to d1
#   d2 (markers)    16G      1:00    sbatch --mem=16G -t 0-1:00 pipeline/run_pipeline.sh --from d2 --to d2
#   d2gs (signaling) 32G     1:30    use pipeline/run_d2gs.sh (right-sized header)
#   d2i-d2k (PITX2) 32G      1:30    sbatch --mem=32G -t 0-1:30 pipeline/run_pipeline.sh --from d2i --to d2k
#   h1-h2 (LR DEG)  48G      1:30    sbatch --mem=48G -t 0-1:30 pipeline/run_pipeline.sh --from h1 --to h2
#   h3 (detection)  24G      1:00    use pipeline/run_h3.sh (right-sized header; O2 RC over-request notice 2026-08-25)
#   e1-fig          32G      0:15    sbatch --mem=32G -t 0-0:15 pipeline/run_pipeline.sh --from e1 --to fig
#   full (a0-fig)   64G      6:00    (defaults)

set -euo pipefail
FIG_FAILS=()   # figure-stage failures (collected; see the figure stage)
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


# ── Canonical stage list ──────────────────────────
# Both forms encode the same 4 stage buckets.  CANONICAL_STAGES_LIST is
# a Bash array consumed by d2g (subclustering + signaling export); the
# embedded "HH17/16 HH16" element groups raw HH16 cells under HH17/16
# per d2g_stage_restricted_subclustering.py's pre-TES split principle
# (HH17 and HH16 share a Leiden clustering, emitted under the merged
# "HH17/16" label).  CANONICAL_STAGES_CSV is the merged-label CSV
# consumed by h2 (DESeq2 sweep).  Adding/removing a stage requires
# editing both.
CANONICAL_STAGES_LIST=(HH13 "HH17/16 HH16" HH21 HH25)
CANONICAL_STAGES_CSV="HH13,HH17/16,HH21,HH25"

SCRATCH="${DUCK_SCRATCH_ROOT}"
PHASE2="$SCRATCH/phase2"
cd "$SCRATCH"

# ── Phase definitions ──────────────────────────────────────────────────
VALID_PHASES="a0 a1 b1 b1b b2 bv c1 c2 c3 c4 d0 d1 d2 d2g d2gs d2i d2j d2k e1 h1 h2 h3 fig"

validate_phase() {
    local phase="$1" label="$2"
    local found=false
    for p in $VALID_PHASES; do [ "$p" = "$phase" ] && found=true; done
    if ! $found; then
        echo "Error: Invalid $label phase '$phase'. Valid: $VALID_PHASES" >&2
        exit 1
    fi
}

# ── Parse arguments ────────────────────────────────────────────────────
START="a0"
END="fig"

while [ $# -gt 0 ]; do
    case "$1" in
        --from)
            [ -z "${2:-}" ] && echo "Error: --from requires a phase" >&2 && exit 1
            START="$2"; shift 2 ;;
        --from=*)
            START="${1#--from=}"; shift ;;
        --to)
            [ -z "${2:-}" ] && echo "Error: --to requires a phase" >&2 && exit 1
            END="$2"; shift 2 ;;
        --to=*)
            END="${1#--to=}"; shift ;;
        *)
            # Bare argument = start phase (backward compat)
            START="$1"; shift ;;
    esac
done

validate_phase "$START" "start"
validate_phase "$END" "end"

should_run() {
    local phase="$1"
    local in_range=false
    for p in $VALID_PHASES; do
        [ "$p" = "$START" ] && in_range=true
        if [ "$p" = "$phase" ] && $in_range; then return 0; fi
        [ "$p" = "$END" ] && in_range=false
    done
    return 1
}

# Verify --from comes before --to
_start_ok=false
for _p in $VALID_PHASES; do
    [ "$_p" = "$START" ] && _start_ok=true
    [ "$_p" = "$END" ] && break
done
if ! $_start_ok; then
    echo "Error: --from '$START' must come before --to '$END'" >&2
    exit 1
fi

# ── Environment setup ──────────────────────────────────────────────────
module purge
module load gcc/14.2.0 samtools conda/miniforge3 R/4.4.2
conda activate vireo

# NUMBA_CACHE_DIR — required on o2 compute nodes.  scanpy is installed
# in user-site (/home/cl266/.local/lib/python3.9/site-packages/scanpy),
# and numba's default "write next to the source" cache strategy fails
# with `cannot cache function 'sparse_mean_var_minor_axis': no locator
# available` because the user-site's __pycache__ isn't writable from
# compute.  Point numba at a scratch-side cache instead.  (Verified
# 2026-04-17: without this, scanpy imports crash immediately — see
# job 37302988.)
export NUMBA_CACHE_DIR="$SCRATCH/.numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"

echo "=== Installing dependencies ==="
# Pin harmonypy only. Do NOT pip install numpy or scikit-misc —
# they are ABI-incompatible with conda-forge packages and cause segfaults.
# HVG uses seurat_v3_paper (pure scipy) instead of seurat_v3 (skmisc).
#
# The comment said "pin" while the command was a bare unpinned install, so every rerun could
# silently take a new Harmony implementation. the manuscript Methods report v0.2.0, which
# makes that the version of record; pin to it and report what actually got installed.
pip install --quiet "harmonypy==0.2.0"
python3 -c "import harmonypy; print('harmonypy', harmonypy.__version__)"

echo ""
echo "=== Production Pipeline (from: $START, to: $END) ==="
echo "Date: $(date)"

# ── Phase A: Genotype Map + Pre-CB Filter ─────────────────────────────
if should_run a0; then
    echo ""
    echo "--- Phase A0: Produce CB Genotype Map ---"
    python3 FC_08656/20260326_35_produce_cb_genotype_map_v2.py
else
    echo "--- Phase A0: SKIPPED ---"
fi

if should_run a1; then
    echo ""
    echo "--- Phase A1: Pre-CellBender Genotype Filter ---"
    python3 FC_08656/20260326_step02e_preCB_genotype.py
else
    echo "--- Phase A1: SKIPPED ---"
fi

# ── Phase B: Genotype Consolidation ────────────────────────────────────
if should_run b1; then
    echo ""
    echo "--- Phase B1: Donor Consolidation ---"
    python3 pipeline/b1_donor_consolidation.py
else
    echo "--- Phase B1: SKIPPED ---"
fi

if should_run b1b; then
    echo ""
    echo "--- Phase B1b: Over-merge Reclassification ---"
    python3 pipeline/b1b_overmerge_reclassify.py
else
    echo "--- Phase B1b: SKIPPED ---"
fi

if should_run b2; then
    echo ""
    echo "--- Phase B2: Clean Genotypes (apply b1+b1b to h5ad) ---"
    python3 FC_08656/20260330_step02f_clean_genotypes.py
else
    echo "--- Phase B2: SKIPPED ---"
fi

if should_run bv; then
    echo ""
    echo "--- Phase BV: Validate Genotype-Design Consistency ---"
    python3 pipeline/validate_genotype_design.py
else
    echo "--- Phase BV: SKIPPED ---"
fi

# ── Phase C: Cell QC + Genotype + Sex ──────────────────────────────────
if should_run c1; then
    echo ""
    echo "--- Phase C1: QC Filter ---"
    python3 pipeline/c1_qc_filter.py
else
    echo "--- Phase C1: SKIPPED ---"
fi

if should_run c2; then
    echo ""
    echo "--- Phase C2: Genotype Assignment ---"
    python3 pipeline/c2_genotype_assign.py
else
    echo "--- Phase C2: SKIPPED ---"
fi

if should_run c3; then
    echo ""
    echo "--- Phase C3: Sex Determination ---"
    python3 pipeline/c3_sex_determination.py
else
    echo "--- Phase C3: SKIPPED ---"
fi

if should_run c4; then
    echo ""
    echo "--- Phase C4: Build Genotype-Embryo Table ---"
    python3 pipeline/build_genotype_embryo_table.py
else
    echo "--- Phase C4: SKIPPED ---"
fi

# ── Phase D0: Patch genotype metadata (skip Harmony/clustering) ────────
# D0 is a lightweight alternative to D1 — only patches genotype labels
# without re-running Harmony. Skip D0 if D1 will also run (D1 overwrites step04).
if should_run d0 && ! should_run d1; then
    echo ""
    echo "--- Phase D0: Patch genotype metadata into step04 ---"
    echo "    (skips Harmony/clustering when only genotype labels changed)"
    python3 pipeline/d0_patch_genotype_metadata.py
else
    echo "--- Phase D0: SKIPPED (D1 will rebuild step04) ---"
fi

# ── Phase D: Clustering ────────────────────────────────────────────────
if should_run d1; then
    echo ""
    echo "--- Phase D1: Normalize + PCA + Harmony ---"
    python3 pipeline/d1_normalize_pca.py
else
    echo "--- Phase D1: SKIPPED ---"
fi

if should_run d2; then
    echo ""
    echo "--- Phase D2: Clustering + Markers ---"
    python3 pipeline/d2_cluster_markers.py
    # Publish the UMAP panel's single declared input on EVERY d2 run (not only in `fig`), so a
    # partial --from d2 --to d2 rerun can never leave stale coords under fresh annotations.
    cp -v "$PHASE2/figures/step04a_umap_clusters.tsv.gz" "figures/data/d2_umap_boundary_data.tsv.gz"
else
    echo "--- Phase D2: SKIPPED ---"
fi

# ── Phase D2b: Cluster validation (auto-assign cell types from markers) ──
if should_run d2; then
    echo ""
    echo "--- Phase D2b: Validate cluster annotations ---"
    python3 pipeline/d2b_validate_clusters.py
fi

# ── Phase D2e: Silhouette export (per-cluster + focused-mes silhouette) ──
# Runs after d2b because it reads the freshly-written
# references/cluster_annotations.tsv to look up current cluster IDs
# dynamically (no hardcoded leiden IDs).  Writes
# d2e_silhouette_global_res02.tsv + d2e_silhouette_focused_mes.tsv
# which both appear in pipeline/derived_data_manifest.tsv.
if should_run d2; then
    echo ""
    echo "--- Phase D2e: Silhouette export ---"
    python3 pipeline/d2e_silhouette_export.py
fi

# ── Phase D2g: Stage-restricted subclustering + signaling export ──────
# Three steps chained in this phase:
#
#   1. d2g_stage_restricted_subclustering.py — 4 stages × 2 variants
#      (default + --regress-cc).  Writes d2g_lpm_subclustering_*.tsv.gz
#      to figures/data/.  All four stages must run together — leaving
#      any out risks subcluster files holding cells bucketed under the
#      wrong stage (the 2026-04-17 regression where HH17/16 + HH21
#      shared g05/g09/g10).
#   2. d2g_validate_subclusters.py — invokes the compendium scorer
#      (scripts/annotate_d2g_via_compendium.py), applies the HH17/16
#      Dorsal-biased LPM override, writes
#      references/d2g_subcluster_annotations.tsv at canonical resolutions
#      (HH13=0.2, HH17/16=HH21=HH25=0.3).  Hard-fails on a stale or
#      mismatched annotation table.
#   3. d2g_signaling_export.py — 4 stages, per-cell signaling-gene
#      CP10K sidecar for the figS02E extended dotplot.
#
# Wired between d2 and e1/h1 because h1 pseudobulk scripts load the
# d2g_lpm_subclustering_*.tsv.gz sidecars; without this phase those
# sidecars stay stale relative to a freshly-rebuilt step04 h5ad
if should_run d2g; then
    echo ""
    echo "--- Phase D2g: Stage-restricted subclustering ---"
    for stage_args in "${CANONICAL_STAGES_LIST[@]}"; do
        # shellcheck disable=SC2086
        python3 pipeline/d2g_stage_restricted_subclustering.py --stage $stage_args
        # shellcheck disable=SC2086
        python3 pipeline/d2g_stage_restricted_subclustering.py --stage $stage_args --regress-cc
    done

    echo ""
    echo "--- Phase D2g: Compendium scorer + subcluster annotations ---"
    python3 pipeline/d2g_validate_subclusters.py

    echo ""
    echo "--- Phase D2g: Signaling-gene CP10K sidecar ---"
    for stage_args in "${CANONICAL_STAGES_LIST[@]}"; do
        # shellcheck disable=SC2086
        python3 pipeline/d2g_signaling_export.py --stage $stage_args
    done
else
    echo "--- Phase D2g: SKIPPED ---"
fi

# ── Phase D2gs: signaling-gene sidecar ONLY (no re-subclustering) ─────
# Standalone re-export of d2g_signaling_<stage>.tsv.gz — for gene-list
# extensions (2026-08-27: +BMP7/WNT5B/SFRP2/SHISA2 for the S15 onset
# panel) without touching leiden or the annotation chain (the export
# deliberately avoids the full d2g rerun; see d2g_signaling_export.py).
# Skipped when phase d2g is in range (d2g already runs the export).
if should_run d2gs && ! should_run d2g; then
    echo ""
    echo "--- Phase D2gs: Signaling-gene CP10K sidecar (standalone) ---"
    for stage_args in "${CANONICAL_STAGES_LIST[@]}"; do
        # shellcheck disable=SC2086
        python3 pipeline/d2g_signaling_export.py --stage $stage_args
    done
elif should_run d2gs; then
    echo "--- Phase D2gs: covered by d2g ---"
else
    echo "--- Phase D2gs: SKIPPED ---"
fi

# ── Phase D2i: PITX2 locus snRNA-seq read pileup ───────────────────────
# Pseudobulks per-base read depth at the PITX2 locus across the 9
# stage-restricted populations consumed by fig02 panel D.  Reads BAMs
# directly via samtools (hence the samtools module load above) and the
# population × library barcode manifest at
# figures/data/d2i_barcodes/manifest.tsv.  Strand-filters to PITX2-sense
# so antisense ncRNAs (LOC113843538/LOC113843539) don't contaminate.
if should_run d2i; then
    echo ""
    echo "--- Phase D2i: PITX2 locus RNA pileup ---"
    python3 pipeline/d2i_pitx2_rna_pileup.py
else
    echo "--- Phase D2i: SKIPPED ---"
fi

# ── Phase D2j: PITX2 isoform pileups (per-genotype + per-cell) ─────────
# Two passes share the cell-type annotation map; chained in one phase
# because they read the same h5ad (the 32G allocation in the resource
# guide above covers d2j_per_genotype's h5ad-into-memory peak).  Both
# apply the sense-of-mRNA strand filter (-F 16); a Tier 2 baseline
# (--no-strand) is available via FC_08656/run_d2j_baseline_nostrand.sh
# but is NOT part of the production pipeline.
if should_run d2j; then
    echo ""
    echo "--- Phase D2j: PITX2 isoform per-genotype ---"
    python3 pipeline/d2j_pitx2_isoform_per_genotype.py
    echo ""
    echo "--- Phase D2j: PITX2 isoform per-cell ---"
    python3 pipeline/d2j_pitx2_isoform_percell.py
else
    echo "--- Phase D2j: SKIPPED ---"
fi

# ── Phase D2k: PITX2c laterality QC ────────────────────────────────────
# Validates per-genotype L/R assignments against d2j PITX2c isoform
# enrichment.  Reads existing TSVs only (no BAM I/O); kept on cluster
# for sync uniformity rather than because cluster compute is needed.
# Exits non-zero if any FAIL_LR_FLIP verdicts surface (data-integrity
# fail-fast).
if should_run d2k; then
    echo ""
    echo "--- Phase D2k: PITX2c laterality QC ---"
    python3 pipeline/d2k_pitx2c_laterality_qc.py
else
    echo "--- Phase D2k: SKIPPED ---"
fi

# ── Phase E: Embryo-level ──────────────────────────────────────────────
if should_run e1; then
    echo ""
    echo "--- Phase E1: Pseudobulk PCA ---"
    python3 pipeline/e1_pseudobulk_pca.py
    echo ""
    echo "--- Phase E1b: Stage Inference ---"
    python3 pipeline/e1b_stage_inference.py
else
    echo "--- Phase E1: SKIPPED ---"
fi

# ── Phase H1: snRNA pseudobulk for L vs R DEG ──────────────────────────
# Pseudobulk variants by lineage scope: LPM (h1_pseudobulk_snrna),
# excluded (h1_excluded_pseudobulk), dorsal (h1_dorsal_pseudobulk).
# All read step04_preCB_clustered.h5ad + lr_inference_multiome.tsv +
# e1b_stage_assignments.tsv; running them in the same job (immediately
# after e1) ensures identical input snapshots (cross-job
# reuse of the same h5ad is avoided by design).
#
# All three h1 scripts also load ccreg d2g_lpm_subclustering_HH*_ccreg.tsv.gz +
# references/d2g_subcluster_annotations.tsv (from figures/data/) as the
# subcluster-population SoT.  When d2g is
# NOT in the run range (e.g. --from h1 --to h2), those sidecars are
# committed state and may be stale relative to a freshly rebuilt step04
# h5ad.  Fail-fast on missing sidecars; warn loudly if d2g was skipped.
if should_run h1; then
    for f in d2g_lpm_subclustering_HH13_ccreg.tsv.gz d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz \
             d2g_lpm_subclustering_HH21_ccreg.tsv.gz d2g_lpm_subclustering_HH25_ccreg.tsv.gz; do
        [[ -f "figures/data/$f" ]] || { echo "ERROR: h1 input missing: figures/data/$f (did d2g fail or get skipped?)"; exit 1; }
    done
    if ! should_run d2g; then
        echo ""
        echo "WARNING: h1 is using COMMITTED d2g sidecars (d2g not in --from..--to range)."
        echo "  If step04 h5ad was rebuilt since the committed d2g sidecars, h1 outputs"
        echo "  will combine fresh h5ad with stale subcluster labels — biologically wrong."
        echo "  Re-run with --from d2g (or earlier) to regenerate sidecars."
    fi
    echo ""
    echo "--- Phase H1: snRNA pseudobulk (LPM + excluded + dorsal) ---"
    python3 pipeline/h1_pseudobulk_snrna.py
    python3 pipeline/h1_excluded_pseudobulk.py
    python3 pipeline/h1_dorsal_pseudobulk.py
else
    echo "--- Phase H1: SKIPPED ---"
fi

# ── Phase H2: paired + unpaired+sex DESeq2 across the scopes ──────────
# paired model × {LPM, excluded, dorsal}, then
# unpaired+sex model × {LPM, excluded}.  Output naming
# matches the manifest's produced_by labels:
#   paired:   h2_deg_<scope>_*.tsv          (h2_deseq2_<scope>)
#   unpaired: h2_deg_unpaired_<scope>_*.tsv (h2_deseq2_unpaired_<scope>)
# Reads h1 outputs from figures/data/ within the same job — chained
# with H1 to avoid cross-job h5ad reuse race (per
# feedback_same_job_chain memory).  DESeq2 is supplied by ~/R/4.4.2/
# (see pipeline/install_o2_r_packages.R for the install record).
if should_run h2; then
    # H3: DESeq2 must be installed in ~/R/4.4.2/.
    # Without this guard, a fresh O2 setup that skipped
    # run_install_o2_r_packages.sh wastes cluster time before
    # h2_deseq2_snrna.R fails with "could not find function".
    Rscript -e '.libPaths(c("~/R/4.4.2", .libPaths())); if (!requireNamespace("DESeq2", quietly=TRUE)) stop("DESeq2 missing in ~/R/4.4.2/. Run: sbatch pipeline/run_install_o2_r_packages.sh", call.=FALSE)'

    # assert h1 produced every expected output.
    # A `--from h2 --to h2` run reuses h1 outputs from a prior job; if
    # h1 silently produced an empty/corrupted file, h2 would fail with
    # a cryptic missing-file error that doesn't identify the scope.
    for f in h1_pseudobulk_snrna_counts.tsv.gz h1_pseudobulk_snrna_metadata.tsv \
             h1_excluded_pseudobulk_counts.tsv.gz h1_excluded_pseudobulk_metadata.tsv \
             h1_dorsal_pseudobulk_counts.tsv.gz h1_dorsal_pseudobulk_metadata.tsv; do
        [[ -f "figures/data/$f" ]] || { echo "ERROR: h2 input missing: figures/data/$f (did h1 fail?)"; exit 1; }
    done

    STAGES="$CANONICAL_STAGES_CSV"
    echo ""
    echo "--- Phase H2: DESeq2 L vs R (paired model: LPM, excluded, dorsal) ---"
    Rscript pipeline/h2_deseq2_snrna.R \
        --stages "$STAGES" \
        --counts-name h1_pseudobulk_snrna_counts.tsv.gz \
        --meta-name   h1_pseudobulk_snrna_metadata.tsv \
        --out-prefix  h2_deg
    Rscript pipeline/h2_deseq2_snrna.R \
        --stages "$STAGES" \
        --counts-name h1_excluded_pseudobulk_counts.tsv.gz \
        --meta-name   h1_excluded_pseudobulk_metadata.tsv \
        --out-prefix  h2_deg_excluded
    # Dorsal = bilateral control (fig03 panel B).  No HH13 dorsal population, so
    # h2 skips HH13 for this scope (no samples).
    Rscript pipeline/h2_deseq2_snrna.R \
        --stages "$STAGES" \
        --counts-name h1_dorsal_pseudobulk_counts.tsv.gz \
        --meta-name   h1_dorsal_pseudobulk_metadata.tsv \
        --out-prefix  h2_dorsal

    echo ""
    echo "--- Phase H2: DESeq2 L vs R (unpaired+sex model: LPM, excluded) ---"
    Rscript pipeline/h2_deseq2_snrna.R \
        --stages "$STAGES" --unpaired --control-sex \
        --counts-name h1_pseudobulk_snrna_counts.tsv.gz \
        --meta-name   h1_pseudobulk_snrna_metadata.tsv \
        --out-prefix  h2_deg_unpaired
    Rscript pipeline/h2_deseq2_snrna.R \
        --stages "$STAGES" --unpaired --control-sex \
        --counts-name h1_excluded_pseudobulk_counts.tsv.gz \
        --meta-name   h1_excluded_pseudobulk_metadata.tsv \
        --out-prefix  h2_deg_unpaired_excluded
else
    echo "--- Phase H2: SKIPPED ---"
fi

# ── Phase H3: per-nucleus gene detection-fraction export ──────────────
# Tiny producer of NODAL / CFC1 / PITX2 per-nucleus detection fractions (
# all-nuclei + inference-pool scopes). Reads the same step02 substrate +
# d2g sidecars as h1; safe standalone (--from h3 --to h3).
if should_run h3; then
    for f in d2g_lpm_subclustering_HH13_ccreg.tsv.gz d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz \
             d2g_lpm_subclustering_HH21_ccreg.tsv.gz d2g_lpm_subclustering_HH25_ccreg.tsv.gz; do
        [[ -f "figures/data/$f" ]] || { echo "ERROR: h3 input missing: figures/data/$f (did d2g fail or get skipped?)"; exit 1; }
    done
    echo ""
    echo "--- Phase H3: gene detection-fraction export ---"
    python3 pipeline/h3_gene_detection_export.py
else
    echo "--- Phase H3: SKIPPED ---"
fi

# ── Figures ────────────────────────────────────────────────────────────
if should_run fig; then
    echo ""
    echo "--- Dotplot + composition export ---"
    python3 FC_08656/20260330_generate_dotplot.py

    echo ""
    echo "--- Export figure data TSVs to repo ---"
    FIGDATA="figures/data"
    cp -v "$PHASE2/figures/step03_sex_per_genotype.tsv"  "$FIGDATA/"
    cp -v "$PHASE2/figures/step03_sex_per_cell.tsv.gz"   "$FIGDATA/"
    cp -v "$PHASE2/figures/e1b_stage_assignments.tsv"    "$FIGDATA/"
    # The UMAP panel reads ONLY the declared committed-path input; refresh it from the step04a
    # file d2_cluster_markers just wrote (manifest row), else a d2 rerun would render stale
    # cells/coords against fresh annotations. No `|| true`: set -e aborts if step04a is missing.
    cp -v "$PHASE2/figures/step04a_umap_clusters.tsv.gz" "$FIGDATA/d2_umap_boundary_data.tsv.gz"
    # Update PIPELINE_STATE.json with current counts
    N_GENO=$(tail -n +2 "$FIGDATA/step03_sex_per_genotype.tsv" | wc -l | tr -d ' ')
    N_CELLS=$(awk -F'\t' 'NR>1 {s+=$2} END {print s}' "$FIGDATA/step03_sex_per_genotype.tsv")
    COMMIT=$(git -C . rev-parse --short HEAD 2>/dev/null || echo unknown)
    GENDATE=$(date -I)
    python3 -c "
import json, os
p = '$FIGDATA/PIPELINE_STATE.json'
s = json.load(open(p)) if os.path.exists(p) else {}
s['n_genotypes'] = int('$N_GENO')
s['n_cells'] = int('$N_CELLS')
s['generated_date'] = '$GENDATE'
s['pipeline_commit'] = '$COMMIT'
s['stale_warning'] = None
s['known_inconsistencies'] = []
json.dump(s, open(p, 'w'), indent=2)
print(f'  PIPELINE_STATE: {s[\"n_genotypes\"]} genotypes, {s[\"n_cells\"]} cells, {s[\"generated_date\"]}')
"

    echo ""
    echo "--- Figures (lightweight: TSV-based panels first) ---"
    # Every generator is attempted; failures are COLLECTED and fail the run at the
    # end (never a success banner over missing panels). ALLOW_FIGURE_FAILURES=1
    # downgrades them to a loud warning, e.g. when only the derived TSVs matter.
    fig_fail() { echo "  ERROR: $1 failed"; FIG_FAILS+=("$1"); }
    # UMAP panel reads figures/data/d2_umap_boundary_data.tsv.gz (refreshed from step04a above)
    python3 figures/fig01/generate_umap_res02.py || fig_fail "UMAP panel"
    # R figures (need R/4.4.2 loaded in SBATCH header)
    Rscript figures/fig02/generate_paired_laterality.R || fig_fail "paired laterality"
    Rscript figures/fig02/generate_paired_laterality_isoform.R || fig_fail "paired isoform"
    Rscript figures/fig02/generate_pitx2_pileup_multistage.R || fig_fail "pileup"
    # Compositors intentionally omitted — no active compositor after the
    # figure reorganization; assemble PNGs manually when ready to publish.

    echo ""
    echo "--- Figures (heavyweight: h5ad or data-dependent) ---"
    python3 figures/fig01_supp/generate_libqc.py || fig_fail "libqc"
    python3 figures/fig01_supp/generate_genotype_qc.py || fig_fail "genotype_qc"
    python3 figures/fig01_supp/generate_clustering.py || fig_fail "clustering"
    # SVG diagrams (generate_pca_staging.py, generate_spatial_svg.py) run LOCALLY — they
    # need svgwrite/cairosvg not available in vireo env. Skip on cluster.
else
    echo "--- Figures: SKIPPED ---"
fi

if [ "${#FIG_FAILS[@]}" -gt 0 ]; then
    echo ""
    echo "=== FIGURE STAGE FAILED: ${#FIG_FAILS[@]} generator(s): ${FIG_FAILS[*]} ==="
    if [ "${ALLOW_FIGURE_FAILURES:-0}" != "1" ]; then
        echo "(derived TSVs above were written; set ALLOW_FIGURE_FAILURES=1 to downgrade)"
        exit 1
    fi
    echo "ALLOW_FIGURE_FAILURES=1 — continuing with missing panels"
fi

echo ""
echo "=== Pipeline Complete (from: $START, to: $END) ==="
echo "Finished: $(date)"
