# Chick syrinx control — genotype demux pipeline (Group A: `ARC01-02-0S-chick`)

Adapts the settled duck pipeline (`pipeline/a2_*`, `a3_sub_vireo_g3.sh`,
`a4_vireo_seed_qc.py`) to the chick syrinx control. Experimental design: 12 WT
embryos (`chick-13/15/25 ×4`) in Group A, 12 embryos (6 GFP⁺ + 6 WT) in Group B.

**Inputs (galGal6a / `crarc_galgal6a`; deposited in GEO GSE336327):**
`ARC01` (2,794 cells), `ARC02` (3,579), `ARC0S-chick` (2,956) — multiome;
each `outs/` has `gex_possorted_bam.bam` + `filtered_feature_bc_matrix/`.
`ARC01`/`ARC02` are the **same 12 embryos** (cross-lane concordance check);
`ARC0S-chick` is L-only HH15+HH25 (adds depth for 8 left genotypes).

**Env:** `conda activate vireo` (has `cellsnp-lite`, `vireo`, `bgzip`, `tabix`).
**Workspace:** `$DUCK_SCRATCH_USER/chick_demux/` (Group A) and `$DUCK_SCRATCH_USER/chick_demux_B/`
(Group B); see `docs/audit_hardcoded_paths.md` for the path roots.

## Steps
| Step | Script | Where | What |
|---|---|---|---|
| c0 | `c0_stage_inputs.sh` | transfer/login node | stage GEX BAMs + barcodes + reference into the workspace; build per-lib barcodes + chroms |
| c1 | `c1_cellsnp_discovery.sh` | SLURM | cellSNP **mode-2b bulk** discovery on the 3 GEX BAMs → shared variant base |
| c2 | `c2_cellsnp_genotype.sh` | SLURM array 0-2 | cellSNP **mode-1** per library at the discovered sites (CellRanger barcodes) |
| c2b | `c2b_merge_cellsnp.py` | login/compute | merge 3 per-lib matrices (rows align: same base) → combined cellData |
| c3 | `c3_vireo_multiseed.sh` | SLURM array | 20-seed vireo, swept over N∈{10,12,14,16} (expected **12**) |
| c4 | `c4_seed_qc.py` | login/compute | logLik rank + stability + **cross-lane (ARC01 vs ARC02) concordance** |
| c5* | `c5_validate.py`, `c5c_sex_wgene.py`, `c5b_gfp.py`, `c5d`–`c5h` | login/compute | post-demux validation gates: W/Z-gene sex, GFP⁺ label, over-split / split-concordance / merge-confirm |
| c6b | `c6b_assign_embryo.py` | login/compute | genotype → embryo assignment (library membership + stage + design intersection) |
| cAE | `cAE_full_atlas.py` (`cAE_atlas.sh`) | SLURM | whole-embryo A+B atlas (HVG→PCA→Harmony→Leiden→d2b-marker argmax); CellBender-corrected redraw |
| cAL | `cAL_stage_subclustering.py` | SLURM | LPM-lineage, stage-restricted subclustering (chick `d2g` analog) + UEG/Wilcoxon/lineage-score triad |
| cAM | `cAM_compendium_annotate.py` | login/compute | compendium annotation of the cAL subclusters → `cAL_subcluster_annotations.tsv` SSoT + per-cell populations |
| cAN | `cAN_stage_laterality.py` | local | PITX2 laterality fold: orientation-free Panel F + design-oriented directional control |

(Group B = `GEX09/GEX10` in the separate `chick_demux_B` workspace, with the `ARC0S` left-half lane;
staged + demuxed by `c0b_stage_groupB.sh` and the `g0s*` scripts following the same c1–c4 recipe.)

## Run
```bash
# c0 on a node that sees the raw data (writes only to the workspace)
bash pipeline/chick/c0_stage_inputs.sh
# c1 discovery
sbatch pipeline/chick/c1_cellsnp_discovery.sh           # -> WORK/cellsnp/discovery/
# c2 per-library genotyping (after c1)
sbatch --array=0-2 pipeline/chick/c2_cellsnp_genotype.sh # -> WORK/cellsnp/<LIB>/
# c2b merge (after c2)
conda run -n vireo python3 pipeline/chick/c2b_merge_cellsnp.py
# c3 multi-seed vireo, one submit per N (after c2b)
for N in 10 12 14 16; do sbatch --array=0-19 --export=ALL,N_DONOR=$N pipeline/chick/c3_vireo_multiseed.sh; done
# c4 seed QC + cross-lane concordance (after c3)
conda run -n vireo python3 pipeline/chick/c4_seed_qc.py
```

## Canonical path vs one-shots
The **canonical** figure-feeding chain is `cAE → cAL → cAM → cAN` (laterality, figS04 Panel F), plus
`cAE → cAK` / `cb_libqc_export` / `cA_crosslane_3lane` (figS04-supp QC panels). It has upstream
**prerequisites** that produce the staging/sex/split inputs the chain consumes — **do NOT skip
these**: the demux front-half (`c0`–`c4`), the sex/validation gates (`c5*`, which write the per-donor
sex tables), Group-A staging (`cT`), embryo assignment (`c6b`), and the targeted Group-B donor
splits (`cS`). Concretely: `cAE` reads `cT_groupA_stage_sex.tsv` (`cT`) + `c5h_corrected_donor_table.tsv`
(`c5h`) for stage/sex, and `cAL` applies the `cS` splits via `resolve_donor` — omitting any of these
makes the run fail. The per-step lineage is in `figures/data/README.md` (chick rows) and in
`figures/figS04/PROVENANCE.md` / `figures/figS04_supp/PROVENANCE.md`.

**Other figure-feeding steps** (not read by `cAE`/`cAL`/`cAM`/`cAN`): `cAD` (Group B PITX2
pileup, figS04 panel E), `cAF` (duck-panel dotplot), `cAH` (PITX2 3′UTR recovery), `cAJ`
(LPM pseudobulk PCA), `cA_crosslane_corr`, `cA_sex_scan`, `cb_floor_check` /
`cb_lane_floor_diag` (QC).
