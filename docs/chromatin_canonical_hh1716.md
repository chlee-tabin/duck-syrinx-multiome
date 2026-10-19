# Chromatin-canonical L/R inference at HH17/16

Methods note: how the left/right side of each HH17/16 genotype is called when the snRNA and chromatin readouts disagree.

---

## 1. The claim in one paragraph

At HH17/16, snRNA-seq PITX2c expression is at the limit of detection — the median per-genotype PITX2c-positive cell fraction in the splanchnic-LPM inference cohort is exactly 0%. snATAC-seq chromatin accessibility at the PITX2c promoter does not have this floor: the same per-(genotype, library) cohort that produces 0% PITX2c expression yields a measurable, asymmetric chromatin signal. We adopt **chromatin as the canonical L/R caller at HH17/16** when the two modalities disagree, justified by (i) chromatin's higher sensitivity at the inter-wave transcriptional gap and (ii) chromatin's faithful corroboration of snRNA-based L/R calls at the three flanking stages (HH13, HH21, HH25, paired Wilcoxon p=0.035 at HH13 and p=0.014 each at HH21/HH25). At HH17/16 the chromatin paired Wilcoxon remains significant (paired p=0.013, n=9 paired embryos; unpaired p=7.7e-5 at n_L=19, n_R=14 — `figures/data/step_f5_pitx2_region_stats.tsv` after the 2026-10-14 lane-loading-record correction; it was paired p=0.00166 at n=13, unpaired p=2.08e-5 at n_L=n_R=18 before it), consistent with chromatin retaining the laterality signal that NODAL-driven Wave 1 transcription deposited and through-which BMP-driven Wave 2 transcription will later reactivate.

This rule underlies fig03 (PITX2 chromatin memory) and the figS03 L/R-resolution matrix.

---

## 2. Empirical premises (with verifiable numbers)

### 2a. snRNA PITX2c is at the noise floor at HH17/16

Per-genotype PITX2c-positive cell fraction in the d2l splanchnic-LPM inference cohort, by stage:

| stage | mean | median | n_obs | non-direct firmness rate |
|---|---|---|---|---|
| HH13 | 1.45% | 0.50% | 16 | 25.0% (4/16) |
| **HH17/16** | **0.56%** | **0.00%** | **42** | **19.0% (8/42)** |
| HH21 | 1.07% | 0.46% | 20 | 20.0% (4/20) |
| HH25 | 1.24% | 0.63% | 16 | 0.0% (0/16) |

HH17/16 is the only stage where the median is exactly 0 — at least half the genotype × library rows have *no* PITX2c-positive cells in their inference cohort — and its mean PITX2c+ fraction is about half that of the flanking stages. It is also the only stage with an unresolved (`weak`) d2l call. This is the manuscript's "snRNA blind spot" concretised.

Source: `figures/data/d2l_lr_inference.tsv`, columns `pitx2c_pos`, `firmness`, `design_stage`.

### 2b. Chromatin faithfully corroborates snRNA at flanking stages

Stage-by-stage chromatin paired Wilcoxon signed-rank L>R at PITX2c promoter:

| stage | snRNA L>R (population) | chromatin L>R paired Wilcoxon | d2l vs d2n disagreements |
|---|---|---|---|
| HH13 | strong (Wave 1, NODAL-driven) | p=0.035 (n=7) | 0 |
| HH17/16 | weak / absent | **paired p=0.013 (n=9); unpaired p=7.7e-5 (n_L=19, n_R=14)** | 4 of 38 independently inferred rows (g23 and g24, each × GEX07/GEX08) = 34/38 agree (89.5%), plus 4 design-fixed same-side-duplicate rows (G3_unfiltered_donor0, G3_unfiltered_donor3 × {GEX04, GEX0S-duck}) that are not a comparison; see the note below |
| HH21 | strong (Wave 2, BMP-driven) | p=0.014 (n=8) | 0 |
| HH25 | strong (Wave 2 maintenance) | p=0.014 (n=8) | 0 |

**Note.** The HH17/16 counts above are recounted from `figures/data/validation_d2n_vs_d2l_hh1716.tsv` (`transition_kind`: 34 agree, 4 flip, 4 same_side_dup). G3_unfiltered_donor0 and G3_unfiltered_donor3 are each one left half sequenced in two libraries (GEX04 and GEX0S-duck, per the lane-loading record), so both callers carry the design side for those rows and they are not an independent comparison. d2l (direct) and d2n (placement) both put G3_unfiltered_donor5 at GEX04 R / GEX0S-duck L, and both call G3_unfiltered_donor8 GEX04 L (singleton rank).

Disagreement is concentrated entirely at HH17/16 — exactly where snRNA loses signal and chromatin retains it. At all three stages where snRNA gives strong PITX2c L/R signal, chromatin agrees on every (genotype, library) pair. This pattern is the basis for trusting chromatin at HH17/16 specifically.

Source: per-stage statistics in `figures/data/step_f5_pitx2_region_stats.tsv` (`PITX2c_promoter` rows); the disagreement counts compare `figures/data/d2l_lr_inference.tsv` with `figures/data/lr_inference_multiome.tsv`.

### 2c. Chromatin retains a per-genotype signal at HH17/16

The MIN_CELLS=30 per-genotype HH17/16 panel (`figures/exploratory/f4_lr_pitx2_pileup/f4_pitx2_lr_pileup_per_genotype.png` *(not distributed)*) shows 12 paired + 7 singleton chromatin tracks at PITX2c. Asymmetry is visible per-genotype across all 12 paired tracks. The g52 case study confirmed that even unusual high-baseline embryos preserve L>R chromatin at the PITX2c TSS; chromatin asymmetry is a per-embryo phenomenon, not a population-marginal artifact.

---

## 3. The decision tree, walked

### 3a. Scope of "what chromatin resolves"

Three categories of HH17/16 ambiguity could in principle benefit from a second modality:

| category | description | chromatin's leverage |
|---|---|---|
| (α) Chromatin memory hypothesis | population-level chromatin-memory claim | **only chromatin can test it** |
| (β) Per-genotype L/R sign at HH17/16 | snRNA call quality is low; chromatin can corroborate | **direct contribution** |
| (γ) Single-library asymmetric multi-candidate (e.g., G2_genotype_11) | design has multiple candidate embryos for one library | none — chromatin shares the library coordinate |
| (δ) Symmetric peer-group multi-candidate (interchangeable) | design symmetric between e.g. {17-5, 17-6} | none — same library coordinate |

We adopt chromatin for (α) + (β). (γ) and (δ) remain irreducible from any single (library, side) measurement.

### 3b. Is chromatin a tie-breaker, or canonical?

The disagreement decision: at HH17/16, when d2l (snRNA) and d2n (chromatin) disagree on a (genotype, library) L/R sign, is chromatin's call **canonical** (always wins) or a **tie-breaker** (wins only when snRNA is internally inconsistent, e.g., pair_constraint=False)?

We adopt **canonical at HH17/16**. The justification rests on §2: snRNA expression is at noise floor at this stage (median 0% PITX2c+), so a snRNA "direct firmness" call comparing 0 vs 1 cell is methodologically weaker than a chromatin "placement firmness" call comparing accessibility profiles across hundreds of cells. We do not extend this to other stages — at HH13/HH21/HH25 snRNA expression is well above background and the canonical-snRNA framework applies.

The 4 (genotype, library) HH17/16 rows where the canonical side differs from the d2l call:

| genotype × library | d2l side | d2l firmness | chromatin side | chromatin firmness | change |
|---|---|---|---|---|---|
| g23 GEX07 | R | direct | L | direct | flip |
| g23 GEX08 | L | direct | R | direct | flip |
| g24 GEX07 | R | direct | L | direct | flip |
| g24 GEX08 | L | direct | R | direct | flip |

Every flipped row has a d2l PITX2c-positive fraction of 0 or near 0 (g23: 0 of 32 and 0 of 88 nuclei; g24: 0 of 59 nuclei in GEX07 vs 6 of 345 in GEX08). These "direct" calls rest on the total-PITX2 metrics with PITX2c itself at the floor, which is the kind of snRNA call the chromatin-canonical rule is designed to override.

### 3c. Block-level vs per-genotype `pair_constraint_ok`

**Scope note.** §3c–§5b walk the calls on the design as recorded before the lane-loading-record correction, when 17-5/17-6 were listed as an L/R pair across GEX04 and GEX0S-duck. On the corrected design, G3_unfiltered_donor0 and G3_unfiltered_donor3 are same-side duplicates with a design-fixed L side; d2l and d2n both call G3_unfiltered_donor8 GEX04 L (singleton rank; uniquely 16-4); and every HH17/16 G3 genotype reconciles with a design candidate (`pc_fail` False in `figures/data/d2m_multiome_genotype_to_embryo.tsv`), while the snRNA-only rollup still flags the GEX04 × HH17/16 block (`pair_constraint_ok` False in `figures/data/d2l_lr_inference.tsv`).

`d2l` defines `pair_constraint_ok` per (library × design_stage) block: the count of d2l-resolved L's and R's must equal the design's expected counts. This is a **conservative block-level rollup** that propagates a single weak/missing call across all genotypes in the block. Concretely: if any genotype in a HH17/16 block is `?` (weak), the block is `pair_constraint_ok=False`, and ALL genotypes in the block inherit that flag in `d2m`'s `_pc_fail` rollup.

For figS02D (snRNA-only), this conservative rollup is correct: with only one modality, you cannot tell whether a block's failure is genotype-specific or pervasive. The red rows honestly mean "this block is broken; without a second modality we don't know exactly who."

For the figS03 L/R-resolution matrix (chromatin-arbitrated), we adopt a **genotype-specific** rule: a genotype is flagged only if its OWN sides cannot be reconciled with at least one design-candidate embryo, ignoring block-mates' weak calls. After chromatin arbitration, donor8 remains the sole unresolvable case — its single GEX04-only fingerprint with `?` side has no chromatin partner to compare against. Under the genotype-specific rule, donor0/3/5/9 (whose individual sides DO match design candidates after chromatin arbitration) are no longer dragged into red by donor8's failure.

The two rules show the **methodological progression**: figS02 panel D shows the snRNA-only state with conservative block rollup (the problem statement); the figS03 L/R-resolution matrix shows the chromatin-arbitrated state with per-genotype check (the resolution).

**Implementation note — convergence-loop placement.** `pipeline/d2l_lr_inference.py::apply_placement` iterates to convergence: placement is re-run until no new sides are filled, so a partner update at one library (e.g., GEX08 placing g20=L) can unlock placement by elimination at a sibling block (e.g., GEX03-1's G2_genotype_03). With this, the GEX03-1 and GEX03-2 HH17/16 blocks reconcile via cascade (G2_genotype_03 → R/L placement; G2_genotype_11 → L placement matching design), and the snRNA-only conservative rollup isolates the residual failure to the GEX04 × HH17/16 d2l block (1L/3R inferred + 1 weak donor8): 5 recovered G3 genotypes for 6 design embryos at HH17/16, with one HH16-unpaired genotype undetected in demultiplexing. The four paired G3 donors (donor0/3/5/9) inherit pc_fail via the conservative rollup from the GEX04 d2l block. d2n's `apply_placement` has the same convergence-loop structure. This is the irreducible single-modality + design-count limit; chromatin's per-genotype rule reduces it to 1 (donor8 alone) per §4.

**Implementation note — 1-of-many singleton rule.** `pipeline/d2l_lr_inference.py::resolve_singletons` and `d2n_chromatin_lr_inference_hh1716.py::resolve_singletons` rank singletons when the singletons in a (library × design_stage) block exactly fill the unfilled design slots (`n_L_remaining + n_R_remaining == n_singletons`). A lone singleton facing both an L and an R slot (donor8 at GEX04 HH17/16: 1L + 1R unfilled after paired-direct calls) cannot be split that way, so a 1-of-many branch applies: if the singleton's directional metric is unambiguously in one tail of the within-library firm-direct cohort distribution (above max(firm_R) for L, below min(firm_L) for R), the singleton is assigned firmly; otherwise the rule abstains. d2n uses the single metric VOTE_METRIC (`pitx2c_frags_per_cell`); d2l uses a majority vote (≥2/3 of `pitx2_pos`/`pitx2c_pos`/`pitx2_mean`). Low-n firm rows are excluded from the comparison cohort. Donor8's chromatin reading exceeds every firm-R reading at GEX04 HH17/16 (pitx2c_frags=0.278 > max(donor3=0.171, donor9=0.020)), so the chromatin call resolves donor8 → L firm; d2l's snRNA call for donor8 stays weak (`?`). See §5b for the cohort-internal rationale.

**Implementation note — cohort-indirect override sidecar (inactive).** A hand-curated `references/genotype_overrides.tsv` can carry editorial soft-calls (schema: `genotype, override_type, design_stage, inferred_side, justification`). `pipeline/d2o_merge_multiome_inference.py::apply_overrides` stamps each override row into a separate `cohort_indirect_side` column on `lr_inference_multiome.tsv` *without* touching `inferred_side`; the renderer's `pivot_by_genotype` translates the column into per-slot `_cohort_indirect` flags that paint the panel as orange `=L?` soft-calls. The mechanism is for cases where firm-tier inference cannot fire (e.g., two genotypes with cardinality 2 and no within-library firm-direct cohort to compare against). `apply_overrides` raises if any override row matches no candidate rows. The shipped TSV is header-only (no active overrides).

### 3d. Full-library-set match

`d2m`'s candidate-matching logic checks (library, side) per row but does not enforce that a candidate embryo's FULL library set match the genotype's library set. This causes a single-library genotype like donor8 (cells in GEX04 only) to inherit candidates that require both GEX04 and GEX0S-duck cells (e.g., 17-5 has GEX04 L + GEX0S-duck R per design — donor8 has no GEX0S-duck cells and cannot be 17-5).

The rule: enforce `set(genotype_libraries) == set(candidate_libraries)`. This narrows donor8 from 6 candidates {16-4, 16-5, 17-5, 17-6, 17-7, 17-8} to 2 {16-4, 16-5} — both of which are explicitly recorded in `design/embryo_lr_design.tsv` as "G3; unpaired (GEX04 only)" HH16 single-half embryos. The narrowing is pure design enforcement (no chromatin needed), but compounds with chromatin arbitration.

Applied uniformly to BOTH `d2m_genotype_to_embryo.py` (figS02D) and `d2m_multiome_genotype_to_embryo.py` (the figS03 L/R-resolution matrix). The narrowing alone leaves donor8 a 2-candidate case; §5b resolves it.

---

## 4. Net effect — figS02 panel D vs the figS03 L/R-resolution matrix

The same matrix layout, two data sources, two rules. The difference IS the figure's argument.

| genotype | figS02 panel D (snRNA-only, block-rollup) | figS03 matrix (chromatin-arb. + library-set match + genotype-specific) | mechanism |
|---|---|---|---|
| G3_donor0 | red (pc_fail block) | peer-asterisk with donor5 → {17-5, 17-6} | chromatin agrees, peer group forms |
| G3_donor3 | red (pc_fail block) | peer-asterisk with donor9 → {17-7, 17-8} | chromatin agrees, peer group forms |
| G3_donor5 | red (pc_fail block) | peer-asterisk with donor0 (after chrom. flip) | chromatin arbitrates the snRNA call (d2l R/L on 4 GEX04 nuclei; chromatin L/R) |
| G3_donor9 | red (pc_fail block) | peer-asterisk with donor3 | chromatin agrees, peer group forms |
| G3_donor8 | red (pc_fail; cardinality 2 after §3d) | **firm L** → uniquely 16-4 | 1-of-many singleton rule (§3c implementation note) — donor8's chromatin reading exceeds every firm-R reference at GEX04 HH17/16 (d2l leaves it weak). Narrows cardinality 2 → 1 |
| G2_genotype_11 | uniquely 17-KC-1 (post §3d) | uniquely 17-KC-1 | strict library-set match → only KC HH17/16 embryo with `lib_m={GEX03-2}` (§3d's literal rule applied symmetrically) |

**figS02 panel D**: 5 reds + 0 oranges (the snRNA-only problem statement
under §3d's library-set match).

**figS03 matrix**: 0 reds + 0 oranges + 4 peer-asterisks (the chromatin
+ library-set match + 1-of-many singleton resolution).  donor8 is firm-resolved
to 16-4 (L) via the 1-of-many singleton rule (§3c implementation note);
no residual oranges remain.

---

## 5. Design-intrinsic edge-case resolutions (no residual oranges)

After all refinements (§3a chromatin-canonical, §3b conservative pc_fail rollup, §3c convergence-loop placement + 1-of-many singleton rule, §3d full-library-set match), **no genotypes remain orange in the figS03 matrix**.  Two genotypes that would otherwise be orange — G2_genotype_11 and G3_donor8 — are both firm-resolved.  They are described below because both illustrate **design-intrinsic ambiguity classes**: any single-library single-half genotype can take the same shape, and the resolution mechanism (§3d for G2_11; the 1-of-many singleton rule for donor8) is what a reader needs to follow the calls.

### 5a. G2_genotype_11 (resolved by §3d)

Single library (GEX03-2), single side (L), 3 KC HH17 design candidates pre-§3d: {16-KC-1, 17-KC-1, 17-Md-14}. No paired library and no symmetric peer genotype with the same fingerprint. Chromatin shares the (library, side) coordinate so cannot break the embryo-identity ambiguity. Disambiguation via design-side library-set match (§3d) uniquely picks 17-KC-1 (the only HH17/16 KC embryo with the observed single-library fingerprint), so the matrix shows G2_genotype_11 as **uniquely resolved (no orange tint)**.  The *kind* of design ambiguity it illustrates is general: it applies to any single-library asymmetric multi-candidate genotype that is not uniquely resolved by library-set match.

### 5b. G3_donor8 (resolved by the 1-of-many singleton rule)

Single library (GEX04), 2 HH16 design candidates: {16-4 (L), 16-5 (R)}.  The
wetlab design records 16-4 and 16-5 as **separate embryos** of which
only one half each was retained (16-4-L, 16-5-R) — no L/R partner for
either embryo exists in any library.  Without a within-embryo partner,
no within-embryo chromatin contrast can establish "this is 16-4-L" vs
"this is 16-5-R" at the paired-direct firmness tier.  The exact-fill
firm-tier `singleton_rank` branch also cannot fire because it
requires `n_singletons` to exactly fill the unfilled (L+R) slots —
GEX04 had 1 singleton (donor8) but 2 unfilled slots (1L + 1R after
direct/placement claimed donor0+donor5=L and donor3+donor9=R), so the
within-singleton rank logic couldn't split 1 thing into 2 groups.

**Resolution via the 1-of-many singleton rule.**  The 1-of-many
firm-tier branch (§3c implementation note) handles this case by
comparing the lone singleton's directional metric to the within-library
firm-direct cohort distribution.  For donor8 in GEX04 HH17/16:

```
within-library firm-direct cohort (low_n donor5 excluded):
                              donor0 (L)  donor3 (R)  donor9 (R)  donor8 (?)
  pitx2c_frags/cell (chrom.)  0.113       0.171       0.020       0.278
  pitx2c_open_frac (chrom.)   0.085       0.155       0.020       0.185
  pitx2_pos (snRNA)           0.183       0.203       0.000       0.212
  pitx2_mean (snRNA)          0.332       0.338       0.000       0.350
  pitx2c_pos (snRNA)          0.000       0.000       0.000       0.000
```

donor8's chromatin VOTE_METRIC (pitx2c_frags/cell = 0.278) exceeds
**every firm-R reference** in this library (max(firm_R) = 0.171),
so the d2n rule assigns chromatin firm L.  On snRNA alone d2l leaves
donor8 weak (`?`): pitx2c_pos is at the noise floor for every genotype in
this block.  The cardinality narrows from 2 → 1 ({16-4}); design_conflict drops
to False; the panel renders firm L (no orange tint).

**Why this is the principled resolution, not a relaxation.**  The
"high pitx2c → L" assumption is the same directional bias the
multi-singleton rank logic already encodes (top of metric rank → L).
Within a single library, that within-library firm cohort IS the
controlled reference for the directional bias — same chemistry,
same batch, same per-cell normalization.  The 1-of-many rule formalizes
that intuition with a strict tail-test (above firm-R max → L; below
firm-L min → R); singletons that fall between the two firm tails
stay weak (the rule abstains when evidence is ambiguous).  donor8
clears the chromatin test.

---

## 6. Where the rules are implemented

| rule | file |
|---|---|
| Full-library-set match | `pipeline/d2m_genotype_to_embryo.py` (shared narrow-genotype function, used by d2m and d2m_multiome) |
| Genotype-specific pc_fail mode | `pipeline/d2m_genotype_to_embryo.py` (flag), `pipeline/d2m_multiome_genotype_to_embryo.py` (consumes it); block rollup is the default |
| 3-tier colour rule (red = pc_fail; orange = design ambiguity; asterisk = peer group) | `figures/fig02_supp/generate_lr_inference.py`, `figures/_shared/lr_resolution_panel.py` |
| Chromatin-arbitrated L/R-resolution matrix | `figures/figS03/generate_lr_resolution.py` (consumes d2m_multiome) |
| Convergence-loop placement | `pipeline/d2l_lr_inference.py::apply_placement`, `pipeline/d2n_chromatin_lr_inference_hh1716.py::apply_placement` |
| 1-of-many singleton rule | `pipeline/d2l_lr_inference.py::resolve_singletons`, `pipeline/d2n_chromatin_lr_inference_hh1716.py::resolve_singletons` (d2n single metric VOTE_METRIC; d2l 2/3 majority across METRICS) |
| Cohort-indirect override sidecar (inactive) | `references/genotype_overrides.tsv`, `pipeline/d2o_merge_multiome_inference.py::load_overrides` + `apply_overrides` |
| Input schema checks | `figures/_shared/lr_resolution_panel.py::assert_lr_input_schemas` + the `pivot_by_genotype` row-set check (raise if d2m lacks `pc_fail`/`design_conflict` or any d2l/multiome genotype is missing from d2m) |

---

## 7. Why this matters scientifically

The chromatin-canonical-at-HH17/16 rule is not a methodological convenience — it is a **direct test of the chromatin memory hypothesis**. The manuscript's central claim is that PITX2c chromatin asymmetry persists across the inter-wave transcriptional gap (HH17/16) when expression itself is silent. If snRNA disagreed with chromatin at HH17/16 *and we trusted snRNA*, we would be discarding the very signal the hypothesis predicts. By accepting chromatin's call where modalities disagree at HH17/16, we are operationalizing the chromatin memory hypothesis in the L/R inference pipeline itself: the data products downstream of d2o (chromatin-arbitrated L/R) reflect the world the hypothesis predicts.

The faithful corroboration at HH13/HH21/HH25 (zero disagreements; paired Wilcoxon p=0.035 at HH13, p=0.014 each at HH21/HH25) provides the external validation: chromatin's disagreement with snRNA at HH17/16 is not a bug of either modality; it is the predicted consequence of biological state — transcription dampened, chromatin retained.
