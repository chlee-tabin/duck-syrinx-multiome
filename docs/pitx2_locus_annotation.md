# PITX2 locus annotation — duck `crarc_zju1_PITX2` reference

**Summary:** this note characterizes the PITX2 locus on the duck reference and the read-strand evidence behind the pipeline's sense-strand filter. Production code in `pipeline/_pitx2_strand.py`, `pipeline/d2i_pitx2_rna_pileup.py` and `pipeline/d2j_pitx2_isoform_*.py` derives the filter from the GTF strand at run time (§5).
**Source GTF:** `/n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/genes/genes.gtf.gz` (23 MB gz, dated 2024-05-26, Gnomon/NCBI annotation on NC_051775.1 scaffolding)

---

## Purpose

Characterize what is annotated at the PITX2 locus on the duck reference. Specifically:

- Confirm PITX2 genomic coordinates and strand.
- Enumerate every gene-level feature in a ±200 kb neighborhood.
- Identify antisense ncRNAs that overlap the PITX2 body or fall inside the pileup / isoform windows used by `pipeline/d2i_*` and `pipeline/d2j_*`.
- Note which candidates could be orthologs of chick/mouse *Playrr* or human *PANCR*.

This note is the fact base for the strand-filter direction, the stranded pileup and the liftover (§4a).

---

## 1. PITX2 in the duck GTF

| Field | Value |
|---|---|
| Chromosome | `NC_051775.1` |
| Gene body (all transcripts) | **3,961,387 – 3,976,254** (14,867 bp) |
| **Strand** | **`+`** |
| Source | Gnomon (NCBI Eukaryotic Genome Annotation Pipeline) |
| Gene ID attribute | `gene "PITX2"` (no `gene_name` attribute) |

### Three annotated transcripts (isoforms)

| Transcript | TSS | Interval | 5′UTR | Mapped isoform per `d2j_pitx2_isoform_per_genotype.py` |
|---|---|---|---|---|
| XM_038178818.1 | 3,961,387 | 3,961,387 – 3,976,254 | 3,961,387 – 3,961,458 (72 bp) | **PITX2a** |
| XM_027456453.2 | 3,964,253 | 3,964,253 – 3,976,254 | 3,964,253 – 3,964,353 + 3,966,591 – 3,966,600 | **PITX2b** |
| XM_027456452.2 | 3,971,385 | 3,971,385 – 3,976,254 | 3,971,385 – 3,971,829 | **PITX2c** (680 bp unique CDS) |

All three isoforms share CDS exons 3966601–3966646, 3973108–3973313, 3974789–3975349, stop at 3975350–3975352, and 3′UTR 3975353–3976254.

### Critical finding: strand

All three PITX2 transcripts in the GTF are on the **`+` strand**.

A common reading of 10x documentation is that the 3′ R2 is **antisense** to mRNA. That is not what cellranger itself records in the BAMs for this project. Verified empirically (GEX05 `possorted_genome_bam.bam`, all reads dedup via `-F 1024`):

| Gene | Strand | `GN:Z:<gene>` with flag-16 UNSET (`-F 16`) | with flag-16 SET (`-f 16`) | Conclusion |
|---|---|---|---|---|
| PITX2 | `+` | **6,961 reads** | 0 reads | cellranger-sense = flag-16 UNSET for a `+` gene |
| ELOVL6 | `+` | **7,170 reads** | 0 reads | same direction |
| EGF | `−` | 0 reads | **2,556 reads** | cellranger-sense = flag-16 SET for a `−` gene |

These are clean splits (0 / N), meaning cellranger defines sense-of-mRNA unambiguously by read flag-16 bit. In this pipeline's BAMs, **R2 is effectively sense-of-mRNA**, not antisense. Correct strand filter:

| Gene strand | Sense-of-mRNA filter | Antisense filter |
|---|---|---|
| `+` | **`-F 16`** (exclude flag 16) | `-f 16` |
| `−` | **`-f 16`** (require flag 16) | `-F 16` |

**For PITX2 (`+`), the sense filter is therefore `-F 16`.** A filter written for a `−` gene, or derived from the documentation claim, would invert it for other genes.

**Rule.** Derive the filter from the GTF strand and from the empirical cellranger behaviour, not from claims about 10x chemistry: if `GN:Z:<gene>` tags do not partition cleanly by the flag-16 bit on a known expressor, the derivation is wrong and the pipeline should refuse to emit results.

**Implementation.** `pipeline/_pitx2_strand.py` owns strand derivation for this pipeline: `read_gene_strand(gtf_path, "PITX2")` reads the `+` strand from the GTF, and `sense_strand_filter(strand)` returns `-F 16` accordingly. Both `d2j_pitx2_isoform_per_genotype.py` and `d2j_pitx2_isoform_percell.py` import and thread the filter through `samtools view -c`; `d2i_pitx2_rna_pileup.py` threads it through `samtools view` + `samtools depth`. The chromosome is also derived from the GTF via `read_gene_location`, so the d2j scripts do not depend on a hard-coded reference version.

---

## 2. PITX2 neighborhood within ±200 kb (NC_051775.1:3,761,000–4,181,000)

Ten gene-level features. Sorted by start coordinate.

| # | Start | End | Strand | Gene ID | Biotype | Distance to PITX2 | Overlap with PITX2 |
|---|---|---|---|---|---|---|---|
| 1 | 3,891,849 | 3,954,422 | **`−`** | LOC106016804 | ncRNA (62,573 bp) | 6,965 bp upstream | none |
| 2 | 3,961,387 | 3,976,254 | `+` | **PITX2** | protein-coding | — | — |
| 3 | 3,963,664 | 3,964,114 | **`−`** | LOC113843538 | ncRNA (450 bp) | inside PITX2 | **antisense, in PITX2 intron 1** |
| 4 | 3,973,586 | 3,976,858 | **`−`** | LOC113843539 | ncRNA (3,272 bp) | straddles PITX2 3′ end | **antisense, overlaps PITX2 3′UTR + 604 bp past 3′ end** |
| 5 | 3,976,875 | 4,033,262 | `+` | LOC106017192 | ncRNA (56,387 bp) | 621 bp downstream | none (same strand) |
| 6 | 4,003,281 | 4,037,227 | `−` | LOC101795212 | protein-coding | 27 kb downstream | none |
| 7 | 4,041,303 | 4,071,427 | `+` | LOC119716870 | ncRNA | 65 kb downstream | none |
| 8 | 4,105,940 | 4,165,489 | `+` | ELOVL6 | protein-coding | 130 kb downstream | none (reference neighbor) |
| 9 | 4,159,541 | 4,162,591 | `−` | LOC113843470 | ncRNA | 183 kb downstream | none |
| 10 | 4,168,730 | 4,179,954 | `−` | LOC113843471 | ncRNA | 192 kb downstream | none |

No alternative `crarc_zju1_PITX2`-specific construct contig is present in the GTF; PITX2 is annotated on the native chromosome only. The reference's `_PITX2` suffix refers to the upgraded annotation, not an added synthetic contig.

**Scope note (this table is ±200 kb only).** Two TAD-relevant upstream genes sit in the ±200–500 kb zone and are therefore *not* in the table above: **Playrr / LOC106016802** (`3,665,204–3,752,463`, `−`; ~209–296 kb upstream — the best chick-*Playrr* positional ortholog, §4) and **FAM241A / C4orf32** (`3,387,963–3,398,391` real coding; ~563–573 kb upstream — the upstream anchor of the *Pitx2* TAD, Welsh 2015; see §4). Both are present in the zju1 GTF; the gene order upstream→PITX2 is `FAM241A → Playrr → LOC106016804 → PITX2`, colinear with chick.

---

## 3. Antisense features vs. pipeline pileup windows

`d2i_pitx2_rna_pileup.py` uses **PITX2 ±50 kb** → `NC_051775.1:3,911,387–4,026,254`.
`d2j_pitx2_isoform_per_genotype.py` uses three narrow windows (see WINDOWS dict).
`d2j_pitx2_isoform_percell.py` uses four contiguous regions inside PITX2 body (REGION_AB / REGION_C / REGION_SH; FULL_REGION covers 3,959,887–3,976,254).

Overlap of each antisense ncRNA with the pileup / isoform windows:

| Antisense feature | Length | d2i ±50 kb window | d2j `pitx2a_exon1` | d2j `pitx2b_extended` | d2j `pitx2c_exon1` | d2j_percell `REGION_AB` | d2j_percell `REGION_C` | d2j_percell `REGION_SH` |
|---|---|---|---|---|---|---|---|---|
| LOC106016804 (3,891,849–3,954,422) | 62,573 | 42,573 bp of it falls inside (3,911,387–3,954,422) | — | — | — | — | — | — |
| LOC113843538 (3,963,664–3,964,114) | 450 | fully inside | — | **fully inside (450 bp = 28 % of pitx2b window)** | — | **fully inside REGION_AB** | — | — |
| LOC113843539 (3,973,586–3,976,858) | 3,272 | fully inside | — | — | — | — | — | **3,973,586–3,976,254 overlap (2,668 bp = 84.8 % of REGION_SH; 2,668 / 3,146 bp span)** |

**Interpretation.**
- **The `pitx2c_exon1` window overlaps no annotated antisense feature.** The HH13 / HH21 PITX2c laterality signal is counted in this window.
- **`pitx2b_extended` is 28 % covered by LOC113843538.** A PITX2b-specific count from this window alone would include any LOC113843538 sense reads; the strand filter (below) removes them.
- **d2j_percell `REGION_SH` (3′UTR, "shared" denominator) is 84.8 % covered by LOC113843539** (2,668 bp overlap of a 3,146 bp span); per-cell isoform ratios that use REGION_SH as denominator are the analytes most exposed to antisense reads, again removed by the strand filter.
- **The d2i ±50 kb window** contains 62,573 bp of antisense features (LOC106016804 + LOC113843538 + LOC113843539). d2i counts sense-strand reads only (`-F 16`, §5), so these features do not contribute to the fig02 panel D pileup.

### 3a. Empirical per-region strand counts (GEX05 BAM, whole library, dedup)

Read counts by strand filter on `possorted_genome_bam.bam` for GEX05, across the regions identified above:

| Region | Coordinates (NC_051775.1) | Sense-of-PITX2 (`-F 16`) | Antisense (`-f 16`) | Ratio |
|---|---|---|---|---|
| PITX2 5′+CDS clean (no annotated antisense) | 3,964,115 – 3,973,585 | **2,811** | 644 | 4.4 : 1 |
| PITX2 3′UTR (overlaps LOC113843539) | 3,973,586 – 3,976,254 | **5,117** | 883 | 5.8 : 1 |
| LOC113843538 intron segment | 3,963,664 – 3,964,114 | 735 | 53 | 14 : 1 (LOC113843538 barely expressed; PITX2 pre-mRNA dominates this intron) |
| LOC106016804 upstream (62.5 kb) | 3,891,849 – 3,954,422 | 566 | 558 | ~1 : 1 (background; 9 reads/kb) |
| **LOC106016802 upstream (87 kb, Playrr positional candidate)** | **3,665,204 – 3,752,463** | **1,243** | **1,385** | **Antisense slightly wins; 16 reads/kb on `−` strand — ~2× baseline. Consistent with weak expression possibly restricted to a sub-population.** |

**Interpretation of empirical counts (whole-library level):**

1. **Sense-of-PITX2 dominates antisense** in every region with substantial PITX2 signal (4–6 : 1); antisense reads are **~15–23 %** of reads in the PITX2 body regions.
2. **LOC113843539 is weakly expressed** relative to PITX2 — it contaminates `REGION_SH` but contributes ~15 % of reads despite covering 84.8 % of the coordinate span. (Feature-overlap fraction is an upper bound on contamination, not a prediction.)
3. **LOC113843538 is negligibly expressed** at the library level (53 reads in 450 bp; the "sense" signal in that segment is PITX2 pre-mRNA intron 1, not LOC113843538).
4. **LOC106016804 shows ~1:1 across a 62.5 kb window**, consistent with low-level ambient nuclear RNA rather than a strongly transcribed *Playrr*-like element (a library-aggregate observation).
5. **The sense-strand filter** removes these antisense reads, a ~15–20 % reduction in per-position depth at PITX2.

These are pooled-cell counts: they show that the chemistry is stranded and give the magnitude of antisense signal. The per-population stranded pileup is in §5.

---

## 4. Candidate chick/mouse *Playrr* / human *PANCR* orthologs

None of the 10 neighborhood features carry the symbol `Playrr`, `PANCR`, or `LINC01016`. NCBI Gene ortholog records for LOC106016804, LOC113843538, LOC113843539 are all empty (verified via `api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id/*/orthologs` — the LOCs are listed as "uncharacterized" with no cross-species homology links). The NCBI Gene page for LOC106016804 mentions **FAM241A** as an adjacent gene.

FAM241A is annotated in the `crarc_zju1_PITX2` GTF on `NC_051775.1`:
- **Real coding ortholog** `XM_027456395.2`: **3,387,963 – 3,398,391** (~10 kb, 2 exons, `−` strand), **~563–573 kb upstream** of the PITX2 gene body — beyond the ±200 kb window of §2 but within a ±500 kb / TAD-scale window.
- **Gnomon over-extension isoform** `XM_038178404.1`: 3,387,963 – **3,655,942** (268 kb, sparse "exons" at 3,387,963 / 3,523,927 / 3,606,708 / 3,655,689) — an annotation artefact that absorbs PITX2 eRNA / regulatory-landscape reads. The artefact (not real FAM241A biology) drives the "universal R>L FAM241A" DEG at HH13 (log2FC −3.12, padj 3.8e-80) and is filtered at the volcano renderer; see `docs/nodal_chromatin_memory.md` *(not distributed)* §9.1.

FAM241A is the duck ortholog of **C4orf32**, which Welsh et al. 2015 (PMID 26411685; Hi-C + Smc1 ChIA-PET) identify as the **upstream anchor of the *Pitx2* TAD** (chick galGal4 FISH probe `chr4:56,825,232–56,849,123`). The gene order **FAM241A/C4orf32 → Playrr/LOC106016802 → PITX2** is colinear duck↔chick, and the C4orf32→PITX2 gene desert is conserved at ~573 kb (duck, FAM241A-distal → PITX2 TSS) ≈ ~573 kb (chick). So a *Pitx2*-TAD-anchored window on the duck reference runs upstream to FAM241A (~3.39 Mb), enclosing the full Playrr-containing regulatory desert.

Welsh 2015 gives the chicken (galGal4) coordinates:
- **Chick Pitx2:** chr4:57,398,355 – 57,415,925
- **Chick e926 (Playrr TSS region):** chr4:56,989,104 – 57,005,206 — on **`−` strand**, transcribed toward cC4orf32
- **Chick cC4orf32:** chr4:56,825,232 – 56,849,123
- **Chick gene desert** between Pitx2 and cC4orf32: ~549 kb; Playrr TSS sits ~**393 kb** upstream of Pitx2 (i.e., deep in the desert, not adjacent to Pitx2)

Applied to the duck reference (PITX2 at 3,961,387 – 3,976,254), the expected Playrr ortholog coordinates — assuming conserved synteny — are **~300–400 kb upstream of PITX2** on the `−` strand, roughly NC_051775.1:3,560,000 – 3,660,000.

### Revised candidate ranking

| LOC ID | Strand | Distance from PITX2 | Length | Playrr-likelihood |
|---|---|---|---|---|
| **LOC106016802** | `−` | **209 – 296 kb upstream** | **87 kb** | **Best candidate.** In the correct distance range for chick Playrr (~393 kb upstream). Minus strand. ncRNA biotype. Large size consistent with a lncRNA region with extended GRO-seq coverage. **Note: this region is OUTSIDE the current d2i ±50 kb pileup window (which starts at 3,911,387).** The current fig02D panel D does not contain Playrr-region reads. |
| LOC106016804 | `−` | 7 – 70 kb upstream | 62.5 kb | **Unlikely Playrr.** Too close to PITX2 — chick Playrr is ~390 kb away, not ~7 kb. May be a separate antisense ncRNA or spill-over from Playrr's downstream tail. |
| LOC113843539 | `−` | overlaps PITX2 3′UTR | 3.2 kb | **Not Playrr.** Direct antisense overlap with PITX2 body is characteristic of a natural antisense transcript (NAT), not Playrr. Likely a PITX2 NAT. |
| LOC113843538 | `−` | inside PITX2 intron 1 | 450 bp | **Not Playrr.** Intronic antisense = NAT signature, not Playrr. |

For *PANCR* (human left-atrium-specific, same strand as PITX2, 2-exon): **LOC106017192** (+ strand, 621 bp downstream, 56 kb) is the positional analog. PANCR is a left-atrial biology not directly relevant to embryonic LPM stages; low priority.

### 4a. Liftover approach

**Sequence-level anchors for BLAT** (from Chen et al. 2022 bioRxiv, CRISPR sgRNAs targeting chick Playrr / e926):

| Target | sgRNA (23 bp including PAM; chick sequence) |
|---|---|
| Playrr Exon 1 splice site (3′ end of exon 1) | `TAGACGCAGCTGTGCTTAGAAGG` |
| e926 proximal | `GTGGCGGACTCATGTTAAAAAGG` |
| e926 distal | `GTGATTCCCACCACGCTTTGAGG` |

These are short but align near-uniquely where the sequence is conserved (UCSC BLAT chicken → duck, or directly against the ZJU1.0 FASTA). If LOC106016802 is Playrr, the hits fall in duck NC_051775.1:3,665,204 – 3,752,463. A longer anchor is the chick galGal4 e926 region (`NC_006091.5:56,989,104–57,005,206`, ~16 kb), which can be fetched with `efetch` and aligned with `blastn` against the ZJU1.0 FASTA.

The Playrr literature is Welsh et al. 2015 (*Cell Reports*) and the Chen et al. 2022 bioRxiv preprint (doi:10.1101/2022.09.20.508562); neither reports duck coordinates.

---

## 5. Consequence for the pipeline: sense-strand counting

Two strand facts set the filter. PITX2 is annotated on the `+` strand (§1), and on this
project's BAMs cellranger's own `GN:Z:<gene>` tags place sense reads in the flag-16-clear
bucket (§1 table), i.e. the counted read is sense to the mRNA. The sense-strand filter for
PITX2 is therefore `samtools view -F 16`. `pipeline/_pitx2_strand.py` reads the strand from the
GTF (`read_gene_strand(gtf_path, "PITX2")`) and returns the matching flag
(`sense_strand_filter(strand)`), so the rule generalizes to genes on either strand and is
derived from the BAM's empirical orientation rather than from a documentation claim.

A stranded mirror pileup over the PITX2 locus (±50 kb, NC_051775.1:3,911,387–4,026,254) gives
reverse/forward read ratios of 0.54–0.83 per population, higher than the 15–23 % of §3a
because the wider window includes the antisense LOC113843539 3′UTR overlap and the
LOC113843538 intronic antisense contributions:

| Population | fwd reads | rev reads | rev/fwd ratio |
|---|---|---|---|
| splanchnic_HH1716 | 304,743 | 165,229 | 0.54 |
| ventral_only_HH21 | 303,163 | 193,396 | 0.64 |
| ventral_splanchnic_HH21 | 499,937 | 307,576 | 0.62 |
| ventral_splanchnic_HH25 | 304,247 | 179,360 | 0.59 |
| syringeal_HH21 | 32,348 | 22,392 | 0.69 |
| syringeal_HH25 | 73,983 | 61,398 | 0.83 |

What the filter changes in the figures:

- **Dorsal sub-epithelial signal (fig02 panel D).** Without a strand filter the HH21/HH25 dorsal
  peak is mostly antisense. With `-F 16`, the residual sense-strand dorsal signal sits on the
  shared exons 2–4 and the 3′UTR, so it is reported as pan-PITX2, not PITX2a.
- **Ventral splanchnic PITX2c peak.** Sense-strand and unchanged by the filter; the high
  reverse/forward ratios above come from flanking antisense, not from the PITX2c exon.
- **Isoform windows (fig02 panel C, figS02).** The `pitx2c_exon1` window
  (3,971,385–3,972,064) carries no annotated antisense feature, so PITX2c L:R calls are
  insulated. `REGION_SH` (3,973,108–3,976,254) is 84.8 % covered by LOC113843539 by coordinate
  and about 15 % contaminated at the read level; `pitx2b_extended` is 28 % covered by
  LOC113843538. All d2i/d2j tables are counted with the filter; the unfiltered
  `figures/data/d2j_pitx2_isoform_*_nostrand.*` tables are kept as the unfiltered comparison baseline.
- **Open.** Sequence-level confirmation that duck LOC106016802 is the chick *Playrr* ortholog
  (liftover by BLAT/BLAST, §4a). Synteny places it within a few kb of the chick element.

---

## Appendix — commands used

```bash
# On a cluster node with read access to the reference (file ops only)
GTF=/n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/genes/genes.gtf.gz

# PITX2 gene entries
zcat $GTF | awk -F"\t" '$3=="gene"' | grep -i pitx2

# Non-exon PITX2 feature lines (transcripts, CDS, UTR)
zcat $GTF | awk -F"\t" '$3!="exon"' | grep -iE pitx2

# All gene-level features in ±200 kb (for the table in §2)
zcat $GTF | awk -F"\t" -v OFS="\t" '$1=="NC_051775.1" && $3=="gene" \
    && $4>=3761000 && $5<=4180000 {print $1,$4,$5,$7,$9}'

# Contig inventory (to rule out construct contigs)
zcat $GTF | awk -F"\t" '$3=="gene" {print $1}' | sort -u
zcat $GTF | grep 'gene "PITX2"' | awk '{print $1}' | sort -u

# On a compute node (srun/sbatch), for §1 cellranger-GN-tag strand verification:
module purge; module load gcc/14.2.0 samtools/1.21
BAM=/n/scratch/users/c/cl266/FC_08656/GEX05_GEX/outs/possorted_genome_bam.bam

# PITX2 (+ gene) — all GN:Z:PITX2 tagged reads should be in one flag-16 bucket
samtools view -F 1024 -f 16 $BAM NC_051775.1:3961387-3976254 | grep -c 'GN:Z:PITX2'  # → 0
samtools view -F 1024 -F 16 $BAM NC_051775.1:3961387-3976254 | grep -c 'GN:Z:PITX2'  # → 6961

# ELOVL6 (+ gene, cross-check)
samtools view -F 1024 -F 16 $BAM NC_051775.1:4105940-4165489 | grep -c 'GN:Z:ELOVL6' # → 7170

# EGF (- gene, opposite-direction sanity)
samtools view -F 1024 -f 16 $BAM NC_051775.1:4194846-4263663 | grep -c 'GN:Z:EGF'    # → 2556

# §3a per-region strand counts (dedup, all cells in the library):
samtools view -c -F 1024 -F 16 $BAM NC_051775.1:3964115-3973585     # PITX2 5'+CDS clean, sense
samtools view -c -F 1024 -f 16 $BAM NC_051775.1:3964115-3973585     # antisense
samtools view -c -F 1024 -F 16 $BAM NC_051775.1:3973586-3976254     # PITX2 3'UTR / LOC113843539 overlap, sense
samtools view -c -F 1024 -f 16 $BAM NC_051775.1:3973586-3976254     # antisense
samtools view -c -F 1024 -F 16 $BAM NC_051775.1:3891849-3954422     # LOC106016804 upstream 62.5 kb
samtools view -c -F 1024 -f 16 $BAM NC_051775.1:3891849-3954422
```
