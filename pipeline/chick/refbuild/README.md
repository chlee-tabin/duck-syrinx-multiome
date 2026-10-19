# cr7.reference.gg6a.ncbi.EGFP — build provenance + public re-derivation recipe

The chick Group-B snRNA libraries (GEX09, GEX10) were counted with **cellranger
7.1.0** against a custom reference `cr7.reference.gg6a.ncbi.EGFP` = NCBI RefSeq
GRCg6a (galGal6a) + an appended **EGFP transgene contig** (the pan-GFP
demultiplexing label — see "Chick Control Conventions"; no biological role).
That reference is NOT deposited in GEO (only the cellranger-arc `crarc_galgal6a`
fasta/GTF are). This directory makes it re-derivable.

## Contents

| file | what it is |
|---|---|
| `20231014_01_build_gg6a_ncbi.sh` | **VERBATIM archival record** of the mkref run (SLURM job 19590176, 2023-10-14, cellranger 7.1.0). Historical paths intentionally preserved — exempted in `scripts/verify_portable_paths.py`, do not "fix" them. |
| `EGFP.fa` | the EGFP transgene contig (2,903 bp) appended to the base genome |
| `EGFP.gtf` | its 3-line annotation (`EGFP` gene/transcript/exon, biotype protein_coding) |
| `reference.json` | the built reference's own manifest — **hash-pins the exact concatenated inputs**: `fasta_hash debd13707bcbcf87dd9e22d7b421b160260469d0` (`genome.EGFP.fa`), `gtf_hash.gz 022c03eb8e1d2651743f96adf29a5de3f1f0de36` (`genes.EGFP.gtf`), mkref cellranger-7.1.0 |

## Re-derivation recipe (public, from GEO)

1. Base genome + annotation: NCBI RefSeq **GRCg6a** (GCF_000002315.6) genome
   FASTA and GTF (the same assembly as the deposited `galGal6a_genome.fa`; the
   GEX build used the NCBI RefSeq annotation, not the cellranger-arc GTF).
2. `cat genome.fa EGFP.fa > genome.EGFP.fa` and `cat genes.gtf EGFP.gtf >
   genes.EGFP.gtf` (EGFP appended last).
3. Verify against `reference.json`: the two hashes above are cellranger mkref's
   input digests — a matching build proves input identity.
4. `cellranger mkref` exactly as in `20231014_01_build_gg6a_ncbi.sh`
   (cellranger 7.1.0).
5. Reads: `cellranger-arc bamtofastq` on the deposited
   `ARC09_gex_possorted_bam.bam` / `ARC10_gex_possorted_bam.bam` (GSE336327;
   same GEM wells as GEX09/GEX10 — BAMs retain all reads including the
   unmapped EGFP reads) → `cellranger count` (7.1.0) vs the rebuilt reference.

Caveat: step 3 is the ground truth for step 1's annotation choice — if the
NCBI GTF you fetch hashes differently after concatenation, diff it against the
canonical `star.reference.gg6a.ncbi/genes/genes.gtf` on the lab reference share
(`/n/groups/.../singlecell/reference`, alias `cdw`) before trusting counts.
