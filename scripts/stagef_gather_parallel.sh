#!/bin/bash
# PARALLEL relaunch of the Stage-F genome gather (replaces the sequential stagef_gather.sh).
# Each genome runs in its own subshell (no shared-global races); concurrency capped to be
# courteous on a shared transfer node. Idempotent: skips any genome already complete
# (completion marker = the cdw genome.fa.gz, written last). Reuses the static tools the
# sequential run already fetched under stagef_genomes/bin/. Portable (no `wait -n`).
set -u
. "$(dirname "${BASH_SOURCE[0]}")/_paths.sh"   # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

SCRATCH_ROOT="${DUCK_SCRATCH_ROOT}"
WORK="${SCRATCH_ROOT}/stagef_genomes"
CDW="${DUCK_REFERENCE_ROOT}/stagef"
BIN="${WORK}/bin"
MANIFEST="${WORK}/MANIFEST.tsv"
DATASETS="${BIN}/datasets"; FATOTWOBIT="${BIN}/faToTwoBit"; TWOBITINFO="${BIN}/twoBitInfo"
PAR=4
cd "$WORK" || exit 1
mkdir -p "$CDW"

GENOMES=(
  "Anas_acuta_pintail|GCF_963932015.1"
  "Aythya_fuligula_tufted_duck|GCF_009819795.1"
  "Cygnus_olor_mute_swan|GCF_009769625.2"
  "Cygnus_atratus_black_swan|GCF_013377495.2"
  "Anser_cygnoides_swan_goose|GCF_040182565.1"
  "Anseranas_semipalmata_magpie_goose|GCA_013399115.1"
  "Chauna_torquata_screamer|GCA_013399475.1"
  "Meleagris_gallopavo_turkey|GCF_000146605.3"
  "Struthio_camelus_ostrich|GCF_040807025.1"
)

do_genome() {
  local sp="$1" acc="$2"
  local outdir="${WORK}/${sp}" cdwdir="${CDW}/${sp}"
  exec > "${WORK}/gather_${sp}.log" 2>&1   # per-genome log, no interleave
  mkdir -p "$outdir" "$cdwdir"
  if [[ -s "${outdir}/genome.2bit" && -s "${cdwdir}/genome.fa.gz" ]]; then
    echo "SKIP $sp (already complete)"; return 0
  fi
  local zip="${acc}.par.zip" dldir="${acc}_par_dl" fna="" i
  for i in 1 2 3; do
    rm -f "$zip"
    if "$DATASETS" download genome accession "$acc" --include genome --filename "$zip"; then
      rm -rf "$dldir"
      if unzip -o -q "$zip" -d "$dldir"; then
        fna=$(find "$dldir/ncbi_dataset/data/$acc" -maxdepth 1 -name '*_genomic.fna' 2>/dev/null | head -1)
        [[ -n "$fna" && -s "$fna" ]] && break
      fi
    fi
    fna=""; sleep 10
  done
  if [[ -z "$fna" ]]; then                       # FTP fallback (GCA drafts)
    local pfx="${acc:0:3}" digits="${acc:4:9}"
    local d1="${digits:0:3}" d2="${digits:3:3}" d3="${digits:6:3}"
    local basedir="https://ftp.ncbi.nlm.nih.gov/genomes/all/${pfx}/${d1}/${d2}/${d3}/"
    local asmdir; asmdir=$(curl -fsL "$basedir" 2>/dev/null | grep -oE "${acc}[^/\"<> ]*" | head -1)
    if [[ -n "$asmdir" ]]; then
      for i in 1 2 3; do
        rm -rf "$dldir"; mkdir -p "$dldir"
        if curl -fSL "${basedir}${asmdir}/${asmdir}_genomic.fna.gz" -o "${dldir}/g.fna.gz" && gunzip -f "${dldir}/g.fna.gz"; then
          fna="${dldir}/g.fna"; [[ -s "$fna" ]] && break
        fi
        fna=""; sleep 10
      done
    fi
  fi
  if [[ -z "$fna" ]]; then
    echo "FAIL $sp : download"; printf "%s\t%s\tNA\tNA\tNA\t%s\tNA\tFAIL\n" "$sp" "$acc" "$outdir" >> "$MANIFEST"; return 1
  fi
  mv -f "$fna" "${outdir}/genome.fa"          || { echo "FAIL $sp mv"; return 1; }
  "$FATOTWOBIT" "${outdir}/genome.fa" "${outdir}/genome.2bit" || { echo "FAIL $sp faToTwoBit"; return 1; }
  "$TWOBITINFO" "${outdir}/genome.2bit" "${outdir}/genome.sizes" || { echo "FAIL $sp twoBitInfo"; return 1; }
  local n_seqs total_bp n50
  n_seqs=$(wc -l < "${outdir}/genome.sizes")
  total_bp=$(awk '{s+=$2} END{print s}' "${outdir}/genome.sizes")
  n50=$(awk '{print $2}' "${outdir}/genome.sizes" | sort -nr | awk -v tot="$total_bp" '{c+=$1; if(c>=tot/2){print $1; exit}}')
  gzip -c "${outdir}/genome.fa" > "${cdwdir}/genome.fa.gz" || { echo "FAIL $sp gzip cdw"; return 1; }
  rm -f "$zip"; rm -rf "$dldir"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n" "$sp" "$acc" "$n_seqs" "$total_bp" "$n50" "$outdir" "${cdwdir}/genome.fa.gz" >> "$MANIFEST"
  echo "DONE $sp n_seqs=$n_seqs total_bp=$total_bp N50=$n50"
}

echo "=== Stage-F PARALLEL gather START $(date) cap=$PAR ==="
[[ -s "$MANIFEST" ]] || printf "species\taccession\tn_seqs\ttotal_bp\tscaffold_N50\tscratch_dir\tcdw_fasta_gz\tstatus\n" > "$MANIFEST"
i=0
for entry in "${GENOMES[@]}"; do
  do_genome "${entry%%|*}" "${entry##*|}" &
  i=$((i+1)); (( i % PAR == 0 )) && wait
done
wait
echo "=== Stage-F PARALLEL gather DONE $(date) ==="
echo "--- complete genomes (have cdw gz) ---"; ls -1 "$CDW"/*/genome.fa.gz 2>/dev/null | wc -l
