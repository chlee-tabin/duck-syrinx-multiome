#!/bin/bash
# Stage-F bird-genome gathering — runs DETACHED on a cluster transfer node.
# Transfer node has internet + /n/scratch + /n/groups, but NO lmod/module system
# and only a minimal /n/app, so UCSC tools are fetched as STATIC binaries (no install).
#
# For each of 9 assemblies: download genome FASTA via NCBI `datasets` (FTP fallback),
# build genome.2bit + genome.sizes, stage to scratch, copy gzipped FASTA to persistent cdw.
# One bad genome does not kill the rest (|| continue), each download retried 3x.

set -u

# ---- paths ----
. "$(dirname "${BASH_SOURCE[0]}")/_paths.sh"   # portable roots; docs/audit_hardcoded_paths.md
SCRATCH_ROOT="$DUCK_SCRATCH_ROOT"
WORK="${SCRATCH_ROOT}/stagef_genomes"
CDW="${DUCK_REFERENCE_ROOT}/stagef"
BIN="${WORK}/bin"
MANIFEST="${WORK}/MANIFEST.tsv"

mkdir -p "$WORK" "$BIN" "$CDW"
cd "$WORK" || { echo "FATAL: cannot cd $WORK"; exit 1; }

echo "=== Stage-F gather START $(date) on $(hostname) ==="
echo "WORK=$WORK"
echo "CDW=$CDW"

# ---- fetch static tools once ----
fetch_tool() {
  # $1 = url, $2 = dest path
  local url="$1" dest="$2" i
  if [[ -x "$dest" ]]; then echo "tool exists: $dest"; return 0; fi
  for i in 1 2 3; do
    echo "fetch tool (try $i): $url"
    if curl -fSL "$url" -o "$dest" && [[ -s "$dest" ]]; then
      chmod +x "$dest"
      echo "tool OK: $dest"
      return 0
    fi
    sleep 5
  done
  echo "FATAL: could not fetch tool $url"
  return 1
}

DATASETS="${BIN}/datasets"
FATOTWOBIT="${BIN}/faToTwoBit"
TWOBITINFO="${BIN}/twoBitInfo"

fetch_tool "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets" "$DATASETS" || exit 1
fetch_tool "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit" "$FATOTWOBIT" || exit 1
fetch_tool "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo" "$TWOBITINFO" || exit 1

echo "--- tool versions ---"
"$DATASETS" --version 2>&1 | head -1 || true
"$FATOTWOBIT" 2>&1 | head -1 || true
echo "---------------------"

# ---- manifest header (overwrite fresh) ----
printf "species\taccession\tn_seqs\ttotal_bp\tscaffold_N50\tscratch_dir\tcdw_fasta_gz\tstatus\n" > "$MANIFEST"

# ---- genome list: "species_label<TAB>accession" ----
# RefSeq GCF preferred; magpie goose + screamer are GCA drafts (datasets handles GCA too).
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

# ---- download one accession's FASTA via datasets (3x retry), with FTP fallback ----
# echoes the path to the decompressed genomic .fna on success (to stdout via a marker file)
download_fasta() {
  # $1 = accession ; writes decompressed fna path to global FNA_PATH
  local acc="$1" i zip dldir fna
  zip="${acc}.zip"
  dldir="${acc}_dl"
  FNA_PATH=""

  # try datasets first
  for i in 1 2 3; do
    echo "  [download] datasets try $i for $acc"
    rm -f "$zip"
    if "$DATASETS" download genome accession "$acc" --include genome --filename "$zip" 2>&1; then
      rm -rf "$dldir"
      if unzip -o -q "$zip" -d "$dldir" 2>&1; then
        fna=$(find "$dldir/ncbi_dataset/data/$acc" -maxdepth 1 -name '*_genomic.fna' 2>/dev/null | head -1)
        if [[ -n "$fna" && -s "$fna" ]]; then
          FNA_PATH="$fna"
          echo "  [download] OK via datasets: $fna"
          return 0
        fi
      fi
    fi
    echo "  [download] datasets attempt $i failed for $acc"
    sleep 10
  done

  # FTP fallback: resolve genomic.fna.gz URL from assembly_summary via datasets summary
  echo "  [download] FALLBACK to FTP for $acc"
  local url=""
  # datasets summary gives the ftp_path-like info; build the NCBI FTP genomic.fna.gz path
  # GC[AF]_XXXXXXXXX.V -> /genomes/all/GCX/XXX/XXX/XXX/<full_asm_name>/<full_asm_name>_genomic.fna.gz
  local pfx="${acc:0:3}"           # GCA or GCF
  local digits="${acc:4:9}"        # 9 digits
  local d1="${digits:0:3}" d2="${digits:3:3}" d3="${digits:6:3}"
  # need full assembly name (accession + asm name); get via the all/ dir listing
  local basedir="https://ftp.ncbi.nlm.nih.gov/genomes/all/${pfx}/${d1}/${d2}/${d3}/"
  echo "  [download] listing $basedir"
  local asmdir
  asmdir=$(curl -fsL "$basedir" 2>/dev/null | grep -oE "${acc}[^/\"<> ]*" | head -1)
  if [[ -n "$asmdir" ]]; then
    url="${basedir}${asmdir}/${asmdir}_genomic.fna.gz"
    echo "  [download] FTP url: $url"
    for i in 1 2 3; do
      rm -rf "$dldir"; mkdir -p "$dldir/ncbi_dataset/data/$acc"
      if curl -fSL "$url" -o "${dldir}/${acc}_genomic.fna.gz" 2>&1; then
        if gunzip -f "${dldir}/${acc}_genomic.fna.gz" 2>&1; then
          fna="${dldir}/${acc}_genomic.fna"
          if [[ -s "$fna" ]]; then
            FNA_PATH="$fna"
            echo "  [download] OK via FTP: $fna"
            return 0
          fi
        fi
      fi
      echo "  [download] FTP attempt $i failed for $acc"
      sleep 10
    done
  else
    echo "  [download] could not resolve FTP asm dir for $acc"
  fi

  return 1
}

# ---- per-genome processing ----
N_OK=0
N_FAIL=0
for entry in "${GENOMES[@]}"; do
  sp="${entry%%|*}"
  acc="${entry##*|}"
  echo ""
  echo "==================================================================="
  echo "### GENOME: $sp ($acc)  $(date)"
  echo "==================================================================="

  process_one() {
    local outdir="${WORK}/${sp}"
    local cdwdir="${CDW}/${sp}"
    mkdir -p "$outdir" "$cdwdir"

    # 1. download
    if ! download_fasta "$acc"; then
      echo "FAIL $sp : download failed"
      return 1
    fi
    local src_fna="$FNA_PATH"

    # 2. place genome.fa
    echo "  [stage] moving FASTA -> ${outdir}/genome.fa"
    mv -f "$src_fna" "${outdir}/genome.fa" || { echo "FAIL $sp : mv genome.fa"; return 1; }

    # 3. build 2bit + sizes
    echo "  [build] faToTwoBit"
    "$FATOTWOBIT" "${outdir}/genome.fa" "${outdir}/genome.2bit" || { echo "FAIL $sp : faToTwoBit"; return 1; }
    echo "  [build] twoBitInfo"
    "$TWOBITINFO" "${outdir}/genome.2bit" "${outdir}/genome.sizes" || { echo "FAIL $sp : twoBitInfo"; return 1; }

    # 4. stats from sizes
    local n_seqs total_bp n50
    n_seqs=$(wc -l < "${outdir}/genome.sizes")
    total_bp=$(awk '{s+=$2} END{print s}' "${outdir}/genome.sizes")
    # scaffold N50: sort sizes desc, cumulative >= half total
    n50=$(awk '{print $2}' "${outdir}/genome.sizes" | sort -nr | awk -v tot="$total_bp" '
      { c+=$1; if (c >= tot/2) { print $1; exit } }')
    echo "  [stats] n_seqs=$n_seqs total_bp=$total_bp scaffold_N50=$n50"

    # 5. gzipped FASTA to persistent cdw
    echo "  [cdw] writing ${cdwdir}/genome.fa.gz"
    gzip -c "${outdir}/genome.fa" > "${cdwdir}/genome.fa.gz" || { echo "FAIL $sp : gzip to cdw"; return 1; }

    # 6. cleanup zip + dl dir to save scratch
    rm -f "${WORK}/${acc}.zip"
    rm -rf "${WORK}/${acc}_dl"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sp" "$acc" "$n_seqs" "$total_bp" "$n50" "$outdir" "${cdwdir}/genome.fa.gz" "OK" >> "$MANIFEST"
    echo "DONE $sp : n_seqs=$n_seqs total_bp=$total_bp N50=$n50"
    return 0
  }

  if process_one; then
    N_OK=$((N_OK+1))
  else
    N_FAIL=$((N_FAIL+1))
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sp" "$acc" "NA" "NA" "NA" "${WORK}/${sp}" "NA" "FAIL" >> "$MANIFEST"
    echo "FAIL $sp ($acc) — continuing to next genome"
    continue
  fi
done

echo ""
echo "=== Stage-F gather DONE $(date) ==="
echo "OK=$N_OK FAIL=$N_FAIL of ${#GENOMES[@]}"
echo "--- MANIFEST ---"
cat "$MANIFEST"
echo "=== END ==="
