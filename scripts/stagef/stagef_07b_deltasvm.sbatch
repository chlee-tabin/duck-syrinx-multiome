#!/bin/bash
#SBATCH -J stagef07b_dsvm
#SBATCH -p short
#SBATCH -t 0-00:30:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/07b_deltasvm_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/07b_deltasvm_%j.err
# Stage-F step 7b: deltaSVM on clade-FIXED substitutions within the duck L-peaks, scored
# with the duck h2 (left-vs-bilateral) gkm-SVM. Alignment-free corroboration; no
# accessibility data needed in the added species. Depends on step-6 per-element FASTAs.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
ROOT=${DUCK_SCRATCH_ROOT}/stagef
SC=$ROOT/scripts
FADIR=$ROOT/phyloacc/elements_fasta
DUCK=$ROOT/ref/duck_locus.fa
MAP=$ROOT/elements/elements_lpeaks.map.tsv
MODEL=${DUCK_SCRATCH_ROOT}/phase2/atac/gkmsvm/h2_model_strict_gc.model.txt
OUT=$ROOT/deltasvm
mkdir -p "$OUT"
[[ -d "$FADIR" ]] && ls "$FADIR"/*.fa >/dev/null 2>&1 || { echo "FATAL: no per-element FASTAs in $FADIR (run step 6 first)"; exit 1; }
[[ -s "$MODEL" ]] || { echo "FATAL: gkm-SVM model missing: $MODEL"; exit 1; }

# 1. clade-fixed substitutions -> ref/alt duck peak sequences
conda activate stagef
python3 "$SC/deltasvm_clade_fixed.py" "$FADIR" "$DUCK" "$MAP" "$OUT"

# 2. gkm-SVM scores for ref (duck actual) and alt (symmetric-reverted) peaks
conda activate gkmqc
echo "[07b] gkmpredict ref + alt"
gkmpredict "$OUT/ref_peaks.fa" "$MODEL" "$OUT/ref_scores.txt" 2> "$OUT/gkmpredict_ref.log"
gkmpredict "$OUT/alt_peaks.fa" "$MODEL" "$OUT/alt_scores.txt" 2> "$OUT/gkmpredict_alt.log"

# 3. deltaSVM = score(ref) - score(alt); join with sub counts
conda activate stagef
python3 - "$OUT" <<'PY'
import sys, csv, os
out = sys.argv[1]
def load(p):
    d = {}
    for ln in open(p):
        a = ln.split()
        if len(a) >= 2:
            try: d[a[0]] = float(a[1])
            except ValueError: pass
    return d
ref = load(os.path.join(out, "ref_scores.txt"))
alt = load(os.path.join(out, "alt_scores.txt"))
nsub = {}
with open(os.path.join(out, "elem_subs_summary.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        nsub[r["elem_id"]] = r["n_clade_fixed_subs"]
rows = []
for e in sorted(ref):
    if e in alt:
        rows.append((e, nsub.get(e, "0"), ref[e], alt[e], ref[e]-alt[e]))
with open(os.path.join(out, "deltasvm_results.tsv"), "w") as fh:
    fh.write("elem_id\tn_clade_fixed_subs\tscore_ref_duck\tscore_alt_symm\tdeltaSVM\n")
    for e, n, r, a, d in rows:
        fh.write(f"{e}\t{n}\t{r:.4f}\t{a:.4f}\t{d:.4f}\n")
print(f"deltaSVM rows: {len(rows)}  -> {out}/deltasvm_results.tsv")
# headline: positive deltaSVM with >=1 sub
hits = [r for r in rows if int(r[1]) > 0 and r[4] > 0]
print(f"peaks with clade-fixed subs AND positive deltaSVM (bulla bases raise L-access): {len(hits)}")
for e, n, r, a, d in sorted(hits, key=lambda x: -x[4])[:10]:
    print(f"  {e}  nsub={n}  deltaSVM={d:+.3f}")
PY
echo "STAGEF_07B_DONE __OK__"
