#!/bin/bash
# Is the reported f5 PITX2-region L/R DESeq2 result an artifact of the
# MIN_CELLS=20 pseudobulk floor (chosen, per generate_beeswarm.R's own comment, after seeing
# that 30 pushed a Wilcoxon toward marginal)?  Re-run the IDENTICAL model
# (figures/figS03/beeswarm_deseq2_stats.R) at 10 / 20 / 30 / 50 and stack the results.
# Local, deps = conda duckfig R + DESeq2.  The reported table (threshold 20) is untouched by
# this script; the R script writes `min_cells` and `note` columns itself, so the stack is a
# plain concatenation with one header.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=figures/data/step_f5_pitx2_region_deseq2_sweep.tsv
FLOORS=(10 20 30 50)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
: > "$TMP/stack.tsv"                                # never append to a stale partial file
for t in "${FLOORS[@]}"; do
    MIN_CELLS=$t OUT_PATH="$TMP/t$t.tsv" conda run -n duckfig Rscript figures/figS03/beeswarm_deseq2_stats.R \
        > "$TMP/t$t.log" 2>&1 \
        || { echo "FAIL at MIN_CELLS=$t — R log follows:" >&2; tail -30 "$TMP/t$t.log" >&2; exit 1; }
    if [ ! -s "$TMP/stack.tsv" ]; then head -1 "$TMP/t$t.tsv" > "$TMP/stack.tsv"; fi   # header once, from the first floor that ran
    tail -n +2 "$TMP/t$t.tsv" >> "$TMP/stack.tsv"
done
# ponytail: one runnable check BEFORE installing — the reported (t=20) rows must reproduce
# the committed step_f5_pitx2_region_deseq2.tsv. Compared at 6 significant figures because
# cluster-vs-local R differ in BLAS rounding at the 13th-15th digit, which is not a discrepancy.
# Column positions are resolved from the HEADER, not hardcoded: the R script owns the schema
# and has gained columns twice (min_cells/note, then n_L_fit/n_R_fit).
if ! diff <(head -1 "$TMP/stack.tsv") <(head -1 figures/data/step_f5_pitx2_region_deseq2.tsv); then
    echo "MISMATCH: sweep and reported table have different headers" >&2
    exit 1
fi
# round the three float columns by NAME (BLAS rounding differs at the 13th-15th digit between
# cluster and local R, which is not a discrepancy)
fmt6() { awk -F'\t' -v OFS='\t' '
    NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    { for (c in col) if (c ~ /^deseq2_/ && $col[c] != "NA") $col[c] = sprintf("%.6g", $col[c]); print }' "$1"; }
sel20() { awk -F'\t' '
    NR==1 { for (i = 1; i <= NF; i++) if ($i == "min_cells") m = i;
            if (!m) { print "no min_cells column" > "/dev/stderr"; exit 3 } print; next }
    $m == 20' "$1"; }
if ! diff <(sel20 "$TMP/stack.tsv" | fmt6 /dev/stdin) <(fmt6 figures/data/step_f5_pitx2_region_deseq2.tsv) > "$TMP/t20.diff"; then
    echo "MISMATCH: sweep t=20 rows do not reproduce the committed step_f5_pitx2_region_deseq2.tsv (6 s.f.):" >&2
    cat "$TMP/t20.diff" >&2
    exit 1
fi
mv "$TMP/stack.tsv" "$OUT"
echo "wrote $OUT (floors: ${FLOORS[*]}); OK: t=20 rows reproduce the committed table (6 s.f.)"
