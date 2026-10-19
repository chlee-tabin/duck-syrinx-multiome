#!/usr/bin/env python3
"""cAH — surgical recovery of chick PITX2 3'UTR snRNA reads beyond the galGal6a annotation.

galGal6a RefSeq PITX2 ends at chr4:57,798,939 (shared last-exon/3'UTR end for all 3 isoforms).
snRNA-seq is 3'-biased; cAD's local pileup showed ~6.6% of sense PITX2 reads fall BEYOND
57,798,939, peaking ~57,799,078. cellranger (unextended GTF) therefore undercounts PITX2.

This script does NOT re-run cellranger. It mines the existing GEX09/GEX10 possorted BAMs by
CB barcode + UB UMI for reads in the 3'UTR extension, deduplicates (CB,UB) per cell (a
cellranger-UMI approximation — see ASSUMPTIONS), and compares to the current matrix PITX2 UMIs.

Mirrors cAD's depth_for() (samtools view -F 16 / -f 16, sense/antisense). Group B BAMs +
donor map per cAD. Restricts to assigned Group B barcodes (donor_ids.tsv, excl doublet/unassigned).

Steps:
  1. polyA endpoint: sense (-F 16) + antisense (-f 16) pileup over a wide window; endpoint =
     where sense depth returns to baseline (clean continuation, before any secondary feature).
  2. feasibility: for reads in [57,798,940, endpoint] sense, MAPQ>=255: confirm CB:Z + UB:Z
     present, +strand, report gene tag (GX/GN) profile, antisense fraction at the peak.
  3. per-cell: dedup (CB,UB) UMIs in the extension, sense, MAPQ>=255, per cell, per lane.
  4. recovery table: barcode<TAB>lane<TAB>recovered_pitx2_umis (cells with >0).
  5. matrix compare: per-cell current PITX2 UMIs from matrix; boost stats; L/R (lane) split.

ASSUMPTIONS / CAVEATS (honest):
  * UMI dedup here = unique (CB,UB) pairs among reads OVERLAPPING the extension window. cellranger
    dedups per (CB,UB,gene) over the whole transcript with UMI-collapse (1-mismatch) + a more
    refined multi-mapper/confidence model. So a (CB,UB) we count as "new" could already be counted
    by cellranger via an in-annotation read of the SAME molecule (then it is NOT additive). We
    bound this by also reporting (CB,UB) pairs whose read is EXCLUSIVELY in the extension (no GX
    PITX2 tag) — the conservative additive set. Reported both ways.
  * Antisense reads at the same coordinates (-f 16) are a potential confound (overlapping
    convergent transcript / internal-priming artifact). We quantify the antisense fraction and
    EXCLUDE antisense from the recovery tally (sense-only).
  * MAPQ>=255 (cellranger 10x convention for uniquely-mapped) gates out multi-mappers.

Outputs (scratch seed_qc/):
  cAH_pitx2_3utr_endpoint.tsv     per-base sense+antisense depth over the scan window (+ chosen endpoint)
  cAH_pitx2_3utr_feasibility.tsv  read-level QC: n reads, CB/UB-tagged frac, MAPQ profile, gene-tag profile, antisense frac
  cAH_pitx2_3utr_recovery.tsv     per-cell barcode/lane/recovered_pitx2_umis (>0 only) — gzipped
  cAH_pitx2_3utr_impact.tsv       per-cell join: recovered vs current matrix PITX2 UMIs, % boost, lane

  sbatch pipeline/chick/cAH_pitx2_3utr_recovery.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip, subprocess, traceback
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

WORK = f"{SCRATCH_USER}/chick_demux_B"
BDIR = f"{WORK}/vireo_seeded/N14/seed_105"
OUT = f"{WORK}/vireo_seeded/seed_qc"
LIBS = [("GEX09_GEX", "9", f"{WORK}/bams/GEX09_GEX/possorted_genome_bam.bam"),
        ("GEX10_GEX", "10", f"{WORK}/bams/GEX10_GEX/possorted_genome_bam.bam")]
PITX2_GENE_END = 57798939          # galGal6a RefSeq PITX2 last-exon / shared 3'UTR end
SCAN_START = 57798940              # one past the annotated end
SCAN_END = 57801500               # wide window to see the full decay + any secondary feature
PITX2_LASTEXON_START = 57797673    # shared last exon start (for in-annotation context)
CHROM_CANDIDATES = ["chr4", "NC_052535.1", "4"]
MAPQ_MIN = 255                     # 10x/cellranger uniquely-mapped convention
BASELINE_FRAC = 0.02               # endpoint = first sustained return to <=2% of peak depth
BASELINE_WIN = 50                  # sustained over this many bp


def log(m): print(m, flush=True); sys.stdout.flush()


def sh(cmd, timeout=1200):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)


def detect_chrom_suffix(bam):
    """Mirror cAD: find the PITX2-locus chrom name + the CB:Z suffix appended in this BAM."""
    chrom = None
    for c in CHROM_CANDIDATES:
        if sh(f"samtools view {bam} {c}:{PITX2_LASTEXON_START}-{PITX2_GENE_END} 2>/dev/null | head -50").stdout.strip():
            chrom = c; break
    if chrom is None:
        raise RuntimeError(f"no reads at PITX2 locus in {bam}")
    sfx = ""
    for line in sh(f"samtools view {bam} {chrom}:{PITX2_LASTEXON_START}-{PITX2_GENE_END} | head -3000").stdout.split("\n"):
        if "CB:Z:" in line:
            for tok in line.split("\t"):
                if tok.startswith("CB:Z:"):
                    cb = tok[5:]; sfx = cb[16:] if len(cb) > 16 else ""; break
            break
    return chrom, sfx


def write_bc(path, raw16_list, sfx):
    with open(path, "w") as fh:
        fh.write("\n".join(f"{b}{sfx}" for b in raw16_list) + "\n")


def depth_for(bam, bcfile, region, flag="-F 16"):
    """Barcode-filtered, strand-specific samtools depth. flag='-F 16' sense, '-f 16' antisense.
    -q MAPQ_MIN restricts to uniquely-mapped (matches recovery gate)."""
    tmp = f"/tmp/cah_{os.getpid()}_{abs(hash(region+flag+bcfile))%10**8}.bam"
    try:
        r1 = sh(f"samtools view -b {flag} -q {MAPQ_MIN} -D CB:'{bcfile}' '{bam}' {region} > '{tmp}'")
        if r1.returncode != 0:
            raise RuntimeError(f"samtools view failed (rc={r1.returncode}): {r1.stderr[:300]}")
        sh(f"samtools index '{tmp}'")
        out = sh(f"samtools depth -a -r {region} '{tmp}'").stdout
        d = {}
        for line in out.strip().split("\n"):
            if line:
                pp = line.split("\t")
                if len(pp) >= 3:
                    d[int(pp[1])] = int(pp[2])
        return d
    finally:
        for fn in (tmp, tmp + ".bai"):
            try: os.unlink(fn)
            except OSError: pass


# ── load assigned Group B barcodes (donor_ids.tsv, excl doublet/unassigned) ───
dm = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
dm = dm[~dm.donor_id.isin(["doublet", "unassigned"])].copy()
# cell = <16bp>-1-<lanenum>; lane suffix is the trailing -9 / -10
dm["raw16"] = dm.cell.str.extract(r"^([ACGT]+)")
dm["lanenum"] = dm.cell.str.extract(r"-(\d+)$")[0]
dm["lane"] = dm.lanenum.map({"9": "GEX09", "10": "GEX10"})
assert dm.lane.notna().all(), "unmapped lane suffix in donor_ids cell names"
lane_bcs = {ln: g for ln, g in dm.groupby("lane")}
log(f"[barcodes] {len(dm)} assigned Group B cells: " +
    ", ".join(f"{ln}={len(g)}" for ln, g in lane_bcs.items()))

bam_meta = {}
for lib, suf, bam in LIBS:
    c, sfx = detect_chrom_suffix(bam)
    bam_meta[lib] = (c, sfx)
    log(f"[{lib}] chrom={c} CB-suffix='{sfx}'")

BCDIR = f"{OUT}/cAH_bc"; os.makedirs(BCDIR, exist_ok=True)
chrom0 = bam_meta[LIBS[0][0]][0]


# ── STEP 1: polyA endpoint — sense + antisense pileup over the scan window ────
log("\n=== STEP 1: 3'UTR-extension endpoint (sense -F 16 + antisense -f 16) ===")
sense_tot = defaultdict(int); anti_tot = defaultdict(int)
# also pileup the in-annotation 3'-proximal region to define the local "peak" reference
ref_start = PITX2_LASTEXON_START
for lib, suf, bam in LIBS:
    lane = "GEX09" if suf == "9" else "GEX10"
    g = lane_bcs.get(lane)
    if g is None or not len(g):
        continue
    c, sfx = bam_meta[lib]
    bcfile = f"{BCDIR}/ALL_{lib}.txt"
    write_bc(bcfile, list(g.raw16), sfx)
    region = f"{c}:{ref_start}-{SCAN_END}"
    ds = depth_for(bam, bcfile, region, flag="-F 16")
    da = depth_for(bam, bcfile, region, flag="-f 16")
    for p, v in ds.items(): sense_tot[p] += v
    for p, v in da.items(): anti_tot[p] += v

positions = list(range(ref_start, SCAN_END + 1))
# peak depth in the in-annotation 3'-proximal window (last exon ... gene end) = reference level
inann = [sense_tot.get(p, 0) for p in range(PITX2_LASTEXON_START, PITX2_GENE_END + 1)]
peak_inann = max(inann) if inann else 0
# peak depth in the EXTENSION
ext_vals = {p: sense_tot.get(p, 0) for p in range(SCAN_START, SCAN_END + 1)}
ext_peak_pos = max(ext_vals, key=ext_vals.get) if ext_vals else SCAN_START
ext_peak = ext_vals.get(ext_peak_pos, 0)
# endpoint: first position p>SCAN_START where sense depth stays <= BASELINE_FRAC*ext_peak
# for a sustained BASELINE_WIN-bp run (so a secondary downstream bump does NOT get included)
thr = max(1, BASELINE_FRAC * ext_peak)
endpoint = SCAN_END
p = SCAN_START
while p <= SCAN_END - BASELINE_WIN:
    window = [sense_tot.get(q, 0) for q in range(p, p + BASELINE_WIN)]
    if max(window) <= thr:
        endpoint = p
        break
    p += 1
log(f"[endpoint] in-annotation 3' peak depth={peak_inann}; extension peak={ext_peak} @ {ext_peak_pos}; "
    f"baseline thr={thr:.1f}; endpoint(first {BASELINE_WIN}bp run <=thr)={endpoint} "
    f"(~{endpoint - PITX2_GENE_END} bp past annotated end)")
ep_rows = [{"position": p, "sense_depth": sense_tot.get(p, 0), "antisense_depth": anti_tot.get(p, 0)}
           for p in positions]
ep_df = pd.DataFrame(ep_rows)
ep_df["is_extension"] = ep_df.position > PITX2_GENE_END
ep_df["chosen_endpoint"] = endpoint
ep_df["chrom"] = chrom0
ep_df.to_csv(f"{OUT}/cAH_pitx2_3utr_endpoint.tsv", sep="\t", index=False)
log(f"[wrote] cAH_pitx2_3utr_endpoint.tsv ({len(ep_df)} positions)")


# ── STEP 2+3: read-level feasibility + per-cell UMI recovery in [SCAN_START, endpoint] ──
log(f"\n=== STEP 2/3: per-read feasibility + per-cell (CB,UB) recovery in "
    f"{chrom0}:{SCAN_START}-{endpoint} (sense, MAPQ>={MAPQ_MIN}) ===")
allowed = {ln: set(g.raw16 + bam_meta["GEX09_GEX" if ln == "GEX09" else "GEX10_GEX"][1])
           for ln, g in lane_bcs.items()}
# build allowed CB-with-suffix sets per lane (BAM CB form)
allowed_cb = {}
for ln, g in lane_bcs.items():
    lib = "GEX09_GEX" if ln == "GEX09" else "GEX10_GEX"
    sfx = bam_meta[lib][1]
    allowed_cb[ln] = set(f"{b}{sfx}" for b in g.raw16)

# read-level QC accumulators
qc = {"n_reads_window": 0, "n_sense": 0, "n_antisense": 0, "n_sense_mapq_pass": 0,
      "n_sense_mapq_pass_cbub": 0, "n_no_cb": 0, "n_no_ub": 0, "n_cb_not_in_set": 0,
      "n_sense_with_gx_pitx2": 0, "n_sense_with_gx_other": 0, "n_sense_no_gx": 0}
gx_counter = defaultdict(int)   # gene-tag (GX) profile of sense reads in the window
# per-cell UMI sets
umis_overlap = defaultdict(set)        # (lane,cb)->{ub}: any sense read overlapping the window (broad)
umis_exclusive = defaultdict(set)      # (lane,cb)->{ub}: sense read with no PITX2 GX tag (conservative additive)
mapq_hist = defaultdict(int)

for lib, suf, bam in LIBS:
    lane = "GEX09" if suf == "9" else "GEX10"
    if lane not in allowed_cb:
        continue
    c, sfx = bam_meta[lib]
    cbset = allowed_cb[lane]
    region = f"{c}:{SCAN_START}-{endpoint}"
    # stream all reads in the window (both strands) for QC; filter in python
    proc = subprocess.Popen(f"samtools view '{bam}' {region}", shell=True,
                            stdout=subprocess.PIPE, text=True, bufsize=1)
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qc["n_reads_window"] += 1
        flag = int(f[1])
        mapq = int(f[4])
        is_rev = bool(flag & 16)
        # tags
        cb = ub = gx = None
        for t in f[11:]:
            if t.startswith("CB:Z:"): cb = t[5:]
            elif t.startswith("UB:Z:"): ub = t[5:]
            elif t.startswith("GX:Z:"): gx = t[5:]
        if is_rev:
            qc["n_antisense"] += 1
            continue
        # sense (forward) reads only from here
        qc["n_sense"] += 1
        mapq_hist[mapq] += 1
        # gene-tag profile (sense)
        if gx is None:
            qc["n_sense_no_gx"] += 1
        else:
            gx_counter[gx] += 1
            if "PITX2" in gx.upper():
                qc["n_sense_with_gx_pitx2"] += 1
            else:
                qc["n_sense_with_gx_other"] += 1
        if mapq < MAPQ_MIN:
            continue
        qc["n_sense_mapq_pass"] += 1
        if cb is None:
            qc["n_no_cb"] += 1; continue
        if cb not in cbset:
            qc["n_cb_not_in_set"] += 1; continue
        if ub is None:
            qc["n_no_ub"] += 1; continue
        qc["n_sense_mapq_pass_cbub"] += 1
        key = (lane, cb)
        umis_overlap[key].add(ub)
        # conservative: read carries no PITX2 GX (i.e. cellranger would NOT have assigned it to PITX2)
        if gx is None or "PITX2" not in gx.upper():
            umis_exclusive[key].add(ub)
    proc.wait()

# antisense fraction at the sense peak coordinate (confound metric)
peak_sense = sense_tot.get(ext_peak_pos, 0)
peak_anti = anti_tot.get(ext_peak_pos, 0)
anti_frac_peak = peak_anti / max(1, peak_sense + peak_anti)
total_anti_frac = qc["n_antisense"] / max(1, qc["n_sense"] + qc["n_antisense"])

mapq_top = sorted(mapq_hist.items(), key=lambda x: -x[1])[:6]
gx_top = sorted(gx_counter.items(), key=lambda x: -x[1])[:10]
feas_rows = [{"metric": k, "value": v} for k, v in qc.items()]
feas_rows += [
    {"metric": "endpoint", "value": endpoint},
    {"metric": "bp_past_annotated_end", "value": endpoint - PITX2_GENE_END},
    {"metric": "ext_peak_pos", "value": ext_peak_pos},
    {"metric": "ext_peak_sense_depth", "value": ext_peak},
    {"metric": "antisense_frac_at_peak", "value": round(anti_frac_peak, 4)},
    {"metric": "antisense_frac_window_overall", "value": round(total_anti_frac, 4)},
    {"metric": "mapq_top", "value": str(mapq_top)},
    {"metric": "gx_top", "value": str(gx_top)},
    {"metric": "n_cells_with_recovery_overlap", "value": sum(1 for v in umis_overlap.values() if v)},
    {"metric": "n_cells_with_recovery_exclusive", "value": sum(1 for v in umis_exclusive.values() if v)},
    {"metric": "total_umis_overlap", "value": sum(len(v) for v in umis_overlap.values())},
    {"metric": "total_umis_exclusive", "value": sum(len(v) for v in umis_exclusive.values())},
]
pd.DataFrame(feas_rows).to_csv(f"{OUT}/cAH_pitx2_3utr_feasibility.tsv", sep="\t", index=False)
log("[feasibility] " + "; ".join(f"{k}={v}" for k, v in qc.items()))
log(f"[feasibility] antisense_frac_at_peak={anti_frac_peak:.3f} overall={total_anti_frac:.3f}; "
    f"mapq_top={mapq_top}; gx_top={gx_top}")
log(f"[recovery] cells(overlap)={sum(1 for v in umis_overlap.values() if v)} "
    f"umis(overlap)={sum(len(v) for v in umis_overlap.values())}; "
    f"cells(excl)={sum(1 for v in umis_exclusive.values() if v)} "
    f"umis(excl)={sum(len(v) for v in umis_exclusive.values())}")


# ── per-cell recovery table (overlap = broad; exclusive = conservative additive) ──
rec_rows = []
all_keys = set(umis_overlap) | set(umis_exclusive)
for (lane, cb) in all_keys:
    ov = len(umis_overlap.get((lane, cb), ()))
    ex = len(umis_exclusive.get((lane, cb), ()))
    # strip the BAM CB suffix back to the -1 matrix barcode form for the deliverable
    sfx = bam_meta["GEX09_GEX" if lane == "GEX09" else "GEX10_GEX"][1]
    bc16 = cb[:-len(sfx)] if sfx and cb.endswith(sfx) else cb
    rec_rows.append({"barcode": bc16, "lane": lane,
                     "recovered_pitx2_umis": ov, "recovered_exclusive_umis": ex})
rec = pd.DataFrame(rec_rows)
rec = rec[rec.recovered_pitx2_umis > 0].sort_values(["lane", "recovered_pitx2_umis"],
                                                    ascending=[True, False])
rec.to_csv(f"{OUT}/cAH_pitx2_3utr_recovery.tsv.gz", sep="\t", index=False, compression="gzip")
log(f"[wrote] cAH_pitx2_3utr_recovery.tsv.gz ({len(rec)} cells with >0 recovery)")


# ── STEP 5: matrix compare — current per-cell PITX2 UMIs vs recovered ─────────
log("\n=== STEP 5: impact vs current matrix PITX2 UMIs ===")
cur_rows = []
for lib, suf, bam in LIBS:
    lane = "GEX09" if suf == "9" else "GEX10"
    md = f"{WORK}/matrix/{lib}"
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
    sym = [f[1] if len(f) > 1 else f[0] for f in feats]
    pidx = [i for i, s in enumerate(sym) if s.strip().upper() == "PITX2"]
    if not pidx:
        log(f"[{lib}] PITX2 not in features?!"); continue
    pidx = pidx[0]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    prow = np.asarray(M[pidx].todense()).ravel()
    for b, v in zip(bcs, prow):
        cur_rows.append({"barcode": b, "lane": lane, "current_pitx2_umis": int(v)})
cur = pd.DataFrame(cur_rows)
# join recovery (overlap + exclusive) onto current matrix counts (restricted to assigned cells)
assigned_set = set(zip(dm.raw16 + "-1", dm.lane))   # matrix barcodes are <16bp>-1
cur["bc16"] = cur.barcode.str.extract(r"^([ACGT]+)")
cur["is_assigned"] = list(zip(cur.bc16 + "-1", cur.lane))
cur["is_assigned"] = cur["is_assigned"].isin(assigned_set)
cur_a = cur[cur.is_assigned].copy()
rec_j = rec.rename(columns={"barcode": "bc16"})[["bc16", "lane", "recovered_pitx2_umis", "recovered_exclusive_umis"]]
imp = cur_a.merge(rec_j, on=["bc16", "lane"], how="left")
imp["recovered_pitx2_umis"] = imp.recovered_pitx2_umis.fillna(0).astype(int)
imp["recovered_exclusive_umis"] = imp.recovered_exclusive_umis.fillna(0).astype(int)
imp["new_total_overlap"] = imp.current_pitx2_umis + imp.recovered_pitx2_umis
imp["new_total_exclusive"] = imp.current_pitx2_umis + imp.recovered_exclusive_umis
imp.to_csv(f"{OUT}/cAH_pitx2_3utr_impact.tsv", sep="\t", index=False)

def summarize(df, lbl):
    cur_pitx2 = df.current_pitx2_umis
    rec_ov = df.recovered_pitx2_umis
    rec_ex = df.recovered_exclusive_umis
    n_cells = len(df)
    n_cur_pos = int((cur_pitx2 > 0).sum())
    cur_total = int(cur_pitx2.sum())
    ov_total = int(rec_ov.sum())
    ex_total = int(rec_ex.sum())
    # per-cell % boost among cells that currently have PITX2 (avoid div0)
    pos = df[cur_pitx2 > 0]
    pct_ov = (pos.recovered_pitx2_umis / pos.current_pitx2_umis * 100) if len(pos) else pd.Series(dtype=float)
    pct_ex = (pos.recovered_exclusive_umis / pos.current_pitx2_umis * 100) if len(pos) else pd.Series(dtype=float)
    # cells newly gaining PITX2 (currently 0, recovery >0)
    new_cells_ov = int(((cur_pitx2 == 0) & (rec_ov > 0)).sum())
    new_cells_ex = int(((cur_pitx2 == 0) & (rec_ex > 0)).sum())
    log(f"\n[{lbl}] n_cells={n_cells}; PITX2+ cells(current)={n_cur_pos}; "
        f"current PITX2 UMIs total={cur_total}")
    log(f"[{lbl}] recovered UMIs — overlap={ov_total} ({100*ov_total/max(1,cur_total):.1f}% of current total), "
        f"exclusive={ex_total} ({100*ex_total/max(1,cur_total):.1f}% of current total)")
    log(f"[{lbl}] per-cell %boost among PITX2+ cells (overlap): "
        f"median={pct_ov.median():.1f}% mean={pct_ov.mean():.1f}%; "
        f"(exclusive): median={pct_ex.median():.1f}% mean={pct_ex.mean():.1f}%")
    log(f"[{lbl}] cells newly gaining PITX2 (was 0): overlap={new_cells_ov} exclusive={new_cells_ex}")
    return dict(label=lbl, n_cells=n_cells, pitx2_pos_cells_current=n_cur_pos,
                current_total=cur_total, recovered_overlap=ov_total, recovered_exclusive=ex_total,
                pct_total_overlap=round(100*ov_total/max(1,cur_total),2),
                pct_total_exclusive=round(100*ex_total/max(1,cur_total),2),
                median_pct_boost_overlap=round(float(pct_ov.median()),2) if len(pct_ov) else 0.0,
                mean_pct_boost_overlap=round(float(pct_ov.mean()),2) if len(pct_ov) else 0.0,
                median_pct_boost_exclusive=round(float(pct_ex.median()),2) if len(pct_ex) else 0.0,
                mean_pct_boost_exclusive=round(float(pct_ex.mean()),2) if len(pct_ex) else 0.0,
                new_cells_overlap=new_cells_ov, new_cells_exclusive=new_cells_ex)

summaries = [summarize(imp, "ALL")]
for ln, g in imp.groupby("lane"):
    summaries.append(summarize(g, ln))
pd.DataFrame(summaries).to_csv(f"{OUT}/cAH_pitx2_3utr_impact_summary.tsv", sep="\t", index=False)
log(f"\n[wrote] cAH_pitx2_3utr_impact.tsv + cAH_pitx2_3utr_impact_summary.tsv")

# L/R-proportionality check: does recovery preserve the GEX09 vs GEX10 ratio of current PITX2?
g09 = imp[imp.lane == "GEX09"]; g10 = imp[imp.lane == "GEX10"]
cur09, cur10 = g09.current_pitx2_umis.sum(), g10.current_pitx2_umis.sum()
rec09, rec10 = g09.recovered_pitx2_umis.sum(), g10.recovered_pitx2_umis.sum()
log(f"\n[L/R-proportionality] current PITX2 UMIs GEX09={cur09} GEX10={cur10} "
    f"(GEX09 frac={cur09/max(1,cur09+cur10):.3f})")
log(f"[L/R-proportionality] recovered UMIs GEX09={rec09} GEX10={rec10} "
    f"(GEX09 frac={rec09/max(1,rec09+rec10):.3f})")
log(f"[L/R-proportionality] => recovery {'PRESERVES' if abs(cur09/max(1,cur09+cur10) - rec09/max(1,rec09+rec10)) < 0.10 else 'SHIFTS'} "
    f"the lane ratio (|Δfrac| {abs(cur09/max(1,cur09+cur10) - rec09/max(1,rec09+rec10)):.3f}); "
    f"adds power, direction-change risk {'LOW' if abs(cur09/max(1,cur09+cur10) - rec09/max(1,rec09+rec10)) < 0.10 else 'CHECK'}")
log("\n=== cAH done ===")
