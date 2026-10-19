#!/usr/bin/env python3
"""cS — targeted split of a suspect donor with select variants.

donor11 (Group B) / donor12 (Group A) are wide on the mesenchyme stage axis = candidate
same-sex heterogeneous donors the sex gate can't flag. This isolates the donor's cells,
keeps only INFORMATIVE variants (covered in >=K cells AND within-subset AF in [0.1,0.9] —
the sites that vary among its cells, i.e. could discriminate two sub-genotypes), writes a
subset cellSNP dir, and re-runs vireo N=2 and N=3. If the donor is really two embryos, vireo
splits it into GT-distinct sub-donors; we then sex each sub-donor (W/Z) to see if it's an
M+F or same-sex pair.

  sbatch pipeline/chick/run_python.sh pipeline/chick/cS_targeted_split.py <WORK> <BDIR> <DONOR>
"""
import os, sys, gzip, subprocess
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

WORK, BDIR, DONOR = sys.argv[1], sys.argv[2], sys.argv[3]
CELL = f"{WORK}/cellsnp/combined"; OUT = f"{WORK}/vireo_seeded/seed_qc/split_{DONOR}"
PCELL = f"{WORK}/vireo_seeded/seed_qc/per_cell_sex_wgene.tsv.gz"
os.makedirs(OUT, exist_ok=True)
KCELL = 15  # min cells covering a variant within the subset


def loadvcf(p):
    op = gzip.open if p.endswith(".gz") else open
    head, body = [], []
    with op(p, "rt") as f:
        for ln in f:
            (head if ln.startswith("#") else body).append(ln)
    return head, body


# ── isolate donor cells ──────────────────────────────────────────────────────
samples = [l.strip() for l in open(f"{CELL}/cellSNP.samples.tsv") if l.strip()]
sidx = {bc: i for i, bc in enumerate(samples)}
did = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
cells = [c for c in did.loc[did.donor_id == DONOR, "cell"] if c in sidx]
cols = [sidx[c] for c in cells]
print(f"[{DONOR}] {len(cells)} cells")
AD = sio.mmread(f"{CELL}/cellSNP.tag.AD.mtx").tocsc()[:, cols]
DP = sio.mmread(f"{CELL}/cellSNP.tag.DP.mtx").tocsc()[:, cols]
OTH = sio.mmread(f"{CELL}/cellSNP.tag.OTH.mtx").tocsc()[:, cols]

# ── select informative variants ──────────────────────────────────────────────
ncov = np.asarray((DP > 0).sum(1)).ravel()
ad_s = np.asarray(AD.sum(1)).ravel(); dp_s = np.asarray(DP.sum(1)).ravel()
af = np.divide(ad_s, dp_s, out=np.zeros_like(ad_s, float), where=dp_s > 0)
keep = (ncov >= KCELL) & (af > 0.1) & (af < 0.9)
print(f"[{DONOR}] informative variants: {int(keep.sum())} / {len(keep)} "
      f"(>= {KCELL} cells, 0.1<AF<0.9)")
ridx = np.where(keep)[0]

# ── write subset cellSNP dir ─────────────────────────────────────────────────
# Format MUST match the combined cellSNP dir EXACTLY (which vireo reads fine):
#   cellSNP.tag.{AD,DP,OTH}.mtx  -> UNCOMPRESSED (vireoSNP read_cellSNP reads the
#       non-.gz path via scipy.mmread; a .mtx.gz makes mmread choke on the gzip
#       magic bytes -> "Line 1: Not a Matrix Market file. Missing banner").
#   cellSNP.base.vcf.gz          -> gzipped   (load_VCF reads the .gz).
#   cellSNP.samples.tsv          -> uncompressed.
# remove any stale .mtx.gz left by a prior (failed) gzip-the-mtx run
subprocess.run(["rm", "-f"] + [f"{OUT}/cellSNP.tag.{e}.mtx.gz" for e in ("AD", "DP", "OTH")])
sio.mmwrite(f"{OUT}/cellSNP.tag.AD.mtx", sp.csr_matrix(AD[ridx]), field="integer")
sio.mmwrite(f"{OUT}/cellSNP.tag.DP.mtx", sp.csr_matrix(DP[ridx]), field="integer")
sio.mmwrite(f"{OUT}/cellSNP.tag.OTH.mtx", sp.csr_matrix(OTH[ridx]), field="integer")
with open(f"{OUT}/cellSNP.samples.tsv", "w") as f:
    f.write("\n".join(cells) + "\n")
head, body = loadvcf(next(f"{CELL}/{e}" for e in ("cellSNP.base.vcf.gz", "cellSNP.base.vcf")
                          if os.path.exists(f"{CELL}/{e}")))
with open(f"{OUT}/cellSNP.base.vcf", "w") as f:
    f.writelines(head); f.writelines(body[i] for i in ridx)
subprocess.run(["gzip", "-f", f"{OUT}/cellSNP.base.vcf"], check=True)

# ── per-cell sex for sub-donor labelling ─────────────────────────────────────
pc = pd.read_csv(PCELL, sep="\t").set_index("cell")


def gtcorr(vcf):
    head, body = loadvcf(vcf); donors = head[-1].rstrip("\n").split("\t")[9:]
    G = []
    for ln in body:
        x = ln.rstrip("\n").split("\t"); gi = x[8].split(":").index("GT") if "GT" in x[8] else 0
        row = []
        for s in x[9:]:
            g = s.split(":")[gi].replace("|", "/").split("/")
            try: row.append(sum(int(a) for a in g))
            except ValueError: row.append(np.nan)
        G.append(row)
    G = np.array(G, float).T
    if len(donors) < 2: return donors, np.nan
    m = ~(np.isnan(G[0]) | np.isnan(G[1]))
    return donors, (np.corrcoef(G[0][m], G[1][m])[0, 1] if m.sum() >= 50 else np.nan)


# ── run vireo N=2, N=3 ───────────────────────────────────────────────────────
failures = []
for N in (2, 3):
    od = f"{OUT}/N{N}"; os.makedirs(od, exist_ok=True)
    res = subprocess.run(["vireo", "--cellData", OUT, "--nDonor", str(N), "--nInit", "20",
                          "--extraDonor", "0", "--randSeed", "1", "-p", "2", "--outDir", od],
                         capture_output=True, text=True)
    if res.returncode or not os.path.exists(f"{od}/donor_ids.tsv"):
        print(f"VIREO FAILED N={N} (rc={res.returncode}):\n{(res.stderr or res.stdout)[-2500:]}", flush=True)
        failures.append(N)
        continue
    d = pd.read_csv(f"{od}/donor_ids.tsv", sep="\t")
    d = d[~d.donor_id.isin(["doublet", "unassigned"])]
    print(f"\n=== {DONOR} re-vireo N={N} ===")
    for sub, g in d.groupby("donor_id"):
        w = pc.reindex(g.cell)["percent_w"].median(); z = pc.reindex(g.cell)["percent_z"].median()
        print(f"  {sub}: {len(g)} cells  median_pct_w={w:.4f} pct_z={z:.4f}  "
              f"sex={'M' if w < 0.0025 else 'F'}")
    gt = next((f"{od}/{e}" for e in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf")
               if os.path.exists(f"{od}/{e}")), None)
    if gt and N == 2:
        donors, r = gtcorr(gt)
        print(f"  GT-distance between the 2 sub-donors: r={r:.3f} "
              f"({'distinct (real split)' if r < 0.9 else 'same genotype (NOT a split)'})")
print(f"\nwrote {OUT}/N2, {OUT}/N3")
if failures:
    sys.exit(f"cS: vireo failed for N={failures} (see traceback above) — non-zero exit so SLURM flags it")
