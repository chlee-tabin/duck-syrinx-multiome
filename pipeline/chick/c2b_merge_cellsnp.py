#!/usr/bin/env python3
"""c2b — merge per-library cellSNP mode-1 outputs into one cellData dir for vireo.

All three libraries were genotyped against the SAME discovered variant base
(c1 -> c2 -R), so the matrices share identical SNP rows; merging is an hstack of
the cell columns. Cell IDs are made unique by appending a per-library suffix.

Run (after c2):  conda run -n vireo python3 pipeline/chick/c2b_merge_cellsnp.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip, shutil
import scipy.io as sio
import scipy.sparse as sp

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux"
LIBS = [("ARC01", "1"), ("ARC02", "2")]  # Group A switched to v2.1.0 ARC01+ARC02 (both lanes = 12 embryos); ARC0S-chick dropped (no v2.1.0 rerun)
OUT = f"{WORK}/cellsnp/combined"
os.makedirs(OUT, exist_ok=True)


def read_mtx(stem):
    """Read a possibly-gzipped MatrixMarket file (stem has no .gz)."""
    for p in (stem + ".gz", stem):
        if os.path.exists(p):
            if p.endswith(".gz"):
                with gzip.open(p, "rb") as fh:
                    return sio.mmread(fh).tocsc()
            return sio.mmread(p).tocsc()
    raise FileNotFoundError(f"missing {stem}[.gz]")


AD, DP, OTH, samples = [], [], [], []
base_src = None
for L, suf in LIBS:
    d = f"{WORK}/cellsnp/{L}"
    ad = read_mtx(f"{d}/cellSNP.tag.AD.mtx")
    dp = read_mtx(f"{d}/cellSNP.tag.DP.mtx")
    AD.append(ad); DP.append(dp)
    try:
        OTH.append(read_mtx(f"{d}/cellSNP.tag.OTH.mtx"))
    except FileNotFoundError:
        pass
    sfile = f"{d}/cellSNP.samples.tsv"
    with open(sfile) as fh:
        bcs = [ln.strip() for ln in fh if ln.strip()]
    samples += [f"{bc}-{suf}" for bc in bcs]
    print(f"{L}: {ad.shape[0]} variants x {ad.shape[1]} cells")
    assert ad.shape[1] == len(bcs), f"{L}: AD cols ({ad.shape[1]}) != barcodes ({len(bcs)})"
    if base_src is None:
        for ext in ("cellSNP.base.vcf.gz", "cellSNP.base.vcf"):
            if os.path.exists(f"{d}/{ext}"):
                base_src = (f"{d}/{ext}", ext); break

# Rows (variants) must be identical across libraries (shared -R base)
nrows = {m.shape[0] for m in AD}
assert len(nrows) == 1, f"variant-row mismatch across libraries: {nrows}"

AD_c = sp.hstack(AD).tocoo()
DP_c = sp.hstack(DP).tocoo()
sio.mmwrite(f"{OUT}/cellSNP.tag.AD.mtx", AD_c, field="integer")
sio.mmwrite(f"{OUT}/cellSNP.tag.DP.mtx", DP_c, field="integer")
if OTH:
    sio.mmwrite(f"{OUT}/cellSNP.tag.OTH.mtx", sp.hstack(OTH).tocoo(), field="integer")

with open(f"{OUT}/cellSNP.samples.tsv", "w") as fh:
    fh.write("\n".join(samples) + "\n")

if base_src:
    dst = f"{OUT}/cellSNP.base.vcf.gz"
    if base_src[1].endswith(".gz"):
        shutil.copy(base_src[0], dst)
    else:
        with open(base_src[0], "rb") as fi, gzip.open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
else:
    sys.exit("ERROR: no cellSNP.base.vcf[.gz] found in any library dir")

print(f"\nMERGED -> {OUT}")
print(f"  {AD_c.shape[0]} variants x {AD_c.shape[1]} cells "
      f"({'+'.join(str(m.shape[1]) for m in AD)})")
