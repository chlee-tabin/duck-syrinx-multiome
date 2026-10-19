#!/usr/bin/env python3
"""c2c — merge Group B per-library cellSNP mode-1 outputs into one cellData for vireo.

Same shared variant base across libraries => hstack cell columns; suffix cell IDs by
lane (9 / 10), which encodes anatomical side. Run after c2b.
  conda run -n vireo python3 pipeline/chick/c2c_merge_cellsnp.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip, shutil
import scipy.io as sio
import scipy.sparse as sp

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
LIBS = [("GEX09_GEX", "9"), ("GEX10_GEX", "10")]
OUT = f"{WORK}/cellsnp/combined"; os.makedirs(OUT, exist_ok=True)


def read_mtx(stem):
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
    ad = read_mtx(f"{d}/cellSNP.tag.AD.mtx"); dp = read_mtx(f"{d}/cellSNP.tag.DP.mtx")
    AD.append(ad); DP.append(dp)
    try:
        OTH.append(read_mtx(f"{d}/cellSNP.tag.OTH.mtx"))
    except FileNotFoundError:
        pass
    with open(f"{d}/cellSNP.samples.tsv") as fh:
        bcs = [ln.strip() for ln in fh if ln.strip()]
    samples += [f"{bc}-{suf}" for bc in bcs]
    assert ad.shape[1] == len(bcs), f"{L}: AD cols != barcodes"
    print(f"{L}: {ad.shape[0]} variants x {ad.shape[1]} cells")
    if base_src is None:
        for ext in ("cellSNP.base.vcf.gz", "cellSNP.base.vcf"):
            if os.path.exists(f"{d}/{ext}"):
                base_src = (f"{d}/{ext}", ext); break

assert len({m.shape[0] for m in AD}) == 1, "variant-row mismatch across libraries"
sio.mmwrite(f"{OUT}/cellSNP.tag.AD.mtx", sp.hstack(AD).tocoo(), field="integer")
sio.mmwrite(f"{OUT}/cellSNP.tag.DP.mtx", sp.hstack(DP).tocoo(), field="integer")
if OTH:
    sio.mmwrite(f"{OUT}/cellSNP.tag.OTH.mtx", sp.hstack(OTH).tocoo(), field="integer")
with open(f"{OUT}/cellSNP.samples.tsv", "w") as fh:
    fh.write("\n".join(samples) + "\n")
if not base_src:
    sys.exit("ERROR: no cellSNP.base.vcf[.gz] found")
dst = f"{OUT}/cellSNP.base.vcf.gz"
if base_src[1].endswith(".gz"):
    shutil.copy(base_src[0], dst)
else:
    with open(base_src[0], "rb") as fi, gzip.open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)
print(f"\nMERGED -> {OUT}: {AD[0].shape[0]} variants x {sum(m.shape[1] for m in AD)} cells")
