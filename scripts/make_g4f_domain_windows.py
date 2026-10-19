#!/usr/bin/env python3
"""make_g4f_domain_windows.py — multi-stage PITX2c cis-block motif-window BED (fig04 panel C).

Supersedes make_g4f_hh21_windows.py.  Builds the windows for fig04 panel C (the HH25 wave-2
snapshot; panel C was consolidated to a single stage since HH17/16's L-sig focal peaks are a
strict subset of HH25's).  The focal block is the CONTIGUOUS L-biased domain = the HH17/16 ∪ HH25
ventral L-sig union, with per-stage flags written to the sidecar; plus a within-locus BILATERAL
background (replacing the 4 cherry-picked "ns ctrl") and the distal R-desert artifact landmark.

Windows (name = <label>|<class>):
  focal     : contiguous PITX2c cis-block L-biased peaks (HH17/16 ∪ HH25 ventral Lsig,
              padj_locus<0.1, dist −5..+20 kb), EXCLUDING the broken-out c-prom / +5.4 kb.
  cprom     : PITX2c promoter (dist≈0) — the CpG-island motif DESERT (broken out).
  p5p4      : +5.4 kb putative enhancer (broken out).
  aprom     : PITX2a promoter (−10 kb) — the SMAD5 / BMP-SBE anchor (+ve control).
  bprom     : PITX2b promoter (−7 kb).
  Rdesert   : distal R-biased desert (HH25, −400..−443 kb) — FAM241A/eRNA artifact landmark
              (DUCK-only; 460+ kb distal, outside the proximal chick lift window → chick ×).
  bilateral : ALL non-laterally-biased peaks in the proximal locus (±BIL_KB, non-sig at BOTH
              stages) = the within-locus GC/locale-matched background motif rate.

Per-peak stage significance (which stage each focal/bilateral peak is L-sig in) is written to
the sidecar g4f_domain_window_stageflags.tsv so the plot can aggregate the L-block column
per stage (n grows HH17/16 → HH25) WITHOUT re-scanning (motif content is genomic/stage-invariant).

Consistent with panel A (same g4f floor-20 census, ventral inference pool).

  conda run -n duckfig python3 scripts/make_g4f_domain_windows.py
"""
import collections
import os

import pandas as pd

# fig04 panel C's scan-window builder; the BED output lives in figures/data/.
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, "figures", "data")
TSS = 3971385
FOCAL_LO, FOCAL_HI = -5000, 20000     # contiguous cis-block dist window
BIL_KB = 50000                        # bilateral background half-width = the fig03/panel-A focal
                                      # ±50 kb scope (GC/locale-matched, liftable to chick)
PADJ = 0.1
CPROM_LO, CPROM_HI = -200, 300        # dist window that tags the c-promoter peak
P54_LO, P54_HI = 5000, 5800           # dist window that tags the +5.4 kb peak


def coords(pid):
    c, rng = pid.split(":")
    s, e = rng.split("-")
    return c, int(s), int(e)


def mid(pid):
    _, s, e = coords(pid)
    return (s + e) // 2


def main():
    g = pd.read_csv(os.path.join(DATA, "g4f_lr_long_floor20.tsv"), sep="\t")
    v = g[(g.design == "unpaired") & (g.compartment == "ventral")].copy()
    v["dist"] = v.peak_id.map(lambda p: mid(p) - TSS)

    def lsig(stage):
        h = v[v.stage == stage]
        return set(h[(h.padj_locus < PADJ) & (h.log2FoldChange > 0)].peak_id)

    def anysig(stage):
        h = v[v.stage == stage]
        return set(h[h.padj_locus < PADJ].peak_id)

    L1716, L25 = lsig("HH17/16"), lsig("HH25")
    sig_any = anysig("HH17/16") | anysig("HH25")

    # one row per distinct peak_id with its dist (ventral table may repeat across stages)
    peak_dist = v.drop_duplicates("peak_id").set_index("peak_id")["dist"]

    # ── focal contiguous cis-block (HH17/16 ∪ HH25 Lsig, in the cis window) ──────
    focal = sorted({p for p in (L1716 | L25)
                    if FOCAL_LO <= peak_dist[p] <= FOCAL_HI}, key=lambda p: peak_dist[p])

    rows, flags = [], []   # rows = (chrom,start,end,name); flags = stage-sig sidecar

    def add(pid, label, cls):
        c, s, e = coords(pid)
        d = int(peak_dist[pid])
        rows.append((c, s, e, f"{label}|{cls}"))
        flags.append((label, cls, d, int(pid in L1716), int(pid in L25)))

    for pid in focal:
        d = peak_dist[pid]
        if CPROM_LO <= d <= CPROM_HI:
            add(pid, "PITX2c_prom", "cprom")
        elif P54_LO <= d <= P54_HI:
            add(pid, "Peak_p5.4", "p5p4")
        else:
            add(pid, f"L{d:+d}", "focal")

    # ── a / b promoters (from the promoter SSOT bed) ────────────────────────────
    pb = pd.read_csv(os.path.join(REPO, "references", "pitx2_promoter_regions.bed"),
                     sep="\t", comment="#", header=None)
    pb.columns = (["chrom", "start", "end", "name", "score", "strand"])[:pb.shape[1]]
    prom_ctrl_iv = []   # named promoter-control intervals to keep OUT of the bilateral null
    for iso, label, cls in [("promoter|PITX2a", "PITX2a_prom", "aprom"),
                            ("promoter|PITX2b", "PITX2b_prom", "bprom")]:
        r = pb[pb.name == iso].iloc[0]
        d = int((r.start + r.end) / 2 - TSS)
        rows.append((r.chrom, int(r.start), int(r.end), f"{label}|{cls}"))
        flags.append((label, cls, d, 0, 0))
        prom_ctrl_iv.append((r.chrom, int(r.start), int(r.end)))

    # ── R-desert: HH25 ventral R-sig, distal (−350 kb+), 4 most-significant ──────
    h25 = v[v.stage == "HH25"]
    Rdes = h25[(h25.padj_locus < PADJ) & (h25.log2FoldChange < 0)
               & (h25.dist < -350000)].nsmallest(4, "padj_locus")
    for _, r in Rdes.iterrows():
        c, s, e = coords(r.peak_id)
        rows.append((c, s, e, f"R{int(r.dist // 1000)}kb|Rdesert"))
        flags.append((f"R{int(r.dist // 1000)}kb", "Rdesert", int(r.dist), 0, 0))

    # ── bilateral background: ALL proximal (±BIL_KB) peaks non-sig at BOTH stages ─
    # EXCLUDING any peak overlapping a named promoter control (a/b prom): otherwise the
    # within-locus null double-counts the very promoter element plotted as a positive column
    # (the a-prom g4f peak both anchors the SMAD5/SBE column AND would inflate the bilateral
    # SMAD5 rate it is compared against).  c-prom / +5.4 kb are L-sig → in focal_set.
    def hits_prom_ctrl(pid):
        c, s, e = coords(pid)
        return any(c == pc and s < pe and e > ps for pc, ps, pe in prom_ctrl_iv)

    focal_set = set(focal)
    bil = sorted({p for p in v.peak_id
                  if abs(peak_dist[p]) <= BIL_KB and p not in sig_any
                  and p not in focal_set and not hits_prom_ctrl(p)}, key=lambda p: peak_dist[p])
    for pid in bil:
        d = int(peak_dist[pid])
        c, s, e = coords(pid)
        rows.append((c, s, e, f"bil{d:+d}|bilateral"))
        flags.append((f"bil{d:+d}", "bilateral", d, 0, 0))

    out_bed = os.path.join(DATA, "g4f_domain_motif_windows.bed")
    with open(out_bed, "w") as f:
        f.write("# g4f multi-stage PITX2c cis-block motif windows: contiguous focal L-block "
                "(c-prom/+5.4 broken out) + a/b prom + R-desert (duck-only) + proximal bilateral "
                "background. name=<label>|<class>. 0-based. See make_g4f_domain_windows.py.\n")
        for c, s, e, n in rows:
            f.write(f"{c}\t{s}\t{e}\t{n}\t0\t+\n")

    fl = pd.DataFrame(flags, columns=["window", "class", "dist", "lsig_HH1716", "lsig_HH25"])
    fl.to_csv(os.path.join(DATA, "g4f_domain_window_stageflags.tsv"), sep="\t", index=False)

    print(f"wrote {os.path.relpath(out_bed, REPO)}: {len(rows)} windows")
    print(dict(collections.Counter(n.split('|')[1] for *_, n in rows)))
    blk = fl[fl["class"].isin(["focal", "cprom", "p5p4"])]
    print(f"focal cis-block ({len(blk)} peaks): HH17/16-sig={int(blk.lsig_HH1716.sum())}  "
          f"HH25-sig={int(blk.lsig_HH25.sum())}")
    print(f"bilateral background: {len(bil)} peaks (±{BIL_KB//1000} kb, non-sig both stages)")


if __name__ == "__main__":
    main()
