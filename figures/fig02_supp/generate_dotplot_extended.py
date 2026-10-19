#!/usr/bin/env python3
"""FigS02 Panel D: lineage + signaling + PITX2 + cell-cycle dotplot,
L/R-split for HH13/HH21/HH25, collapsed for HH17/16.

(Was panel E until the 2026-08-24 recompose, which dropped the
per-genotype L/R-inference matrix from the figure — it is now
supplementary Table S4 — and relettered E→D, F→E.)

Folds the old figS02 cell-cycle laterality panel into this one by
adding G1 / S / G2M phase-fraction columns at the right.  The old
L-vs-R total cell count panel is dropped — that information is already
in the per-library cell counts.

Reads
  figures/data/d2g_lpm_subclustering_HH*.tsv.gz   (per-cell markers + leiden)
  figures/data/d2g_signaling_HH*.tsv.gz           (per-cell CP10K, signaling)
  figures/data/d2j_pitx2_isoform_percell.tsv.gz   (per-cell PITX2 isoform UMIs)
  figures/data/d2l_lr_inference.tsv               (per (genotype, library_m) side)
Output
  figures/fig02_supp/panel_D_dotplot_extended.png
"""
import sys
import os

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.dotplot import render_dotplot_staged
from _shared.layout import W_FULL, apply_layout_fonts
from _shared.subcluster_map import (
    assert_display_keys_emitted, display_population, get_leiden_col,
    get_population_clusters)
from _shared.library_naming import TECH_REP, collapse_lib

apply_layout_fonts()

# Fail loud if a publication hedge names a population d2g never emits:
# the import-time assert only proves the key is a POOLS label, not that it fires.
assert_display_keys_emitted()

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Must match stage_gap arg passed to render_dotplot_staged below.  Hoisted
# to module scope so the post-render L/R pair-banding loop and the
# renderer call share a single source of truth instead of two 0.7s drift-
# ing apart.
STAGE_GAP = 0.7

# Minimum per-row cell count required for a row to contribute to the
# per-column max-normalization scale of PITX2 isoform and cell-cycle
# phase columns.  Without this guard, a 2-cell row where both happen to
# be in S phase yields frac=1.0 for S and halves the entire S column's
# biological scale (H5 / M5).
MIN_N_FOR_NORM = 30

# ── Lineage marker columns (existing fig02-A set, BMP4/ALDH1A2 moved out) ──
GENES_LINEAGE = [
    "NKX6_1", "IRX3",                        # dorsal
    "FOXF1", "FOXF2", "GATA4", "ZFPM2",      # splanchnic
    "WT1",                                    # splanchnic / IM boundary
    "TBX5", "TBX4",                           # airway / ventral
    "BARX1",                                  # T-E septum
    "ISL1", "HAND1", "HAND2", "NR2F2",       # LPM / cardiopharyngeal
    "PRRX1", "TWIST1",                        # somatic / pan-mes
]

# ── Signaling marker columns (sidecar from d2g_signaling_export.py) ───────
# WNT2 peaks at HH17/16; WNT2B surfaces at HH21 + HH25 splanchnic (U-method
# UEG rank 1 at HH21 cluster 2, u = 0.31).  Including both keeps the "WNT"
# column header coherent end-to-end across all four stages.
GENES_SIGNALING = [
    "BMP2", "BMP4", "BMPR1B",                # BMP
    "FGF9", "FGF10", "FGFR2",                # FGF
    "WNT2", "WNT2B", "WNT5A", "WNT11",       # WNT
    "RDH10", "ALDH1A2", "CRABP1",            # RA
    "TGFB2",                                  # TGF-β
    "PTCH1",                                  # Hedgehog
    "ESR1", "ESRRG",                         # Estrogen
]

# lncRNA candidate for chick Playrr.  LOC106016802 (NC_051775.1:
# 3,665,204-3,752,463, - strand, 87 kb ncRNA) maps syntenically to the
# chick Playrr / e926 region on NC_052535.1 ~57.0 Mb (GRCg7b; LASTZ
# synteny evidence not distributed); see docs/pitx2_locus_annotation.md §4 for the
# Welsh et al. 2015 positional rationale.  Rendered as PLAYRR* (the
# asterisk flags the positional-candidate status pending sequence
# verification — if verification fails the footnote is downgraded).
GENES_LNCRNA = [
    "LOC106016802",
]

# Display rename (NKX6-1 canonical symbol + Playrr positional candidate)
GENE_DISPLAY = {
    "NKX6_1": "NKX6-1",
    "LOC106016802": "PLAYRR*",
}

PITX2_COLS = ["PITX2", "PITX2a/b", "PITX2c"]
PHASE_COLS = ["G1", "S", "G2M"]


def disp(g):
    return GENE_DISPLAY.get(g, g)


# Final column order (left → right) and bracket grouping
GENE_ORDER = (
    [disp(g) for g in GENES_LINEAGE]
    + [disp(g) for g in GENES_SIGNALING]
    + [disp(g) for g in GENES_LNCRNA]
    + PITX2_COLS
    + PHASE_COLS
)

MARKER_GROUPS = {}
for g in ["NKX6-1", "IRX3"]:
    MARKER_GROUPS[g] = "Dorsal"
for g in ["FOXF1", "FOXF2", "GATA4", "ZFPM2"]:
    MARKER_GROUPS[g] = "Splanchnic"
for g in ["WT1"]:
    MARKER_GROUPS[g] = "Spl./IM"
for g in ["TBX5", "TBX4"]:
    MARKER_GROUPS[g] = "Airway"
for g in ["BARX1"]:
    MARKER_GROUPS[g] = "T-E"
for g in ["ISL1", "HAND1", "HAND2", "NR2F2"]:
    MARKER_GROUPS[g] = "LPM"
for g in ["PRRX1", "TWIST1"]:
    MARKER_GROUPS[g] = "Skeletog."   # markers of the Skeletogenic LPM population (formerly "Somatic")
for g in ["BMP2", "BMP4", "BMPR1B"]:
    MARKER_GROUPS[g] = "BMP"
for g in ["FGF9", "FGF10", "FGFR2"]:
    MARKER_GROUPS[g] = "FGF"
for g in ["WNT2", "WNT2B", "WNT5A", "WNT11"]:
    MARKER_GROUPS[g] = "WNT"
for g in ["RDH10", "ALDH1A2", "CRABP1"]:
    MARKER_GROUPS[g] = "RA"
for g in ["TGFB2", "PTCH1"]:
    MARKER_GROUPS[g] = "TGF-β / Hh"
for g in ["ESR1", "ESRRG"]:
    MARKER_GROUPS[g] = "Estrogen"
for g in [disp(g) for g in GENES_LNCRNA]:
    MARKER_GROUPS[g] = "lncRNA"
for g in PITX2_COLS:
    MARKER_GROUPS[g] = "PITX2"
for g in PHASE_COLS:
    MARKER_GROUPS[g] = "Cell cycle"

# ── Stage / population layout ─────────────────────────────────────────────
STAGE_TSV = {
    "HH13":    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}

SIGNALING_TSV = {
    "HH13":    "d2g_signaling_HH13.tsv.gz",
    "HH17/16": "d2g_signaling_HH1716_HH16.tsv.gz",
    "HH21":    "d2g_signaling_HH21.tsv.gz",
    "HH25":    "d2g_signaling_HH25.tsv.gz",
}

# figS02 extended dotplot — full d2g compendium, including populations
# omitted from fig02 main per docs/lpm_mesenchyme_subcluster_biology.md §4.
# 2026-05-04: "Splanchnic" → "Venous-pole spl." rename (off-axis
# dissection-edge introgression at the cardiac venous pole; exposed here
# for transparency, omitted from fig02 main).
# "Epicardial mes." → "Mesothelium" rename + display added.
STAGE_ROWS = [
    ("HH13",    "Venous-pole spl.",     "Venous-pole spl."),
    ("HH13",    "Early LPM",            "Early LPM"),
    ("HH13",    "Cardiopharyngeal LPM", "Cardiopharyngeal"),
    ("HH13",    "PRRX1+/FOXF1-low mes.", "PRRX1+/FOXF1-low"),  # annotation vet: HH13-only soften (excluded; not skeletal pre-HH17)
    ("HH17/16", "Dorsal-biased LPM",    "Dorsal-biased"),
    ("HH17/16", "Early LPM",            "Early LPM"),
    ("HH17/16", "Venous-pole spl.",     "Venous-pole spl."),
    ("HH17/16", "Cardiopharyngeal LPM", "Cardiopharyngeal"),
    ("HH17/16", "BARX1+",               "BARX1+"),
    ("HH17/16", "Mesothelium",          "Mesothelium"),
    ("HH17/16", "Skeletogenic LPM",     "Skeletogenic"),        # cl2+cl9 (EYA1/EBF/TBX1 earned)
    ("HH17/16", "PRRX1+/FOXF1-low mes.", "PRRX1+/FOXF1-low"),   # cl4 fibroblast (branchiomeric gate)
    ("HH21",    "Dorsal sub-epi.",      "Dorsal sub-epi."),
    ("HH21",    "Ventral spl.",         "Ventral spl."),
    ("HH21",    "Venous-pole spl.",     "Venous-pole spl."),
    ("HH21",    "BARX1+",               "BARX1+"),
    ("HH21",    "Syringeal",            "Syringeal"),
    ("HH21",    "Mesothelium",          "Mesothelium"),
    ("HH21",    "Skeletogenic LPM",     "Skeletogenic"),
    ("HH25",    "Dorsal sub-epi.",      "Dorsal sub-epi."),
    ("HH25",    "Ventral spl.",         "Ventral spl."),
    ("HH25",    "Syringeal",            "Syringeal"),           # cl3 IRX1+ (was "Chondrogenic")
    ("HH25",    "Venous-pole spl.",     "Venous-pole spl."),    # cl10 (was mislabeled Mesothelium)
    ("HH25",    "Skeletogenic LPM",     "Skeletogenic"),
    # 2026-07-02 vet: HH25 "Chondrogenic LPM" (cl3 → Syringeal) + "Mesothelium"
    # (cl8 → Ventral spl., cl10 → Venous-pole spl.) rows dropped.
]

# Stages where L/R is firm enough to split rows (per d2l + figS02D).
# HH17/16 inference is weaker and is collapsed.
SPLIT_LR = {"HH13", "HH21", "HH25"}

# Tech-rep collapse: d2g cell-level files use raw library names (GEX03-1A);
# d2l_lr_inference.tsv uses tech-rep-collapsed library_m (GEX03-1).
# Single source of truth: figures/_shared/library_naming.py (TECH_REP
# + collapse_lib). collapse_lib is what this script calls on line ~347.


def load_side_map():
    """{(genotype, library_m): inferred_side ∈ {L, R, ?}} from d2l."""
    df = pd.read_csv(os.path.join(DATA, "d2l_lr_inference.tsv"), sep="\t")
    return {(r["genotype"], r["library_m"]): r["inferred_side"]
            for _, r in df.iterrows()}


# ── Stat computation ─────────────────────────────────────────────────────
def compute_marker_stats(sub, gene_cols):
    """frac_expressing + mean_raw per gene over a cell subset.

    Raises KeyError if any `gene_cols` entry is missing from the
    dataframe — silently filling (0.0, 0.0) would render a valid-looking
    "not expressed" dot for a gene that simply never got joined into
    the sidecar (e.g., a typo in GENES_LNCRNA, or a signaling-export
    rerun that dropped a column).
    """
    missing = [c for c in gene_cols if c not in sub.columns]
    if missing:
        raise KeyError(
            f"compute_marker_stats: columns {missing} absent from the "
            f"merged cell-level frame.  A missing column silently "
            f"renders as zero expression — rebuild the d2g signaling "
            f"sidecar to match the gene list in this figure script."
        )
    out = []
    n = len(sub)
    for col in gene_cols:
        vals = sub[col].astype(float).values
        if n == 0:
            out.append((col, 0.0, 0.0))
        else:
            out.append((col, float((vals > 0).mean()),
                        float(vals.mean())))
    return out


def compute_isoform_stats(sub, d2j):
    """frac_expressing + mean for PITX2 / PITX2a-b / PITX2c."""
    n = len(sub)
    out = []
    if n > 0 and "PITX2_cp10k" in sub.columns:
        v = sub["PITX2_cp10k"].astype(float).values
        out.append(("PITX2", float((v > 0).mean()), float(v.mean())))
    else:
        out.append(("PITX2", 0.0, 0.0))

    merged = sub[["barcode", "library"]].merge(d2j, on=["barcode", "library"],
                                                how="left")
    merged["pitx2ab_umi"] = merged["pitx2ab_umi"].fillna(0).astype(int)
    merged["pitx2c_umi"] = merged["pitx2c_umi"].fillna(0).astype(int)
    if len(merged) > 0:
        out.append(("PITX2a/b",
                    float((merged["pitx2ab_umi"] > 0).mean()),
                    float(merged["pitx2ab_umi"].mean())))
        out.append(("PITX2c",
                    float((merged["pitx2c_umi"] > 0).mean()),
                    float(merged["pitx2c_umi"].mean())))
    else:
        out.append(("PITX2a/b", 0.0, 0.0))
        out.append(("PITX2c", 0.0, 0.0))
    return out


def compute_phase_stats(sub):
    """Fraction of cells in G1 / S / G2M phases."""
    n = len(sub)
    if n == 0 or "phase" not in sub.columns:
        return [(p, 0.0, 0.0) for p in PHASE_COLS]
    counts = sub["phase"].value_counts()
    out = []
    for p in PHASE_COLS:
        frac = float(counts.get(p, 0)) / n
        out.append((p, frac, frac))   # color encodes same value (Greys)
    return out


def main():
    print("Panel E: lineage + signaling + PITX2 + cell-cycle, L/R-split")

    side_map = load_side_map()
    d2j = pd.read_csv(os.path.join(DATA, "d2j_pitx2_isoform_percell.tsv.gz"),
                      sep="\t")
    print(f"  d2j: {len(d2j)} rows; side_map: {len(side_map)} entries")
    # M6: the failure mode is "d2j regen missed one or more libraries";
    # not "most cells lack PITX2 reads" (d2j is sparse by nature — only
    # PITX2-expressing cells carry UMIs).  A library-level check
    # catches the real failure without firing on normal sparsity.
    d2j_libs = set(d2j["library"].unique())

    # Load per-stage tables.  Merge signaling sidecar onto the lineage table
    # so all gene columns coexist as <gene>_cp10k.
    stage_cells = {}
    for stage, fname in STAGE_TSV.items():
        path = os.path.join(DATA, fname)
        df = pd.read_csv(path, sep="\t")
        sig_path = os.path.join(DATA, SIGNALING_TSV[stage])

        # H3: missing sidecar used to render a silently-empty signaling
        # block.  Fail loudly so a stale pipeline state cannot produce a
        # half-empty panel that looks legitimate.
        if not os.path.exists(sig_path):
            sys.exit(
                f"ERROR: {sig_path} missing — run "
                f"`sbatch pipeline/run_d2gs.sh` "
                f"(d2g_signaling_export.py) first."
            )

        sig = pd.read_csv(sig_path, sep="\t")
        sig_cols = [c for c in sig.columns if c.endswith("_cp10k")]
        keep = ["genotype", "library", "barcode"] + sig_cols
        sig = sig[keep].rename(columns={"genotype": "genotype_sig"})
        # join on (library, barcode); d2g 'genotype' col == d2g sig 'genotype_sig'
        df = df.merge(
            sig, left_on=["library", "barcode"],
            right_on=["library", "barcode"], how="left",
            suffixes=("", "_sig"),
        )

        # H2 canary: (library, barcode) uniquely identifies a physical
        # droplet, so the genotype on both sides of the join must agree.
        # Disagreement ⇒ sidecar regenerated against a different h5ad
        # than the lineage TSV — the silent-wrong-results failure mode
        # Track B is designed to catch.
        mism = (df["genotype"] != df["genotype_sig"]).sum()
        if mism:
            raise RuntimeError(
                f"signaling sidecar genotype mismatch on {mism} "
                f"rows for stage {stage} — sidecar stale vs lineage TSV"
            )
        df = df.drop(columns=["genotype_sig"])

        # H3 NaN check: after M3's small-libs drop the sidecar is a
        # strict superset of the lineage cells, so a NaN in any
        # signaling column means a regen mismatch — fail rather than
        # silently propagate NaN → mean() → 0 "zero expression".
        miss = df[sig_cols].isna().any(axis=1).sum()
        if miss:
            raise RuntimeError(
                f"signaling sidecar misses {miss} d2g cells at {stage} "
                f"— regenerate the sidecar to match current d2g cell set"
            )

        # If a gene exists in both (BMP4, ALDH1A2), prefer the sidecar
        for col in sig_cols:
            if col in df.columns and (col + "_sig") in df.columns:
                df[col] = df[col + "_sig"].combine_first(df[col])
                df = df.drop(columns=[col + "_sig"])
        print(f"  Loaded {fname} + signaling sidecar: {len(df)} cells")

        # M6: every library present in this stage's d2g must appear in
        # d2j.  A missing library silently renders as "zero PITX2
        # isoform expression" for every cell — an easy-to-miss regen
        # mismatch (d2j covers 12 libraries today; a stale d2j would
        # zero-out whole panels of rows).
        stage_libs = set(df["library"].unique())
        missing_in_d2j = stage_libs - d2j_libs
        if missing_in_d2j:
            raise RuntimeError(
                f"{stage}: libraries present in d2g but missing from "
                f"d2j: {sorted(missing_in_d2j)} — rerun "
                f"pipeline/d2j_pitx2_isoform_percell.py"
            )

        # Per-cell side from (genotype, collapsed library)
        df["lib_m"] = df["library"].map(collapse_lib)
        df["side"] = [side_map.get((g, l), "?")
                      for g, l in zip(df["genotype"], df["lib_m"])]

        # H4a: every library_m after collapse should exist in d2l's
        # inference output.  A new tech-rep not yet added to TECH_REP
        # (e.g. GEX09-1A alongside GEX09-1B) would fail to collapse and
        # the raw name would be missing from d2l's library_m universe,
        # silently routing every cell of that library into side=?.
        known_libs_m = {lib_m for (_, lib_m) in side_map.keys()}
        present_libs_m = set(df["lib_m"].unique())
        unknown = present_libs_m - known_libs_m
        if unknown:
            raise RuntimeError(
                f"{stage}: libraries present in d2g but missing from "
                f"d2l: {sorted(unknown)} — rerun "
                f"pipeline/d2l_lr_inference.py or add the new raw "
                f"name to TECH_REP in figures/_shared/library_naming.py"
            )

        # H4b: surface any ?-side dropouts on L/R-split stages.  Today
        # only 2 cells total; tighten to hard error if this grows into
        # material low-N slot loss.
        if stage in SPLIT_LR:
            n_dropped = (df["side"] == "?").sum()
            if n_dropped:
                print(f"  {stage}: {n_dropped} cells with side=? "
                      f"(dropped from L/R rows)")

        stage_cells[stage] = df

    # ── Build the row layout (interleaved L/R per population) ────────
    all_records = []
    group_order = []
    stage_groups = []   # (bracket_label, n_rows)
    pair_ids = []       # per-emitted-row: (stage_idx, pop_idx_within_stage)
                        # — used to band each (population × {L,R}) pair
    n_by_row = {}       # row_label → n_cells, used for M5 norm guard

    # Map gene name → cp10k column (for the marker-stat helper).  The
    # d2g exporter rewrites LOC113843492 (duck NKX6-1) to NKX6_1_cp10k,
    # so a single `{gene}_cp10k` lookup now covers every gene in the
    # lineage + signaling panels (L5: the special-case branch for
    # NKX6_1 was dead code).
    def gene_col(g):
        return f"{g}_cp10k"

    # Per-stage: which sides to emit (interleaved L,R for L/R-split stages;
    # single None for HH17/16).  Per-population the rows appear adjacent so
    # the eye can compare L vs R directly within a population.
    def sides_for(stage):
        if stage in SPLIT_LR:
            return [("L", " L"), ("R", " R")]
        return [(None, "")]

    # Walk stage by stage so brackets stay contiguous
    for stage_idx, stage in enumerate(["HH13", "HH17/16", "HH21", "HH25"]):
        df = stage_cells[stage]
        leiden_col = get_leiden_col(stage)
        rows_for_stage = [r for r in STAGE_ROWS if r[0] == stage]

        n_rows_this_stage = 0
        for pop_idx, (_, pop, label) in enumerate(rows_for_stage):
            # Hedged SSOT names override the figure-local shortening;
            # clusters still resolve from the SSOT `pop` string.
            label = display_population(pop, short=True, unhedged=label)
            cluster_ids = get_population_clusters(stage, pop)
            cluster_str = [str(c) for c in cluster_ids]

            for side, side_suffix in sides_for(stage):
                if side is not None:
                    stage_df = df[df["side"] == side]
                else:
                    stage_df = df
                sub = stage_df[stage_df[leiden_col].astype(str).isin(cluster_str)]
                n = len(sub)
                if n == 0:
                    continue   # skip empty rows entirely

                row_label = f"{label}{side_suffix} ({n:,})"

                # Lineage + signaling + lncRNA markers (all per-cell
                # CP10K, same dot encoding).  lncRNA block carries
                # LOC106016802 (PLAYRR*) from the d2g signaling sidecar.
                stats = compute_marker_stats(
                    sub,
                    [gene_col(g) for g in GENES_LINEAGE
                                        + GENES_SIGNALING
                                        + GENES_LNCRNA],
                )
                # PITX2 isoforms
                iso = compute_isoform_stats(sub, d2j)
                # Cell-cycle phases
                phase = compute_phase_stats(sub)

                # Map (col_name, frac, mean) → gene-display name + records
                gene_records = []
                for (col, frac, mean) in stats:
                    raw_gene = col.replace("_cp10k", "")
                    gname = disp(raw_gene)
                    gene_records.append({
                        "group": row_label, "gene": gname,
                        "frac_expressing": round(frac, 4),
                        "mean_raw": round(mean, 4),
                    })
                for (gname, frac, mean) in iso:
                    gene_records.append({
                        "group": row_label, "gene": gname,
                        "frac_expressing": round(frac, 4),
                        "mean_raw": round(mean, 4),
                    })
                for (gname, frac, mean) in phase:
                    gene_records.append({
                        "group": row_label, "gene": gname,
                        "frac_expressing": round(frac, 4),
                        "mean_raw": round(mean, 4),
                    })
                all_records.extend(gene_records)
                group_order.append(row_label)
                pair_ids.append((stage_idx, pop_idx))
                n_by_row[row_label] = n
                n_rows_this_stage += 1

        if n_rows_this_stage > 0:
            stage_groups.append((stage, n_rows_this_stage))

    dotplot_df = pd.DataFrame(all_records)
    print(f"  Built {len(group_order)} rows across "
          f"{len(stage_groups)} bracket blocks")

    # Max-normalize per-column for fractions that exceed dot_max=0.5:
    # PITX2 isoforms (long-standing) and the cell-cycle phase columns
    # (G1 in particular routinely > 50%).
    #
    # M5 guard: rows with fewer than MIN_N_FOR_NORM cells are excluded
    # from the scale-computing max() so a 2-cell outlier that lands
    # entirely in one phase cannot halve the biological scale of the
    # whole column.  The rows themselves are still rendered — they just
    # don't drive normalization.
    NORM_GENES = ["PITX2a/b", "PITX2c"] + PHASE_COLS
    for g in NORM_GENES:
        mask = dotplot_df["gene"] == g
        eligible = dotplot_df["group"].map(
            lambda r: n_by_row.get(r, 0) >= MIN_N_FOR_NORM
        )
        max_frac = dotplot_df.loc[mask & eligible, "frac_expressing"].max()
        if max_frac and max_frac > 0:
            scale = 0.5 / max_frac
            dotplot_df.loc[mask, "frac_expressing"] = (
                dotplot_df.loc[mask, "frac_expressing"] * scale
            )
            dotplot_df.loc[mask, "mean_raw"] = (
                dotplot_df.loc[mask, "mean_raw"] * scale
            )
            contributors = (dotplot_df.loc[mask & eligible]
                            .nlargest(3, "frac_expressing")["group"].tolist())
            print(f"  {g}: max frac {max_frac:.4f} (from ≥{MIN_N_FOR_NORM}-cell "
                  f"rows: {contributors}) → max-normalised")

    tsv_path = os.path.join(DATA, "d2g_dotplot_spatial_extended.tsv")
    dotplot_df.to_csv(tsv_path, sep="\t", index=False)
    print(f"  Saved intermediate: {tsv_path} ({len(dotplot_df)} rows)")

    # Render — height grows with row count
    n_rows_total = sum(n for _, n in stage_groups)
    height = max(7.0, 0.32 * n_rows_total + 2.0)
    out = os.path.join(OUTDIR, "panel_D_dotplot_extended.png")
    fig, ax = render_dotplot_staged(
        tsv_path,
        group_col="group", gene_col="gene",
        frac_col="frac_expressing", mean_col="mean_raw",
        group_order=group_order,
        marker_order=GENE_ORDER,
        marker_groups=MARKER_GROUPS,
        row_stage_groups=stage_groups,
        stage_gap=STAGE_GAP,
        cmap="Greys",
        dot_max=0.5,
        figsize=(W_FULL, height),
        save_path=None,
    )

    # Vertical separators between the five major column blocks
    # (Lineage | Signaling | lncRNA | PITX2 | Cell cycle).  Slightly
    # thicker + darker than the within-block grid.
    for sep_gene in ["BMP2", "PLAYRR*", "PITX2", "G1"]:
        if sep_gene in GENE_ORDER:
            idx = GENE_ORDER.index(sep_gene)
            ax.axvline(x=idx - 0.5, color="#444444", linestyle="--",
                       linewidth=0.9, alpha=0.85, zorder=2)

    # ── L/R pair banding ────────────────────────────────────────────
    # Recompute y_positions using the same logic as render_dotplot_staged
    # so the bands sit exactly under each pair's two rows.  Then shade
    # alternating populations within each stage with a faint gray band so
    # each (population × {L, R}) pair reads as one visual unit (user
    # request 2026-04-22).  STAGE_GAP hoisted to module scope (L6) so
    # renderer call + banding loop share one source of truth.
    y_positions = []
    y = 0
    for si, (_, n_rows) in enumerate(stage_groups):
        for _ in range(n_rows):
            y_positions.append(y)
            y += 1
        if si < len(stage_groups) - 1:
            y += STAGE_GAP

    # Group consecutive rows by (stage_idx, pop_idx) → a "pair" unit
    pair_runs = []   # list of (start_y, end_y, stage_idx, pop_idx)
    if pair_ids:
        run_start = 0
        for i in range(1, len(pair_ids) + 1):
            if i == len(pair_ids) or pair_ids[i] != pair_ids[run_start]:
                stage_idx, pop_idx = pair_ids[run_start]
                pair_runs.append((y_positions[run_start],
                                  y_positions[i - 1],
                                  stage_idx, pop_idx))
                run_start = i

    # Alternate band fill per stage: every other pop_idx within a stage
    # gets a faint gray band.  Skip the no-band side (white background).
    for top_y, bot_y, stage_idx, pop_idx in pair_runs:
        if pop_idx % 2 == 1:
            # Half-row pad so the band hugs the dot row tightly without
            # bleeding into the next pair.
            ax.axhspan(top_y - 0.5, bot_y + 0.5,
                       facecolor="#eaeaea", alpha=0.55,
                       linewidth=0, zorder=0)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
