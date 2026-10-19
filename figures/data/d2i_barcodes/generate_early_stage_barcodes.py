#!/usr/bin/env python3
"""Generate per-library barcode files for d2i PITX2 pileup populations.

Covers ALL d2g-subclustered stages (HH13 / HH17/16 / HH21 / HH25). Cluster
IDs are resolved from `references/d2g_subcluster_annotations.tsv` via the
shared subcluster_map helper — never hardcoded here, because leiden numbers
are not seed-stable across d2g reruns.

Ventral-splanchnic populations combine `Ventral spl.` + `BARX1+` T-E
septum because both sit ventrally and are folded together in the paired
laterality analysis (Fig02-C).

Output: per-library barcode files under this directory + manifest.tsv.
"""

import os
import sys
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
D2G = os.path.join(REPO, "figures", "data")
OUTDIR = os.path.join(REPO, "figures", "data", "d2i_barcodes")

sys.path.insert(0, os.path.join(REPO, "figures"))
from _shared.subcluster_map import get_leiden_col, get_population_clusters


def _clusters(stage, populations, allow_empty=frozenset()):
    """Union of cluster IDs across one or more population labels.

    Fail-fast: raises if any individual label resolves to []. Pass
    `allow_empty={"BARX1+"}` (etc.) when a label is known-absent at
    a stage (BARX1+ peaks HH17–HH21, absent at HH25 where the
    population has dispersed).

    The fail-fast guard caught the `emerging_dorsal_HH1716` regression
    on 2026-05-04 when post-rename the lookup silently returned 0
    cells.  Same regression class.
    """
    out = []
    for pop in populations:
        cls = get_population_clusters(stage, pop)
        if not cls and pop not in allow_empty:
            raise RuntimeError(
                f"_clusters({stage!r}, {populations!r}): label {pop!r} "
                f"resolves to no clusters in references/d2g_subcluster_annotations.tsv. "
                f"Pass allow_empty={{{pop!r}}} if the label is known-absent at this stage, "
                f"or update d2g_validate_subclusters.py / POPULATIONS dict."
            )
        out.extend(cls)
    return sorted(set(out))


# Shape the pipeline: (tsv, {pop_name: (leiden_col, cluster_ids_or_None)})
# Covers ALL d2g stages so barcode files track current clustering. Every
# pop_name here must also have a matching row in the d2i pileup manifest
# (names are stable once committed — cluster IDs behind them are not).
# Variant must match `pipeline/d2g_validate_subclusters.py::VARIANT` ("ccreg"),
# so cluster IDs resolved from the SSOT TSV align with the underlying cells.
# Non-ccreg variants share only 1.6%-49% of leiden_0.3 assignments at HH13/HH17/16/HH21/HH25.
POPULATIONS = {
    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz": {
        # 2026-05-04: HH13 had no Ventral spl. (D/V unresolved pre-TES).
        # Pileup track now uses Early LPM + Cardiopharyngeal LPM (the HH13
        # inference pool members per POOLS).  Venous-pole spl. (cl2) is
        # excluded as off-axis introgression.
        "splanchnic_HH13":        (get_leiden_col("HH13"),
                                   _clusters("HH13", ["Early LPM",
                                                       "Cardiopharyngeal LPM"])),
        # HH13 has no "Skeletogenic LPM" label post-split: the PRRX1-marked,
        # FOXF1-low non-splanchnic mesenchyme keeps its own label there
        # (skeletal TFs absent until HH17/16). pop_suffix stays "somatic_HH13"
        # for column-name stability in the d2i manifest.
        "somatic_HH13":           (get_leiden_col("HH13"),
                                   _clusters("HH13", ["PRRX1+/FOXF1-low mes."])),
        "all_lpm_HH13":           (get_leiden_col("HH13"), None),
    },
    "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz": {
        # 2026-05-04: HH17/16 inference pool (Early LPM + Cardiopharyngeal
        # LPM + BARX1+) — pre-TES, D/V unresolved.  Venous-pole spl. (cl8)
        # is excluded as off-axis introgression.  Old definition included
        # "Ventral spl." which doesn't exist at HH17/16 (resolved post-TES).
        # 2026-08-09 (seeded d2g rerun): HH17/16 emits no "Early LPM"
        # cluster; the inference pool is Cardiopharyngeal + BARX1+ only (see
        # POOLS in figures/_shared/subcluster_map.py).
        "splanchnic_HH1716":      (get_leiden_col("HH17/16"),
                                   _clusters("HH17/16", ["Cardiopharyngeal LPM",
                                                          "BARX1+"])),
        "emerging_dorsal_HH1716": (get_leiden_col("HH17/16"),
                                   _clusters("HH17/16", ["Dorsal-biased LPM"])),
        "barx1_pos_HH1716":       (get_leiden_col("HH17/16"),
                                   _clusters("HH17/16", ["BARX1+"])),
        "somatic_HH1716":         (get_leiden_col("HH17/16"),
                                   _clusters("HH17/16", ["Skeletogenic LPM"])),
        "all_lpm_HH1716":         (get_leiden_col("HH17/16"), None),
    },
    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz": {
        "dorsal_subepi_HH21":     (get_leiden_col("HH21"),
                                   _clusters("HH21", ["Dorsal sub-epi."])),
        # ventral_splanchnic = ventral splanchnic core + BARX1+ (ventral
        # position, folded together in paired-laterality)
        "ventral_splanchnic_HH21":(get_leiden_col("HH21"),
                                   _clusters("HH21", ["Ventral spl.", "BARX1+"])),
        "ventral_only_HH21":      (get_leiden_col("HH21"),
                                   _clusters("HH21", ["Ventral spl."])),
        "barx1_pos_HH21":         (get_leiden_col("HH21"),
                                   _clusters("HH21", ["BARX1+"])),
        "syringeal_HH21":         (get_leiden_col("HH21"),
                                   _clusters("HH21", ["Syringeal"])),
        "somatic_HH21":           (get_leiden_col("HH21"),
                                   _clusters("HH21", ["Skeletogenic LPM"])),
    },
    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz": {
        # BARX1+ is HH17/16-HH21 only — declines by HH25 (Woo 2011 mouse window).
        # ventral_splanchnic_HH25 union with BARX1+ resolves to just Ventral spl.
        # at this stage; explicit allow_empty={"BARX1+"} so the union doesn't
        # silently mask future regressions.
        "dorsal_subepi_HH25":     (get_leiden_col("HH25"),
                                   _clusters("HH25", ["Dorsal sub-epi."])),
        "ventral_splanchnic_HH25":(get_leiden_col("HH25"),
                                   _clusters("HH25", ["Ventral spl.", "BARX1+"],
                                             allow_empty={"BARX1+"})),
        "syringeal_HH25":         (get_leiden_col("HH25"),
                                   _clusters("HH25", ["Syringeal"])),
        "somatic_HH25":           (get_leiden_col("HH25"),
                                   _clusters("HH25", ["Skeletogenic LPM"])),
    },
}


def main():
    manifest_lines = []

    for d2g_file, pops in POPULATIONS.items():
        path = os.path.join(D2G, d2g_file)
        print(f"\nLoading {d2g_file}...")
        df = pd.read_csv(path, sep="\t")
        print(f"  {len(df):,} cells, {df['library'].nunique()} libraries")

        for pop_name, (leiden_col, cluster_ids) in pops.items():
            # Filter cells for this population
            if cluster_ids is None:
                subset = df
            else:
                subset = df[df[leiden_col].isin(cluster_ids)]

            print(f"\n  Population: {pop_name}")
            print(f"    Cells: {len(subset):,}")

            # Split by library and write barcode files
            for lib in sorted(subset["library"].unique()):
                lib_cells = subset[subset["library"] == lib]
                n = len(lib_cells)
                if n == 0:
                    continue

                # Barcode format: 16-mer + "-1" suffix
                barcodes = lib_cells["barcode"].tolist()
                # Ensure -1 suffix
                barcodes = [
                    bc if "-" in bc else f"{bc}-1"
                    for bc in barcodes
                ]

                fname = f"{pop_name}_{lib}.txt"
                fpath = os.path.join(OUTDIR, fname)
                with open(fpath, "w") as f:
                    for bc in sorted(barcodes):
                        f.write(f"{bc}\n")

                manifest_lines.append(
                    f"{pop_name}\t{lib}\t{n}\t{fname}"
                )
                print(f"    {lib}: {n} barcodes -> {fname}")

    # Write the complete manifest (covers every population generated here).
    manifest_path = os.path.join(OUTDIR, "manifest.tsv")
    with open(manifest_path, "w") as f:
        f.write("population\tlibrary\tn_barcodes\tfilename\n")
        for line in manifest_lines:
            f.write(line + "\n")
    print(f"\nWrote {manifest_path} ({len(manifest_lines)} rows)")

    # Remove orphan .txt files no longer produced by the current POPULATIONS dict
    # (e.g., a stage drops a population, or a library has zero cells now).
    # Directory must match the manifest 1:1 so consumers globbing *.txt agree
    # with consumers reading the manifest.
    current = {line.split("\t")[3] for line in manifest_lines}
    orphans = sorted(
        f for f in os.listdir(OUTDIR)
        if f.endswith(".txt") and f not in current
    )
    for fname in orphans:
        os.remove(os.path.join(OUTDIR, fname))
    if orphans:
        print(f"Removed {len(orphans)} orphan barcode file(s):")
        for f in orphans:
            print(f"  - {f}")


if __name__ == "__main__":
    main()
