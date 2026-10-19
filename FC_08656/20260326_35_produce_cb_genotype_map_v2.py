#!/usr/bin/env python3
"""Produce barcode->genotype group mapping for all 12 GEX libraries (v2).

Updates from v1:
- Adds G3 (GEX04 + GEX0S-duck) from vireo N=11 best seed (multi-seed, 100% consistency)
- Adds dp500 donor29 (25-8) rescue for cells fragmented in libpair
- Updated dedup: libpair preferred for paired donors, dp500 for 25-8
- Fixed: dp500 only wins in GEX05 (25-8's only library by design).
  Previously dp500 had priority 0 in PILOTGEX3, stealing ~839 barcodes
  from correct libpair genotypes (25-5..7, 13-5..8).

v2.1 update:
- G2 now uses OLD G2 vireo (GEX-biased variant discovery, n_donor_16, 4 chunks)
  instead of tri_dp50.  The old vireo correctly separates 17-KC-1 (g11) which
  was merged into donor10 in the tri_dp50 run due to the universal variant set
  lacking KC-discriminative positions.
- Multi-chunk assembly: reads all 4 chunk best-seeds, maps chunk-level donors
  to unified genotype names via genotype_groups.json from the QC directory.

Output: cb_genotype_map_v2.tsv.gz with columns:
  cb          - 16bp cell barcode (no suffix)
  library     - library identity (e.g. GEX05, GEX03-1AB, GEX04)
  suffix      - barcode suffix in the cellSNP/vireo context
  demux_group - which vireo run assigned this cell
  donor_id    - vireo donor ID within that run
  embryo_group - suffix-based embryo group
  status      - singlet, doublet, or unassigned

Run on O2:
    module load gcc/14.2.0 conda/miniforge3 && conda activate vireo
    cd /n/scratch/users/c/cl266/FC_08656
    python3 20260326_35_produce_cb_genotype_map_v2.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import pandas as pd
import json
import csv
from collections import defaultdict
import os

SCRATCH = SCRATCH_ROOT
OUTBASE = f"{SCRATCH}/universal_genotyping"
OUTDIR = f"{OUTBASE}/identity_table"
os.makedirs(OUTDIR, exist_ok=True)

# ── Canonical vireo result paths (must match vireo_provenance.tsv) ─────
# When a better seed is found, update BOTH this script AND vireo_provenance.tsv.
# The g9_g6_mapping.tsv builder (b2_g9_concordance.py) must use the same G9 path.
VIREO_PATHS = {
    "G6_libpair": f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv",
    "G6_dp500":   f"{OUTBASE}/group6/vireo_filt_dp500/donor_ids.tsv",
    "G9_dp100":   f"{SCRATCH}/vireo_seeded/g9_dp100/seed_126746/donor_ids.tsv",
    "G3":         f"{SCRATCH}/vireo_seeded/g3_denovo/seed_46/donor_ids.tsv",
}

# Old G2 vireo paths (GEX-biased variant discovery, chunked)
G2_OLD_VIREO_BASE = (
    f"{SCRATCH}/GEX03-1A_GEX/outs/pseudobulk_group2_merged/"
    "vireo_output/n_donor_16_40000"
)
G2_OLD_QC_DIR = f"{G2_OLD_VIREO_BASE}/demux_qc"

# ============================================================
# Suffix -> library mapping per demux group
# ============================================================
SUFFIX_TO_LIBRARY_G6 = {
    '-5': 'GEX05', '-6': 'GEX06', '-7': 'GEX07', '-8': 'GEX08',
    '-2': 'PILOTGEX4B', '-3': 'PILOTGEX3',
}
SUFFIX_TO_LIBRARY_G9 = {
    '-1': 'PILOTGEX4A', '-2': 'PILOTGEX4B', '-3': 'PILOTGEX3',
}
SUFFIX_TO_LIBRARY_G2 = {
    '-1': 'GEX03-1AB', '-2': 'GEX03-1AB', '-3': 'GEX03-2',
}
SUFFIX_TO_LIBRARY_G3 = {
    '-4': 'GEX04', '-S': 'GEX0S-duck',
}

# Embryo group classification from suffix pattern (G6/G9 context)
SUFFIX_TO_EMBRYO = {
    frozenset({'-5', '-6'}): "21-5..8 or 16-6/7",
    frozenset({'-5', '-3'}): "25-5..7",
    frozenset({'-5'}): "25-5..7 or 25-8",
    frozenset({'-6', '-2'}): "25-11..14",
    frozenset({'-7', '-8'}): "17-9..12 or 16-8/9",
    frozenset({'-7', '-3'}): "13-5..8",
    frozenset({'-7'}): "13-5..8",
    frozenset({'-8', '-2'}): "13-11..14",
    frozenset({'-8'}): "17-Md-13",
    frozenset({'-2'}): "17-A..D",
    frozenset({'-3'}): "PILOT3",
}

# G9 cross-group embryo links (from AF correlation r>0.93)
G9_TO_EMBRYO = {
    'donor10': '13-5..8', 'donor6': '13-5..8', 'donor16': '13-5..8',
    'donor12': '25-5..7', 'donor13': '25-5..7', 'donor9': '25-5..7',
    'donor15': '25-11..14', 'donor14': '25-11..14', 'donor17': '25-11..14',
    'donor5': '17-A..D',
    'donor3': '13-11..14', 'donor8': '13-11..14', 'donor7': '13-11..14',
}

# G3 exclusions: chick contaminant (confirmed by W/Z ratio: pct_z=0.030 vs ~0.065 duck)
# De novo variant run (vireo_N11_denovo_best): chick = donor6
# Previous 665K universal run (vireo_N11_best): chick = donor7
G3_EXCLUDE = {'donor6'}


def get_embryo_group(suffix_counts):
    """Classify donor by significant suffixes."""
    sig = {s for s, c in suffix_counts.items() if c > 5}
    if not sig:
        sig = set(suffix_counts.keys())
    return SUFFIX_TO_EMBRYO.get(frozenset(sig), "unknown")


def strip_suffix(barcode):
    """Extract 16bp CB from suffixed barcode."""
    # Handles both numeric (-5) and letter (-S) suffixes
    return barcode.rsplit('-', 1)[0]


def get_suffix(barcode):
    """Get suffix (e.g., -5, -S) from barcode."""
    return '-' + barcode.rsplit('-', 1)[-1]


def classify_status(donor_id):
    if donor_id == 'doublet':
        return 'doublet'
    elif donor_id == 'unassigned':
        return 'unassigned'
    return 'singlet'


def load_g2_old_vireo(base_dir, qc_dir):
    """Assemble old G2 vireo result from multi-chunk best seeds.

    Reads genotype_groups.json + seed_ranking.tsv from the QC directory,
    then loads donor_ids.tsv from each chunk's best seed.  Translates
    chunk-level donor names to unified genotype names (genotype_01, etc.)
    and deduplicates cells that appear in overlapping chunks.

    Returns a DataFrame with columns: cell, donor_id (=genotype name or
    doublet/unassigned).
    """
    # 1. Load genotype groups: maps genotype name -> merged_pattern
    #    merged_pattern[i] = chunk-level donor name for chunk i (or null)
    gg_path = os.path.join(qc_dir, 'genotype_groups.json')
    if not os.path.isfile(gg_path):
        raise FileNotFoundError(f"genotype_groups.json not found: {gg_path}")
    with open(gg_path) as f:
        genotype_groups = json.load(f)
    print(f"  Loaded {len(genotype_groups)} genotype groups from QC")

    # 2. Load seed ranking: get best seed per chunk
    sr_path = os.path.join(qc_dir, 'seed_ranking.tsv')
    if not os.path.isfile(sr_path):
        raise FileNotFoundError(f"seed_ranking.tsv not found: {sr_path}")
    best_seeds = {}
    with open(sr_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chunk_raw = row['chunk_raw']
            if chunk_raw not in best_seeds:
                best_seeds[chunk_raw] = row['seed']
    print(f"  Best seeds: { {k: v for k, v in best_seeds.items()} }")

    # 3. Determine chunk ordering from genotype_groups merged_pattern length
    first_group = next(iter(genotype_groups.values()))
    n_chunks = len(first_group['merged_pattern'])

    # Chunk names in genotype_groups order (must match seed_ranking chunk_raw)
    # Convention: chunk_0_overlap50, chunk_1_overlap50, ..., chunk_edge_NNNNN
    chunk_names = sorted(best_seeds.keys())
    # Reorder: numeric chunks first, edge last
    numeric_chunks = sorted([c for c in chunk_names if 'edge' not in c])
    edge_chunks = [c for c in chunk_names if 'edge' in c]
    chunk_order = numeric_chunks + edge_chunks
    assert len(chunk_order) == n_chunks, (
        f"Mismatch: {len(chunk_order)} chunks vs "
        f"{n_chunks} positions in merged_pattern"
    )
    print(f"  Chunk order ({n_chunks} chunks): {chunk_order}")

    # 4. Build reverse map: (chunk_index, donor_name) -> genotype_name
    chunk_donor_to_geno = {}
    for geno_name, info in genotype_groups.items():
        for ci, donor_name in enumerate(info['merged_pattern']):
            if donor_name is not None and not pd.isna(donor_name):
                chunk_donor_to_geno[(ci, str(donor_name))] = geno_name

    # 5. Load donor_ids from each chunk's best seed and translate
    all_cells = {}  # barcode -> (genotype_or_status, chunk_index)
    for ci, chunk_name in enumerate(chunk_order):
        seed = best_seeds[chunk_name]
        donor_ids_path = os.path.join(
            base_dir, chunk_name, f"seed_{seed}", "donor_ids.tsv"
        )
        if not os.path.isfile(donor_ids_path):
            print(f"  WARNING: missing {donor_ids_path}, skipping chunk")
            continue

        chunk_df = pd.read_csv(donor_ids_path, sep='\t')
        n_chunk = len(chunk_df)
        n_mapped = 0

        for _, row in chunk_df.iterrows():
            barcode = row['cell']
            raw_donor = str(row['donor_id'])

            if raw_donor in ('doublet', 'unassigned'):
                mapped = raw_donor
            else:
                key = (ci, raw_donor)
                mapped = chunk_donor_to_geno.get(key, f"unmapped_{raw_donor}")
                if not mapped.startswith('unmapped'):
                    n_mapped += 1

            # Dedup: prefer mapped genotype over doublet/unassigned;
            # for cells in multiple chunks, keep first mapped assignment
            existing = all_cells.get(barcode)
            if existing is None:
                all_cells[barcode] = mapped
            elif existing in ('doublet', 'unassigned') or existing.startswith('unmapped'):
                # Upgrade to a real assignment if available
                if mapped not in ('doublet', 'unassigned') and not mapped.startswith('unmapped'):
                    all_cells[barcode] = mapped

        print(f"  {chunk_name}/seed_{seed}: {n_chunk:,} cells, "
              f"{n_mapped:,} mapped to genotypes")

    # 6. Build output DataFrame
    records = []
    for barcode, donor in all_cells.items():
        records.append({'cell': barcode, 'donor_id': donor})
    result = pd.DataFrame(records)
    print(f"  Old G2 vireo assembled: {len(result):,} unique cells, "
          f"{result['donor_id'].nunique()} distinct donor values")
    return result


# ============================================================
# Load vireo results
# ============================================================
print("Loading vireo results...\n")

rows = []

# --- G6 libpair reference (PRIMARY for GEX05-08 + PILOT4B + PILOT3) ---
print("G6 libpair reference...")
g6_lp = pd.read_csv(f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv",
                     sep='\t')

g6_donor_sfx = defaultdict(lambda: defaultdict(int))
for _, row in g6_lp.iterrows():
    if row['donor_id'] not in ('doublet', 'unassigned'):
        g6_donor_sfx[row['donor_id']][get_suffix(row['cell'])] += 1

g6_embryo = {d: get_embryo_group(sfx) for d, sfx in g6_donor_sfx.items()}

for _, row in g6_lp.iterrows():
    bc = row['cell']
    donor = row['donor_id']
    sfx = get_suffix(bc)
    lib = SUFFIX_TO_LIBRARY_G6.get(sfx, f"unknown{sfx}")
    status = classify_status(donor)
    rows.append({
        'cb': strip_suffix(bc), 'library': lib, 'suffix': sfx,
        'demux_group': 'G6_libpair', 'donor_id': donor,
        'embryo_group': g6_embryo.get(donor, '') if status == 'singlet' else '',
        'status': status,
    })
n_g6 = len(rows)
print(f"  {n_g6:,} cells")

# --- G6 combined HH13 recovery (PG3+GEX07 vireo, best seed) ---
# Overrides libpair for GEX07 and PILOTGEX3 cells.
# The combined vireo merges PG3+GEX07 in a single run, resolving
# clutch-mates (g32/g35) that the libpair runB couldn't separate.
HH13_RECOVERY_DIR = f"{SCRATCH}/hh13_recovery/13-5..8/vireo_N15_seeds"
HH13_DONOR_MAPPING = f"{SCRATCH}/hh13_recovery/13-5..8/recovered_barcode_map.tsv"

# Auto-detect best seed by ELBO
import re, glob
_best_seed_dir = None
_best_elbo = -float("inf")
for _s in range(50):
    _sdir = f"{HH13_RECOVERY_DIR}/seed_{_s}"
    if not os.path.isfile(f"{_sdir}/donor_ids.tsv"):
        continue
    _logs = glob.glob(f"{HH13_RECOVERY_DIR}/seed_{_s}_*.log")
    for _log in _logs:
        with open(_log) as _f:
            for _line in _f:
                _m = re.search(r"lower bound ranges \[([^\]]+)\]", _line)
                if _m:
                    _elbo = max(float(x.strip()) for x in _m.group(1).split(","))
                    if _elbo > _best_elbo:
                        _best_elbo = _elbo
                        _best_seed_dir = _sdir

if _best_seed_dir and os.path.isfile(HH13_DONOR_MAPPING):
    print(f"G6 combined HH13 recovery (best seed: {os.path.basename(_best_seed_dir)}, "
          f"ELBO={_best_elbo:.1f})...")

    # Load donor mapping: recovery_donor -> mapped_genotype
    _rec_map_df = pd.read_csv(HH13_DONOR_MAPPING, sep="\t")
    _rec_map_df["suffix"] = _rec_map_df["suffix"].astype(str)
    _rec_map_df.loc[~_rec_map_df["suffix"].str.startswith("-"), "suffix"] = \
        "-" + _rec_map_df.loc[~_rec_map_df["suffix"].str.startswith("-"), "suffix"]
    # Build recovery_donor -> genotype from the mapping
    _donor_to_geno = _rec_map_df.groupby("recovery_donor")["mapped_genotype"].first().to_dict()
    print(f"  Donor->genotype: {_donor_to_geno}")

    # Load best seed donor_ids.tsv
    _hh13_df = pd.read_csv(f"{_best_seed_dir}/donor_ids.tsv", sep="\t")

    # Suffix mapping for combined vireo barcodes
    HH13_SUFFIX_TO_LIB = {"-7": "GEX07", "-3": "PILOTGEX3"}

    n_hh13 = 0
    for _, row in _hh13_df.iterrows():
        bc = row["cell"]
        donor = row["donor_id"]
        sfx = get_suffix(bc)
        lib = HH13_SUFFIX_TO_LIB.get(sfx)
        if lib is None:
            continue  # skip non-target suffixes

        status = classify_status(donor)
        # Map recovery donor to existing genotype ID
        mapped_geno = _donor_to_geno.get(donor) if status == "singlet" else None

        # Only emit rows for donors with known genotype mappings (HH13).
        # Non-HH13 donors (HH17/16, HH25) fall back to libpair via dedup.
        if mapped_geno is None:
            continue

        rows.append({
            "cb": strip_suffix(bc), "library": lib, "suffix": sfx,
            "demux_group": "G6_combined_hh13", "donor_id": mapped_geno,
            "embryo_group": "13-5..8" if mapped_geno in {"g16", "g25", "g32", "g35"} else "",
            "status": status,
        })
        n_hh13 += 1
    print(f"  {n_hh13:,} cells from combined HH13 vireo")
else:
    print("G6 combined HH13 recovery: not available (skipping)")
    if not _best_seed_dir:
        print("  No completed seeds found")
    if not os.path.isfile(HH13_DONOR_MAPPING):
        print(f"  Donor mapping not found: {HH13_DONOR_MAPPING}")

# --- G6 dp500 donor29 ONLY (25-8 rescue) ---
# dp500 is authoritative for GEX05 (where libpair fragmented 25-8).
# 25-8's PG3 cells are handled separately: libpair assigns them to g36,
# and b1 merges g36 → dp500_donor29 via cross_pair_merge.
# Emit dp500 donor29 for GEX05 only — no suffix hardcoding for PG3.
print("G6 dp500 (donor29 = 25-8 rescue, GEX05)...")
g6_dp500 = pd.read_csv(f"{OUTBASE}/group6/vireo_filt_dp500/donor_ids.tsv",
                        sep='\t')
n_dp500_25_8 = 0
for _, row in g6_dp500.iterrows():
    if row['donor_id'] != 'donor29':
        continue
    bc = row['cell']
    sfx = get_suffix(bc)
    if sfx != '-5':
        continue  # Only GEX05; PG3 cells come via g36 → cross_pair_merge in b1
    lib = SUFFIX_TO_LIBRARY_G6.get(sfx, f"unknown{sfx}")
    rows.append({
        'cb': strip_suffix(bc), 'library': lib, 'suffix': sfx,
        'demux_group': 'G6_dp500', 'donor_id': 'donor29',
        'embryo_group': '25-8', 'status': 'singlet',
    })
    n_dp500_25_8 += 1
print(f"  {n_dp500_25_8:,} cells in GEX05 (donor29 = 25-8)")

# --- G9 dp500 (PILOTGEX4A primary, also PILOT4B + PILOT3) ---
# CRITICAL: must use the SAME G9 result as g9_g6_mapping.tsv (seeded best seed).
# The single-seed vireo_filt_dp500 has DIFFERENT donor numbering!
# See vireo_provenance.tsv for canonical paths.
G9_CANONICAL = f"{SCRATCH}/vireo_seeded/g9_dp100/seed_126746/donor_ids.tsv"
print(f"G9 dp100 (seeded best seed: seed_126746)...")
assert os.path.isfile(G9_CANONICAL), f"G9 canonical result not found: {G9_CANONICAL}"
g9 = pd.read_csv(G9_CANONICAL, sep='\t')

g9_donor_sfx = defaultdict(lambda: defaultdict(int))
for _, row in g9.iterrows():
    if row['donor_id'] not in ('doublet', 'unassigned'):
        g9_donor_sfx[row['donor_id']][get_suffix(row['cell'])] += 1

for _, row in g9.iterrows():
    bc = row['cell']
    donor = row['donor_id']
    sfx = get_suffix(bc)
    lib = SUFFIX_TO_LIBRARY_G9.get(sfx, f"unknown{sfx}")
    status = classify_status(donor)

    embryo = ''
    if status == 'singlet':
        embryo = G9_TO_EMBRYO.get(donor, '')
        if not embryo:
            sfx_d = g9_donor_sfx[donor]
            total = sum(sfx_d.values())
            embryo = 'PILOT4_unlinked' if total >= 100 else 'tiny_unlinked'

    rows.append({
        'cb': strip_suffix(bc), 'library': lib, 'suffix': sfx,
        'demux_group': 'G9_dp500', 'donor_id': donor,
        'embryo_group': embryo, 'status': status,
    })
n_g9 = sum(1 for r in rows if r['demux_group'] == 'G9_dp500')
print(f"  {n_g9:,} cells")

# --- G2 old vireo (GEX-biased variants, multi-chunk) ---
print("G2 old vireo (GEX-biased variant discovery)...")
g2 = load_g2_old_vireo(G2_OLD_VIREO_BASE, G2_OLD_QC_DIR)

g2_donor_sfx = defaultdict(lambda: defaultdict(int))
for _, row in g2.iterrows():
    if row['donor_id'] not in ('doublet', 'unassigned') and not row['donor_id'].startswith('unmapped'):
        g2_donor_sfx[row['donor_id']][get_suffix(row['cell'])] += 1

for _, row in g2.iterrows():
    bc = row['cell']
    donor = row['donor_id']
    sfx = get_suffix(bc)
    lib = SUFFIX_TO_LIBRARY_G2.get(sfx, f"unknown{sfx}")
    status = classify_status(donor)

    # Unmapped donors (below QC threshold) treated as unassigned
    if donor.startswith('unmapped'):
        status = 'unassigned'

    embryo = ''
    if status == 'singlet':
        sfx_d = g2_donor_sfx[donor]
        has_3 = sfx_d.get('-3', 0)
        has_12 = sfx_d.get('-1', 0) + sfx_d.get('-2', 0)
        if has_3 == 0 and has_12 > 0:
            embryo = '17-Md-13'
        elif has_12 == 0 and has_3 > 0:
            embryo = '17-KC-1'
        else:
            embryo = 'G2'

    rows.append({
        'cb': strip_suffix(bc), 'library': lib, 'suffix': sfx,
        'demux_group': 'G2_old', 'donor_id': donor,
        'embryo_group': embryo, 'status': status,
    })
n_g2 = sum(1 for r in rows if r['demux_group'] == 'G2_old')
print(f"  {n_g2:,} cells")

# --- G3 vireo N=11 best seed (de novo variants, multi-seed) ---
# Uses 2024 de novo variant discovery (173K pool-specific positions)
# instead of 665K universal set which caused donor13 shattering.
print("G3 vireo N=11 (de novo variants, best seed from 100-seed run)...")
g3 = pd.read_csv(f"{OUTBASE}/group3/vireo_N11_denovo_best/donor_ids.tsv", sep='\t')

n_chick = 0
for _, row in g3.iterrows():
    bc = row['cell']
    donor = row['donor_id']

    # Skip chick contaminant
    if donor in G3_EXCLUDE:
        n_chick += 1
        continue

    sfx = get_suffix(bc)
    lib = SUFFIX_TO_LIBRARY_G3.get(sfx, f"unknown{sfx}")
    status = classify_status(donor)

    rows.append({
        'cb': strip_suffix(bc), 'library': lib, 'suffix': sfx,
        'demux_group': 'G3_unfiltered', 'donor_id': donor,
        'embryo_group': 'G3',
        'status': status,
    })
if n_chick:
    print(f"  Excluded {n_chick:,} chick cells (G3_EXCLUDE: {G3_EXCLUDE})")
n_g3 = sum(1 for r in rows if r['demux_group'] == 'G3_unfiltered')
print(f"  {n_g3:,} cells")

# ============================================================
# Build DataFrame
# ============================================================
df = pd.DataFrame(rows)
print(f"\nTotal rows (before dedup): {len(df):,}")

# Summary by demux_group
for grp in ['G6_libpair', 'G6_dp500', 'G9_dp500', 'G2_old', 'G3_unfiltered']:
    gdf = df[df['demux_group'] == grp]
    if len(gdf) == 0:
        continue
    singlet = (gdf['status'] == 'singlet').sum()
    print(f"  {grp}: {len(gdf):,} total, {singlet:,} singlets")

# ============================================================
# Deduplicate: cells appearing in multiple demux groups
# ============================================================
print("\nDeduplicating (cells in multiple demux groups)...")
n_before = len(df)

# Status priority: singlet wins over doublet wins over unassigned
status_priority = {"singlet": 0, "doublet": 1, "unassigned": 2}
df["_priority"] = df["status"].map(status_priority).fillna(3).astype(int)

# Demux group preference per library (lower = preferred)
# Key design: dp500 authoritative for 25-8 ONLY in GEX05 (its only library by design).
# dp500 ran on the full G6 pool (all 6 libraries), so donor29 picks up barcodes from
# PILOTGEX3 and other libraries where 25-8 has NO cells by experimental design.
# Those barcodes belong to their correct libpair genotypes (25-5..7, 13-5..8, etc.).
# libpair is preferred for all paired libraries; dp500 only wins in GEX05.
# Sort by (_pref, _priority) so group preference dominates, then status breaks ties.
LIB_DEMUX_PREF = {
    "GEX05": {"G6_dp500": 0, "G6_libpair": 1},
    "GEX06": {"G6_libpair": 0},
    # GEX07 + PILOTGEX3: combined HH13 vireo preferred (resolves
    # clutch-mates g32/g35 that libpair runB couldn't separate)
    "GEX07": {"G6_combined_hh13": 0, "G6_libpair": 1},
    "GEX08": {"G6_libpair": 0},
    "PILOTGEX3": {"G6_combined_hh13": 0, "G6_libpair": 1, "G9_dp500": 2},
    "PILOTGEX4A": {"G9_dp500": 0},
    "PILOTGEX4B": {"G6_libpair": 0, "G9_dp500": 1},
    "GEX03-1AB": {"G2_old": 0},
    "GEX03-2": {"G2_old": 0},
    "GEX04": {"G3_unfiltered": 0},
    "GEX0S-duck": {"G3_unfiltered": 0},
}

df["_pref"] = df.apply(
    lambda r: LIB_DEMUX_PREF.get(r["library"], {}).get(r["demux_group"], 99),
    axis=1,
)

df = (
    df.sort_values(["_pref", "_priority"])
    .drop_duplicates(subset=["cb", "library"], keep="first")
    .drop(columns=["_priority", "_pref"])
)
n_after = len(df)
print(f"  {n_before:,} -> {n_after:,} (removed {n_before - n_after:,} duplicates)")

# ============================================================
# Summary
# ============================================================
print(f"\n--- Cells per library ---")
lib_summary = df.groupby(['library', 'demux_group', 'status']).size().reset_index(name='n')
for lib in sorted(df['library'].unique()):
    ldf = lib_summary[lib_summary['library'] == lib]
    total = ldf['n'].sum()
    singlet = ldf[ldf['status'] == 'singlet']['n'].sum()
    groups = sorted(ldf['demux_group'].unique())
    print(f"  {lib:<16} total:{total:>7,}  singlet:{singlet:>7,} "
          f"({singlet/total*100:.1f}%)  [{', '.join(groups)}]")

print(f"\n--- Singlets per embryo group ---")
singlets = df[df['status'] == 'singlet']
eg_counts = singlets.groupby('embryo_group').size().sort_values(ascending=False)
for eg, n in eg_counts.items():
    if eg:
        print(f"  {eg:<30} {n:>7,}")

# ============================================================
# Save
# ============================================================
outfile = f"{OUTDIR}/cb_genotype_map_v2.tsv.gz"
df.to_csv(outfile, sep='\t', index=False, compression='gzip')
print(f"\nSaved: {outfile}")
print(f"  {len(df):,} rows, {df['cb'].nunique():,} unique barcodes")

# Library-level summary
summary_file = f"{OUTDIR}/library_summary_v2.tsv"
lib_detail = df.groupby(['library', 'demux_group']).agg(
    total=('cb', 'size'),
    singlet=('status', lambda x: (x == 'singlet').sum()),
    doublet=('status', lambda x: (x == 'doublet').sum()),
    unassigned=('status', lambda x: (x == 'unassigned').sum()),
    n_donors=('donor_id', lambda x: x[~x.isin(['doublet', 'unassigned'])].nunique()),
).reset_index()
lib_detail['singlet_pct'] = (lib_detail['singlet'] / lib_detail['total'] * 100).round(1)
lib_detail.to_csv(summary_file, sep='\t', index=False)
print(f"Saved: {summary_file}")

print("\n=== Done ===")
