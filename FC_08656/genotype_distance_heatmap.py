#!/usr/bin/env python3
"""
Pairwise genotype distance heatmap within each group.

Computes pairwise Pearson correlation (and 1-r distance) between all
genotypes within G2, G5, and G6, then generates heatmaps with
hierarchical clustering to reveal breed/population structure.

For G2: annotates known breed labels (Md vs KC).
For G5/G6: all Mallard, but shows sub-structure from over-splitting.

Usage:
    python genotype_distance_heatmap.py \
        --g2_base <g2_vireo_base> --g2_qc <g2_qc> \
        --g5_base <g5_vireo_base> --g5_qc <g5_qc> \
        --g6_base <g6_vireo_base> --g6_qc <g6_qc> \
        --output_dir genotype_distance_output
"""

import argparse
import csv
import gzip
import json
import os
import sys
from datetime import datetime

import numpy as np


def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# --- GT loading (shared) ---

def parse_PL_to_prob(pl_str):
    try:
        pls = [int(x) for x in pl_str.split(',')]
        if len(pls) != 3:
            return np.array([0.33, 0.34, 0.33])
        probs = np.array([10.0 ** (-min(pl, 100) / 10.0) for pl in pls])
        total = probs.sum()
        if total == 0:
            return np.array([0.33, 0.34, 0.33])
        return probs / total
    except (ValueError, AttributeError):
        return np.array([0.33, 0.34, 0.33])


def parse_GT_to_prob(gt_str):
    if gt_str in ['./.', '.']:
        return np.array([0.33, 0.34, 0.33])
    sep = '/' if '/' in gt_str else ('|' if '|' in gt_str else None)
    if sep is None:
        return np.array([0.33, 0.34, 0.33])
    try:
        a1, a2 = int(gt_str.split(sep)[0]), int(gt_str.split(sep)[1])
        if a1 == 0 and a2 == 0:
            return np.array([1.0, 0.0, 0.0])
        elif (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
            return np.array([0.0, 1.0, 0.0])
        elif a1 == 1 and a2 == 1:
            return np.array([0.0, 0.0, 1.0])
    except (ValueError, IndexError):
        pass
    return np.array([0.33, 0.34, 0.33])


def load_GT_from_vcf(vcf_file):
    donor_names, GT_list, variant_ids = [], [], []
    field_used = None
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                donor_names = line.strip().split('\t')[9:]
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            variant_ids.append(f"{parts[0]}:{parts[1]}")
            format_fields = parts[8].split(':')
            if field_used is None:
                field_used = 'PL' if 'PL' in format_fields else ('GT' if 'GT' in format_fields else None)
                if field_used is None:
                    continue
            try:
                field_idx = format_fields.index(field_used)
            except ValueError:
                continue
            variant_probs = []
            for sample_data in parts[9:]:
                fields = sample_data.split(':')
                if len(fields) > field_idx:
                    val = fields[field_idx]
                    variant_probs.append(parse_PL_to_prob(val) if field_used == 'PL' else parse_GT_to_prob(val))
                else:
                    variant_probs.append(np.array([0.33, 0.34, 0.33]))
            if len(variant_probs) == len(donor_names):
                GT_list.append(np.array(variant_probs))
    if not GT_list:
        return None, donor_names, variant_ids
    return np.stack(GT_list, axis=0), donor_names, variant_ids


def load_genotype_groups(json_path):
    with open(json_path) as f:
        groups = json.load(f)
    n_chunks = len(next(iter(groups.values()))['merged_pattern'])
    chunk_mapping = [{} for _ in range(n_chunks)]
    for geno, info in groups.items():
        cells = info.get('total_cells', 0)
        for ci, donor in enumerate(info['merged_pattern']):
            if donor is None:
                continue
            existing = chunk_mapping[ci].get(donor)
            if existing is None or cells > existing[1]:
                chunk_mapping[ci][donor] = (geno, cells)
    return [{k: v[0] for k, v in d.items()} for d in chunk_mapping]


def load_best_seeds(seed_ranking_path):
    best = {}
    with open(seed_ranking_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            c = row['chunk_raw']
            if c not in best:
                best[c] = row['seed']
    return best


def load_donor_group_assignment(tsv_path):
    assignment = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            assignment[row['genotype']] = {
                'assigned_group': row.get('assigned_group', ''),
                'candidate_donors': row.get('candidate_donors', ''),
                'total_cells': int(float(row.get('total_cells', 0))),
            }
    return assignment


def load_group_GT(base_dir, qc_dir, group_label):
    gg_path = os.path.join(qc_dir, 'genotype_groups.json')
    sr_path = os.path.join(qc_dir, 'seed_ranking.tsv')
    chunk_to_geno = load_genotype_groups(gg_path)
    best_seeds = load_best_seeds(sr_path)
    log(f"  {group_label}: {len(best_seeds)} chunks, {len(chunk_to_geno)} chunk positions")
    chunk_raws_sorted = sorted(best_seeds.keys())
    chunk_raw_to_idx = {name: i for i, name in enumerate(chunk_raws_sorted)}
    geno_data = {}
    for chunk_raw in chunk_raws_sorted:
        seed = best_seeds[chunk_raw]
        vcf_path = os.path.join(base_dir, chunk_raw, f"seed_{seed}", "GT_donors.vireo.vcf.gz")
        if not os.path.exists(vcf_path):
            vcf_path = os.path.join(base_dir, chunk_raw, f"seed_{seed}", "GT_donors.vireo.vcf")
            if not os.path.exists(vcf_path):
                log(f"  WARNING: VCF not found: {vcf_path}")
                continue
        log(f"  Loading {chunk_raw}/seed_{seed}")
        GT, donors, variants = load_GT_from_vcf(vcf_path)
        if GT is None:
            continue
        log(f"    {len(variants)} variants, {len(donors)} donors")
        chunk_idx = chunk_raw_to_idx[chunk_raw]
        if chunk_idx >= len(chunk_to_geno):
            continue
        mapping = chunk_to_geno[chunk_idx]
        for j, donor in enumerate(donors):
            geno_name = mapping.get(donor)
            if geno_name is None:
                continue
            entry = (len(variants), GT[:, j, :], variants)
            if geno_name not in geno_data:
                geno_data[geno_name] = entry
            elif len(variants) > geno_data[geno_name][0]:
                geno_data[geno_name] = entry
    log(f"  {group_label}: {len(geno_data)} genotypes loaded")
    return geno_data


# --- Pairwise correlation ---

def compute_within_group_corr(geno_data, names):
    """Compute pairwise Pearson correlation on flattened GT prob vectors."""
    n = len(names)
    # All genotypes within a group share variants (same cellSNP run),
    # so we can use the full variant set from each genotype directly.
    # But different chunks may have different variant subsets.
    # Use intersection of all variant sets for this group.
    all_var_sets = [set(geno_data[name][2]) for name in names if name in geno_data]
    if not all_var_sets:
        return np.zeros((n, n)), 0
    shared = sorted(set.intersection(*all_var_sets))
    log(f"    Shared variants within group: {len(shared)}")

    shared_map = {v: i for i, v in enumerate(shared)}
    ns = len(shared)

    # Build vectors
    vectors = {}
    for name in names:
        if name not in geno_data:
            continue
        _, gt_prob, var_ids = geno_data[name]
        vec = np.full((ns, 3), 1.0 / 3)
        for k, vid in enumerate(var_ids):
            if vid in shared_map:
                vec[shared_map[vid]] = gt_prob[k]
        vectors[name] = vec.flatten()

    corr_mat = np.zeros((n, n))
    for i, ni in enumerate(names):
        if ni not in vectors:
            continue
        vi = vectors[ni]
        vi_c = vi - vi.mean()
        norm_i = np.sqrt(np.sum(vi_c ** 2))
        if norm_i == 0:
            continue
        for j, nj in enumerate(names):
            if j <= i:
                continue
            if nj not in vectors:
                continue
            vj = vectors[nj]
            vj_c = vj - vj.mean()
            norm_j = np.sqrt(np.sum(vj_c ** 2))
            if norm_j == 0:
                continue
            r = np.sum(vi_c * vj_c) / (norm_i * norm_j)
            corr_mat[i, j] = r
            corr_mat[j, i] = r
    np.fill_diagonal(corr_mat, 1.0)
    return corr_mat, len(shared)


# --- Hierarchical clustering for reordering ---

def cluster_order(dist_mat):
    """Simple agglomerative clustering to reorder matrix. Returns index order."""
    n = dist_mat.shape[0]
    if n <= 2:
        return list(range(n))
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        # Ensure valid distance matrix
        d = np.clip(dist_mat, 0, None)
        np.fill_diagonal(d, 0)
        condensed = squareform(d, checks=False)
        Z = linkage(condensed, method='average')
        return list(leaves_list(Z))
    except ImportError as e:
        raise SystemExit(f"scipy is required for the dendrogram ordering ({e})") from e


def make_heatmap(corr_mat, names, labels, title, output_path, breed_info=None):
    """Generate and save a correlation heatmap with clustering."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    dist_mat = 1.0 - corr_mat
    order = cluster_order(dist_mat)
    ordered_corr = corr_mat[np.ix_(order, order)]
    ordered_names = [names[i] for i in order]
    ordered_labels = [labels[i] for i in order]

    n = len(names)
    fig_size = max(8, n * 0.25)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    im = ax.imshow(ordered_corr, cmap='RdBu_r', vmin=-0.2, vmax=1.0, aspect='equal')
    plt.colorbar(im, ax=ax, label='Pearson r', shrink=0.7)

    # Tick labels: genotype name + annotation
    tick_labels = []
    tick_colors = []
    for i, idx in enumerate(order):
        name = names[idx]
        label = labels[idx]
        short = name.replace('genotype_', 'g')
        if breed_info and idx in breed_info:
            breed = breed_info[idx]
            tick_labels.append(f"{short} [{breed}]")
            if 'KC' in breed:
                tick_colors.append('#e41a1c')
            elif 'Md' in breed or 'Mallard' in breed:
                tick_colors.append('#377eb8')
            else:
                tick_colors.append('black')
        else:
            tick_labels.append(f"{short} ({label})")
            tick_colors.append('black')

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=max(5, min(9, 200 // n)))
    ax.set_yticklabels(tick_labels, fontsize=max(5, min(9, 200 // n)))

    for i, color in enumerate(tick_colors):
        ax.get_xticklabels()[i].set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    log(f"  Saved: {output_path}")

    return order


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--g2_base', required=True)
    parser.add_argument('--g2_qc', required=True)
    parser.add_argument('--g5_base', required=True)
    parser.add_argument('--g5_qc', required=True)
    parser.add_argument('--g6_base', required=True)
    parser.add_argument('--g6_qc', required=True)
    parser.add_argument('--output_dir', default='genotype_distance_output')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load G2 ---
    log("Loading G2...")
    g2_data = load_group_GT(args.g2_base, args.g2_qc, "G2")
    g2_assign = load_donor_group_assignment(os.path.join(args.g2_qc, 'donor_group_assignment.tsv'))
    g2_names = sorted(g2_data.keys(), key=lambda x: int(x.replace('genotype_', '')))

    # G2 breed labels
    g2_breed = {}
    for i, name in enumerate(g2_names):
        info = g2_assign.get(name, {})
        candidates = info.get('candidate_donors', '')
        group = info.get('assigned_group', '')
        if group == '{-1, -2}' and '17-Md-13' in candidates:
            g2_breed[i] = 'Md (17-Md-13)'
        elif group == '{-3}' and '17-KC-1' in candidates:
            g2_breed[i] = 'KC (17-KC-1)'
        else:
            g2_breed[i] = '?'

    g2_labels = [g2_assign.get(n, {}).get('assigned_group', '') for n in g2_names]

    # --- Load G5 ---
    log("Loading G5...")
    g5_data = load_group_GT(args.g5_base, args.g5_qc, "G5")
    g5_assign = load_donor_group_assignment(os.path.join(args.g5_qc, 'donor_group_assignment.tsv'))
    g5_names = sorted(g5_data.keys(), key=lambda x: int(x.replace('genotype_', '')))
    g5_labels = [g5_assign.get(n, {}).get('assigned_group', '') for n in g5_names]

    # --- Load G6 ---
    log("Loading G6...")
    g6_data = load_group_GT(args.g6_base, args.g6_qc, "G6")
    g6_assign = load_donor_group_assignment(os.path.join(args.g6_qc, 'donor_group_assignment.tsv'))
    g6_names = sorted(g6_data.keys(), key=lambda x: int(x.replace('genotype_', '')))
    g6_labels = [g6_assign.get(n, {}).get('assigned_group', '') for n in g6_names]

    # --- Compute pairwise correlations ---
    log("Computing G2 pairwise correlations...")
    g2_corr, g2_nvar = compute_within_group_corr(g2_data, g2_names)

    log("Computing G5 pairwise correlations...")
    g5_corr, g5_nvar = compute_within_group_corr(g5_data, g5_names)

    log("Computing G6 pairwise correlations...")
    g6_corr, g6_nvar = compute_within_group_corr(g6_data, g6_names)

    # --- Save correlation matrices ---
    for label, corr, names in [('g2', g2_corr, g2_names), ('g5', g5_corr, g5_names), ('g6', g6_corr, g6_names)]:
        path = os.path.join(args.output_dir, f'{label}_pairwise_corr.tsv')
        with open(path, 'w') as f:
            f.write('genotype\t' + '\t'.join(names) + '\n')
            for i, name in enumerate(names):
                vals = '\t'.join(f"{corr[i,j]:.4f}" for j in range(len(names)))
                f.write(f"{name}\t{vals}\n")
        log(f"  Saved: {path}")

    # --- Generate heatmaps ---
    try:
        import matplotlib
        matplotlib.use('Agg')

        log("Generating G2 heatmap...")
        g2_order = make_heatmap(g2_corr, g2_names, g2_labels,
                                f'Group 2 Genotype Correlation (n={len(g2_names)}, {g2_nvar} variants)\n'
                                'Blue=Mallard, Red=KC',
                                os.path.join(args.output_dir, 'g2_corr_heatmap.png'),
                                breed_info=g2_breed)

        log("Generating G5 heatmap...")
        make_heatmap(g5_corr, g5_names, g5_labels,
                     f'Group 5 Genotype Correlation (n={len(g5_names)}, {g5_nvar} variants)',
                     os.path.join(args.output_dir, 'g5_corr_heatmap.png'))

        log("Generating G6 heatmap...")
        make_heatmap(g6_corr, g6_names, g6_labels,
                     f'Group 6 Genotype Correlation (n={len(g6_names)}, {g6_nvar} variants)',
                     os.path.join(args.output_dir, 'g6_corr_heatmap.png'))

    except ImportError as e:
        raise SystemExit(f"matplotlib/scipy are required for the heatmaps ({e})") from e

    # --- G2 breed analysis ---
    log("\n" + "=" * 80)
    log("G2 BREED ANALYSIS")
    log("=" * 80)

    # Compute mean within-breed and between-breed correlations
    md_idx = [i for i in range(len(g2_names)) if 'Md' in g2_breed.get(i, '')]
    kc_idx = [i for i in range(len(g2_names)) if 'KC' in g2_breed.get(i, '')]
    unk_idx = [i for i in range(len(g2_names)) if g2_breed.get(i, '?') == '?']

    log(f"  Known Mallard: {len(md_idx)} ({[g2_names[i] for i in md_idx]})")
    log(f"  Known KC: {len(kc_idx)} ({[g2_names[i] for i in kc_idx]})")
    log(f"  Unknown: {len(unk_idx)} ({[g2_names[i] for i in unk_idx]})")

    # For each unknown, compute mean correlation to Md reference vs KC reference
    if md_idx and kc_idx:
        log("\n  Per-donor breed inference (mean r to Md vs KC reference):")
        log(f"  {'Genotype':>13} {'Cells':>6} {'r_to_Md':>8} {'r_to_KC':>8} {'Breed_call':>12}")
        log(f"  {'-'*13} {'-'*6} {'-'*8} {'-'*8} {'-'*12}")

        # First show reference donors
        for i in md_idx:
            r_md = np.mean([g2_corr[i, j] for j in md_idx if j != i]) if len(md_idx) > 1 else float('nan')
            r_kc = np.mean([g2_corr[i, j] for j in kc_idx])
            cells = g2_assign.get(g2_names[i], {}).get('total_cells', 0)
            log(f"  {g2_names[i]:>13} {cells:>6} {r_md:>8.4f} {r_kc:>8.4f} {'Md (ref)':>12}")

        for i in kc_idx:
            r_md = np.mean([g2_corr[i, j] for j in md_idx])
            r_kc = np.mean([g2_corr[i, j] for j in kc_idx if j != i]) if len(kc_idx) > 1 else float('nan')
            cells = g2_assign.get(g2_names[i], {}).get('total_cells', 0)
            log(f"  {g2_names[i]:>13} {cells:>6} {r_md:>8.4f} {r_kc:>8.4f} {'KC (ref)':>12}")

        log(f"  {'-'*50}")

        # Now unknown donors
        for i in unk_idx:
            r_md = np.mean([g2_corr[i, j] for j in md_idx])
            r_kc = np.mean([g2_corr[i, j] for j in kc_idx])
            cells = g2_assign.get(g2_names[i], {}).get('total_cells', 0)
            breed = 'Mallard' if r_md > r_kc else 'KC'
            diff = abs(r_md - r_kc)
            conf = 'strong' if diff > 0.02 else 'weak'
            log(f"  {g2_names[i]:>13} {cells:>6} {r_md:>8.4f} {r_kc:>8.4f} {breed:>12} ({conf}, d={diff:.4f})")

    # Overall within-breed vs between-breed summary
    if md_idx and kc_idx:
        within_md = [g2_corr[i, j] for i in md_idx for j in md_idx if i < j]
        within_kc = [g2_corr[i, j] for i in kc_idx for j in kc_idx if i < j]
        between = [g2_corr[i, j] for i in md_idx for j in kc_idx]
        all_unknown_pairs = [(i, j) for i in unk_idx for j in unk_idx if i < j]

        log(f"\n  Within-Md correlation:  {np.mean(within_md):.4f} +/- {np.std(within_md):.4f} (n={len(within_md)} pairs)" if within_md else "")
        log(f"  Within-KC correlation:  {np.mean(within_kc):.4f} +/- {np.std(within_kc):.4f} (n={len(within_kc)} pairs)" if within_kc else "")
        log(f"  Between Md-KC:          {np.mean(between):.4f} +/- {np.std(between):.4f} (n={len(between)} pairs)")

    log("\nDone.")


if __name__ == '__main__':
    main()
