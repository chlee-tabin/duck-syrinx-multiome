#!/usr/bin/env python3
"""Compare vireo donor assignments between two seeds."""
import csv, sys
from collections import Counter, defaultdict

base = sys.argv[1]  # vireo output base dir
s1, s2 = sys.argv[2], sys.argv[3]  # seed names

def load(path):
    d = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            d[row['cell']] = row['donor_id']
    return d

a = load(f'{base}/{s1}/donor_ids.tsv')
b = load(f'{base}/{s2}/donor_ids.tsv')
common = set(a) & set(b)
print(f'Shared cells: {len(common)}')

# Build donor mapping from a -> b
a2b = defaultdict(Counter)
for c in common:
    if a[c] not in ('doublet','unassigned') and b[c] not in ('doublet','unassigned'):
        a2b[a[c]][b[c]] += 1

mapping = {}
for d in sorted(a2b):
    best, n = a2b[d].most_common(1)[0]
    mapping[d] = (best, n, sum(a2b[d].values()))

agree = disagree = both_dbl = both_una = mixed = 0
for c in common:
    x, y = a[c], b[c]
    if x == 'doublet' and y == 'doublet': both_dbl += 1
    elif x == 'unassigned' and y == 'unassigned': both_una += 1
    elif x in ('doublet','unassigned') or y in ('doublet','unassigned'): mixed += 1
    elif x in mapping and mapping[x][0] == y: agree += 1
    else: disagree += 1

tot = agree + disagree
print(f'Singlet-singlet: {tot} (agree={agree}, disagree={disagree})')
if tot: print(f'Concordance: {agree/tot*100:.1f}%')
print(f'Both doublet: {both_dbl}')
print(f'Both unassigned: {both_una}')
print(f'Mixed: {mixed}')
print(f'\nDonor mapping ({s1} -> {s2}):')
for d, (best, n, t) in sorted(mapping.items(), key=lambda x: -x[1][1]):
    print(f'  {d:>10} -> {best:<10}: {n:>5}/{t:<5} ({n/t*100:.1f}%)')
