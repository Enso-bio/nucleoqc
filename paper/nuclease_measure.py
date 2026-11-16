#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Natural experiment, step 2: measures the fraction of MS2 spectra with nucleobase (atlas metric) in projects comparing
(A) the same purification with and without RNase and (B) nuclear vs cytoplasmic fractions. Files are assigned to arms by
name with hand-reviewed rules (below). Output: nuclease_experiment.tsv"""
import re, json, os
from concurrent.futures import ThreadPoolExecutor
import atlas

# project -> (type, [(arm, regex on the file name)], exclusion regex)   -- first arm = "treated"
DESIGN = {
    'PXD011948': ('RNase', [('+RNase', r'rnase'), ('-RNase', r'.')], None),
    'PXD023516': ('RNase', [('+RNase', r'rnase'), ('-RNase', r'_exp')], None),
    'PXD042819': ('RNase', [('+RNase', r'plusrnase'), ('-RNase', r'minusrnase')], r'N2a_cells'),
    'PXD024388': ('RNase', [('+RNase', r'_rnase'), ('-RNase', r'_nouv')], None),
    'PXD013343': ('fraction', [('nuclear', r'nuc'), ('cytoplasm', r'cyto')], None),
    'PXD045677': ('fraction', [('nuclear', r'_nuc_'), ('cytoplasm', r'_cyto_')], None),
}
OUTPUT = 'nuclease_experiment.tsv'
COLS = ['project', 'type', 'arm', 'file', 'file_MB', 'n_ms2', 'valid', 'reason', 'coverage_136', 'resolution_score', 'nb10', 'nb2', 'A2', 'G2', 'C2', 'ribose2']


def tasks():
    d = json.load(open('nuclease_candidates.json'))['candidates']
    out = []
    for p, (kind, arms, excl) in DESIGN.items():
        fs = [f for g in d[p].values() for f in g]
        for n, size, fmt, url in fs:
            if excl and re.search(excl, n, re.I):
                continue
            b = next((name for name, rx in arms if re.search(rx, n, re.I)), None)
            if b:
                out.append((p, kind, b, n, size, fmt, url))
    return out


def one(t):
    p, kind, b, n, size, fmt, url = t
    try:
        r = atlas.measure(url, size, fmt)
    except Exception as e:
        r = {'n_ms2': 0, 'valid': False, 'reason': 'error: %s' % str(e)[:60]}
    row = [p, kind, b, n, '%.0f' % (size / 1e6)] + [('%.5f' % r[k]) if isinstance(r.get(k), float) else str(r.get(k, '')) for k in COLS[5:]]
    with open(OUTPUT, 'a', encoding='utf-8') as f:
        f.write('\t'.join(row) + '\n')
    print('%s %-11s %-60s n=%-5s nb10=%s %s' % (p, b, n[:60], r.get('n_ms2'), ('%.4f' % r['nb10']) if 'nb10' in r else '-', r.get('reason', '')), flush=True)


if __name__ == '__main__':
    T = tasks()
    done = set()
    if os.path.exists(OUTPUT):
        done = {l.split('\t')[3] for l in open(OUTPUT, encoding='utf-8').read().split('\n')[1:] if l}
    else:
        open(OUTPUT, 'w', encoding='utf-8').write('\t'.join(COLS) + '\n')
    T = [t for t in T if t[3] not in done]
    print('%d files to measure' % len(T), flush=True)
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(one, T))
    print('DONE')
