#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Natural experiment: within-project comparisons (with/without RNase; nuclear/cytoplasmic). Reads
nuclease_experiment.tsv; paired test for PXD042819; Figure S4."""
import collections
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nuclease_measure import DESIGN

L = [l.split('\t') for l in open('nuclease_experiment.tsv', encoding='utf-8').read().split('\n') if l]
header = L[0]; F = [dict(zip(header, r)) for r in L[1:] if len(r) == len(header)]
print('Files measured: %d; valid: %d' % (len(F), sum(f['valid'] == 'True' for f in F)))
for m, n in collections.Counter(f['reason'] for f in F if f['valid'] != 'True').items():
    print('   discarded because of %s: %d' % (m, n))
V = [f for f in F if f['valid'] == 'True']
EPS = 1e-4
res = []
for p, (kind, arms, _) in DESIGN.items():
    a, b = arms[0][0], arms[1][0]
    xa = np.array([float(f['nb10']) for f in V if f['project'] == p and f['arm'] == a])
    xb = np.array([float(f['nb10']) for f in V if f['project'] == p and f['arm'] == b])
    if len(xa) < 2 or len(xb) < 2:
        print('%s: not enough valid files (%d / %d)' % (p, len(xa), len(xb)))
        continue
    pv = stats.mannwhitneyu(xa, xb, alternative='two-sided').pvalue
    fc = (np.median(xa) + EPS) / (np.median(xb) + EPS)
    z = np.sign(np.log(fc)) * stats.norm.isf(pv / 2)
    res.append((p, kind, a, b, xa, xb, fc, pv, z))
    print('%s %-8s %-11s n=%2d median %6.3f %% | %-11s n=%2d median %6.3f %% | ratio %6.2f | Mann-Whitney p = %.2g'
          % (p, kind, a, len(xa), 100 * np.median(xa), b, len(xb), 100 * np.median(xb), fc, pv))

print()
for kind, name in (('RNase', 'with RNase / without RNase'), ('fraction', 'nuclear / cytoplasm')):
    r = [x for x in res if x[1] == kind]
    if r:
        zs = np.array([x[8] for x in r]); Z = zs.sum() / np.sqrt(len(zs))
        print('%s: %d projects; ratios %s; in the expected direction %d of %d; Stouffer Z = %.2f, p = %.2g'
              % (name, len(r), ', '.join('%.2f' % x[6] for x in r),
                 sum((x[6] < 1) if kind == 'RNase' else (x[6] > 1) for x in r), len(r), Z, 2 * stats.norm.sf(abs(Z))))

# PXD042819: paired design (same co-IP and replicate, with and without RNase)
import re
pair = {}
for f in F:
    m = re.search(r'IP_(Ctrl|T3)(plus|minus)RNA[sS]e_(\d)', f['file']) if f['project'] == 'PXD042819' else None
    if m and f['valid'] == 'True':
        pair[(m.group(1), m.group(3), m.group(2))] = (float(f['nb10']), float(f['nb2']))
keys = sorted({k[:2] for k in pair if k[:2] + ('plus',) in pair and k[:2] + ('minus',) in pair})
if keys:
    print('\nPXD042819 paired (co-IP of ADAR3 or control, same replicate, without -> with RNase):')
    for k in keys:
        s, c = pair[k + ('minus',)], pair[k + ('plus',)]
        print('   %-4s rep %s  nb10 %5.2f %% -> %5.2f %%   nb2 %5.1f %% -> %5.1f %%' % (k[0], k[1], 100 * s[0], 100 * c[0], 100 * s[1], 100 * c[1]))
    a = [pair[k + ('minus',)][0] for k in keys]; b = [pair[k + ('plus',)][0] for k in keys]
    print('   lower with RNase in %d of %d pairs; paired Wilcoxon p = %.3g' % (sum(y < x for x, y in zip(a, b)), len(a), stats.wilcoxon(a, b).pvalue))

if res:
    fig, ax = plt.subplots(1, len(res), figsize=(2.1 * len(res), 3.6), sharey=True)
    ax = np.atleast_1d(ax)
    for k, (p, kind, a, b, xa, xb, fc, pv, z) in enumerate(res):
        for j, (x, c) in enumerate(((xb, '#9bbbd9'), (xa, '#c0504d'))):
            ax[k].scatter(j + np.random.default_rng(k + j).uniform(-.12, .12, len(x)), 100 * x + 100 * EPS, s=14, c=c)
            ax[k].plot([j - .25, j + .25], [100 * np.median(x) + 100 * EPS] * 2, 'k-', lw=1)
        LABELS = {'+RNase': '+RNase', '-RNase': '−RNase', 'nuclear': 'nuclear', 'cytoplasm': 'cytoplasm'}
        ax[k].set_xticks([0, 1]); ax[k].set_xticklabels([LABELS.get(b, b), LABELS.get(a, a)], fontsize=7, rotation=20)
        ax[k].set_title('%s\np = %.2g (unpaired)' % (p, pv), fontsize=8)
    ax[0].set_yscale('log'); ax[0].set_ylabel('% of MS2 with adenine/guanine $\\geq$ 10 %')
    plt.tight_layout(); plt.savefig('nuclease_figure.png', dpi=150); print('\nfigure: nuclease_figure.png')
