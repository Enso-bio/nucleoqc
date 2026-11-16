#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analysis of identifications*.tsv / ids_spectra*.tsv (Sage, 1 % FDR): within-file ratio, Mantel-Haenszel OR,
dose-response, loss per file, Figure S3.  Usage: python identifications_analysis.py _1p"""
import collections, math
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
SUF = sys.argv[1] if len(sys.argv) > 1 else ''      # '_1p' = one-step search against the whole proteome
A = [l.split('\t') for l in open('identifications%s.tsv' % SUF, encoding='utf-8').read().split('\n') if l]
ca = A[0]; A = [dict(zip(ca, r)) for r in A[1:]]
S = [l.split('\t') for l in open('ids_spectra%s.tsv' % SUF, encoding='utf-8').read().split('\n')[1:] if l]
by_file = collections.defaultdict(list)
for p, fn, i, pmz, z, nb, ide, hs in S:
    by_file[(p, fn)].append((float(pmz), float(nb), int(ide), float(hs)))
group = {(a['project'], a['file']): a['group'] for a in A if float(a['id_rate_all']) >= 0.02}
failed = [a for a in A if float(a['id_rate_all']) < 0.02]
print('Failed searches (<2 %% identified, excluded): %s' % ', '.join('%s %s' % (a['project'], a['file'][:30]) for a in failed))
A = [a for a in A if float(a['id_rate_all']) >= 0.02]
cont = [a for a in A if a['group'] == 'contaminated']; clean = [a for a in A if a['group'] == 'clean']
print('Files searched: %d contaminated (%d projects), %d clean (%d projects)'
      % (len(cont), len({a['project'] for a in cont}), len(clean), len({a['project'] for a in clean})))

# 1) per file: observed vs expected (same file, precursor +-10 Th)
print('\n1) Identification rate of spectra with nucleobase (>=10 % of the base peak), per file')
rows = [(a, float(a['id_rate_nucleobase']), float(a['id_rate_expected_nucleobase'])) for a in cont
        if a['id_rate_nucleobase'] and a['id_rate_expected_nucleobase'] and int(a['n_nucleobase_spectra']) >= 20]
for a, o, e in sorted(rows, key=lambda t: -int(t[0]['n_nucleobase_spectra']) / int(t[0]['n_ms2'])):
    print('   %s %-34s %-7s phospho=%s  nb %4d/%-5d  ID with nb %5.1f %%  expected %5.1f %%  (without nb %5.1f %%)'
          % (a['project'], a['file'][:34], a['labelling'], a['phospho_enrichment'], int(a['n_nucleobase_spectra']), int(a['n_ms2']),
             100 * o, 100 * e, 100 * float(a['id_rate_no_nucleobase'])))
o = np.array([f[1] for f in rows]); e = np.array([f[2] for f in rows])
print('   %d files: median observed %.1f %% vs expected %.1f %%; lower in %d of %d; Wilcoxon p = %.2g; median ratio %.2f'
      % (len(o), 100 * np.median(o), 100 * np.median(e), int((o < e).sum()), len(o), stats.wilcoxon(o, e).pvalue, np.median(o / e)))

# 2) Mantel-Haenszel: strata = file x 20 Th precursor window
num = den = 0.0
for k, v in by_file.items():
    if group.get(k) != 'contaminated':
        continue
    v = np.array(v)
    for b in np.unique(np.floor(v[:, 0] / 20)):
        s = v[np.floor(v[:, 0] / 20) == b]; n = len(s)
        x = s[:, 1] >= 0.10; y = s[:, 2] == 1
        a_, b_, c_, d_ = (x & y).sum(), (x & ~y).sum(), (~x & y).sum(), (~x & ~y).sum()
        num += a_ * d_ / n; den += b_ * c_ / n
print('\n2) Mantel-Haenszel odds ratio (identified ~ nucleobase; strata file x 20 Th): %.2f' % (num / den))

# 3) dose-response: relative intensity of the nucleobase
print('\n3) Dose-response (contaminated files; observed / expected rate per stratum file x 20 Th)')
cuts = [(0, 0.001, '<0.1 %'), (0.001, 0.02, '0.1-2 %'), (0.02, 0.10, '2-10 %'), (0.10, 0.30, '10-30 %'), (0.30, 1.01, '>=30 %')]
for lo, hi, name in cuts:
    obs = exp_ = n = 0
    for k, v in by_file.items():
        if group.get(k) != 'contaminated':
            continue
        v = np.array(v); bins = np.floor(v[:, 0] / 20)
        rate = {b: (v[(bins == b) & (v[:, 1] < 0.001), 2].mean() if ((bins == b) & (v[:, 1] < 0.001)).sum() >= 5 else np.nan) for b in np.unique(bins)}
        sel = (v[:, 1] >= lo) & (v[:, 1] < hi) & np.array([not np.isnan(rate[b]) for b in bins])
        obs += v[sel, 2].sum(); exp_ += sum(rate[b] for b in bins[sel]); n += sel.sum()
    print('   nucleobase %-15s n=%6d  identified %5.1f %%  observed/expected %.2f' % (name, n, 100 * obs / max(n, 1), obs / max(exp_, 1e-9)))

# 4) identifications lost per file
print('\n4) Identifications lost per file = n_nb x (expected - observed) / (identified + lost)')
loss = []
for a, o_, e_ in rows:
    n_id = float(a['id_rate_all']) * int(a['n_ms2'])
    lost = int(a['n_nucleobase_spectra']) * (e_ - o_)
    loss.append((lost / max(n_id + lost, 1), int(a['n_nucleobase_spectra']) / int(a['n_ms2']), a['project']))
for x, f, p in sorted(loss, key=lambda t: -t[1])[:10]:
    print('   %s  nb in %5.1f %% of MS2 -> %+.1f %% identifications' % (p, 100 * f, -100 * x))
print('   median %+.1f %%; in files with >20 %% of spectra with nucleobase: %s'
      % (-100 * np.median([x for x, _, _ in loss]), ', '.join('%+.1f %%' % (-100 * x) for x, f, _ in loss if f > 0.2)))

# 5) between files (confounded: different laboratory, sample and enrichment)
print('\n5) Overall identification rate: contaminated median %.1f %% vs clean %.1f %% (Mann-Whitney p = %.2g; confounded by project)'
      % (100 * np.median([float(a['id_rate_all']) for a in cont]), 100 * np.median([float(a['id_rate_all']) for a in clean]),
         stats.mannwhitneyu([float(a['id_rate_all']) for a in cont], [float(a['id_rate_all']) for a in clean]).pvalue))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].scatter(100 * e, 100 * o, s=18, c='#c0504d'); m = max(100 * e.max(), 100 * o.max()) * 1.05
ax[0].plot([0, m], [0, m], 'k--', lw=.7); ax[0].set_xlabel('% identified, expected (neighbours without nucleobase)')
ax[0].set_ylabel('% identified, spectra with nucleobase'); ax[0].set_title('Sage, 1 % FDR, per file', fontsize=10)
xs, ys = [], []
for lo, hi, name in cuts:
    obs = exp_ = 0
    for k, v in by_file.items():
        if group.get(k) != 'contaminated':
            continue
        v = np.array(v); bins = np.floor(v[:, 0] / 20); rate = {b: (v[(bins == b) & (v[:, 1] < 0.001), 2].mean() if ((bins == b) & (v[:, 1] < 0.001)).sum() >= 5 else np.nan) for b in np.unique(bins)}
        sel = (v[:, 1] >= lo) & (v[:, 1] < hi) & np.array([not np.isnan(rate[b]) for b in bins]); obs += v[sel, 2].sum(); exp_ += sum(rate[b] for b in bins[sel])
    xs.append(name); ys.append(obs / max(exp_, 1e-9))
ax[1].bar(range(len(xs)), ys, color='#3a6ea5'); ax[1].axhline(1, c='k', lw=.6); ax[1].set_xticks(range(len(xs))); ax[1].set_xticklabels(xs, fontsize=8)
ax[1].set_xlabel('adenine/guanine intensity (% of base peak)'); ax[1].set_ylabel('identified, observed / expected')
ax[1].set_title('Dose-response', fontsize=10)
plt.tight_layout(); plt.savefig('identifications_figure%s.png' % SUF, dpi=150); print('\nfigure: identifications_figure%s.png' % SUF)
