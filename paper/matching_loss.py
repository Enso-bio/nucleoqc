#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Identification loss per file under the adjusted model (charge, m/z, RT, precursor intensity): identifications
predicted with the nucleobase indicators set to zero minus those predicted with the observed values, divided by
identified + lost. -> appended to matching_results.txt"""
import numpy as np
import statsmodels.api as sm
import matching as E

D = E.load()
X, y, g = [], [], []
for k, d in enumerate(D):
    has_rt = np.ptp(d['rt']) > 60; has_i = np.isfinite(d['logi']).mean() > 0.5
    rt = (d['rt'] - d['rt'].min()) / np.ptp(d['rt']) if has_rt else np.full(len(d['rt']), 0.5)
    li = d['logi'] - np.nanmedian(d['logi']) if has_i else np.zeros(len(d['rt']))
    missing = ~np.isfinite(li) if has_i else np.ones(len(li), bool); li = np.where(np.isfinite(li), li, 0.0)
    mz = (d['pmz'] - 800) / 400
    X.append(np.column_stack([(d['nb'] >= 0.10) & (d['nb'] < 0.30), d['nb'] >= 0.30, d['z'] == 3, d['z'] >= 4, d['z'] == 0, mz, mz ** 2,
                              rt, rt ** 2, li, li ** 2, missing]).astype(float))
    y.append(d['id'].astype(float)); g.append(np.full(len(d['id']), k))
Xs, y, g = np.vstack(X), np.concatenate(y), np.concatenate(g)
F = (g[:, None] == np.arange(g.max() + 1)[None, :]).astype(float)
keep = np.where(Xs.std(0) > 0)[0]
Z = np.column_stack([Xs[:, keep], F]); r = sm.Logit(y, Z).fit(disp=0, maxiter=200)
Z0 = Z.copy(); Z0[:, :2] = 0          # the first two columns are the nucleobase indicators
p1, p0 = r.predict(Z), r.predict(Z0)
lines = ['\nD) Identification loss per file with the adjusted model (vs the original unadjusted measure)']
orig = {(a['project'], a['file']): a for a in E.tab('identifications_1p.tsv')}
res = []
for k, d in enumerate(D):
    s = g == k; m = s & (Xs[:, 0] + Xs[:, 1] > 0)
    if (d['nb'] >= 0.10).sum() < 20 or d['group'] != 'contaminated':
        continue
    lost = (p0[m] - p1[m]).sum(); adjusted = lost / max(d['id'].sum() + lost, 1)
    a = orig[(d['p'], d['fn'])]; n_nb = int(a['n_nucleobase_spectra'])
    pe = n_nb * (float(a['id_rate_expected_nucleobase']) - float(a['id_rate_nucleobase']))
    unadjusted = pe / max(float(a['id_rate_all']) * int(a['n_ms2']) + pe, 1)
    res.append((d['p'], d['fn'], (d['nb'] >= 0.10).mean(), unadjusted, adjusted))
for p, fn, f, s0, s1 in sorted(res, key=lambda t: -t[2]):
    lines.append('   %s %-34s nb in %5.1f %% of MS2 | unadjusted loss %5.1f %% | adjusted %5.1f %%' % (p, fn[:34], 100 * f, 100 * s0, 100 * s1))
s0 = np.array([x[3] for x in res]); s1 = np.array([x[4] for x in res]); f = np.array([x[2] for x in res])
lines.append('   median: unadjusted %.1f %%, adjusted %.1f %%; in files with >50 %% of MS2 with nb: unadjusted %s | adjusted %s'
             % (100 * np.median(s0), 100 * np.median(s1), ', '.join('%.0f' % (100 * x) for x in s0[f > 0.5]), ', '.join('%.0f' % (100 * x) for x in s1[f > 0.5])))
print('\n'.join(lines))
open('matching_results.txt', 'a', encoding='utf-8').write('\n'.join(lines) + '\n')
