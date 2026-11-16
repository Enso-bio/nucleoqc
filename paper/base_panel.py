#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Why adenine and guanine? Comparison with an extended panel of nucleobase ions (manuscript: Methods and Results,
Table S10).
 a) Atlas (692 valid files): A, G and C (>=2 % of base peak) already measured per file.
 b) Stored spectra of the 38 files searched with Sage: A, G, C, U, T and the ribose ion (97.0284), each with a chance-match
    control (same nominal mass shifted +0.035 Da).
 c) Type of nucleic acid: RNA- vs DNA-specific signals in contaminated files.
Usage (from this folder): python base_panel.py  -> base_panel_results.txt"""
import os, zlib, pickle
import numpy as np
from scipy import stats
import nucleic_census as NC

H = 1.007276
IONS = {'A': 135.05450 + H, 'G': 151.04941 + H, 'C': 111.04326 + H, 'U': 112.02728 + H, 'T': 126.04293 + H, 'ribose': 97.02841}
TOL, SHIFT = 0.005, 0.035
lines = []
out = lambda s='': (print(s, flush=True), lines.append(s))
tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])

# ---------------------------------------------------------------- a) atlas
V = [d for d in tab('atlas.tsv') if d['valid'] == 'True']
A2 = np.array([float(d['A2']) for d in V]); G2 = np.array([float(d['G2']) for d in V]); C2 = np.array([float(d['C2']) for d in V])
nb2 = np.array([float(d['nb2']) for d in V]); nb10 = np.array([float(d['nb10']) for d in V])
out('a) Atlas, %d valid files: fraction of MS2 with each ion >= 2 %% of the base peak' % len(V))
out('   median per file: A %.3f %%, G %.3f %%, C %.3f %%' % (100 * np.median(A2), 100 * np.median(G2), 100 * np.median(C2)))
out('   in contaminated files (nb10 > 1 %%, n = %d): A %.1f %%, G %.1f %%, C %.1f %% (medians)'
    % ((nb10 > 0.01).sum(), 100 * np.median(A2[nb10 > 0.01]), 100 * np.median(G2[nb10 > 0.01]), 100 * np.median(C2[nb10 > 0.01])))
out('   Spearman between files: C2 vs nb2 %.2f; A2 vs G2 %.2f' % (stats.spearmanr(C2, nb2)[0], stats.spearmanr(A2, G2)[0]))
for u in (0.01, 0.05):
    ag = nb2 > u; c = C2 > u
    out('   threshold %.0f %%: positive by A/G %d, by C %d; C positive and A/G negative: %d; A/G positive and C negative: %d'
        % (100 * u, ag.sum(), c.sum(), (c & ~ag).sum(), (ag & ~c).sum()))

# ---------------------------------------------------------------- b) full spectra of the searched files


def rel(s, m):
    j = np.searchsorted(s['mz'], m); b = 0.0
    for k in (j - 1, j):
        if 0 <= k < len(s['mz']) and abs(s['mz'][k] - m) <= TOL:
            b = max(b, s['it'][k])
    return b / s['it'].max()


rows = []
for a in tab('identifications_1p.tsv'):
    f = os.path.join(NC.WORK, '%s_%d.pkl' % (a['project'], zlib.crc32(a['file'].encode()) % 100000))
    if not os.path.exists(f):
        continue
    spectra = pickle.load(open(f, 'rb'))
    # does the scan range reach cytosine (112.05)? 5th percentile of the first recorded peak, below the cytosine m/z
    low = np.percentile([s['mz'][0] for s in spectra], 5) < IONS['C'] - TOL
    r = {'p': a['project'], 'group': a['group'], 'low': low}
    for k, m in IONS.items():
        x = np.array([rel(s, m) for s in spectra]); y = np.array([rel(s, m + SHIFT) for s in spectra])
        r[k + '2'] = np.mean(x >= 0.02); r[k + '10'] = np.mean(x >= 0.10); r[k + '_ctrl2'] = np.mean(y >= 0.02)
        r[k + '_vec'] = x
    ag = np.maximum(r['A_vec'], r['G_vec'])
    r['nb10'] = np.mean(ag >= 0.10)
    r['only_others10'] = np.mean((ag < 0.10) & (np.maximum.reduce([r['C_vec'], r['U_vec'], r['T_vec']]) >= 0.10))
    for k in IONS:
        del r[k + '_vec']
    rows.append(r)
cont = [r for r in rows if r['nb10'] > 0.01]; clean = [r for r in rows if r['nb10'] <= 0.001]
out('\nb) Full spectra of %d files (%d contaminated, %d clean); %% of MS2 with the ion >= 2 %% of the base peak (median per file)'
    % (len(rows), len(cont), len(clean)))
out('   %-7s %12s %12s %16s %16s' % ('ion', 'contaminated', 'clean', 'control +0.035', 'r with nb10'))
nbv = np.array([r['nb10'] for r in rows])
for k in IONS:
    out('   %-7s %11.2f %% %11.2f %% %15.2f %% %16.2f' % (k, 100 * np.median([r[k + '2'] for r in cont]), 100 * np.median([r[k + '2'] for r in clean]),
        100 * np.median([r[k + '_ctrl2'] for r in cont]), stats.spearmanr([r[k + '2'] for r in rows], nbv)[0]))
out('   intensity >= 10 %% (contaminated, median): ' + ', '.join('%s %.2f %%' % (k, 100 * np.median([r[k + '10'] for r in cont])) for k in IONS))
out('   spectra with C, U or T >= 10 %% but neither A nor G >= 10 %%: median %.3f %% of MS2 (maximum %.2f %%)'
    % (100 * np.median([r['only_others10'] for r in cont]), 100 * max(r['only_others10'] for r in cont)))
out('   scan range reaching cytosine (5th percentile of the first peak < m/z %.3f): %d of %d files; cannot record it: %d'
    % (IONS['C'] - TOL, sum(r['low'] for r in rows), len(rows), sum(not r['low'] for r in rows)))

# ---------------------------------------------------------------- c) RNA or DNA
out('\nc) RNA (uracil, ribose) versus DNA (thymine) in contaminated files')
for r in sorted(cont, key=lambda r: -r['nb10']):
    out('   %s nb10 %5.1f %% | U %5.2f %% (control %.2f) | T %5.2f %% (control %.2f) | ribose %5.2f %%'
        % (r['p'], 100 * r['nb10'], 100 * r['U2'], 100 * r['U_ctrl2'], 100 * r['T2'], 100 * r['T_ctrl2'], 100 * r['ribose2']))
open('base_panel_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
