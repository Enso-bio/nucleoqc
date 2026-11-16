#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sensitivity of the conclusions to the thresholds chosen (manuscript: Results, robustness).
 a) frequency in PRIDE with other project cut-offs (nb10 and nb2; maximum or mean of files)
 b) protocol features with other contamination cut-offs and with the continuous measure
 c) identification effect with other per-spectrum thresholds
No network: uses atlas.tsv, atlas_texts.json, identifications_1p.tsv and ids_spectra_1p.tsv.
-> sensitivity_results.txt"""
import re, json, math, collections
import numpy as np
from scipy import stats
import atlas

lines = []
out = lambda s='': (print(s, flush=True), lines.append(s))
tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def wilson(k, n, z=1.96):
    p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0, c - h), min(1, c + h)


# ---------------------------------------------------------------- projects (same construction as atlas_analysis.py)
val = [d for d in tab('atlas.tsv') if d['valid'] == 'True']
proj = collections.OrderedDict()
for d in val:
    proj.setdefault(d['project'], []).append(d)
texts = json.load(open('atlas_texts.json', encoding='utf-8'))
P = []
for acc, ds in proj.items():
    m = dict(ds[0])
    m['nb10_max'] = max(float(d['nb10']) for d in ds); m['nb10_mean'] = float(np.mean([float(d['nb10']) for d in ds]))
    m['nb2_max'] = max(float(d['nb2']) for d in ds)
    m['dia'] = '1' if re.search(r'data[- ]independent|\bdia\b|swath|diann|dia-nn|spectronaut', texts.get(acc, '')) else '0'
    if re.search(r'dia|swath', ' '.join(d['file'] for d in ds), re.I):
        m['dia'] = '1'
    m['by_design'] = acc in atlas.MANUAL_BY_DESIGN or bool(re.search(atlas.BY_DESIGN, texts.get(acc, '')))
    P.append(m)
Pn = [p for p in P if not p['by_design']]
out('a) Frequency of unexplained contamination (n = %d projects) by measure and cut-off' % len(Pn))
for meas, cuts in (('nb10_max', (0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20)), ('nb10_mean', (0.01, 0.05)), ('nb2_max', (0.01, 0.05, 0.10, 0.20))):
    x = np.array([p[meas] for p in Pn])
    for u in cuts:
        k = int((x > u).sum()); lo, hi = wilson(k, len(x))
        out('   %-9s > %5.1f %%: %3d projects = %5.1f %% (95%% CI %.1f-%.1f)' % (meas, 100 * u, k, 100 * k / len(x), 100 * lo, 100 * hi))
x = np.array([p['nb10_max'] for p in Pn]); y = np.array([p['nb2_max'] for p in Pn])
out('   Spearman nb10 vs nb2 per project: %.2f' % stats.spearmanr(x, y)[0])


# ---------------------------------------------------------------- b) features
def logistic(X, y, lam=0.5):
    from scipy.optimize import minimize
    X = np.column_stack([np.ones(len(y)), X]); pen = np.r_[0, np.full(X.shape[1] - 1, lam)]
    f = lambda b: -(y * (X @ b) - np.logaddexp(0, X @ b)).sum() + 0.5 * (pen * b * b).sum()
    g = lambda b: -X.T @ (y - 1 / (1 + np.exp(-X @ b))) + pen * b
    b = minimize(f, np.zeros(X.shape[1]), jac=g, method='BFGS').x
    pr = 1 / (1 + np.exp(-X @ b)); H = X.T @ (X * (pr * (1 - pr))[:, None]) + np.diag(pen)
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return b[1:], se[1:], 2 * stats.norm.sf(np.abs(b[1:] / se[1:]))


FEATURES = [f for f in list(atlas.FLAGS) + ['dia'] if sum(p[f] == '1' for p in Pn) >= 8]
KEY = ['sonication', 'phospho_enrichment', 'tissue', 'bacteria', 'nuclease', 'isobaric_labelling']
X = np.array([[float(p[f] == '1') for f in FEATURES] + [float(p['year'] >= '2018')] for p in Pn])
out('\nb) Protocol features: multivariable OR (penalised logistic regression, the same %d features + year) with different cut-offs' % len(FEATURES))
out('   %-11s' % 'cut-off' + ''.join('%-20s' % f for f in KEY) + 'contaminated')
for u in (0.005, 0.01, 0.02, 0.05):
    y = np.array([p['nb10_max'] > u for p in Pn], float)
    b, se, pv = logistic(X, y)
    cells = []
    for f in KEY:
        j = FEATURES.index(f)
        cells.append('%.2f (p %.3f)' % (math.exp(b[j]), pv[j]))
    out('   nb10 > %3.1f%% ' % (100 * u) + ''.join('%-20s' % c for c in cells) + '%d' % int(y.sum()))
y = np.array([p['nb2_max'] > 0.05 for p in Pn], float); b, se, pv = logistic(X, y)
out('   nb2  > 5.0% ' + ''.join('%-20s' % ('%.2f (p %.3f)' % (math.exp(b[FEATURES.index(f)]), pv[FEATURES.index(f)])) for f in KEY) + '%d' % int(y.sum()))
# continuous measure: linear regression of log10(nb10 + 1e-4) on the same features
yl = np.log10(np.array([p['nb10_max'] for p in Pn]) + 1e-4)
Xl = np.column_stack([np.ones(len(yl)), X]); beta, *_ = np.linalg.lstsq(Xl, yl, rcond=None)
res = yl - Xl @ beta; s2 = res @ res / (len(yl) - Xl.shape[1]); se = np.sqrt(np.diag(s2 * np.linalg.inv(Xl.T @ Xl)))
pv = 2 * stats.t.sf(np.abs(beta / se), len(yl) - Xl.shape[1])
out('   continuous (log10 nb10, OLS): ' + ', '.join('%s x%.2f (p %.3f)' % (f, 10 ** beta[1 + FEATURES.index(f)], pv[1 + FEATURES.index(f)]) for f in KEY))
out('   (x = multiplicative factor on nb10)')

# ---------------------------------------------------------------- c) identification with other per-spectrum thresholds
A = {(a['project'], a['file']): a for a in tab('identifications_1p.tsv') if float(a['id_rate_all']) >= 0.02 and a['group'] == 'contaminated'}
E = collections.defaultdict(list)
for r in tab('ids_spectra_1p.tsv'):
    if (r['project'], r['file']) in A:
        E[(r['project'], r['file'])].append((float(r['pmz']), float(r['nb']), r['id'] == '1'))
out('\nc) Identification: spectra with nucleobase >= threshold versus neighbours in the same file (+-10 Th)')
for ctrl_name, ctrl_fn in (('controls < threshold', lambda nb, t: nb < t), ('controls < 2 % fixed', lambda nb, t: nb < 0.02)):
    out('   %s:' % ctrl_name)
    for t in (0.02, 0.05, 0.10, 0.20, 0.30, 0.50):
        if ctrl_name.endswith('fixed') and t <= 0.02:
            continue
        oo, ee = [], []
        for k, v in E.items():
            pmz = np.array([x[0] for x in v]); nb = np.array([x[1] for x in v]); idn = np.array([x[2] for x in v])
            m = nb >= t; c = np.where(ctrl_fn(nb, t))[0]; o = np.argsort(pmz[c]); c, cp = c[o], pmz[c][o]
            obs, exp_ = [], []
            for i in np.where(m)[0]:
                a, b = np.searchsorted(cp, [pmz[i] - 10, pmz[i] + 10])
                if b - a >= 5:
                    obs.append(idn[i]); exp_.append(idn[c[a:b]].mean())
            if len(obs) >= 20:
                oo.append(np.mean(obs)); ee.append(np.mean(exp_))
        oo, ee = np.array(oo), np.array(ee)
        out('      nb >= %2.0f %%: %2d files | median ratio %.2f | lower in %2d of %2d | Wilcoxon p = %.2g'
            % (100 * t, len(oo), np.median(oo / ee), int((oo < ee).sum()), len(oo), stats.wilcoxon(oo, ee).pvalue))

open('sensitivity_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
