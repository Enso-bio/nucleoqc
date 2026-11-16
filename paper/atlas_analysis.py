#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analysis of atlas.tsv: frequency, associated protocol features, within-project concordance, Figure S1."""
import math, collections, re
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import atlas

THRESHOLD = 0.01       # "contaminated" project: >1 % of spectra with adenine or guanine >=10 % of the base peak
rows = [l.split('\t') for l in open('atlas.tsv', encoding='utf-8').read().split('\n') if l]
header, rows = rows[0], rows[1:]
D = [dict(zip(header, f)) for f in rows if len(f) == len(header)]
val = [d for d in D if d['valid'] == 'True']


def wilson(k, n, z=1.96):
    if n == 0:
        return (0, 0)
    p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0, c - h), min(1, c + h))


proj = collections.OrderedDict()
for d in val:
    proj.setdefault(d['project'], []).append(d)
P = []
for acc, ds in proj.items():
    x = [float(d['nb10']) for d in ds]
    m = dict(ds[0]); m['nb10_max'] = max(x); m['nb10_mean'] = float(np.mean(x)); m['nfiles'] = len(x); m['vals'] = x
    m['cont'] = int(max(x) > THRESHOLD)
    P.append(m)

# extra features from the full project text (cached) and from the file name
import os, json
from concurrent.futures import ThreadPoolExecutor
import nucleic_census as C
CACHE = 'atlas_texts.json'
texts = json.load(open(CACHE, encoding='utf-8')) if os.path.exists(CACHE) else {}
missing = [p['project'] for p in P if p['project'] not in texts]
with ThreadPoolExecutor(8) as ex:
    for acc, j in zip(missing, ex.map(lambda a: C.get_json(C.API + 'projects/' + a), missing)):
        if isinstance(j, dict) and j:
            texts[acc] = ' '.join(str(j.get(k, '')) for k in ('title', 'projectDescription', 'sampleProcessingProtocol',
                                                                 'dataProcessingProtocol', 'keywords', 'experimentTypes')).lower()
json.dump(texts, open(CACHE, 'w', encoding='utf-8'))
EXTRA = {'dia': r'data[- ]independent|\bdia\b|swath|diann|dia-nn|spectronaut',
         'nucleotide_added': r'\batp\b|\badp\b|\bgtp\b|\bnad\+?\b|nucleotide|adenosine|lip-ms|limited proteolysis'}
# studies with nucleic acid or nucleotide by design (not contamination): rule applied to ALL projects + manual
# review of the contaminated ones (atlas.BY_DESIGN, atlas.MANUAL_BY_DESIGN)
for p in P:
    t = texts.get(p['project'], '')
    for k, pat in EXTRA.items():
        p[k] = '1' if re.search(pat, t) else '0'
    p['by_design'] = '1' if (p['project'] in atlas.MANUAL_BY_DESIGN or re.search(atlas.BY_DESIGN, t)) else '0'
    if re.search(r'dia|swath', ' '.join(d['file'] for d in proj[p['project']]), re.I):
        p['dia'] = '1'

print('Atlas: %d projects with an eligible file examined, %d files measured, %d valid files, %d valid projects'
      % (len({d['project'] for d in D}), len(D), len(val), len(P)))
print('Files discarded, by reason: %s' % dict(collections.Counter(d['reason'] for d in D if d['valid'] != 'True')))
x = np.array([p['nb10_max'] for p in P])
print('\nFraction of spectra with nucleobase (maximum per project): median %.3f %%, p90 %.2f %%, max %.1f %%'
      % (100 * np.median(x), 100 * np.percentile(x, 90), 100 * x.max()))
for u in (0.001, 0.01, 0.05):
    k = int((x > u).sum()); lo, hi = wilson(k, len(x))
    print('  projects > %4.1f %%: %3d of %d = %4.1f %%  (95%% CI %.1f-%.1f %%)' % (100 * u, k, len(x), 100 * k / len(x), 100 * lo, 100 * hi))

E = [p for p in P if p['by_design'] == '1']
print('\nStudies with nucleic acid/nucleotide by design: %d projects, %d above 1 %%' % (len(E), sum(p['cont'] for p in E)))
P_all = P
P = [p for p in P if p['by_design'] != '1']
x = np.array([p['nb10_max'] for p in P])
print('Contamination NOT explained by the design (n = %d):' % len(P))
for u in (0.001, 0.01, 0.05):
    k = int((x > u).sum()); lo, hi = wilson(k, len(x))
    print('  projects > %4.1f %%: %3d of %d = %4.1f %%  (95%% CI %.1f-%.1f %%)' % (100 * u, k, len(x), 100 * k / len(x), 100 * lo, 100 * hi))

print('\nFeatures (Fisher, project contaminated >1 %, by-design studies excluded; BH over all features):')
res = []
FEATURES = list(atlas.FLAGS) + ['dia']
for f in FEATURES:
    a = sum(1 for p in P if p[f] == '1' and p['cont']); b = sum(1 for p in P if p[f] == '1' and not p['cont'])
    c = sum(1 for p in P if p[f] != '1' and p['cont']); d = sum(1 for p in P if p[f] != '1' and not p['cont'])
    if a + b < 3:
        continue
    orr, pv = stats.fisher_exact([[a, b], [c, d]])
    # also on the continuous measure (Mann-Whitney on the maximum nb10)
    s1 = [p['nb10_max'] for p in P if p[f] == '1']; s0 = [p['nb10_max'] for p in P if p[f] != '1']
    pmw = stats.mannwhitneyu(s1, s0).pvalue if s1 and s0 else 1
    res.append([f, a, a + b, c, c + d, (a + .5) * (d + .5) / ((b + .5) * (c + .5)), pv, pmw])
pvs = np.array([r[6] for r in res]); o = np.argsort(pvs); q = np.empty(len(pvs))
q[o] = np.minimum.accumulate((pvs[o] * len(pvs) / np.arange(1, len(pvs) + 1))[::-1])[::-1]
for r, qq in sorted(zip(res, q), key=lambda t: t[0][6]):
    print('  %-21s with feature %2d/%-3d (%4.1f %%)  without %2d/%-3d (%4.1f %%)  OR %5.2f  p %.3f  q %.3f  MW p %.3f'
          % (r[0], r[1], r[2], 100 * r[1] / r[2], r[3], r[4], 100 * r[3] / max(r[4], 1), r[5], r[6], qq, r[7]))


def logistic(X, y, lam=0.5):
    """Logistic regression with a mild L2 penalty (stable with few events); OR, 95 % CI and Wald p."""
    from scipy.optimize import minimize
    X = np.column_stack([np.ones(len(y)), X]); pen = np.r_[0, np.full(X.shape[1] - 1, lam)]
    f = lambda b: -(y * (X @ b) - np.logaddexp(0, X @ b)).sum() + 0.5 * (pen * b * b).sum()
    g = lambda b: -X.T @ (y - 1 / (1 + np.exp(-X @ b))) + pen * b
    b = minimize(f, np.zeros(X.shape[1]), jac=g, method='BFGS').x
    pr = 1 / (1 + np.exp(-X @ b)); H = X.T @ (X * (pr * (1 - pr))[:, None]) + np.diag(pen)
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return b[1:], se[1:], 2 * stats.norm.sf(np.abs(b[1:] / se[1:]))


sel = [r[0] for r in res if r[2] >= 8]
Pm = P
X = np.array([[float(p[f] == '1') for f in sel] + [float(p['year'] >= '2018')] for p in Pm]); y = np.array([p['cont'] for p in Pm], float)
b, se, pv = logistic(X, y)
print('\nMultivariable model (penalised logistic regression, %d projects without by-design studies, %d contaminated):' % (len(y), int(y.sum())))
for name, bb, ss, pp in sorted(zip(sel + ['year>=2018'], b, se, pv), key=lambda t: t[3]):
    print('  %-21s OR %5.2f (95%% CI %.2f-%.2f)  p %.3f' % (name, math.exp(bb), math.exp(bb - 1.96 * ss), math.exp(bb + 1.96 * ss), pp))

print('\nBy instrument / year:')
fam = lambda s: next((k for k in ('Astral', 'timsTOF', 'TripleTOF', 'Exploris', 'Eclipse', 'Lumos', 'Fusion', 'Q Exactive', 'Elite', 'Velos') if k.lower() in s.lower()), 'other')
for name, key in (('instrument', lambda p: fam(p['instrument'])), ('period', lambda p: '<=2017' if p['year'] <= '2017' else ('2018-2021' if p['year'] <= '2021' else '>=2022'))):
    g = collections.defaultdict(list)
    for p in P:
        g[key(p)].append(p)
    for k, ps in sorted(g.items()):
        print('  %-11s %-11s n=%3d  contaminated %2d (%4.1f %%)  median %.3f %%' % (name, k, len(ps), sum(p['cont'] for p in ps),
              100 * np.mean([p['cont'] for p in ps]), 100 * np.median([p['nb10_max'] for p in ps])))

pairs = [p['vals'] for p in P if p['nfiles'] == 2]
if len(pairs) >= 10:
    a = np.array(pairs)
    rho = stats.spearmanr(a[:, 0], a[:, 1])
    both = sum(1 for u, v in pairs if u > THRESHOLD and v > THRESHOLD); one = sum(1 for u, v in pairs if (u > THRESHOLD) != (v > THRESHOLD))
    print('\nWithin-project concordance (%d projects with 2 files): Spearman %.2f (p %.1e); both contaminated %d, only one %d'
          % (len(pairs), rho.statistic, rho.pvalue, both, one))
    # variance between projects vs within (log10)
    l = np.log10(a + 1e-4); between = np.var(l.mean(1), ddof=1); within = np.mean(np.var(l, axis=1, ddof=1))
    print('  log10: variance between projects %.3f, within %.3f -> approximate ICC %.2f' % (between, within, (between - within / 2) / (between + within / 2)))

# figure
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
xs = np.sort(100 * x)
ax[0].plot(np.maximum(xs, 1e-3), 1 - np.arange(len(xs)) / len(xs), drawstyle='steps-post', color='#3a6ea5')
LABELS = {'phospho_enrichment': 'phospho enrichment', 'nuclease': 'nuclease', 'sonication': 'sonication', 'ip_ap_ms': 'IP / AP-MS',
          'nucleic_acid_binding': 'nucleic-acid binding', 'cross_linking': 'cross-linking', 'isobaric_labelling': 'isobaric labelling',
          'body_fluids': 'body fluids', 'tissue': 'tissue', 'in_gel_digestion': 'in-gel digestion', 'sp3_strap_fasp': 'SP3/S-Trap/FASP',
          'glyco': 'glyco', 'acetyl_ubiquitin': 'acetyl / ubiquitin', 'bacteria': 'bacteria', 'dia': 'DIA'}
ax[0].set_xscale('log'); ax[0].set_xlabel('% of MS2 spectra with nucleobase (maximum per project)'); ax[0].set_ylabel('fraction of projects $\\geq$ x')
ax[0].axvline(1, ls='--', c='gray', lw=.8); ax[0].set_title('PRIDE atlas, %d random projects' % len(P), fontsize=10)
rs = sorted(res, key=lambda r: r[5])
ax[1].barh(range(len(rs)), [math.log2(r[5]) for r in rs], color=['#c0504d' if r[6] < 0.05 else '#9bbbd9' for r in rs])
ax[1].set_yticks(range(len(rs))); ax[1].set_yticklabels(['%s (%d)' % (LABELS.get(r[0], r[0]), r[2]) for r in rs], fontsize=8)
ax[1].axvline(0, c='k', lw=.6); ax[1].set_xlabel('log2 odds ratio (red: Fisher p < 0.05)'); ax[1].set_title('Protocol feature vs contamination > 1 % (text rules)', fontsize=10)
plt.tight_layout(); plt.savefig('atlas_figure.png', dpi=150)
print('\nfigure: atlas_figure.png')
print('\nContaminated projects (>1 %):')
for p in sorted(P, key=lambda p: -p['nb10_max']):
    if p['cont']:
        print('  %s %s %-22s %5.1f %%  [%s]  %s' % (p['project'], p['year'], p['instrument'][:22], 100 * p['nb10_max'],
              ','.join(f for f in FEATURES if p[f] == '1'), p['title'][:60]))
