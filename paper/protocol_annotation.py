#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protocol annotation (manuscript: protocol features, Figure 2).
Step 1: downloads the sampleProcessingProtocol field of every valid atlas project -> annotation_protocols.json
Step 2: extracts snippets (+-90 characters) mentioning sonication, nuclease or phosphopeptide enrichment for manual
        review -> annotation_snippets.txt
Step 3 (with annotation_manual.json written after the review): re-analyses the features -> protocol_annotation_results.txt
Usage: python protocol_annotation.py [analyse]"""
import sys, os, re, json, math, collections
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import stats
import nucleic_census as C
import atlas

PAT = {'sonication': r'sonicat|sonif|bioruptor|ultrason|covaris|probe tip|tip sonic',
       'nuclease': r'benzonase|nuclease|dnase|rnase(?! ?(inhibitor|out|in\b))|micrococcal|turbonuclease|salt active nuclease|denarase',
       'phospho_enrichment': r'phospho(peptide|-?enrich|proteom)|tio2|titanium|imac|fe-nta|fe\(iii\)|fe3\+|ti4\+|ti-imac|zr-imac|mag resyn'}
NEG = r'(without|no|omitt?ed|not|w/o|instead of)\W+(\w+\W+){0,3}$'

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])
val = [d for d in tab('atlas.tsv') if d['valid'] == 'True']
projects = sorted({d['project'] for d in val})
CACHE = 'annotation_protocols.json'
prot = json.load(open(CACHE, encoding='utf-8')) if os.path.exists(CACHE) else {}
missing = [p for p in projects if p not in prot]
if missing:
    with ThreadPoolExecutor(8) as ex:
        for p, j in zip(missing, ex.map(lambda a: C.get_json(C.API + 'projects/' + a), missing)):
            if isinstance(j, dict) and j:
                prot[p] = {'sample': str(j.get('sampleProcessingProtocol', '')), 'data': str(j.get('dataProcessingProtocol', '')),
                           'title': str(j.get('title', '')), 'description': str(j.get('projectDescription', ''))}
    json.dump(prot, open(CACHE, 'w', encoding='utf-8'))
print('%d valid projects, %d with protocol downloaded, %d with a non-empty sample protocol'
      % (len(projects), len(prot), sum(len(v['sample']) > 30 for v in prot.values())))

if len(sys.argv) < 2:
    with open('annotation_snippets.txt', 'w', encoding='utf-8') as f:
        for k, pat in PAT.items():
            f.write('\n######## %s\n' % k)
            for p in projects:
                t = prot.get(p, {}).get('sample', '')
                for m in re.finditer(pat, t, re.I):
                    a, b = max(0, m.start() - 90), min(len(t), m.end() + 90)
                    neg = bool(re.search(NEG, t[max(0, m.start() - 40):m.start()], re.I))
                    f.write('%s %s| %s\n' % (p, 'NEG? ' if neg else '', t[a:b].replace('\n', ' ')))
    print('snippets: annotation_snippets.txt')
    sys.exit()

# ---------------------------------------------------------------- step 3: re-analysis with the manual annotation
M = json.load(open('annotation_manual.json', encoding='utf-8'))
texts = json.load(open('atlas_texts.json', encoding='utf-8'))
by_project = collections.defaultdict(list)
for d in val:
    by_project[d['project']].append(d)
P = []
for p, ds in by_project.items():
    if p in atlas.MANUAL_BY_DESIGN or re.search(atlas.BY_DESIGN, texts.get(p, '')):
        continue
    t = prot.get(p, {}).get('sample', '')
    r = {'p': p, 'nb': max(float(d['nb10']) for d in ds), 'year': ds[0]['year']}
    r['sonication'] = int(bool(re.search(PAT['sonication'], t, re.I)) and p not in M['sonication_false_positives'])
    r['nuclease'] = int(p in M['nuclease_confirmed'])
    r['phospho_enrichment'] = int(p in M['phospho_confirmed'])
    for k in ('tissue', 'bacteria', 'isobaric_labelling', 'ip_ap_ms', 'nucleic_acid_binding', 'body_fluids', 'in_gel_digestion',
              'sp3_strap_fasp', 'glyco', 'acetyl_ubiquitin', 'cross_linking'):
        r[k] = int(ds[0][k] == '1')
    r['son_without_nuc'] = int(r['sonication'] and not r['nuclease'])
    P.append(r)
lines = []
out = lambda s='': (print(s, flush=True), lines.append(s))
out('Manual annotation of the sample protocol (%d projects without nucleic acid by design)' % len(P))
for k in ('sonication', 'nuclease', 'phospho_enrichment'):
    out('   %-18s %3d projects' % (k, sum(r[k] for r in P)))


def logistic(X, y, lam=0.5):
    from scipy.optimize import minimize
    X = np.column_stack([np.ones(len(y)), X]); pen = np.r_[0, np.full(X.shape[1] - 1, lam)]
    f = lambda b: -(y * (X @ b) - np.logaddexp(0, X @ b)).sum() + 0.5 * (pen * b * b).sum()
    g = lambda b: -X.T @ (y - 1 / (1 + np.exp(-X @ b))) + pen * b
    b = minimize(f, np.zeros(X.shape[1]), jac=g, method='BFGS').x
    pr = 1 / (1 + np.exp(-X @ b)); H = X.T @ (X * (pr * (1 - pr))[:, None]) + np.diag(pen)
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return b[1:], se[1:], 2 * stats.norm.sf(np.abs(b[1:] / se[1:]))


VARS = ['sonication', 'nuclease', 'phospho_enrichment', 'tissue', 'bacteria', 'isobaric_labelling', 'ip_ap_ms', 'nucleic_acid_binding',
        'body_fluids', 'in_gel_digestion', 'sp3_strap_fasp', 'glyco', 'acetyl_ubiquitin']
VARS = [v for v in VARS if sum(r[v] for r in P) >= 8]
X = np.array([[r[v] for v in VARS] + [float(r['year'] >= '2018')] for r in P], float)
KEY = ('sonication', 'nuclease', 'phospho_enrichment', 'bacteria', 'tissue')
for u in (0.01, 0.02, 0.05):
    y = np.array([r['nb'] > u for r in P], float)
    out('\nContaminated = nb10 > %.0f %% (%d projects)' % (100 * u, int(y.sum())))
    for v in KEY:
        a = sum(1 for r, yy in zip(P, y) if r[v] and yy); n1 = sum(r[v] for r in P)
        c = sum(1 for r, yy in zip(P, y) if not r[v] and yy); n0 = len(P) - n1
        orr, pf = stats.fisher_exact([[a, n1 - a], [c, n0 - c]])
        out('   %-18s %2d/%-3d (%4.1f %%) vs %2d/%-3d (%4.1f %%)  Fisher p = %.3f' % (v, a, n1, 100 * a / max(n1, 1), c, n0, 100 * c / n0, pf))
    b, se, pv = logistic(X, y)
    out('   multivariable: ' + ', '.join('%s OR %.2f (95%% CI %.2f-%.2f; p %.3f)' % (v, math.exp(b[j]), math.exp(b[j] - 1.96 * se[j]), math.exp(b[j] + 1.96 * se[j]), pv[j])
                                         for j, v in enumerate(VARS) if v in KEY))
# sonication without nuclease versus sonication with nuclease
y = np.array([r['nb'] > 0.01 for r in P], float)
for name, sel in (('sonication without nuclease', lambda r: r['sonication'] and not r['nuclease']), ('sonication with nuclease', lambda r: r['sonication'] and r['nuclease']),
                  ('nuclease (any)', lambda r: r['nuclease']), ('neither sonication nor nuclease', lambda r: not r['sonication'] and not r['nuclease'])):
    s = [yy for r, yy in zip(P, y) if sel(r)]
    out('   %-31s %2d/%-3d contaminated (%.1f %%)' % (name, int(sum(s)), len(s), 100 * np.mean(s) if s else np.nan))
open('protocol_annotation_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
