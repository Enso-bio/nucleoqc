#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Representativeness of the atlas. The atlas can only measure projects with MGF/mzML; most of PRIDE deposits raw
files only. A new random PRIDE sample (same procedure, another seed) is split into high-resolution projects with peak-list
files and raw-only projects, compared on year, instrument family and protocol features (same rules as the atlas), and the
atlas model (penalised logistic regression) is applied to raw-only projects to estimate their expected contamination.
Usage: python representativeness.py  -> representativeness.tsv, representativeness_results.txt"""
import os, re, json, math, random, collections, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import stats
import nucleic_census as C
import atlas

OUT = 'representativeness.tsv'
TARGET = 700
RNG = random.Random(20261114)
FAM = ('Astral', 'timsTOF', 'TripleTOF', 'ZenoTOF', 'Exploris', 'Eclipse', 'Lumos', 'Fusion', 'Ascend', 'Q Exactive', 'Elite', 'Velos')
COLS = ['project', 'group', 'year', 'instrument_family'] + list(atlas.FLAGS)


def record(p):
    acc = p.get('accession'); ins = atlas.as_text(p.get('instruments'))
    if not re.search(atlas.INSTR, ins, re.I):
        return None
    t = C.get_json(C.API + 'files/getCountOfFilesByType/' + acc)
    if not isinstance(t, dict):
        return None
    peak, raw = t.get('PEAK', 0) or 0, t.get('RAW', 0) or 0
    if not peak and not raw:
        return None
    txt = ' '.join(str(p.get(k, '')) for k in ('title', 'projectDescription', 'sampleProcessingProtocol', 'keywords')).lower()
    r = {'project': acc, 'group': 'peak_lists' if peak else 'raw_only', 'year': str(p.get('submissionDate', ''))[:4],
         'instrument_family': next((f for f in FAM if f.lower() in ins.lower()), 'other')}
    for k, pat in atlas.FLAGS.items():
        r[k] = int(bool(re.search(pat, txt)))
    return r


if __name__ == '__main__':
    done = {l.split('\t')[0] for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
    if not done:
        open(OUT, 'w', encoding='utf-8').write('\t'.join(COLS) + '\n')
    n = len(done)
    pages = list(range(int(C.get_json(C.API + 'projects/count') or 42000) // 100)); RNG.shuffle(pages)
    with ThreadPoolExecutor(8) as ex:
        for pg in pages:
            if n >= TARGET:
                break
            ps = C.get_json(C.API + 'projects?pageSize=100&page=%d' % pg)
            ps = [p for p in (ps if isinstance(ps, list) else []) if p.get('accession') not in done]
            random.Random(pg).shuffle(ps)
            for r in ex.map(record, ps[:15]):
                if r:
                    open(OUT, 'a', encoding='utf-8').write('\t'.join(str(r[k]) for k in COLS) + '\n'); done.add(r['project']); n += 1
            print('page %d: %d projects' % (pg, n), flush=True)

    # ---- analysis
    L = [l.split('\t') for l in open(OUT, encoding='utf-8').read().split('\n') if l]
    R = [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])]
    lines = []
    out = lambda s='': (print(s, flush=True), lines.append(s))
    G = {g: [r for r in R if r['group'] == g] for g in ('peak_lists', 'raw_only')}
    out('New random PRIDE sample: %d high-resolution projects (%d with peak-list files, %d raw only)'
        % (len(R), len(G['peak_lists']), len(G['raw_only'])))
    out('Feature (project text)         peak lists   raw only   Fisher p')
    for k in list(atlas.FLAGS) + ['year>=2018']:
        f = (lambda r: r['year'] >= '2018') if k == 'year>=2018' else (lambda r, k=k: r[k] == '1')
        a = sum(f(r) for r in G['peak_lists']); b = sum(f(r) for r in G['raw_only'])
        p = stats.fisher_exact([[a, len(G['peak_lists']) - a], [b, len(G['raw_only']) - b]])[1]
        out('   %-26s %5.1f %%    %5.1f %%    %.3f' % (k, 100 * a / len(G['peak_lists']), 100 * b / len(G['raw_only']), p))
    out('Instrument family:')
    for f in FAM + ('other',):
        a = sum(r['instrument_family'] == f for r in G['peak_lists']); b = sum(r['instrument_family'] == f for r in G['raw_only'])
        if a + b:
            out('   %-12s %5.1f %%  %5.1f %%' % (f, 100 * a / len(G['peak_lists']), 100 * b / len(G['raw_only'])))

    # atlas model (text features + year) applied to the raw-only projects
    A = [l.split('\t') for l in open('atlas.tsv', encoding='utf-8').read().split('\n') if l]
    header = A[0]; A = [dict(zip(header, r)) for r in A[1:] if len(r) == len(header)]
    by_project = collections.defaultdict(list)
    for d in A:
        if d['valid'] == 'True':
            by_project[d['project']].append(d)
    feats = list(atlas.FLAGS)
    Xa = np.array([[float(ds[0][k] == '1') for k in feats] + [float(ds[0]['year'] >= '2018')] for ds in by_project.values()])
    ya = np.array([max(float(d['nb10']) for d in ds) > 0.01 for ds in by_project.values()], float)
    from scipy.optimize import minimize
    Xc = np.column_stack([np.ones(len(ya)), Xa]); pen = np.r_[0, np.full(Xa.shape[1], 0.5)]
    fobj = lambda b: -(ya * (Xc @ b) - np.logaddexp(0, Xc @ b)).sum() + 0.5 * (pen * b * b).sum()
    b = minimize(fobj, np.zeros(Xc.shape[1]), method='BFGS').x
    pred = lambda g: np.mean(1 / (1 + np.exp(-np.column_stack([np.ones(len(G[g])), np.array([[float(r[k] == '1') for k in feats] + [float(r['year'] >= '2018')] for r in G[g]])]) @ b)))
    out('\nAtlas model (all valid projects, nb10 > 1 %%; observed %.1f %%) applied to the new sample:' % (100 * ya.mean()))
    out('   expected frequency: peak lists %.1f %% | raw only %.1f %%' % (100 * pred('peak_lists'), 100 * pred('raw_only')))
    open('representativeness_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
