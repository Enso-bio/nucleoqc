#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Impact of contamination without peptide identification (manuscript: spectra section, Figure S2).

For each MS2 spectrum: amino-acid ladder (fraction of the 40 most intense peaks above m/z 140 linked to another peak by
a residue mass, tol 0.01, minus the same with masses shifted +0.37 Da) and nucleotide ladder (RNA/DNA residue masses).
Compared between spectra with and without nucleobase in contaminated atlas files and clean files of the same instrument
type.

Usage: python impact.py   -> impact_spectra.npz"""
import re, random, collections
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import nucleic_census as C
import atlas

AA = np.array([57.02146, 71.03711, 87.03203, 97.05276, 99.06841, 101.04768, 160.03065, 113.08406, 114.04293,
               115.02694, 128.05858, 128.09496, 129.04259, 131.04049, 137.05891, 147.06841, 156.10111, 163.06333,
               186.07931, 147.03540])
NUC = np.array([329.05252, 345.04744, 305.04129, 306.02530,      # RNA: pA pG pC pU
                313.05761, 289.04637, 304.04604])                # DNA: dA dC dT (dG = pA)
TOL, NULL_SHIFT = 0.01, 0.37


def ladder(mz, it, masses, top=40):
    sel = mz > 140
    m, i = mz[sel], it[sel]
    if len(m) < 8:
        return np.nan, np.nan
    o = np.argsort(-i)[:top]; m = np.sort(m[o])
    d = m[None, :] - m[:, None]

    def frac(ms):
        hit = np.zeros(d.shape, bool)
        for a in ms:
            hit |= np.abs(d - a) < TOL
        return np.mean(hit.any(0) | hit.any(1))
    return frac(masses), frac(masses + NULL_SHIFT)


def file_spectra(project, fn):
    fs = [x for x in atlas.eligible_files(project) if x[0] == fn]
    if not fs:
        return []
    n, size, fmt, url = fs[0]
    offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * atlas.CHUNK else [0]
    spectra = []
    for a in offs:
        t = C.byte_range(url, a, min(size - 1, a + atlas.CHUNK - 1)).decode('latin-1')
        spectra += (C.mgf_spectra if fmt == 'mgf' else atlas.mzml_spectra)(t) or []
    C.recalibrate(spectra)
    return spectra


def measure(arg):
    group, project, fn = arg
    try:
        spectra = file_spectra(project, fn)
    except Exception as ex:
        return arg, None, str(ex)[:60]
    rows = []
    for s in spectra:
        nb = max(C.rel(s, C.BASES['A']), C.rel(s, C.BASES['G']))
        a, an = ladder(s['mz'], s['it'], AA)
        u, un = ladder(s['mz'], s['it'], NUC)
        rows.append((nb, a, an, u, un, s['pmz']))
    return arg, np.array(rows), ''


if __name__ == '__main__':
    rows = [l.split('\t') for l in open('atlas.tsv', encoding='utf-8').read().split('\n') if l]
    header = rows[0]; D = [dict(zip(header, f)) for f in rows[1:] if len(f) == len(header)]
    D = [d for d in D if d['valid'] == 'True']
    # projects with nucleic acid by design, same rule as the atlas analysis
    import json
    T = json.load(open('atlas_texts.json', encoding='utf-8'))
    by_design = {a for a, t in T.items() if re.search(atlas.BY_DESIGN, t)} | set(atlas.MANUAL_BY_DESIGN)
    fam = lambda s: next((k for k in ('Exploris', 'Lumos', 'Q Exactive', 'TripleTOF') if k.lower() in s.lower()), 'other')
    cont = [d for d in D if float(d['nb10']) > 0.01 and d['project'] not in by_design]
    clean = [d for d in D if float(d['nb10']) < 0.001 and d['project'] not in by_design]
    rng = random.Random(7)
    # controls: same number per instrument family, one file per project
    ctrl = []
    per_family = collections.Counter(fam(d['instrument']) for d in cont)
    for f, k in per_family.items():
        cand = {}
        for d in clean:
            if fam(d['instrument']) == f:
                cand.setdefault(d['project'], d)
        ctrl += rng.sample(list(cand.values()), min(k, len(cand)))
    tasks = [('contaminated', d['project'], d['file']) for d in cont] + [('clean', d['project'], d['file']) for d in ctrl]
    print('%d contaminated files, %d clean' % (len(cont), len(ctrl)), flush=True)
    res = {}
    with ThreadPoolExecutor(6) as ex:
        for (g, p, fn), arr, err in ex.map(measure, tasks):
            if arr is None or len(arr) < 100:
                print('  failed', p, fn[:40], err, flush=True); continue
            res['%s|%s|%s' % (g, p, fn)] = arr
            nbm = arr[:, 0] >= 0.10
            print('%-12s %s %-40s n=%5d nb=%5.1f %%  AA with/without nb %.3f/%.3f  NUC with/without nb %.3f/%.3f'
                  % (g, p, fn[:40], len(arr), 100 * nbm.mean(),
                     np.nanmean(arr[nbm, 1] - arr[nbm, 2]) if nbm.any() else np.nan, np.nanmean(arr[~nbm, 1] - arr[~nbm, 2]),
                     np.nanmean(arr[nbm, 3] - arr[nbm, 4]) if nbm.any() else np.nan, np.nanmean(arr[~nbm, 3] - arr[~nbm, 4])), flush=True)
    np.savez_compressed('impact_spectra.npz', **res)
    print('DONE')
