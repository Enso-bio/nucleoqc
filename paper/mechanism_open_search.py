#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mechanism: are spectra with nucleobase peptides carrying a nucleotide adduct? (manuscript: software rescue section)

Open search with Sage (precursor +-500 Da, no variable modifications except Met oxidation; fixed Cys and labels) on the
most contaminated files. For each identified spectrum (1 % FDR by hyperscore, best hit over all shards) the mass shift
(experimental - theoretical) is taken, and the frequency of shifts matching nucleic-acid components is compared between
spectra with and without nucleobase. Internal control: phosphorylation (+79.966) must appear as a shift.

Usage: python mechanism_open_search.py SAGE_EXE FASTA  -> mechanism_open_search_results.txt"""
import sys, os, re, json, zlib, pickle, subprocess, collections
SAGE, FASTA = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], SAGE, FASTA, '--one-step']
import numpy as np
from scipy import stats
import identifications as I

FILES = [('PXD023674', 'Katy_phospho_3b_R.raw.-1.mgf'), ('PXD023674', 'Katy_phospho__2a_R.raw.-1.mgf'),
         ('PXD034105', 'QE1_Ateeq_Exp1_Sample19_R.raw.-1.mgf'), ('PXD034105', 'QE1_Ateeq_Exp1_Sample18_R.raw.-1.mgf'),
         ('PXD015943', '61caa4abd7fc8222___b1790p079_PP1_NTERA_AZD5363_r1.raw.-1.mgf'),
         ('PXD015943', '65563538acf34a7e___b1790p087_MCF7_DMSO_PP2_r1.raw.-1.mgf'),
         ('PXD021039', 'R16058_15_M2.mgf'), ('PXD072006', 'CE3_2024Jan30_Jelinek_P23193_PO3_Fr.mzML')]
NEUTRALS = {'adenine': 135.05450, 'guanine': 151.04941, 'cytosine': 111.04326, 'uracil': 112.02728, 'thymine': 126.04293,
            'adenosine': 267.09675, 'guanosine': 283.09167, 'cytidine': 243.08552, 'uridine': 244.06954,
            'AMP': 347.06308, 'GMP': 363.05800, 'CMP': 323.05185, 'UMP': 324.03587,
            'dAMP': 331.06817, 'dGMP': 347.06308, 'dCMP': 307.05694, 'dTMP': 322.05660}
H2O = 18.010565
CAND = dict(NEUTRALS)
for k in ('AMP', 'GMP', 'CMP', 'UMP', 'dAMP', 'dCMP', 'dTMP', 'adenosine', 'guanosine', 'cytidine', 'uridine'):
    CAND[k + '-H2O'] = NEUTRALS[k] - H2O          # attachment by condensation
CAND['adenine-H2O'] = NEUTRALS['adenine'] - H2O
TOL = 0.02
PHOSPHO = 79.96633

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def open_search(base, n_spec, label_mass):
    n = 4
    best = np.full(n_spec, -np.inf); target = np.zeros(n_spec, bool); delta = np.full(n_spec, np.nan)
    for k, fa in enumerate(I.shards(n)):
        out_dir = '%s_ab_t%d' % (base, k); res = os.path.join(out_dir, 'results.sage.tsv')
        if not os.path.exists(res):
            c = I.config(base + '.mzML', out_dir, False, label_mass, fasta=fa)
            c['database']['variable_mods'] = {'M': [15.99491]}; c['database']['max_variable_mods'] = 1
            c['precursor_tol'] = {'da': [-500, 500]}; c['isotope_errors'] = [0, 0]
            json.dump(c, open(out_dir + '.json', 'w'))
            r = subprocess.run([SAGE, out_dir + '.json', '--disable-telemetry-i-dont-want-to-improve-sage'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(res):
                print('   shard %d failed: %s' % (k, r.stderr[-300:]), flush=True); return None
        lin = open(res, encoding='utf-8').read().split('\n'); J = {c: i for i, c in enumerate(lin[0].split('\t'))}
        for l in lin[1:]:
            if l:
                c = l.split('\t'); i = int(re.search(r'(\d+)$', c[J['scannr']]).group(1)); h = float(c[J['hyperscore']])
                if h > best[i]:
                    best[i] = h; target[i] = c[J['label']] == '1'; delta[i] = float(c[J['expmass']]) - float(c[J['calcmass']])
    ok = np.isfinite(best); idx = np.where(ok)[0][np.argsort(-best[ok])]
    t = np.cumsum(target[idx]); d = np.cumsum(~target[idx]); q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    ident = np.zeros(n_spec, bool); ident[idx] = target[idx] & (q <= 0.01)
    return ident, delta


def classify(dm):
    for name, m in CAND.items():
        if abs(dm - m) <= TOL:
            return name
    return None


if __name__ == '__main__':
    E = collections.defaultdict(dict)
    for r in tab('ids_spectra_1p.tsv'):
        E[(r['project'], r['file'])][int(r['index'])] = (float(r['nb']), r['id'] == '1')
    lines = []
    out = lambda s: (print(s, flush=True), lines.append(s))
    tot = collections.Counter()
    all_nb, all_without = [], []
    for p, fn in FILES:
        base = os.path.join(I.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        spectra = pickle.load(open(base + '.pkl', 'rb'))
        nb = np.array([E[(p, fn)][i][0] for i in range(len(spectra))]); id_closed = np.array([E[(p, fn)][i][1] for i in range(len(spectra))])
        lab, mass = I.labelling(spectra)
        r = open_search(base, len(spectra), mass)
        if r is None:
            continue
        ident, delta = r
        m = nb >= 0.10
        for g, sel in (('nb', m), ('without', ~m)):
            s = sel & ident
            nuc = sum(classify(x) is not None for x in delta[s]); phos = int((np.abs(delta[s] - PHOSPHO) <= TOL).sum())
            zero = int((np.abs(delta[s]) <= TOL).sum())
            tot[(g, 'n')] += int(s.sum()); tot[(g, 'nuc')] += nuc; tot[(g, 'phos')] += phos; tot[(g, 'zero')] += zero
            for x in delta[s]:
                c = classify(x)
                if c:
                    tot[(g, c)] += 1
            (all_nb if g == 'nb' else all_without).extend(delta[s].tolist())
        a_, b_ = ident[m].mean(), ident[~m].mean()
        out('%s %-34s open: ID with nb %5.1f %% vs without nb %5.1f %% (closed %5.1f / %5.1f) | rescued by the open search among nb: %d'
            % (p, fn[:34], 100 * a_, 100 * b_, 100 * id_closed[m].mean(), 100 * id_closed[~m].mean(), int((ident & m & ~id_closed).sum())))
    out('')
    for g, name in (('nb', 'with nucleobase'), ('without', 'without nucleobase')):
        n = max(tot[(g, 'n')], 1)
        out('%s: %d identified | no shift %.1f %% | phosphorylation (+79.966) %.1f %% | nucleic-acid shift %.2f %% (%d)'
            % (name, tot[(g, 'n')], 100 * tot[(g, 'zero')] / n, 100 * tot[(g, 'phos')] / n, 100 * tot[(g, 'nuc')] / n, tot[(g, 'nuc')]))
    tab2 = [[tot[('nb', 'nuc')], tot[('nb', 'n')] - tot[('nb', 'nuc')]], [tot[('without', 'nuc')], tot[('without', 'n')] - tot[('without', 'nuc')]]]
    orr, pv = stats.fisher_exact(tab2)
    out('Fisher (nucleic-acid shift, with vs without nucleobase): OR %.2f, p = %.2g' % (orr, pv))
    out('Breakdown: ' + ', '.join('%s %d/%d' % (c, tot[('nb', c)], tot[('without', c)]) for c in CAND if tot[('nb', c)] + tot[('without', c)]))
    # shifts most enriched in spectra with nucleobase (0.01 Da bins, |dm| > 0.5)
    bn = collections.Counter(round(x, 2) for x in all_nb if abs(x) > 0.5); bs = collections.Counter(round(x, 2) for x in all_without if abs(x) > 0.5)
    N1, N0 = max(len(all_nb), 1), max(len(all_without), 1)
    enr = []
    for k in set(bn) | set(bs):
        a, b = bn[k] + sum(bn[k + d] for d in (-0.01, 0.01)), bs[k] + sum(bs[k + d] for d in (-0.01, 0.01))
        if a >= 5:
            enr.append(((a + .5) / N1 / ((b + .5) / N0), k, a, b))
    out('Shifts most enriched in spectra with nucleobase (>=5 cases, window +-0.01):')
    for e, k, a, b in sorted(enr, reverse=True)[:12]:
        out('   %+9.3f Da  with nb %d, without nb %d, enrichment %.1fx %s' % (k, a, b, e, ('(' + classify(k) + ')') if classify(k) else ''))
    open('mechanism_open_search_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
