#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Software rescue of spectra with co-isolated nucleic acid (manuscript: software rescue section, Figure S5).
For each contaminated file already searched (identifications_1p.tsv) the sharded search is repeated (same proteome,
parameters and hyperscore FDR) in four variants: 'peaks_removed' (34 nucleic-acid ions removed from all spectra),
'chimeric' (Sage chimeric search), 'removed+chimeric' (both) and 'wide_window' (wide window, +-1 Th), and identification
rates of spectra with and without nucleobase are compared with the original search.
Usage: python rescue.py SAGE_EXE FASTA [PXD...]   -> rescue.tsv"""
import sys, os
SAGE, FASTA = sys.argv[1], sys.argv[2]
ONLY = [a for a in sys.argv[3:] if a.startswith('PXD')]
sys.argv = [sys.argv[0], SAGE, FASTA, '--one-step']
import zlib, pickle, shutil
import numpy as np
import nucleic_census as C
import identifications as I

PROTON = 1.007276
# neutral monoisotopic masses
BASES = {'A': 135.05450, 'G': 151.04941, 'C': 111.04326, 'U': 112.02728, 'T': 126.04293}
NUCLEOSIDES = {'A': 267.09675, 'G': 283.09167, 'C': 243.08552, 'U': 244.06954,
               'dA': 251.10184, 'dG': 267.09675, 'dC': 227.09061, 'dT': 242.09027}
H2O, HPO3 = 18.01056, 79.96633
IONS = [m + PROTON for m in BASES.values()]
for m in NUCLEOSIDES.values():
    IONS += [m + PROTON, m + HPO3 + PROTON, m + HPO3 - H2O + PROTON, m - H2O + PROTON]   # nucleoside, NMP, NMP-H2O, -H2O
IONS += [97.97690 + PROTON]                                                              # H3PO4
IONS = np.unique(np.round(IONS, 5))
TOL = 0.005
VARIANTS = {'peaks_removed': (True, None), 'chimeric': (False, {'chimera': True, 'report_psms': 2}),
            'removed+chimeric': (True, {'chimera': True, 'report_psms': 2}),
            'wide_window': (False, {'wide_window': True})}   # any peptide within the isolation window (+-1 Th)


def remove_peaks(spectra):
    out = []
    for s in spectra:
        k = np.ones(len(s['mz']), bool)
        for m in IONS:
            k &= np.abs(s['mz'] - m) > TOL
        t = dict(s); t['mz'] = s['mz'][k]; t['it'] = s['it'][k]
        out.append(t)
    return out


def rates(ident, m):
    return ident.mean(), (ident[m].mean() if m.sum() else np.nan), ident[~m].mean()


if __name__ == '__main__':
    L = [l.split('\t') for l in open('identifications_1p.tsv', encoding='utf-8').read().split('\n') if l]
    header = L[0]; A = [dict(zip(header, r)) for r in L[1:]]
    A = [a for a in A if a['group'] == 'contaminated' and float(a['id_rate_all']) >= 0.02 and int(a['n_nucleobase_spectra']) >= 20
         and (not ONLY or a['project'] in ONLY)]
    print('%d contaminated files; %d nucleic-acid ions to remove' % (len(A), len(IONS)), flush=True)
    OUT = 'rescue.tsv'
    done = set()
    if os.path.exists(OUT):
        done = {tuple(l.split('\t')[:3]) for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l}
    else:
        open(OUT, 'w', encoding='utf-8').write('project\tfile\tvariant\tn_ms2\tn_nucleobase_spectra\tid_rate_all\tid_rate_nucleobase\t'
                                               'id_rate_no_nucleobase\tpeaks_removed_nucleobase\tpeaks_removed_other\n')
    for a in A:
        p, fn, phospho = a['project'], a['file'], a['phospho_enrichment'] == '1'
        base = os.path.join(I.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        spectra = pickle.load(open(base + '.pkl', 'rb'))
        nb = np.array([max(C.rel(s, C.BASES['A']), C.rel(s, C.BASES['G'])) for s in spectra]); m = nb >= 0.10
        lab, mass = I.labelling(spectra)
        rows = []
        if (p, fn, 'original') not in done:
            r = I.sharded_search(base, len(spectra), phospho, mass)          # reuses the results already computed
            rows.append(('original', r, 0, 0))
        cleaned = None
        for v, (remove, extra) in VARIANTS.items():
            if (p, fn, v) in done:
                continue
            if remove and cleaned is None:
                cleaned = remove_peaks(spectra)
                removed = np.array([len(s['mz']) - len(t['mz']) for s, t in zip(spectra, cleaned)])
                I.write_mzml(cleaned, base + '_clean.mzML')
            b = base + ('_clean' if remove else '') + (('_w' if 'wide_window' in extra else '_q') if extra else '')
            if not remove and extra:
                shutil.copy(base + '.mzML', b + '.mzML')
            elif remove and extra:
                shutil.copy(base + '_clean.mzML', b + '.mzML')
            r = I.sharded_search(b, len(spectra), phospho, mass, extra=extra)
            rows.append((v, r, removed[m].mean() if remove else 0, removed[~m].mean() if remove else 0))
        with open(OUT, 'a', encoding='utf-8') as f:
            for v, r, bn, bs in rows:
                if r is None:
                    print('  failed', p, fn, v, flush=True); continue
                t, tn, ts = rates(r[0], m)
                f.write('%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\n' % (p, fn, v, len(spectra), m.sum(), t, tn, ts, bn, bs))
                print('%s %-34s %-16s ID all %5.1f %% | with nb %5.1f %% | without nb %5.1f %% | peaks removed nb %.1f / other %.1f'
                      % (p, fn[:34], v, 100 * t, 100 * tn, 100 * ts, bn, bs), flush=True)
    print('DONE')
