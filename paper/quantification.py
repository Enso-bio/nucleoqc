#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Effect on quantification (TMT/TMTpro files; manuscript: identification section).

For each spectrum the reporter ions are measured (+-2.5 mDa after recalibration): channels detected and total reporter
intensity. Spectra with nucleobase (>=10 %) are compared with controls from the same file (nucleobase < 10 %, same
charge, precursor within +-10 Th), in all spectra and in identified spectra only. Reads the stored spectra from the
search work folder (nucleic_census.WORK).
-> quantification_results.txt"""
import os, zlib, pickle, collections
import numpy as np
from scipy import stats
import nucleic_census as NC

TMTPRO = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825, 130.141145,
          131.138180, 131.144500, 132.141535, 132.147855, 133.144890, 133.151210, 134.148245]
REP = {'TMTpro': TMTPRO, 'TMT': TMTPRO[:11]}
TOL = 0.0025

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def reporters(s, masses):
    v = []
    for m in masses:
        a, b = np.searchsorted(s['mz'], [m - TOL, m + TOL])
        v.append(s['it'][a:b].max() if b > a else 0.0)
    return np.array(v)


if __name__ == '__main__':
    E = collections.defaultdict(dict)
    for r in tab('ids_spectra_1p.tsv'):
        E[(r['project'], r['file'])][int(r['index'])] = (float(r['nb']), r['id'] == '1')
    lines = []
    out = lambda s: (print(s, flush=True), lines.append(s))
    out('Reporter ions in spectra with nucleobase versus matched controls (same file, charge, +-10 Th)')
    res = []
    for a in tab('identifications_1p.tsv'):
        if a['labelling'] not in REP or a['group'] != 'contaminated':
            continue
        p, fn = a['project'], a['file']
        spectra = pickle.load(open(os.path.join(NC.WORK, '%s_%d.pkl' % (p, zlib.crc32(fn.encode()) % 100000)), 'rb'))
        nb = np.array([E[(p, fn)][i][0] for i in range(len(spectra))]); idn = np.array([E[(p, fn)][i][1] for i in range(len(spectra))])
        R = np.array([reporters(s, REP[a['labelling']]) for s in spectra])
        nchan = (R > 0).sum(1); tot = np.log10(R.sum(1) + 1); nch = R.shape[1]
        pmz = np.array([s['pmz'] for s in spectra]); z = np.array([s['z'] for s in spectra])
        m = nb >= 0.10; ctrl = np.where(~m)[0]
        for label, sel in (('all', np.ones(len(spectra), bool)), ('identified', idn)):
            dc, dt, dfull, n = [], [], [], 0
            for i in np.where(m & sel)[0]:
                c = ctrl[(np.abs(pmz[ctrl] - pmz[i]) <= 10) & (z[ctrl] == z[i]) & sel[ctrl]]
                if len(c) >= 3:
                    dc.append(nchan[i] - nchan[c].mean()); dt.append(tot[i] - tot[c].mean())
                    dfull.append(float(nchan[i] == nch) - np.mean(nchan[c] == nch)); n += 1
            if n >= 10:
                res.append((p, fn, label, np.mean(dc), np.mean(dt), np.mean(dfull), n, nch))
                out('  %s %-36s %-13s %s, %4d matched | channels detected %+.2f of %d | log10 reporter signal %+.2f | '
                    'spectra with all channels %+.1f points'
                    % (p, fn[:36], label, a['labelling'], n, np.mean(dc), nch, np.mean(dt), 100 * np.mean(dfull)))
    for label in ('all', 'identified'):
        r = [x for x in res if x[2] == label]
        if r:
            out('  %-13s %d files: channels %+.2f, log10 signal %+.2f (x%.2f), with all channels %+.1f points (medians per file)'
                % (label, len(r), np.median([x[3] for x in r]), np.median([x[4] for x in r]), 10 ** np.median([x[4] for x in r]),
                   100 * np.median([x[5] for x in r])))
    open('quantification_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
