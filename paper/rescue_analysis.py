#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analysis of rescue.tsv: what each search variant recovers and, for the wide window, whether the identified peptide
is the selected precursor or a co-isolated neighbour; Figure S5. Reads the Sage results of the original and wide-window
searches from the search work folder (nucleic_census.WORK)."""
import os, re, zlib, pickle, glob
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nucleic_census as C

PROTON = 1.007276

L = [l.split('\t') for l in open('rescue.tsv', encoding='utf-8').read().split('\n') if l]
header = L[0]; R = [dict(zip(header, r)) for r in L[1:]]
files = sorted({(r['project'], r['file']) for r in R})
tab = {(r['project'], r['file'], r['variant']): r for r in R}
VAR = ['peaks_removed', 'chimeric', 'removed+chimeric', 'wide_window']
files = [a for a in files if all(a + (v,) in tab for v in ['original'] + VAR)]
print('%d contaminated files with all variants (%d projects)\n' % (len(files), len({a[0] for a in files})))

print('Change in the identification rate relative to the original search (percentage points; median per file)')
res = {}
for v in VAR:
    dn = np.array([float(tab[a + (v,)]['id_rate_nucleobase']) - float(tab[a + ('original',)]['id_rate_nucleobase']) for a in files])
    ds = np.array([float(tab[a + (v,)]['id_rate_no_nucleobase']) - float(tab[a + ('original',)]['id_rate_no_nucleobase']) for a in files])
    # total identifications of the file
    dt = np.array([(float(tab[a + (v,)]['id_rate_all']) - float(tab[a + ('original',)]['id_rate_all'])) * int(tab[a + (v,)]['n_ms2']) for a in files])
    nt = np.array([float(tab[a + ('original',)]['id_rate_all']) * int(tab[a + (v,)]['n_ms2']) for a in files])
    p = stats.wilcoxon(dn).pvalue if np.any(dn != 0) else 1.0
    res[v] = (dn, ds)
    print('  %-16s with nucleobase %+5.1f (better in %2d of %d; Wilcoxon p = %.2g) | without nucleobase %+5.1f | total identifications %+5.1f %%'
          % (v, 100 * np.median(dn), int((dn > 0).sum()), len(dn), p, 100 * np.median(ds), 100 * np.median(dt / np.maximum(nt, 1))))

print('\nWide window, per file (ID with nucleobase original -> wide window; without nucleobase original -> wide window):')
for a in sorted(files, key=lambda a: -int(tab[a + ('original',)]['n_nucleobase_spectra']) / int(tab[a + ('original',)]['n_ms2'])):
    o, w = tab[a + ('original',)], tab[a + ('wide_window',)]
    print('  %s %-34s nb %5.1f %% of MS2 | with nb %5.1f -> %5.1f | without nb %5.1f -> %5.1f'
          % (a[0], a[1][:34], 100 * int(o['n_nucleobase_spectra']) / int(o['n_ms2']), 100 * float(o['id_rate_nucleobase']),
             100 * float(w['id_rate_nucleobase']), 100 * float(o['id_rate_no_nucleobase']), 100 * float(w['id_rate_no_nucleobase'])))


# --- in the wide-window search: is the identified peptide the selected precursor?
def best_hits(base, n):
    best = np.full(n, -np.inf); lab = np.zeros(n, bool); mz = np.zeros(n)
    for d in sorted(glob.glob(base + '_t*')):
        if not os.path.isdir(d) or not re.search(r'_t\d+$', d):
            continue
        lin = open(os.path.join(d, 'results.sage.tsv'), encoding='utf-8').read().split('\n'); J = {c: i for i, c in enumerate(lin[0].split('\t'))}
        for l in lin[1:]:
            if not l:
                continue
            c = l.split('\t'); i = int(re.search(r'(\d+)$', c[J['scannr']]).group(1)); h = float(c[J['hyperscore']])
            if h > best[i]:
                best[i] = h; lab[i] = c[J['label']] == '1'; mz[i] = float(c[J['calcmass']]) / int(c[J['charge']]) + PROTON
    ok = np.isfinite(best); o = np.argsort(-best[ok]); idx = np.where(ok)[0][o]
    t = np.cumsum(lab[idx]); d = np.cumsum(~lab[idx]); q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    ident = np.zeros(n, bool); ident[idx] = lab[idx] & (q <= 0.01)
    return ident, mz


print('\nWide window: of the identified spectra, fraction whose peptide is NOT the selected precursor')
print('(|peptide m/z - selected m/z| > 0.02 Th after allowing for isotope peaks +1/z, +2/z)')
rows = []
for p, fn in files:
    base = os.path.join(C.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
    spectra = pickle.load(open(base + '.pkl', 'rb'))
    nb = np.array([max(C.rel(s, C.BASES['A']), C.rel(s, C.BASES['G'])) for s in spectra]) >= 0.10
    pmz = np.array([s['pmz'] for s in spectra]); z = np.array([s['z'] for s in spectra])
    io, _ = best_hits(base, len(spectra)); iw, mzw = best_hits(base + '_w', len(spectra))
    dif = np.abs(mzw - pmz)
    zz = np.where(z > 0, z, 2)
    same = np.min([np.abs(mzw - pmz - k * 1.00336 / zz) for k in (-2, -1, 0)], axis=0) <= 0.02
    for g, sel in (('with nb', nb), ('without nb', ~nb)):
        s = sel & iw
        if s.sum() >= 10:
            rows.append((p, g, (~same[s]).mean(), s.sum(), (sel & iw & ~io).sum(), (sel & io & ~iw).sum()))
    print('  %s %-30s with nb: %4d ident., %5.1f %% neighbour | without nb: %4d ident., %5.1f %% neighbour | nb gained %d, lost %d'
          % (p, fn[:30], (nb & iw).sum(), 100 * (~same[nb & iw]).mean() if (nb & iw).sum() else np.nan,
             (~nb & iw).sum(), 100 * (~same[~nb & iw]).mean() if (~nb & iw).sum() else np.nan,
             (nb & iw & ~io).sum(), (nb & io & ~iw).sum()))
cn = np.array([f[2] for f in rows if f[1] == 'with nb']); sn = np.array([f[2] for f in rows if f[1] == 'without nb'])
if len(cn) and len(sn):
    print('  median per file: co-isolated neighbour in %.1f %% (with nb) vs %.1f %% (without nb)' % (100 * np.median(cn), 100 * np.median(sn)))

fig, ax = plt.subplots(figsize=(6, 4))
for k, v in enumerate(VAR):
    dn, ds = res[v]
    ax.scatter(np.full(len(ds), k - .12) + np.random.default_rng(k).uniform(-.05, .05, len(ds)), 100 * ds, s=12, c='#9bbbd9', label='without nucleobase' if k == 0 else None)
    ax.scatter(np.full(len(dn), k + .12) + np.random.default_rng(k + 9).uniform(-.05, .05, len(dn)), 100 * dn, s=12, c='#c0504d', label='with nucleobase' if k == 0 else None)
ax.axhline(0, c='k', lw=.6); ax.set_xticks(range(len(VAR))); ax.set_xticklabels(['nucleic-acid\npeaks removed', 'chimeric\nsearch', 'removed +\nchimeric', 'wide\nwindow'], fontsize=8)
ax.set_ylabel('change in % identified vs original search'); ax.legend(fontsize=8); ax.set_title('Software rescue (contaminated files)', fontsize=10)
plt.tight_layout(); plt.savefig('rescue_figure.png', dpi=150); print('\nfigure: rescue_figure.png')
