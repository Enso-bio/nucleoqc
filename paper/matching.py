#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Within-file comparison controlling for charge, retention time and precursor intensity (manuscript: robustness,
Figure 4, Figure S6).

Step 1 (network): re-reads the same chunks and extracts precursor intensity (MGF PEPMASS, mzML MS:1000042), aligned with
the stored spectra -> matching_intensity.npz
Step 2: analysis (increasingly strict matching and per-spectrum logistic regression with file fixed effects)
-> matching_results.txt, matching_figure.png
Usage: python matching.py"""
import os, re, zlib, pickle, collections
import numpy as np
from scipy import stats
import nucleic_census as C
import atlas

CACHE = 'matching_intensity.npz'
RT_WINDOW = 300.0          # +-5 min
LOGI_WINDOW = 0.3          # +-0.3 log10 (factor 2)
MZ_WINDOW = 10.0
MIN_CTRL = 3

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


# ---------------------------------------------------------------- step 1: precursor intensity
_RE_PI = re.compile(r'accession="MS:1000042"[^>]*value="([0-9.eE+-]+)"')


def mgf_intensity(txt):
    """Same filter and order as identifications.mgf_with_charge; returns (pmz, intensity)."""
    i, j = txt.find('BEGIN IONS'), txt.rfind('END IONS')
    out = []
    if i < 0 or j < i:
        return out
    for blk in txt[i:j].split('BEGIN IONS')[1:]:
        npk, pmz, pint = 0, 0.0, np.nan
        for l in blk.split('\n'):
            l = l.strip()
            if not l:
                continue
            if l[0].isdigit():
                p = l.split()
                if len(p) >= 2:
                    try:
                        float(p[0]); float(p[1]); npk += 1
                    except ValueError:
                        pass
            elif l.startswith('PEPMASS='):
                p = l[8:].split()
                try:
                    pmz = float(p[0])
                    if len(p) > 1:
                        pint = float(p[1])
                except (ValueError, IndexError):
                    pass
        if npk >= 10 and pmz > 0:
            out.append((pmz, pint))
    return out


def mzml_intensity(txt):
    out = []
    for blk in atlas._RE_SPEC.findall(txt):
        e = atlas.mzml_spectra(blk)
        if not e:
            continue
        m = _RE_PI.search(blk)
        out.append((e[0]['pmz'], float(m.group(1)) if m else np.nan))
    return out


def intensities(p, fn, spectra):
    fs = [x for x in atlas.eligible_files(p) if x[0] == fn]
    if not fs:
        return None, 'not found'
    n, size, fmt, url = fs[0]
    offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * atlas.CHUNK else [0]
    r = []
    for a in offs:
        t = C.byte_range(url, a, min(size - 1, a + atlas.CHUNK - 1)).decode('latin-1')
        r += (mgf_intensity if fmt == 'mgf' else mzml_intensity)(t)
    if len(r) != len(spectra):
        return None, 'misaligned (%d vs %d)' % (len(r), len(spectra))
    pm = np.array([x[0] for x in r]); pe = np.array([s['pmz'] for s in spectra])
    if np.median(np.abs(pm - pe) / pe * 1e6) > 30:          # recalibration moves m/z by <20 ppm
        return None, 'precursor m/z mismatch'
    v = np.array([x[1] for x in r], float)
    v[~(v > 0)] = np.nan
    return v, 'ok' if np.isfinite(v).mean() > 0.5 else 'no intensity in file'


def load():
    A = [a for a in tab('identifications_1p.tsv') if float(a['id_rate_all']) >= 0.02]
    E = collections.defaultdict(dict)
    for r in tab('ids_spectra_1p.tsv'):
        E[(r['project'], r['file'])][int(r['index'])] = (float(r['nb']), r['id'] == '1')
    cache = dict(np.load(CACHE, allow_pickle=True)) if os.path.exists(CACHE) else {}
    data = []
    for a in A:
        p, fn = a['project'], a['file']
        base = os.path.join(C.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        spectra = pickle.load(open(base + '.pkl', 'rb'))
        e = E[(p, fn)]
        if len(e) != len(spectra):
            print('  %s %s: table and spectra do not match, skipped' % (p, fn)); continue
        key = p + '|' + fn
        if key not in cache:
            v, reason = intensities(p, fn, spectra)
            cache[key] = np.array([reason, v if v is not None else np.full(len(spectra), np.nan)], dtype=object)
            np.savez(CACHE, **cache)
            print('  intensity %s %-40s %s' % (p, fn[:40], reason), flush=True)
        reason, v = cache[key]
        data.append({'p': p, 'fn': fn, 'group': a['group'], 'reason': reason,
                     'nb': np.array([e[i][0] for i in range(len(spectra))]), 'id': np.array([e[i][1] for i in range(len(spectra))]),
                     'pmz': np.array([s['pmz'] for s in spectra]), 'z': np.array([s['z'] for s in spectra]),
                     'rt': np.array([s['rt'] for s in spectra], float), 'logi': np.log10(np.asarray(v, float))})
    return data


# ---------------------------------------------------------------- step 2: analysis
def expected(d, use):
    """For each spectrum with nb >= 10 %, identification rate of its controls without nucleobase.
    use: subset of {'z','rt','int'} added to the base criterion (same file, +-10 Th)."""
    m = d['nb'] >= 0.10; ctrl = np.where(d['nb'] < 0.10)[0]       # controls as in the original analysis
    has_rt = np.ptp(d['rt']) > 60; has_i = np.isfinite(d['logi']).mean() > 0.5
    # a file that lacks a requested criterion (no RT or no precursor intensity) is excluded from that row,
    # never matched on the remaining criteria only
    if ('rt' in use and not has_rt) or ('int' in use and not has_i):
        return np.array([], float), np.array([], float)
    obs, exp_ = [], []
    for i in np.where(m)[0]:
        c = ctrl[np.abs(d['pmz'][ctrl] - d['pmz'][i]) <= MZ_WINDOW]
        if 'z' in use:
            c = c[d['z'][c] == d['z'][i]]
        if 'rt' in use and has_rt:
            c = c[np.abs(d['rt'][c] - d['rt'][i]) <= RT_WINDOW]
        if 'int' in use and has_i:
            if not np.isfinite(d['logi'][i]):
                continue
            c = c[np.abs(d['logi'][c] - d['logi'][i]) <= LOGI_WINDOW]
        if len(c) >= MIN_CTRL:
            obs.append(d['id'][i]); exp_.append(d['id'][c].mean())
    return np.array(obs, float), np.array(exp_, float)


def describe(D, out):
    out('A) How spectra with nucleobase (>=10 %) differ from the other spectra of the same file')
    out('   (median per file of the difference with nb - without nb; 18 contaminated files with >=20 spectra with nb)')
    dz, drt, di = [], [], []
    for d in D:
        m = d['nb'] >= 0.10; s = d['nb'] < 0.10
        z2 = lambda k: np.mean(d['z'][k] == 2) if k.sum() else np.nan
        dz.append(z2(m) - z2(s))
        if np.ptp(d['rt']) > 60:
            r = np.ptp(d['rt'])
            drt.append((np.median(d['rt'][m]) - np.median(d['rt'][s])) / r)
        ok = np.isfinite(d['logi'])
        if ok.mean() > 0.5:
            di.append(np.nanmedian(d['logi'][m & ok]) - np.nanmedian(d['logi'][s & ok]))
    W = lambda x: stats.wilcoxon(x).pvalue if len(x) >= 5 else np.nan
    out('   fraction of charge 2+: %+.3f (Wilcoxon p = %.2g, %d files)' % (np.median(dz), W(np.array(dz)), len(dz)))
    out('   retention time (fraction of the measured gradient): %+.3f (p = %.2g, %d files with RT)' % (np.median(drt), W(np.array(drt)), len(drt)))
    out('   log10 precursor intensity: %+.3f (p = %.2g, %d files with intensity)' % (np.median(di) if di else np.nan, W(np.array(di)), len(di)))


def matched(D, out):
    out('\nB) Observed vs expected rate under increasingly strict matching')
    out('   controls: spectra with nucleobase <10 %% of the same file (as in the original analysis); >=%d per spectrum' % MIN_CTRL)
    criteria = [('m/z +-10 Th (original)', set()), ('+ same charge', {'z'}), ('+ charge + RT +-5 min', {'z', 'rt'}),
                ('+ charge + intensity x/2', {'z', 'int'}), ('+ charge + RT + intensity', {'z', 'rt', 'int'})]
    rows = {}
    for name, u in criteria:
        oo, ee, nn, cov = [], [], 0, []
        for d in D:
            o, e = expected(d, u)
            cov.append(len(o) / max((d['nb'] >= 0.10).sum(), 1))
            if len(o) >= 20:
                oo.append(o.mean()); ee.append(e.mean()); nn += len(o)
        oo, ee = np.array(oo), np.array(ee)
        p = stats.wilcoxon(oo, ee).pvalue
        rows[name] = (oo, ee)
        out('   %-28s %2d files, %5d spectra (%3.0f %% matched) | observed %5.1f %% vs expected %5.1f %% | '
            'median ratio %.2f | lower in %d of %d | Wilcoxon p = %.2g'
            % (name, len(oo), nn, 100 * np.median(cov), 100 * np.median(oo), 100 * np.median(ee), np.median(oo / np.maximum(ee, 1e-9)),
               int((oo < ee).sum()), len(oo), p))
    return rows


def logistic(D, out):
    import statsmodels.api as sm
    out('\nC) Per-spectrum logistic regression with file fixed effects')
    X, y, g = [], [], []
    for k, d in enumerate(D):
        has_rt = np.ptp(d['rt']) > 60; has_i = np.isfinite(d['logi']).mean() > 0.5
        rt = (d['rt'] - d['rt'].min()) / np.ptp(d['rt']) if has_rt else np.full(len(d['rt']), 0.5)
        li = d['logi'] - np.nanmedian(d['logi']) if has_i else np.zeros(len(d['rt']))
        li = np.where(np.isfinite(li), li, 0.0); missing_i = ~np.isfinite(d['logi']) if has_i else np.ones(len(li), bool)
        mz = (d['pmz'] - 800) / 400
        n10 = (d['nb'] >= 0.10) & (d['nb'] < 0.30); n30 = d['nb'] >= 0.30
        X.append(np.column_stack([n10, n30, d['z'] == 3, d['z'] >= 4, d['z'] == 0, mz, mz ** 2, rt, rt ** 2, li, li ** 2, missing_i]).astype(float))
        y.append(d['id'].astype(float)); g.append(np.full(len(d['id']), k))
    X, y, g = np.vstack(X), np.concatenate(y), np.concatenate(g)
    names = ['nb 10-30 %', 'nb >=30 %', 'charge 3', 'charge >=4', 'charge unknown', 'm/z', 'm/z^2', 'RT', 'RT^2', 'log int', 'log int^2', 'no int']
    # files with no identification or with all identified add nothing; fixed effect = one column per file (no constant)
    F = (g[:, None] == np.arange(g.max() + 1)[None, :]).astype(float)
    keep = X.std(0) > 0
    for label, cols in (('nucleobase only', [0, 1]), ('adjusted for charge, m/z, RT and intensity', list(np.where(keep)[0]))):
        r = sm.Logit(y, np.column_stack([X[:, cols], F])).fit(disp=0, maxiter=200)
        ci = np.asarray(r.conf_int()); pr = np.asarray(r.params); pv = np.asarray(r.pvalues)
        out('   %s:' % label)
        for j, c in enumerate(cols[:2]):
            out('      %-12s OR %.2f (95%% CI %.2f-%.2f), p = %.2g' % (names[c], np.exp(pr[j]), np.exp(ci[j, 0]), np.exp(ci[j, 1]), pv[j]))
        if len(cols) > 2:
            out('      covariates: ' + ', '.join('%s OR %.2f' % (names[c], np.exp(pr[j])) for j, c in enumerate(cols) if j >= 2))
    out('   (%d spectra, %d files)' % (len(y), len(D)))


if __name__ == '__main__':
    D = load()
    lines = []
    out = lambda s: (print(s, flush=True), lines.append(s))
    out('Files: %d (%d contaminated); precursor intensity: %s' % (len(D), sum(d['group'] == 'contaminated' for d in D),
        ', '.join('%s %d' % kv for kv in collections.Counter(d['reason'] for d in D).items())))
    out('Files with retention time: %d' % sum(np.ptp(d['rt']) > 60 for d in D))
    Dc = [d for d in D if d['group'] == 'contaminated' and (d['nb'] >= 0.10).sum() >= 20]
    describe(Dc, out)
    rows = matched(Dc, out)
    logistic(D, out)
    open('matching_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for k, (name, (oo, ee)) in enumerate(rows.items()):
        r = oo / np.maximum(ee, 1e-9)
        ax.scatter(np.full(len(r), k) + np.random.default_rng(k).uniform(-.12, .12, len(r)), r, s=14, c='#7f7f7f')
        ax.plot([k - .25, k + .25], [np.median(r)] * 2, c='#2c6fbb', lw=2)
    ax.axhline(1, c='k', lw=.6, ls='--'); ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(['m/z ±10 Th\n(main)', '+ same\ncharge', '+ charge\n+ RT ±5 min', '+ charge\n+ intensity ×/2', '+ charge + RT\n+ intensity'], fontsize=7)
    ax.set_ylabel('identified: observed / expected'); ax.set_title('Effect under increasingly strict matching (per file)', fontsize=10)
    plt.tight_layout(); plt.savefig('matching_figure.png', dpi=150)
    print('figure: matching_figure.png')
