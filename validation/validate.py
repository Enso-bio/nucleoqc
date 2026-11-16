# -*- coding: utf-8 -*-
"""Validation of nucleoqc with real data.
1) Reproduces atlas nb10: re-measures with nucleoqc PRIDE files already measured in atlas.tsv.
2) Estimated vs measured loss (Sage, identifications_1p.tsv) on the same spectra (calibration files; not an
   independent test; the independent validation is paper/validation_independent.py and paper/validation_extended.py).
The second part reads the stored spectra from the search work folder (NUCLEOQC_SEARCH_WORK, or paper/search_work).
Usage (from the paper/ folder): python ../validation/validate.py [n_atlas_files]
"""
import sys, os, zlib, pickle, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import numpy as np
from scipy import stats
from nucleoqc import metrics, reading, pride

N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])

# ---------------------------------------------------------------- 1) atlas
A = [a for a in tab('atlas.tsv') if a['valid'] == 'True']
rng = random.Random(7)
dirty = [a for a in A if float(a['nb10']) > 0.01]; clean = [a for a in A if float(a['nb10']) <= 0.01]
chosen = rng.sample(dirty, N // 2) + rng.sample(clean, N - N // 2)
print('1) nb10: atlas vs nucleoqc (same PRIDE files, same chunks)')
xa, xb, va, vb = [], [], [], []
for a in chosen:
    try:
        f = next(f for f in pride.files(a['project']) if f[0] == a['file'])
        r = metrics.measure(reading.parse(pride.read_remote(f[3], f[1]), f[2]))
    except Exception as ex:
        print('   %s %-40s error: %s' % (a['project'], a['file'][:40], ex)); continue
    xa.append(float(a['nb10'])); xb.append(r.get('nb10', np.nan))
    va.append(metrics.verdict(float(a['nb10']))); vb.append(r['verdict'])
    print('   %s %-40s atlas %6.2f %% | nucleoqc %6.2f %%  %-8s valid=%s' % (a['project'], a['file'][:40], 100 * xa[-1],
          100 * xb[-1], r['verdict'], r['valid']), flush=True)
xa, xb = np.array(xa), np.array(xb); ok = np.isfinite(xb)
print('   %d files; median absolute difference %.3f points, maximum %.3f; Spearman %.3f; same verdict in %d of %d'
      % (ok.sum(), 100 * np.median(np.abs(xa - xb)[ok]), 100 * np.max(np.abs(xa - xb)[ok]),
         stats.spearmanr(xa[ok], xb[ok])[0], sum(x == y for x, y in zip(va, vb)), len(va)))

# ---------------------------------------------------------------- 2) loss (calibration files)
WORK = os.environ.get('NUCLEOQC_SEARCH_WORK') or os.path.join(os.getcwd(), 'search_work')
print('\n2) identification loss: estimated by nucleoqc vs measured with Sage (calibration files)')
print('   (measured = n_nb x (expected - observed) / (identified + lost))')
est, obs = [], []
for a in tab('identifications_1p.tsv'):
    if a['group'] != 'contaminated' or float(a['id_rate_all']) < 0.02 or not a['id_rate_expected_nucleobase']:
        continue
    base = os.path.join(WORK, '%s_%d.pkl' % (a['project'], zlib.crc32(a['file'].encode()) % 100000))
    if not os.path.exists(base):
        continue
    spectra = [{'mz': np.asarray(s['mz'], float), 'it': np.asarray(s['it'], float), 'pmz': float(s['pmz'])} for s in pickle.load(open(base, 'rb'))]
    r = metrics.measure(spectra)
    n_nb, n = int(a['n_nucleobase_spectra']), int(a['n_ms2'])
    lost = n_nb * (float(a['id_rate_expected_nucleobase']) - float(a['id_rate_nucleobase']))
    m = lost / max(float(a['id_rate_all']) * n + lost, 1)
    est.append(r['estimated_loss']); obs.append(m)
    print('   %s %-34s nb10 %5.1f %% | estimated %5.1f %% | measured %5.1f %%' % (a['project'], a['file'][:34], 100 * r['nb10'],
          100 * est[-1], 100 * m))
est, obs = np.array(est), np.array(obs)
print('   %d files; Spearman %.2f (p = %.2g); Pearson %.2f; median absolute error %.1f points'
      % (len(est), *stats.spearmanr(est, obs), stats.pearsonr(est, obs)[0], 100 * np.median(np.abs(est - obs))))
