#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Do the three 4 MB chunks (25/50/75 %) represent the whole file? (manuscript: Methods and Results, partial reading)
Valid atlas files of 17-60 MB (contaminated and clean, random) are downloaded entirely and nb10 of the whole file is
compared with nb10 of the three chunks (same procedure as the atlas).
Usage (from this folder, with nucleoqc installed): python chunks_vs_whole.py [n]  -> chunks_vs_whole.tsv / _results.txt"""
import sys, os, random, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from nucleoqc import reading, metrics, pride

N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
OUT = 'chunks_vs_whole.tsv'
tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])

A = [a for a in tab('atlas.tsv') if a['valid'] == 'True' and a['file_MB'].isdigit() and 17 <= int(a['file_MB']) <= 60]
rng = random.Random(11)
dirty = [a for a in A if float(a['nb10']) > 0.01]; clean = [a for a in A if float(a['nb10']) <= 0.01]
chosen = rng.sample(dirty, min(len(dirty), N // 2)) + rng.sample(clean, N - min(len(dirty), N // 2))
done = {tuple(l.split('\t')[:2]) for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
if not done:
    open(OUT, 'w', encoding='utf-8').write('project\tfile\tMB\tn_chunks\tnb10_chunks\tnb2_chunks\tn_whole\tnb10_whole\tnb2_whole\tverdict_chunks\tverdict_whole\n')
lock = threading.Lock()
print('%d contaminated and %d clean eligible; measuring %d' % (len(dirty), len(clean), len(chosen)), flush=True)


def one(a):
    if (a['project'], a['file']) in done:
        return
    try:
        f = next(f for f in pride.files(a['project']) if f[0] == a['file'])
        n, size, fmt, url = f
        ch = metrics.measure(reading.parse(pride.read_remote(url, size), fmt))
        whole = metrics.measure(reading.parse([pride.byte_range(url, 0, size - 1)], fmt))
    except Exception as ex:
        print('  %s %s: error %s' % (a['project'], a['file'][:40], ex), flush=True); return
    g = lambda r, k: ('%.5f' % r[k]) if k in r else ''
    with lock:
        open(OUT, 'a', encoding='utf-8').write('\t'.join([a['project'], a['file'], '%.0f' % (size / 1e6), str(ch['n']), g(ch, 'nb10'), g(ch, 'nb2'),
                                                          str(whole['n']), g(whole, 'nb10'), g(whole, 'nb2'), ch['verdict'], whole['verdict']]) + '\n')
        print('  %s %-40s %3.0f MB | chunks n=%5d nb10 %6.2f %% | whole n=%6d nb10 %6.2f %% | %s / %s'
              % (a['project'], a['file'][:40], size / 1e6, ch['n'], 100 * ch.get('nb10', np.nan), whole['n'], 100 * whole.get('nb10', np.nan),
                 ch['verdict'], whole['verdict']), flush=True)


with ThreadPoolExecutor(4) as ex:
    list(ex.map(one, chosen))

R = [r for r in tab(OUT) if r['nb10_chunks'] and r['nb10_whole']]
t = np.array([float(r['nb10_chunks']) for r in R]); c = np.array([float(r['nb10_whole']) for r in R])
t2 = np.array([float(r['nb2_chunks']) for r in R]); c2 = np.array([float(r['nb2_whole']) for r in R])
fr = np.array([int(r['n_chunks']) / int(r['n_whole']) for r in R])
lines = ['Chunks (3 x 4 MB) versus whole file: %d files (%d with nb10 > 1 %% in the whole file)' % (len(R), int((c > 0.01).sum())),
         '  fraction of MS2 spectra read with the chunks: median %.0f %% (range %.0f-%.0f %%)' % (100 * np.median(fr), 100 * fr.min(), 100 * fr.max()),
         '  nb10: median absolute difference %.2f points, maximum %.2f; Spearman %.3f; Pearson (log10) %.3f'
         % (100 * np.median(np.abs(t - c)), 100 * np.max(np.abs(t - c)), stats.spearmanr(t, c)[0], np.corrcoef(np.log10(t + 1e-4), np.log10(c + 1e-4))[0, 1]),
         '  nb2:  median absolute difference %.2f points, maximum %.2f' % (100 * np.median(np.abs(t2 - c2)), 100 * np.max(np.abs(t2 - c2))),
         '  same classification (nb10 > 1 %%) in %d of %d; same nucleoqc verdict in %d of %d'
         % (int(((t > 0.01) == (c > 0.01)).sum()), len(R), sum(r['verdict_chunks'] == r['verdict_whole'] for r in R), len(R))]
print('\n'.join(lines))
open('chunks_vs_whole_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
