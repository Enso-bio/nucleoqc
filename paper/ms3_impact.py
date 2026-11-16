#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Impact of MS3 spectra on nb10 (reviewer check). atlas.py and nucleoqc <= 0.2.0 read every spectrum with ms level
>= 2, so SPS-MS3 reporter scans were counted as MS2. Re-reads the same three 4 MB chunks of every valid mzML file
of the atlas and measures nb10 with nucleoqc counting MS2 only and MS2 + MS3 (both with the corrected peak search).
Usage (from this folder, with nucleoqc installed): python ms3_impact.py  -> ms3_impact.tsv, ms3_impact_results.txt"""
import os, re, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from nucleoqc import metrics, reading, pride

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])
_LV3 = re.compile(r'(accession="MS:1000511"[^>]*value=")3"')
_LV = re.compile(r'accession="MS:1000511"[^>]*value="(\d+)"')


def measure(a):
    try:
        f = next(f for f in pride.files(a['project']) if f[0] == a['file'])
        chunks = [c.decode('latin-1') for c in pride.read_remote(f[3], f[1])]
    except Exception as ex:
        return a, None, str(ex)
    levels = [lv for t in chunks for lv in _LV.findall(t)]
    n2, n3 = levels.count('2'), levels.count('3')
    ms2_only = metrics.measure([s for t in chunks for s in reading.mzml_spectra(t)])
    with_ms3 = metrics.measure([s for t in chunks for s in reading.mzml_spectra(_LV3.sub(r'\g<1>2"', t))]) if n3 else ms2_only
    return a, (n2, n3, ms2_only.get('nb10', np.nan), with_ms3.get('nb10', np.nan), ms2_only.get('valid'), with_ms3.get('valid')), ''


if __name__ == '__main__':
    A = [a for a in tab('atlas.tsv') if a['valid'] == 'True' and a['format'] == 'mzml']
    print('%d valid mzML files' % len(A), flush=True)
    rows = []
    with ThreadPoolExecutor(4) as ex:
        for a, r, err in ex.map(measure, A):
            if r is None:
                print('  error %s %s: %s' % (a['project'], a['file'][:40], err), flush=True); continue
            rows.append((a['project'], a['file'], float(a['nb10'])) + r)
            print('  %s %-40s MS2 %5d MS3 %5d | atlas %6.2f %% | MS2 only %6.2f %% | MS2+MS3 %6.2f %%'
                  % (a['project'], a['file'][:40], r[0], r[1], 100 * float(a['nb10']), 100 * r[2], 100 * r[3]), flush=True)
    open('ms3_impact.tsv', 'w', encoding='utf-8').write(
        'project\tfile\tnb10_atlas\tn_ms2\tn_ms3\tnb10_ms2\tnb10_ms2_ms3\tvalid_ms2\tvalid_ms2_ms3\n' +
        ''.join('\t'.join(str(x) for x in f) + '\n' for f in rows))
    F = [f for f in rows if np.isfinite(f[5]) and np.isfinite(f[6])]
    with_ms3 = [f for f in F if f[4] > 0]
    d_rel = np.array([f[5] - f[2] for f in F])
    lin = ['Valid atlas mzML files re-read: %d (errors %d)' % (len(F), len(A) - len(rows)),
           'With MS3 spectra: %d files from %d projects' % (len(with_ms3), len({f[0] for f in with_ms3}))]
    for f in with_ms3:
        lin.append('   %s %-40s MS3 %4.0f %% of MS2+MS3 | nb10 MS2+MS3 %6.2f %% -> MS2 only %6.2f %% | atlas %6.2f %%'
                   % (f[0], f[1][:40], 100 * f[4] / (f[3] + f[4]), 100 * f[6], 100 * f[5], 100 * f[2]))
    changed = [f for f in F if (f[5] > 0.01) != (f[2] > 0.01)]
    lin.append('Files that cross the 1 %% threshold (MS2 only vs atlas): %d %s' % (len(changed), [(f[0], f[1][:30]) for f in changed]))
    lin.append('Difference MS2 only minus atlas (all): median %.3f points, maximum absolute %.2f points'
               % (100 * np.median(d_rel), 100 * np.max(np.abs(d_rel))))
    print('\n'.join(lin))
    open('ms3_impact_results.txt', 'w', encoding='utf-8').write('\n'.join(lin) + '\n')
