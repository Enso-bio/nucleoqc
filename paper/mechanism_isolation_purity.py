#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mechanism: do spectra with nucleobase have a less pure precursor isolation? (manuscript: software rescue section)

In mzML files that keep MS1 spectra, the precursor isolation purity (fraction of MS1 intensity in the isolation window
belonging to the isotope envelope of the selected precursor, peaks at k*1.00336/z, k = -1..3, +-10 ppm) is computed for
each MS2 in the preceding MS1, and compared between spectra with and without nucleobase.

Usage: python mechanism_isolation_purity.py [PXD|file ...]  -> mechanism_isolation_purity_results[_extra].txt"""
import re, os, zlib, base64, collections
import numpy as np
from scipy import stats
import nucleic_census as C
import atlas

FILES = [('PXD055857', '20200115_HD_IN_Phospho2R.mzML'), ('PXD072006', 'CE3_2024Jan30_Jelinek_P23193_PO3_Fr.mzML'),
         ('PXD052300', '1227MT16plexfractions.mzML'), ('PXD052300', '1227MT16plexPhosphofractions.mzML'),
         ('PXD036323', 'Mashid_phospho_deglyco_05034_Mahshid.mzML'), ('PXD011616', 'Run3_Kentrophoros_h_pool3_2500ng.mzML'),
         ('PXD055077', 'SUMOylation_PhoX_TBDSPP_FT_FAIMS_140823.mzML'), ('PXD067736', '20-89PD22_3h-PTM.mzML'),
         ('PXD070403', 'Imp70-1.mzML')]
ISO = 1.003355
_RE_T = re.compile(r'accession="MS:1000827"[^>]*value="([0-9.eE+-]+)"')
_RE_LO = re.compile(r'accession="MS:1000828"[^>]*value="([0-9.eE+-]+)"')
_RE_HI = re.compile(r'accession="MS:1000829"[^>]*value="([0-9.eE+-]+)"')
_RE_Z = re.compile(r'accession="MS:1000041"[^>]*value="(\d+)"')

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def arrays(blk):
    arr = {}
    for b in atlas._RE_BDA.findall(blk):
        t = atlas._RE_BIN.search(b)
        if not t or not t.group(1):
            continue
        raw = base64.b64decode(t.group(1))
        if 'MS:1000574' in b:
            raw = zlib.decompress(raw)
        a = np.frombuffer(raw, dtype='<f8' if 'MS:1000523' in b else '<f4').astype(np.float64)
        if 'MS:1000514' in b:
            arr['mz'] = a
        elif 'MS:1000515' in b:
            arr['it'] = a
    if 'mz' not in arr or 'it' not in arr or len(arr['mz']) != len(arr['it']):
        return None
    o = np.argsort(arr['mz'])
    return arr['mz'][o], arr['it'][o]


def purity(ms1, pmz, z, lo, hi):
    mz, it = ms1
    a, b = np.searchsorted(mz, [pmz - lo, pmz + hi])
    tot = it[a:b].sum()
    if tot <= 0:
        return np.nan
    zz = z if z > 0 else 2
    env = 0.0
    for k in range(-1, 4):
        m = pmz + k * ISO / zz
        j, l = np.searchsorted(mz, [m - m * 10e-6, m + m * 10e-6])
        if l > j and a <= j and l <= b:
            env += it[j:l].max()
    return env / tot


def process(txt):
    """Walks through the spectra of the chunk in order; for each MS2 accepted by the atlas reader, returns its purity."""
    out = []; ms1 = None
    for blk in atlas._RE_SPEC.findall(txt):
        lv = atlas._RE_LV.search(blk)
        if lv and lv.group(1) == '1':
            ms1 = arrays(blk); continue
        e = atlas.mzml_spectra(blk)
        if not e:
            continue
        s = e[0]
        t = _RE_T.search(blk); lo = _RE_LO.search(blk); hi = _RE_HI.search(blk); z = _RE_Z.search(blk)
        lo = float(lo.group(1)) if lo else 0.7; hi = float(hi.group(1)) if hi else 0.7
        z = int(z.group(1)) if z else 0
        pif = purity(ms1, s['pmz'], z, lo if lo > 0 else 0.7, hi if hi > 0 else 0.7) if ms1 is not None else np.nan
        nb = max(C.rel(s, C.BASES['A'], 0.01), C.rel(s, C.BASES['G'], 0.01))
        out.append((pif, nb, lo + hi))
    return out


if __name__ == '__main__':
    import sys
    SUF = ''
    if len(sys.argv) > 1:                      # additional files: PXD|name
        FILES = [tuple(a.split('|', 1)) for a in sys.argv[1:]]; SUF = '_extra'
    IDS = collections.defaultdict(dict)
    for r in tab('ids_spectra_1p.tsv'):
        IDS[(r['project'], r['file'])][int(r['index'])] = (float(r['nb']), r['id'] == '1')
    lines = []
    out = lambda s: (print(s, flush=True), lines.append(s))
    out('Isolation purity (PIF) in the preceding MS1; median per file')
    diff, with_id = [], []
    for p, fn in FILES:
        fs = [x for x in atlas.eligible_files(p) if x[0] == fn]
        if not fs:
            out('  %s %s: not found' % (p, fn)); continue
        n, size, fmt, url = fs[0]
        offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * atlas.CHUNK else [0]
        R = []
        for a in offs:
            R += process(C.byte_range(url, a, min(size - 1, a + atlas.CHUNK - 1)).decode('latin-1'))
        pif = np.array([r[0] for r in R]); nb = np.array([r[1] for r in R]); win = np.median([r[2] for r in R]) if R else np.nan
        ok = np.isfinite(pif); m = nb >= 0.10
        if ok.sum() < 100 or (m & ok).sum() < 10:
            out('  %s %-40s %d MS2, %d with a preceding MS1, %d with nb: insufficient' % (p, fn[:40], len(R), ok.sum(), (m & ok).sum())); continue
        a_, b_ = np.median(pif[m & ok]), np.median(pif[~m & ok])
        diff.append(a_ - b_)
        pmw = stats.mannwhitneyu(pif[m & ok], pif[~m & ok]).pvalue
        line = '  %s %-40s window %.1f Th | PIF with nb %.2f (n=%d) vs without nb %.2f (n=%d) | Mann-Whitney p = %.2g' % (
            p, fn[:40], win, a_, (m & ok).sum(), b_, (~m & ok).sum(), pmw)
        # with identifications (files of the identification analysis): aligned by order with the stored spectra
        table = IDS.get((p, fn))
        if table and len(table) == len(R):
            nb_t = np.array([table[i][0] for i in range(len(R))]); idn = np.array([table[i][1] for i in range(len(R))])
            if np.corrcoef(nb_t, nb)[0, 1] > 0.95:
                for lim in (0.8, 0.9):
                    sel = ok & (pif >= lim)
                    mm = (nb_t >= 0.10) & sel; ss = (nb_t < 0.10) & sel
                    if mm.sum() >= 20:
                        with_id.append((p, fn, lim, idn[mm].mean(), idn[ss].mean(), mm.sum()))
                        line += ' | PIF>=%.1f: ID with nb %.1f %% vs without nb %.1f %%' % (lim, 100 * idn[mm].mean(), 100 * idn[ss].mean())
            else:
                line += ' | (alignment with the identification table doubtful; skipped)'
        out(line)
    diff = np.array(diff)
    if len(diff) >= 5:
        out('\n  %d files: median PIF difference (with nb - without nb) %+.3f; lower with nb in %d of %d; Wilcoxon p = %.2g'
            % (len(diff), np.median(diff), int((diff < 0).sum()), len(diff), stats.wilcoxon(diff).pvalue))
    for lim in (0.8, 0.9):
        c = [x for x in with_id if x[2] == lim]
        if c:
            out('  Only spectra with PIF >= %.1f (%d files with identifications): ID with nb %.1f %% vs without nb %.1f %% (medians); lower in %d of %d'
                % (lim, len(c), 100 * np.median([x[3] for x in c]), 100 * np.median([x[4] for x in c]), sum(x[3] < x[4] for x in c), len(c)))
    open('mechanism_isolation_purity_results%s.txt' % SUF, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
