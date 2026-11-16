#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Shared helpers (PRIDE API access, HTTP range reads, MGF parsing, recalibration, nucleobase ion intensity) and the
first exploratory census of nucleic-acid ions in public proteomics runs (PRIDE), by keyword category.
The census itself was superseded by the random atlas (atlas.py) and is not used for the manuscript estimates.

Search steps cache their spectra and Sage results in a work folder: the NUCLEOQC_SEARCH_WORK environment variable, or
search_work/ next to these scripts.

Usage: python nucleic_census.py [files_per_category]      -> nucleic_census.tsv (resumable)"""
import sys, os, re, json, time, urllib.request, urllib.parse, collections
from concurrent.futures import ThreadPoolExecutor
import numpy as np

API = 'https://www.ebi.ac.uk/pride/ws/archive/v3/'
WORK = os.environ.get('NUCLEOQC_SEARCH_WORK') or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'search_work')
OUTPUT = 'nucleic_census.tsv'
PER_CATEGORY = int(sys.argv[1]) if __name__ == '__main__' and len(sys.argv) > 1 else 7
CHUNK = 6_000_000
CATEGORIES = [  # (category, search keyword, preferred pattern in the file name)
    ('phospho', 'phosphoproteome', r'phos|tio2|imac|_p_|pS|pY'),
    ('phospho', 'phosphopeptide enrichment', r'phos|tio2|imac'),
    ('acetyl', 'acetylome', r'ac|kac'),
    ('ubiquitin', 'ubiquitylome', r'ub|gg|gly'),
    ('glyco', 'glycoproteomics', r'gly|hilic'),
    ('whole_proteome', 'label-free quantitative proteome', r''),
    ('whole_proteome', 'whole cell proteome', r''),
    ('plasma_secretome', 'plasma proteome', r''),
    ('interactome_IP', 'affinity purification interactome', r''),
    ('phospho', 'phosphoproteomics', r'phos|tio2|imac'),
    ('phospho', 'TiO2', r'phos|tio2|imac'),
    ('phospho', 'Fe-IMAC', r'phos|tio2|imac'),
    ('phospho', 'phosphorylation sites', r'phos|tio2|imac'),
    ('phospho', 'kinase signaling phosphoproteome', r'phos|tio2|imac'),
]
ONLY = sys.argv[2] if __name__ == '__main__' and len(sys.argv) > 2 else ''
if ONLY:
    CATEGORIES = [c for c in CATEGORIES if c[0] == ONLY]
INSTR = r'Q Exactive|Exploris|Fusion|Lumos|Eclipse'
BASES = {'A': 136.06177, 'G': 152.05669, 'C': 112.05054}
RIBOSE, Y1R, TOL = 97.02841, 175.11895, 0.005
# immonium and y1 ions present in almost every peptide spectrum, used for the global mass recalibration
REFS = np.array([70.06513, 72.08078, 86.09643, 110.07127, 120.08078, 129.10224, 136.07569, 147.11280, 175.11895])


def get_json(u):
    for _ in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Accept': 'application/json', 'User-Agent': 'python-urllib'}), timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2)
    return []


def byte_range(url, a, b):
    for _ in range(4):
        try:
            rq = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (a, b), 'User-Agent': 'python-urllib'})
            with urllib.request.urlopen(rq, timeout=180) as r:
                return r.read()
        except Exception:
            time.sleep(3)
    return b''


def mgf_spectra(txt):
    """MS2 spectra of an MGF chunk; 'low' = peaks below m/z 380 (used for recalibration); 'esc' is an empty array kept
    for compatibility with the cached spectra."""
    i, j = txt.find('BEGIN IONS'), txt.rfind('END IONS')
    if i < 0 or j < i:
        return []
    out = []
    for blk in txt[i:j].split('BEGIN IONS')[1:]:
        mz, it, pmz, z = [], [], 0.0, 0
        for l in blk.split('\n'):
            l = l.strip()
            if not l:
                continue
            if l[0].isdigit():
                p = l.split()
                if len(p) >= 2:
                    try:
                        mz.append(float(p[0])); it.append(float(p[1]))
                    except ValueError:
                        pass
            elif l.startswith('PEPMASS='):
                try:
                    pmz = float(l[8:].split()[0])
                except (ValueError, IndexError):
                    pass
        if len(mz) >= 10 and pmz > 0:
            m = np.array(mz); o = np.argsort(m); m = m[o]; t = np.array(it)[o]
            out.append({'mz': m, 'it': t, 'low': m[m < 380], 'esc': m[:0], 'pmz': pmz, 'z': z})
    return out


def recalibrate(spectra):
    """Global ppm offset estimated from ubiquitous immonium and y1 ions; corrects the spectra in place."""
    err = []
    for s in spectra:
        m = s['low']
        if len(m):
            j = np.clip(np.searchsorted(m, REFS), 1, len(m) - 1) if len(m) > 1 else np.zeros(len(REFS), int)
            for jj in (j - 1, j):
                e = (m[jj] - REFS) / REFS * 1e6
                err.extend(e[np.abs(e) < 20])
    if len(err) < 200:
        return 0.0
    ppm = float(np.median(err))
    f = 1 - ppm * 1e-6
    for s in spectra:
        s['mz'] = s['mz'] * f; s['low'] = s['low'] * f; s['esc'] = s['esc'] * f; s['pmz'] *= f
    return ppm


def rel(s, m, tol=TOL):
    """Intensity at m (+-tol, the two nearest peaks) relative to the base peak."""
    j = np.searchsorted(s['mz'], m); b = 0.0
    for k in (j - 1, j):
        if 0 <= k < len(s['mz']) and abs(s['mz'][k] - m) <= tol:
            b = max(b, s['it'][k])
    return b / s['it'].max()


def measure(url, size):
    offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * CHUNK else [0]
    with ThreadPoolExecutor(3) as ex:
        chunks = list(ex.map(lambda a: byte_range(url, a, min(size - 1, a + CHUNK - 1)), offs))
    spectra = []
    for t in chunks:
        spectra += mgf_spectra(t.decode('latin-1'))
    if len(spectra) < 300:
        return {'n_ms2': len(spectra), 'valid': False, 'reason': 'too few spectra'}
    ppm = recalibrate(spectra)
    n = len(spectra)
    coverage = np.mean([s['mz'][0] < 130 for s in spectra])
    near = [s for s in spectra if rel(s, Y1R, 0.2) > 0]
    hires = np.mean([rel(s, Y1R) > 0 for s in near]) if len(near) > 50 else 0.0
    r = {'n_ms2': n, 'ppm': ppm, 'low_mass_coverage': coverage, 'high_resolution': hires}
    A = np.array([rel(s, BASES['A']) for s in spectra]); G = np.array([rel(s, BASES['G']) for s in spectra])
    C = np.array([rel(s, BASES['C']) for s in spectra]); R = np.array([rel(s, RIBOSE) for s in spectra])
    r.update(nb10=np.mean(np.maximum(A, G) >= 0.10), nb2=np.mean(np.maximum(A, G) >= 0.02), A2=np.mean(A >= 0.02), G2=np.mean(G >= 0.02),
             C2=np.mean(C >= 0.02), AandG=np.mean((A >= 0.02) & (G >= 0.02)), ribose=np.mean(R >= 0.02))
    r['valid'] = bool(coverage >= 0.5 and hires >= 0.7)
    r['reason'] = '' if r['valid'] else ('low-resolution MS2' if hires < 0.7 else 'no low masses recorded')
    return r


def candidates():
    seen, out = set(done), []
    per_cat = collections.Counter()
    for cat, kw, pat in CATEGORIES:
        ps = get_json(API + 'search/projects?keyword=%s&pageSize=60&sortFields=_score' % urllib.parse.quote(kw))
        reviewed = 0
        for p in ps:
            if per_cat[cat] >= PER_CATEGORY + 3 or reviewed >= 25:
                break
            acc = p.get('accession')
            if acc in seen or not re.search(INSTR, ','.join(p.get('instruments', []))):
                continue
            reviewed += 1; seen.add(acc)
            fs = [f for f in get_json(API + 'projects/%s/files?pageSize=300' % acc)
                  if re.search(r'\.mgf$', f.get('fileName', ''), re.I) and 'pride.mgf' not in f['fileName'].lower()
                  and 2.5e7 < f.get('fileSizeBytes', 0) < 1.5e9]
            if not fs:
                continue
            pref = [f for f in fs if pat and re.search(pat, f['fileName'], re.I)] or fs
            f = sorted(pref, key=lambda f: -f['fileSizeBytes'])[len(pref) // 2]      # median size: avoids extremes
            url = next((l['value'] for l in f['publicFileLocations'] if l['name'].startswith('FTP')), '').replace('ftp://', 'https://')
            if url:
                out.append((cat, acc, p.get('title', '')[:70].replace('\t', ' '), f['fileName'], f['fileSizeBytes'], url))
                per_cat[cat] += 1
                print('candidate %-17s %s %s' % (cat, acc, f['fileName'][:60]), flush=True)
    return out


done = set()
if __name__ == '__main__':
    if os.path.exists(OUTPUT):
        done = {l.split('\t')[1] for l in open(OUTPUT, encoding='utf-8').read().split('\n')[1:] if l}
    else:
        open(OUTPUT, 'w', encoding='utf-8').write('category\tproject\ttitle\tfile\tMB\tn_ms2\tvalid\treason\tppm\t'
                                                 'nucleobase_10pct\tnucleobase_2pct\tadenine\tguanine\tcytosine\tA_and_G\tribose\n')
    cands = [c for c in candidates() if c[1] not in done]
    print('%d files to measure' % len(cands), flush=True)

    def one(c):
        cat, acc, tit, fn, size, url = c
        try:
            r = measure(url, size)
        except Exception as ex:
            r = {'n_ms2': 0, 'valid': False, 'reason': 'error: %s' % str(ex)[:60]}
        return c, r
    with ThreadPoolExecutor(4) as ex:
        for (cat, acc, tit, fn, size, url), r in ex.map(one, cands):
            g = lambda k: ('%.4f' % r[k]) if k in r else ''
            with open(OUTPUT, 'a', encoding='utf-8') as f:
                f.write('\t'.join([cat, acc, tit, fn, '%.0f' % (size / 1e6), str(r['n_ms2']), str(r['valid']), r.get('reason', ''),
                                   ('%.1f' % r['ppm']) if 'ppm' in r else '', g('nb10'), g('nb2'), g('A2'), g('G2'), g('C2'), g('AandG'), g('ribose')]) + '\n')
            print('%-17s %s n=%d valid=%s %s nucleobase>=10%%: %s' % (cat, acc, r['n_ms2'], r['valid'], r.get('reason', ''), g('nb10')), flush=True)
    print('DONE')
