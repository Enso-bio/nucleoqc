#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Atlas of nucleic-acid contamination in PRIDE (random sample of projects).

Projects are drawn at random from all of PRIDE (random pages of the full listing), so the estimated frequency does
not depend on a search term. For each project with uncompressed MGF or mzML, up to two files are measured (to test
whether contamination is a property of the project or of single runs) and protocol metadata are stored.

Usage: python atlas.py [target_number_of_projects]          -> atlas.tsv (resumable)"""
import sys, os, re, json, random, base64, zlib, collections, time, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import nucleic_census as C

OUTPUT = 'atlas.tsv'
TARGET = int(sys.argv[1]) if __name__ == '__main__' and len(sys.argv) > 1 else 200
CHUNK = 4_000_000
INSTR = r'Q Exactive|Exploris|Fusion|Lumos|Eclipse|Ascend|Astral|timsTOF|TripleTOF|ZenoTOF|Orbitrap Elite|Velos|Synapt|maXis|impact'
RNG = random.Random(20261112)

FLAGS = {   # feature -> pattern searched in title + description + sample protocol + keywords
    'phospho_enrichment': r'phospho|tio2|imac|fe-nta|fe3\+|titanium dioxide',
    'nuclease': r'benzonase|nuclease|dnase|rnase|micrococcal|pierce universal',
    'sonication': r'sonicat|bioruptor|ultrason|covaris',
    'ip_ap_ms': r'immunoprecipitat|pull-?down|affinity purif|ap-ms|bioid|turboid|proximity label|co-ip',
    'nucleic_acid_binding': r'rna-binding|rna binding|chromatin|histone|ribosom|interactome capture|dna-binding|nucleosome|spliceosom|nuclear extract',
    'cross_linking': r'cross-?link|xl-ms|dsso|dss\b|formaldehyde',
    'isobaric_labelling': r'\btmt|itraq|tandem mass tag',
    'body_fluids': r'plasma|serum|urine|cerebrospinal|\bcsf\b|saliva',
    'tissue': r'tissue|biops|ffpe|tumou?r',
    'in_gel_digestion': r'in-gel|in gel digest|sds-page',
    'sp3_strap_fasp': r'\bsp3\b|s-trap|strap|fasp',
    'glyco': r'glyco|hilic',
    'acetyl_ubiquitin': r'acetyl-?lysine|acetylome|ubiquitin|k-?gg|diglycine',
    'bacteria': r'bacteri|escherichia|e\. coli|bacillus|staphylococ|mycobacter|pseudomonas|streptococ',
}

# Studies in which nucleic acid or a nucleotide is part of the experiment (signal by design, not contamination):
# text rule applied to ALL projects, plus the projects found by manual review of every contaminated project.
BY_DESIGN = (r'chip-ms|chip-mass|\brime\b|\bchip\b|nucleosome|adp-ribos|protein-rna|rna-protein|protein-dna|'
             r'uv[- ]?(induced |inducible )?cross-?link|rna interactome|interactome capture|\b(par-|i|e)?clip\b|oligonucleotide|aptamer|'
             r'limited proteolysis|lip-ms|in vitro.{0,40}\batp\b')
MANUAL_BY_DESIGN = {'PXD067389': 'ADP-ribosylation', 'PXD050648': 'DNA binding + UV cross-linking',
                    'PXD002656': 'protein-RNA cross-linking', 'PXD042388': 'cross-linked RNA pol II complex',
                    'PXD012818': 'ChIP-MS (RIME)', 'PXD068775': 'nucleosome + UV cross-linking',
                    'PXD047868': 'ChIP-MS', 'PXD044412': 'LiP-MS with metabolites', 'PXD038387': 'proteasome in vitro with ATP',
                    'PXD034348': 'in vitro kinase assay (ATP; adenine only)'}


# ---------------- reading mzML chunks ----------------
_RE_SPEC = re.compile(r'<spectrum\b.*?</spectrum>', re.S)
_RE_BDA = re.compile(r'<binaryDataArray\b.*?</binaryDataArray>', re.S)
_RE_BIN = re.compile(r'<binary>([^<]*)</binary>')
_RE_PMZ = re.compile(r'accession="MS:1000744"[^>]*value="([0-9.eE+-]+)"')
_RE_LV = re.compile(r'accession="MS:1000511"[^>]*value="(\d+)"')


def mzml_spectra(txt):
    out = []
    for blk in _RE_SPEC.findall(txt):
        lv = _RE_LV.search(blk)
        if not lv or lv.group(1) == '1':
            continue
        pm = _RE_PMZ.search(blk)
        if not pm:
            continue
        arr = {}
        for b in _RE_BDA.findall(blk):
            if 'MS:1002312' in b or 'MS:1002313' in b or 'MS:1002314' in b:      # numpress
                return None
            t = _RE_BIN.search(b)
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
        if 'mz' not in arr or 'it' not in arr or len(arr['mz']) != len(arr['it']) or len(arr['mz']) < 10:
            continue
        mz, it = arr['mz'], arr['it']
        if 'MS:1000128' in blk:                                       # profile -> local maxima
            k = np.where((it[1:-1] > it[:-2]) & (it[1:-1] >= it[2:]) & (it[1:-1] > 0))[0] + 1
            mz, it = mz[k], it[k]
            if len(mz) < 10:
                continue
        o = np.argsort(mz); mz, it = mz[o], it[o]
        out.append({'mz': mz, 'it': it, 'low': mz[mz < 380], 'esc': mz[:0], 'pmz': float(pm.group(1)), 'z': 0})
    return out


REFS_RES = np.array([86.09643, 110.07127, 120.08078, 147.11280, 175.11895])   # immonium L/I, H, F; y1 K, R


def resolution(spectra):
    """Fraction of reference ions (most intense peak within +-0.2) that fall within +-0.005 of their exact mass."""
    ok = tot = 0
    for s in spectra:
        for m in REFS_RES:
            a, b = np.searchsorted(s['mz'], [m - 0.2, m + 0.2])
            if b > a:
                k = a + int(np.argmax(s['it'][a:b]))
                tot += 1; ok += abs(s['mz'][k] - m) <= C.TOL
    return ok / tot if tot >= 100 else 0.0


def measure(url, size, fmt):
    offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * CHUNK else [0]
    with ThreadPoolExecutor(3) as ex:
        chunks = list(ex.map(lambda a: C.byte_range(url, a, min(size - 1, a + CHUNK - 1)), offs))
    spectra = []
    for t in chunks:
        e = (C.mgf_spectra if fmt == 'mgf' else mzml_spectra)(t.decode('latin-1'))
        if e is None:
            return {'n_ms2': 0, 'valid': False, 'reason': 'numpress'}
        spectra += e
    mb = sum(len(t) for t in chunks) / 1e6
    if len(spectra) < 200:
        return {'n_ms2': len(spectra), 'valid': False, 'reason': 'too few spectra', 'mb': mb}
    ppm = C.recalibrate(spectra)
    cov = np.mean([s['mz'][0] < 136.0 for s in spectra])
    hires = resolution(spectra)
    A = np.array([C.rel(s, C.BASES['A']) for s in spectra]); G = np.array([C.rel(s, C.BASES['G']) for s in spectra])
    Cc = np.array([C.rel(s, C.BASES['C']) for s in spectra]); R = np.array([C.rel(s, C.RIBOSE) for s in spectra])
    r = {'n_ms2': len(spectra), 'mb': mb, 'ppm': ppm, 'coverage_136': cov, 'resolution_score': hires,
         'nb10': np.mean(np.maximum(A, G) >= 0.10), 'nb2': np.mean(np.maximum(A, G) >= 0.02),
         'A2': np.mean(A >= 0.02), 'G2': np.mean(G >= 0.02), 'C2': np.mean(Cc >= 0.02), 'ribose2': np.mean(R >= 0.02)}
    r['valid'] = bool(cov >= 0.5 and hires >= 0.7)
    r['reason'] = '' if r['valid'] else ('low-resolution MS2' if hires < 0.7 else 'MS2 does not cover m/z 136')
    return r


# ---------------- sampling ----------------
def as_text(v):
    if isinstance(v, str):
        return v
    return '; '.join((x.get('name', '') if isinstance(x, dict) else str(x)) for x in (v or []))


def eligible_files(acc):
    types = C.get_json(C.API + 'files/getCountOfFilesByType/' + acc)      # quick filter: no PEAK files, no MGF/mzML
    if not isinstance(types, dict) or not types.get('PEAK'):
        return []
    fs = C.get_json(C.API + 'projects/%s/files?pageSize=500' % acc)
    out = []
    for f in fs if isinstance(fs, list) else []:
        n = f.get('fileName', ''); t = f.get('fileSizeBytes', 0) or 0
        fmt = 'mgf' if re.search(r'\.mgf$', n, re.I) else ('mzml' if re.search(r'\.mzml$', n, re.I) else '')
        if not fmt or 'pride.mgf' in n.lower() or not (1.5e7 < t < 3e9):
            continue
        url = next((l['value'] for l in f.get('publicFileLocations', []) if l['name'].startswith('FTP')), '')
        if url:
            out.append((n, t, fmt, url.replace('ftp://', 'https://')))
    return out


def project_record(p):
    acc = p.get('accession')
    ins = as_text(p.get('instruments'))
    if not re.search(INSTR, ins, re.I):
        return None
    fs = eligible_files(acc)
    if not fs:
        return None
    RNG2 = random.Random(acc)
    chosen = RNG2.sample(fs, min(2, len(fs)))
    txt = ' '.join(str(p.get(k, '')) for k in ('title', 'projectDescription', 'sampleProcessingProtocol', 'keywords')).lower()
    meta = {'project': acc, 'year': str(p.get('submissionDate', ''))[:4], 'instrument': ins[:80],
            'organism': as_text(p.get('organisms'))[:60], 'country': as_text(p.get('countries'))[:40],
            'experiment_type': as_text(p.get('experimentTypes'))[:60], 'n_eligible_files': len(fs),
            'title': ' '.join(str(p.get('title', '')).split())[:80]}     # no tabs or line breaks in the TSV
    for k, pat in FLAGS.items():
        meta[k] = int(bool(re.search(pat, txt)))
    return meta, chosen


COLS_META = ['project', 'year', 'instrument', 'organism', 'country', 'experiment_type', 'n_eligible_files'] + list(FLAGS) + ['title']
COLS_MEAS = ['file', 'format', 'file_MB', 'n_ms2', 'valid', 'reason', 'ppm', 'coverage_136', 'resolution_score', 'nb10', 'nb2', 'A2', 'G2', 'C2', 'ribose2']

if __name__ == '__main__':
    done = set()
    if os.path.exists(OUTPUT):
        done = {l.split('\t')[0] for l in open(OUTPUT, encoding='utf-8').read().split('\n')[1:] if l}
    else:
        open(OUTPUT, 'w', encoding='utf-8').write('\t'.join(COLS_META + COLS_MEAS) + '\n')
    measured = len(done)
    pages = list(range(int(C.get_json(C.API + 'projects/count') or 42000) // 100))
    RNG.shuffle(pages)
    lock = threading.Lock()
    t0 = time.time(); mbtot = [0.0]

    def process(p):
        r = project_record(p)
        if not r:
            return 0
        meta, chosen = r
        rows = []
        for n, t, fmt, url in chosen:
            try:
                m = measure(url, t, fmt)
            except Exception as ex:
                m = {'n_ms2': 0, 'valid': False, 'reason': 'error: %s' % str(ex)[:50]}
            m.update(file=n, format=fmt, file_MB='%.0f' % (t / 1e6))
            mbtot[0] += m.get('mb', 0)
            rows.append(m)
        with lock:
            with open(OUTPUT, 'a', encoding='utf-8') as f:
                for m in rows:
                    g = lambda k: ('%.4f' % m[k]) if isinstance(m.get(k), float) else str(m.get(k, ''))
                    f.write('\t'.join([str(meta[k]) for k in COLS_META] + [g(k) for k in COLS_MEAS]) + '\n')
            print('%s %-6s %s | %s | %.0f MB in %.0f min' % (meta['project'], meta['year'], meta['instrument'][:25],
                  ', '.join('%s nb10=%s' % (m['reason'] or 'ok', ('%.4f' % m['nb10']) if 'nb10' in m else '-') for m in rows),
                  mbtot[0], (time.time() - t0) / 60), flush=True)
        return 1

    def page(pg):
        ps = C.get_json(C.API + 'projects?pageSize=100&page=%d' % pg)
        ps = [p for p in (ps if isinstance(ps, list) else []) if p.get('accession') not in done]
        random.Random(pg).shuffle(ps)
        return ps[:25]                          # few per page: the sample is spread over all of PRIDE

    with ThreadPoolExecutor(10) as ex:
        for i in range(0, len(pages), 4):
            if measured >= TARGET:
                break
            ps = [p for batch in ex.map(page, pages[i:i + 4]) for p in batch]
            for p in ps:
                done.add(p.get('accession'))
            measured += sum(ex.map(process, ps))
            print('-- pages %s: %d projects measured' % (pages[i:i + 4], measured), flush=True)
    print('DONE')
