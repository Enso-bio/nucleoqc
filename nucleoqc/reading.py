# -*- coding: utf-8 -*-
"""Reading MS2 spectra from MGF and mzML (whole files or byte chunks).

Chunks may start and end in the middle of a spectrum: incomplete blocks are discarded.
Each spectrum is a dict with 'mz' (sorted), 'it' and 'pmz' (precursor m/z).
"""
import os
import re
import base64
import zlib
import numpy as np

MIN_PEAKS = 10


def _spectrum(mz, it, pmz):
    mz = np.asarray(mz, dtype=np.float64); it = np.asarray(it, dtype=np.float64)
    if len(mz) < MIN_PEAKS or len(mz) != len(it) or pmz <= 0 or it.max() <= 0:
        return None
    o = np.argsort(mz, kind='stable')
    return {'mz': mz[o], 'it': it[o], 'pmz': float(pmz)}


# ---------------------------------------------------------------- MGF
def mgf_spectra(txt):
    i, j = txt.find('BEGIN IONS'), txt.rfind('END IONS')
    if i < 0 or j < i:
        return []
    out = []
    for blk in txt[i:j].split('BEGIN IONS')[1:]:
        mz, it, pmz = [], [], 0.0
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
        s = _spectrum(mz, it, pmz)
        if s:
            out.append(s)
    return out


# ---------------------------------------------------------------- mzML
_RE_SPEC = re.compile(r'<spectrum\b.*?</spectrum>', re.S)
_RE_BDA = re.compile(r'<binaryDataArray\b.*?</binaryDataArray>', re.S)
_RE_BIN = re.compile(r'<binary>([^<]*)</binary>')
_RE_PMZ = re.compile(r'accession="MS:1000744"[^>]*value="([0-9.eE+-]+)"')
_RE_LEVEL = re.compile(r'accession="MS:1000511"[^>]*value="(\d+)"')
NUMPRESS = ('MS:1002312', 'MS:1002313', 'MS:1002314')


class NumpressError(ValueError):
    """The mzML uses numpress compression, which this reader does not decode."""


def mzml_spectra(txt):
    out = []
    for blk in _RE_SPEC.findall(txt):
        lv = _RE_LEVEL.search(blk)
        if not lv or lv.group(1) != '2':           # MS2 only (MS1 and SPS-MS3 reporter scans are skipped)
            continue
        pm = _RE_PMZ.search(blk)
        if not pm:
            continue
        arr = {}
        for b in _RE_BDA.findall(blk):
            if any(a in b for a in NUMPRESS):
                raise NumpressError('numpress compression is not supported; convert with msconvert without --numpress')
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
        if 'mz' not in arr or 'it' not in arr or len(arr['mz']) != len(arr['it']):
            continue
        mz, it = arr['mz'], arr['it']
        if 'MS:1000128' in blk:                    # profile spectrum -> local maxima
            k = np.where((it[1:-1] > it[:-2]) & (it[1:-1] >= it[2:]) & (it[1:-1] > 0))[0] + 1
            mz, it = mz[k], it[k]
        s = _spectrum(mz, it, float(pm.group(1)))
        if s:
            out.append(s)
    return out


def file_format(name):
    n = name.lower()
    if n.endswith('.mgf'):
        return 'mgf'
    if n.endswith('.mzml'):
        return 'mzml'
    raise ValueError('unrecognised format (uncompressed .mgf and .mzML are accepted): %s' % name)


def parse(chunks, fmt):
    f = mgf_spectra if fmt == 'mgf' else mzml_spectra
    spectra = []
    for t in chunks:
        spectra += f(t.decode('latin-1'))
    return spectra


# ---------------------------------------------------------------- file sampling
CHUNK = 4_000_000
POSITIONS = (0.25, 0.50, 0.75)


def offsets(size, chunk=CHUNK, whole=False):
    """Start and end (inclusive) of each chunk to read. Small files or whole=True: the entire file."""
    if whole or size <= 4 * chunk:
        return [(0, size - 1)]
    return [(int(size * f), min(size - 1, int(size * f) + chunk - 1)) for f in POSITIONS]


def read_local(path, chunk=CHUNK, whole=False):
    size = os.path.getsize(path)
    chunks = []
    with open(path, 'rb') as fh:
        for a, b in offsets(size, chunk, whole):
            fh.seek(a); chunks.append(fh.read(b - a + 1))
    return chunks, size
