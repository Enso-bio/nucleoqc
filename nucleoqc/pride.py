# -*- coding: utf-8 -*-
"""PRIDE access (API v3): list a project's files and download byte ranges over HTTPS."""
import json
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from . import reading

API = 'https://www.ebi.ac.uk/pride/ws/archive/v3/'
UA = 'nucleoqc'


def _get_json(u):
    for i in range(3):
        try:
            rq = urllib.request.Request(u, headers={'Accept': 'application/json', 'User-Agent': UA})
            with urllib.request.urlopen(rq, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 * (i + 1))
    raise ConnectionError('could not query PRIDE: %s' % u)


def byte_range(url, a, b):
    for i in range(4):
        try:
            rq = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (a, b), 'User-Agent': UA})
            with urllib.request.urlopen(rq, timeout=180) as r:
                return r.read()
        except Exception:
            time.sleep(3 * (i + 1))
    raise ConnectionError('could not download %s' % url)


def files(accession):
    """[(name, size, format, url)] of the project's uncompressed MGF/mzML files."""
    fs = _get_json(API + 'projects/%s/files?pageSize=1000' % accession)
    out = []
    for f in fs if isinstance(fs, list) else []:
        n = f.get('fileName', ''); t = f.get('fileSizeBytes', 0) or 0
        if not re.search(r'\.(mgf|mzml)$', n, re.I) or 'pride.mgf' in n.lower() or t <= 0:
            continue
        url = next((l['value'] for l in f.get('publicFileLocations', []) if l.get('name', '').startswith('FTP')), '')
        if url:
            out.append((n, t, reading.file_format(n), url.replace('ftp://', 'https://')))
    return sorted(out)


def read_remote(url, size, chunk=reading.CHUNK):
    offs = reading.offsets(size, chunk)          # files of <=16 MB are downloaded whole
    with ThreadPoolExecutor(3) as ex:
        return list(ex.map(lambda ab: byte_range(url, *ab), offs))
