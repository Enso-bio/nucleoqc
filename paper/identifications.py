#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Effect of nucleic-acid contamination on peptide identification with Sage (manuscript: identification section,
Figure 3).

For the human files of the impact analysis (contaminated and clean), the same three chunks are re-read and written as
mzML, searched with Sage against the reviewed human UniProt proteome (sharded search over the whole proteome, 1 %
spectrum-level FDR), and the identification rate of spectra with nucleobase is compared with that of spectra without
nucleobase in the same file and within +-10 Th. Spectra, mzML files and Sage results are kept in the search work folder
(nucleic_census.WORK).

Usage: python identifications.py SAGE_EXE FASTA --one-step     -> identifications_1p.tsv, ids_spectra_1p.tsv"""
import sys, os, re, json, base64, subprocess, collections, zlib, pickle
import numpy as np
import nucleic_census as C
import atlas

SAGE, FASTA = sys.argv[1], sys.argv[2]
# --one-step: full search directly against the whole proteome (sharded; the version used in the manuscript)
ONE_STEP = '--one-step' in sys.argv
ONLY = [a for a in sys.argv[3:] if a.startswith('PXD')]       # optional: restrict to these projects
SUF = '_1p' if ONE_STEP else ''
OUT_IDS, OUT_SPEC = 'identifications%s.tsv' % SUF, 'ids_spectra%s.tsv' % SUF
FASTA_TXT = open(FASTA, encoding='utf-8').read()
WORK = C.WORK
os.makedirs(WORK, exist_ok=True)

_RE_Z = re.compile(r'accession="MS:1000041"[^>]*value="(\d+)"')
_RE_RT = re.compile(r'accession="MS:1000016"[^>]*value="([0-9.eE+-]+)"[^>]*unitAccession="(UO:\d+)"')


def mgf_with_charge(txt):
    i, j = txt.find('BEGIN IONS'), txt.rfind('END IONS')
    out = []
    if i < 0 or j < i:
        return out
    for blk in txt[i:j].split('BEGIN IONS')[1:]:
        mz, it, pmz, z, rt = [], [], 0.0, 0, 0.0
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
            elif l.startswith('CHARGE='):
                m = re.match(r'(\d+)', l[7:])
                z = int(m.group(1)) if m else 0
            elif l.startswith('RTINSECONDS='):
                try:
                    rt = float(l[12:])
                except ValueError:
                    pass
        if len(mz) >= 10 and pmz > 0:
            m = np.array(mz); o = np.argsort(m)
            out.append({'mz': m[o], 'it': np.array(it)[o], 'pmz': pmz, 'z': z, 'rt': rt})
    return out


def mzml_with_charge(txt):
    out = []
    for blk in atlas._RE_SPEC.findall(txt):
        e = atlas.mzml_spectra(blk)
        if not e:
            continue
        s = e[0]; z = _RE_Z.search(blk); rt = _RE_RT.search(blk)
        s['z'] = int(z.group(1)) if z else 0
        s['rt'] = (float(rt.group(1)) * (60 if rt.group(2) == 'UO:0000031' else 1)) if rt else 0.0
        out.append(s)
    return out


def download(project, fn):
    fs = [x for x in atlas.eligible_files(project) if x[0] == fn]
    if not fs:
        return []
    n, size, fmt, url = fs[0]
    offs = [int(size * f) for f in (0.25, 0.5, 0.75)] if size > 4 * atlas.CHUNK else [0]
    spectra = []
    for a in offs:
        t = C.byte_range(url, a, min(size - 1, a + atlas.CHUNK - 1)).decode('latin-1')
        spectra += (mgf_with_charge if fmt == 'mgf' else mzml_with_charge)(t)
    for s in spectra:
        s['low'] = s['mz'][s['mz'] < 380]; s['esc'] = s['mz'][:0]
    C.recalibrate(spectra)
    return spectra


def b64(a):
    return base64.b64encode(np.asarray(a, '<f8').tobytes()).decode()


def write_mzml(spectra, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">\n'
                '<run id="r"><spectrumList count="%d">\n' % len(spectra))
        for i, s in enumerate(spectra):
            pre = '<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="%.6f"/>' % s['pmz']
            if s['z']:
                pre += '<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="%d"/>' % s['z']
            f.write('<spectrum index="%d" id="scan=%d" defaultArrayLength="%d">'
                    '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>'
                    '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>'
                    '<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%.3f" unitCvRef="UO" unitAccession="UO:0000010"/></scan></scanList>'
                    '<precursorList count="1"><precursor><isolationWindow>'
                    '<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="%.6f"/>'
                    '<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.0"/>'
                    '<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="1.0"/></isolationWindow>'
                    '<selectedIonList count="1"><selectedIon>%s</selectedIon></selectedIonList>'
                    '<activation><cvParam cvRef="MS" accession="MS:1000422" name="HCD"/></activation></precursor></precursorList>'
                    '<binaryDataArrayList count="2">'
                    '<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
                    '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>'
                    '<binary>%s</binary></binaryDataArray>'
                    '<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
                    '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>'
                    '<binary>%s</binary></binaryDataArray></binaryDataArrayList></spectrum>\n'
                    % (i, i, len(s['mz']), s['rt'] / 60.0, s['pmz'], pre, b64(s['mz']), b64(s['it'])))
        f.write('</spectrumList></run></mzML>\n')


def labelling(spectra):
    """Detects isobaric labelling from its reporter ions (fraction of spectra with the ion)."""
    fr = lambda m: np.mean([C.rel(s, m) > 0 for s in spectra])
    if fr(126.12773) > 0.3 and fr(127.12476) + fr(127.13108) > 0.3:
        return ('TMTpro', 304.20715) if fr(134.14816) > 0.1 or fr(134.15448) > 0.1 else ('TMT', 229.16293)
    if fr(114.11068) > 0.3 and fr(117.11497) > 0.3:
        return ('iTRAQ8', 304.20536) if fr(113.10732) > 0.2 and fr(121.12201) > 0.2 else ('iTRAQ4', 144.10206)
    return ('none', 0.0)


def config(mzml, output_dir, phospho, label_mass, fasta=None, light=False, extra=None):
    c = _config(mzml, output_dir, phospho, label_mass, fasta, light)
    c.update(extra or {})
    return c


def _config(mzml, output_dir, phospho, label_mass, fasta=None, light=False):
    static = {'C': 57.02146}
    var = {'M': [15.99491]} if light else {'M': [15.99491], '[': [42.01057]}
    if label_mass:
        static.update({'K': label_mass, '^': label_mass})
    if phospho and not light:
        var.update({'S': [79.96633], 'T': [79.96633], 'Y': [79.96633]})
    return {'database': {'bucket_size': 8192, 'enzyme': {'missed_cleavages': 1 if light else 2, 'min_len': 7, 'max_len': 30 if light else 40,
                                                         'cleave_at': 'KR', 'restrict': 'P'},
                         'fragment_min_mz': 150.0, 'fragment_max_mz': 2000.0, 'peptide_min_mass': 600.0, 'peptide_max_mass': 5000.0,
                         'ion_kinds': ['b', 'y'], 'min_ion_index': 2, 'static_mods': static, 'variable_mods': var,
                         'max_variable_mods': 1 if light else 3, 'decoy_tag': 'rev_', 'generate_decoys': True, 'fasta': fasta or FASTA},
            'precursor_tol': {'ppm': [-20, 20]}, 'fragment_tol': {'ppm': [-20, 20]}, 'isotope_errors': [0, 1],
            'deisotope': False, 'chimera': False, 'wide_window': False, 'predict_rt': False, 'min_peaks': 8, 'max_peaks': 150,
            'max_fragment_charge': 1, 'report_psms': 1, 'precursor_charge': [2, 4], 'output_directory': output_dir, 'mzml_paths': [mzml]}


def shards(n):
    """Splits the whole FASTA into n shards (alternating proteins); written only once."""
    paths = [os.path.join(WORK, 'shard%d_of%d.fasta' % (k, n)) for k in range(n)]
    if not all(os.path.exists(r) for r in paths):
        bl = FASTA_TXT.split('>')[1:]
        for k, r in enumerate(paths):
            open(r, 'w', encoding='utf-8').write(''.join('>' + b for b in bl[k::n]))
    return paths


def sharded_search(base, n_spec, phospho, label_mass, extra=None):
    """Full search against the whole proteome, in shards so that it fits in RAM.
    For each spectrum the best hyperscore over all shards (target or decoy) is kept and the FDR is computed once over
    those best hits (target-decoy competition by hyperscore)."""
    n = 8 if phospho else 2
    best = np.full(n_spec, -np.inf); target = np.zeros(n_spec, bool)
    for k, fa in enumerate(shards(n)):
        out_dir = '%s_t%d' % (base, k)
        res = os.path.join(out_dir, 'results.sage.tsv')
        if not os.path.exists(res):
            json.dump(config(base + '.mzML', out_dir, phospho, label_mass, fasta=fa, extra=extra), open(out_dir + '.json', 'w'))
            r = subprocess.run([SAGE, out_dir + '.json', '--disable-telemetry-i-dont-want-to-improve-sage'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(res):
                print('   shard %d failed: %s' % (k, r.stderr[-300:]), flush=True); return None
        lin = open(res, encoding='utf-8').read().split('\n'); J = {c: i for i, c in enumerate(lin[0].split('\t'))}
        for l in lin[1:]:
            if not l:
                continue
            c = l.split('\t'); idx = int(re.search(r'(\d+)$', c[J['scannr']]).group(1)); h = float(c[J['hyperscore']])
            if h > best[idx]:
                best[idx] = h; target[idx] = c[J['label']] == '1'
    ok = np.isfinite(best); o = np.argsort(-best[ok]); idxs = np.where(ok)[0][o]
    t = np.cumsum(target[idxs]); d = np.cumsum(~target[idxs])
    q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    ident = np.zeros(n_spec, bool); ident[idxs] = target[idxs] & (q <= 0.01)
    return ident, np.where(ok, best, 0.0)


if __name__ == '__main__':
    rows = [l.split('\t') for l in open('atlas.tsv', encoding='utf-8').read().split('\n') if l]
    header = rows[0]; I = {k: i for i, k in enumerate(header)}
    meta = {(r[0], r[I['file']]): r for r in rows[1:] if len(r) == len(header)}
    Z = np.load('impact_spectra.npz')
    tasks = []
    for k in Z.files:
        g, p, fn = k.split('|')
        r = meta[(p, fn)]
        if r[I['organism']].startswith('Homo sapiens') and ';' not in r[I['organism']]:
            tasks.append((g, p, fn, r[I['phospho_enrichment']] == '1'))
    print('%d human files' % len(tasks), flush=True)
    done = set()
    if os.path.exists(OUT_IDS):
        done = {tuple(l.split('\t')[1:3]) for l in open(OUT_IDS, encoding='utf-8').read().split('\n')[1:] if l}
    else:
        open(OUT_IDS, 'w', encoding='utf-8').write('group\tproject\tfile\tphospho_enrichment\tlabelling\tn_ms2\tn_nucleobase_spectra\tid_rate_all\t'
                                                   'id_rate_nucleobase\tid_rate_no_nucleobase\tid_rate_expected_nucleobase\n')
        open(OUT_SPEC, 'w', encoding='utf-8').write('project\tfile\tindex\tpmz\tz\tnb\tid\thyperscore\n')
    for g, p, fn, phospho in tasks:
        if (p, fn) in done or (ONLY and p not in ONLY):
            continue
        base = os.path.join(WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        if os.path.exists(base + '.pkl'):
            spectra = pickle.load(open(base + '.pkl', 'rb'))
        else:
            spectra = download(p, fn)
            pickle.dump(spectra, open(base + '.pkl', 'wb'))
        if len(spectra) < 200:
            print('  too few spectra', p, fn); continue
        nb = np.array([max(C.rel(s, C.BASES['A']), C.rel(s, C.BASES['G'])) for s in spectra])
        lab, mass = labelling(spectra)
        write_mzml(spectra, base + '.mzML')
        if ONE_STEP:
            r = sharded_search(base, len(spectra), phospho, mass)
            if r is None:
                print('  Sage (sharded) failed', p, fn, flush=True); continue
            ident, hs = r
        else:
            # step 1: light search against the whole proteome -> proteins present
            json.dump(config(base + '.mzML', base + '_p1', phospho, mass, light=True), open(base + '_p1.json', 'w'))
            r = subprocess.run([SAGE, base + '_p1.json', '--disable-telemetry-i-dont-want-to-improve-sage'], capture_output=True, text=True)
            res1 = os.path.join(base + '_p1', 'results.sage.tsv')
            if r.returncode != 0 or not os.path.exists(res1):
                print('  Sage (step 1) failed', p, fn, r.stderr[-300:], flush=True); continue
            l1 = open(res1, encoding='utf-8').read().split('\n'); J1 = {k: i for i, k in enumerate(l1[0].split('\t'))}
            prot = set()
            for l in l1[1:]:
                c = l.split('\t')
                if len(c) > 3 and c[J1['label']] == '1' and float(c[J1['spectrum_q']]) <= 0.05:
                    prot.update(x for x in c[J1['proteins']].split(';') if not x.startswith('rev_'))
            sub = base + '_sub.fasta'; n_sub = 0
            with open(sub, 'w', encoding='utf-8') as out:
                for block in FASTA_TXT.split('>')[1:]:
                    if block.split(None, 1)[0] in prot:
                        out.write('>' + block); n_sub += 1
            if n_sub < 20:
                print('  too few proteins in step 1 (%d)' % n_sub, p, fn, flush=True); continue
            # step 2: full search against the proteins present
            json.dump(config(base + '.mzML', base + '_out', phospho, mass, fasta=sub), open(base + '.json', 'w'))
            r = subprocess.run([SAGE, base + '.json', '--disable-telemetry-i-dont-want-to-improve-sage'], capture_output=True, text=True)
            res = os.path.join(base + '_out', 'results.sage.tsv')
            if r.returncode != 0 or not os.path.exists(res):
                print('  Sage (step 2) failed', p, fn, r.stderr[-300:], flush=True); continue
            lin = open(res, encoding='utf-8').read().split('\n')
            h = lin[0].split('\t'); J = {k: i for i, k in enumerate(h)}
            ident = np.zeros(len(spectra), bool); hs = np.zeros(len(spectra))
            for l in lin[1:]:
                if not l:
                    continue
                c = l.split('\t'); idx = int(re.search(r'(\d+)$', c[J['scannr']]).group(1))
                hs[idx] = float(c[J['hyperscore']])
                if c[J['label']] == '1' and float(c[J['spectrum_q']]) <= 0.01:
                    ident[idx] = True
        m = nb >= 0.10
        pmz = np.array([s['pmz'] for s in spectra])
        # expected rate for spectra with nucleobase: that of spectra without nucleobase within +-10 Th in the same file
        ref = np.where(~m)[0]; o = np.argsort(pmz[ref]); ref, rp = ref[o], pmz[ref][o]
        expected = []
        for i in np.where(m)[0]:
            a, b = np.searchsorted(rp, [pmz[i] - 10, pmz[i] + 10])
            if b - a >= 5:
                expected.append(ident[ref[a:b]].mean())
        e_nb = np.mean(expected) if expected else np.nan
        with open(OUT_IDS, 'a', encoding='utf-8') as f:
            f.write('%s\t%s\t%s\t%d\t%s\t%d\t%d\t%.4f\t%s\t%.4f\t%s\n' % (g, p, fn, phospho, lab, len(spectra), m.sum(), ident.mean(),
                    ('%.4f' % ident[m].mean()) if m.sum() else '', ident[~m].mean(), ('%.4f' % e_nb) if expected else ''))
        with open(OUT_SPEC, 'a', encoding='utf-8') as f:
            for i, s in enumerate(spectra):
                f.write('%s\t%s\t%d\t%.4f\t%d\t%.3f\t%d\t%.2f\n' % (p, fn, i, s['pmz'], s['z'], nb[i], ident[i], hs[i]))
        print('%-12s %s %-36s %s phospho=%d n=%d nb=%d | ID all %.1f %% | with nb %s | without nb %.1f %% | expected nb %s'
              % (g, p, fn[:36], lab, phospho, len(spectra), m.sum(), 100 * ident.mean(),
                 ('%.1f %%' % (100 * ident[m].mean())) if m.sum() else '-', 100 * ident[~m].mean(),
                 ('%.1f %%' % (100 * e_nb)) if expected else '-'), flush=True)
    print('DONE')
