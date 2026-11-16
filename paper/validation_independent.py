#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Independent validation of the loss estimated by nucleoqc (manuscript: nucleoqc section, Figure 5).

Contaminated files from projects NOT used to calibrate the dose-response: other species, and three human files with
nucleic acid present without cross-linking (ChIP-MS, nucleosomes without UV). For each, nucleoqc estimates the loss from
the spectra alone; the file is then searched with Sage (same parameters and sharded search, with the proteome of its
species) and the actual loss is measured with the within-file comparison. Proteomes are downloaded from UniProt into
proteomes/ inside the search work folder.

Usage: python validation_independent.py SAGE_EXE HUMAN_FASTA  -> validation_independent.tsv"""
import sys, os, re, json, zlib, pickle, math, subprocess, urllib.request, urllib.parse
SAGE, FASTA_H = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], SAGE, FASTA_H, '--one-step']
import numpy as np
from scipy import stats
import nucleic_census as C
import identifications as I
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from nucleoqc import metrics as NQ

PROT_DIR = os.path.join(I.WORK, 'proteomes'); os.makedirs(PROT_DIR, exist_ok=True)
PROTEOMES = {   # name -> UniProt query
    'human': None,
    'neurospora': 'organism_id:367110 AND keyword:KW-1185',
    'bcereus_atcc14579': 'organism_id:226900 AND keyword:KW-1185',
    'mouse_sp': 'organism_id:10090 AND reviewed:true',
    'dog': 'organism_id:9615 AND keyword:KW-1185',
    'candida': 'organism_id:237561 AND keyword:KW-1185',
    'saureus_nctc8325': 'organism_id:93061 AND keyword:KW-1185',
    'salmonella_lt2': 'organism_id:99287 AND keyword:KW-1185',
}
FILES = [   # project, file name prefix, proteome, phospho, note
    ('PXD013964', 'Benz_Ncrassa_A_04_01_13_phos', 'neurospora', True, 'phosphoproteome'),
    ('PXD030165', 'Q26761_MS20-013_SN2', 'bcereus_atcc14579', False, 'bacterial secretome'),
    ('PXD030165', 'Q26774_MS20-013_SN15', 'bcereus_atcc14579', False, 'bacterial secretome'),
    ('PXD018783', '022018_BL6_DC_FLT3-t1', 'mouse_sp', False, 'dendritic cells'),
    ('PXD018783', '022018_BL6_DC_FLT3_t2', 'mouse_sp', False, 'dendritic cells'),
    ('PXD028074', 'KapB1_MGF__F066073_', 'dog', False, ''),
    ('PXD067736', '20-89PD22_3h-PTM.mzML', 'candida', True, 'phosphoproteome'),
    ('PXD022858', 'ZN05_AllF063786', 'saureus_nctc8325', False, ''),
    ('PXD011041', 'Salmonella_230614_04', 'salmonella_lt2', False, ''),
    ('PXD011041', 'Salmonella_230614_11', 'salmonella_lt2', False, ''),
    ('PXD012818', 'GPM00500047074', 'human', False, 'ChIP-MS (RIME)'),
    ('PXD047868', '230206_Q1_PR_guide_pool', 'human', False, 'ChIP-MS'),
    ('PXD068775', 'Gireesh_etal_2025_NCPs_CPC_1_190_noUV_Rep2', 'human', False, 'nucleosomes without UV'),
]
OUT = 'validation_independent.tsv'


def proteome(name):
    if PROTEOMES[name] is None:
        return FASTA_H
    path = os.path.join(PROT_DIR, name + '.fasta')
    if not os.path.exists(path):
        import time
        u = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=' + urllib.parse.quote(PROTEOMES[name])
        for attempt in range(6):
            try:
                data = urllib.request.urlopen(urllib.request.Request(u, headers={'User-Agent': 'python-urllib'}), timeout=600).read()
                break
            except Exception as ex:
                print('   UniProt: %s; retrying in %d s' % (ex, 30 * (attempt + 1)), flush=True); time.sleep(30 * (attempt + 1))
        else:
            raise ConnectionError('UniProt does not respond')
        if data.count(b'>') < 100:
            raise ValueError('empty or incomplete proteome: %s' % name)
        open(path, 'wb').write(data)
    return path


def search(base, n_spec, phospho, label_mass, fasta):
    txt = open(fasta, encoding='utf-8').read(); bl = txt.split('>')[1:]
    n = max(1, math.ceil(len(bl) / (2500 if phospho else 10000)))
    tag = os.path.splitext(os.path.basename(fasta))[0]
    best = np.full(n_spec, -np.inf); target = np.zeros(n_spec, bool)
    for k in range(n):
        fa = os.path.join(I.WORK, '%s_shard%d_of%d.fasta' % (tag, k, n))
        if not os.path.exists(fa):
            open(fa, 'w', encoding='utf-8').write(''.join('>' + b for b in bl[k::n]))
        out_dir = '%s_%s_t%d' % (base, tag, k); res = os.path.join(out_dir, 'results.sage.tsv')
        if not os.path.exists(res):
            json.dump(I.config(base + '.mzML', out_dir, phospho, label_mass, fasta=fa), open(out_dir + '.json', 'w'))
            r = subprocess.run([SAGE, out_dir + '.json', '--disable-telemetry-i-dont-want-to-improve-sage'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(res):
                print('   shard %d failed: %s' % (k, r.stderr[-300:]), flush=True); return None
        lin = open(res, encoding='utf-8').read().split('\n'); J = {c: i for i, c in enumerate(lin[0].split('\t'))}
        for l in lin[1:]:
            if l:
                c = l.split('\t'); i = int(re.search(r'(\d+)$', c[J['scannr']]).group(1)); h = float(c[J['hyperscore']])
                if h > best[i]:
                    best[i] = h; target[i] = c[J['label']] == '1'
    ok = np.isfinite(best); idx = np.where(ok)[0][np.argsort(-best[ok])]
    t = np.cumsum(target[idx]); d = np.cumsum(~target[idx]); q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    ident = np.zeros(n_spec, bool); ident[idx] = target[idx] & (q <= 0.01)
    return ident


def measured_loss(nb, pmz, ident):
    m = nb >= 0.10; c = np.where(~m)[0]; o = np.argsort(pmz[c]); c, cp = c[o], pmz[c][o]
    exp_ = []
    for i in np.where(m)[0]:
        a, b = np.searchsorted(cp, [pmz[i] - 10, pmz[i] + 10])
        if b - a >= 5:
            exp_.append(ident[c[a:b]].mean())
    if not exp_:
        return np.nan, np.nan, np.nan
    e = np.mean(exp_); lost = m.sum() * (e - ident[m].mean())
    return lost / max(ident.sum() + lost, 1), ident[m].mean(), e


if __name__ == '__main__':
    import atlas
    done = {tuple(l.split('\t')[:2]) for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
    if not done:
        open(OUT, 'w', encoding='utf-8').write('project\tfile\tproteome\tnote\tn_ms2\tn_nucleobase_spectra\tnb10\tid_rate_all\tid_rate_nucleobase\t'
                                               'id_rate_expected_nucleobase\testimated_loss\tmeasured_loss\n')
    for p, pref, prot, phospho, note in FILES:
        fs = [f for f in atlas.eligible_files(p) if f[0].startswith(pref)]
        if not fs:
            print('%s %s: not found' % (p, pref), flush=True); continue
        fn = fs[0][0]
        if (p, fn) in done:
            continue
        base = os.path.join(I.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        if os.path.exists(base + '.pkl'):
            spectra = pickle.load(open(base + '.pkl', 'rb'))
        else:
            spectra = I.download(p, fn); pickle.dump(spectra, open(base + '.pkl', 'wb'))
        if len(spectra) < 200:
            print('%s %s: too few spectra' % (p, fn), flush=True); continue
        nb = np.array([max(C.rel(s, C.BASES['A']), C.rel(s, C.BASES['G'])) for s in spectra]); pmz = np.array([s['pmz'] for s in spectra])
        est = NQ.estimated_loss(nb)                       # the estimate is fixed BEFORE searching
        lab, mass = I.labelling(spectra)
        if not os.path.exists(base + '.mzML'):
            I.write_mzml(spectra, base + '.mzML')
        try:
            ident = search(base, len(spectra), phospho, mass, proteome(prot))
        except Exception as ex:
            print('%s %s: error %s' % (p, fn, ex), flush=True); continue
        if ident is None:
            continue
        meas, o_, e_ = measured_loss(nb, pmz, ident)
        with open(OUT, 'a', encoding='utf-8') as f:
            f.write('%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (p, fn, prot, note, len(spectra), (nb >= 0.1).sum(),
                    (nb >= 0.1).mean(), ident.mean(), o_, e_, est, meas))
        print('%s %-36s %-18s nb10 %5.1f %% | ID %5.1f %% | with nb %5.1f %% vs expected %5.1f %% | estimated loss %5.1f %% measured %5.1f %%'
              % (p, fn[:36], prot, 100 * (nb >= 0.1).mean(), 100 * ident.mean(), 100 * o_, 100 * e_, 100 * est, 100 * meas), flush=True)
    print('DONE')
