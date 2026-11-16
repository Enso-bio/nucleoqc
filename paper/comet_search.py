#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Second search engine (manuscript: robustness, Figure 4). Repeats the search with Comet (XCorr / E-value) on exactly
the same spectra of the contaminated files and repeats the within-file comparison.

Parameters equivalent to Sage: trypsin, 2 missed cleavages, 7-40 residues, 600-5000 Da, precursor +-20 ppm, isotope
error 0/+1, fragment bin 0.02 Th, fixed carbamidomethyl-Cys, variable Met oxidation and protein N-terminal acetylation,
phospho S/T/Y for phosphoproteomes, fixed TMT/TMTpro/iTRAQ if detected. Concatenated reversed decoys; 1 % spectrum-level
FDR on the best hit per spectrum (target-decoy competition by E-value). Whole proteome.
The stored spectra are read from the search work folder (see nucleic_census.WORK); Comet writes its files to work/ next
to the Comet executable, which must also contain comet.params.new (comet -p).

Usage: python comet_search.py COMET_EXE FASTA   -> comet_ids.tsv, comet_results.txt"""
import sys, os, re, zlib, pickle, subprocess, collections
import numpy as np
from scipy import stats
import nucleic_census as NC

COMET, FASTA = sys.argv[1], sys.argv[2]
DIR = os.path.join(os.path.dirname(COMET), 'work'); os.makedirs(DIR, exist_ok=True)
PARAMS = open(os.path.join(os.path.dirname(COMET), 'comet.params.new'), encoding='utf-8').read()
LABELS = {'TMT': 229.162932, 'TMTpro': 304.207146, 'iTRAQ4': 144.102063, 'iTRAQ8': 304.205360}

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def set_param(txt, key, value):
    t, n = re.subn(r'(?m)^%s\s*=.*$' % re.escape(key), lambda _: '%s = %s' % (key, value), txt)
    assert n == 1, key
    return t


def params(phospho, label):
    t = PARAMS
    for k, v in [('database_name', FASTA), ('decoy_search', 1), ('num_threads', 0), ('peptide_mass_tolerance_upper', 20.0),
                 ('peptide_mass_tolerance_lower', -20.0), ('peptide_mass_units', 2), ('isotope_error', 1),
                 ('allowed_missed_cleavage', 2), ('fragment_bin_tol', 0.02), ('fragment_bin_offset', 0.0),
                 ('theoretical_fragment_ions', 0), ('output_pepxmlfile', 0), ('output_txtfile', 1), ('num_output_lines', 1),
                 ('precursor_charge', '2 4'), ('digest_mass_range', '600.0 5000.0'), ('peptide_length_range', '7 40'),
                 ('max_variable_mods_in_peptide', 3), ('minimum_peaks', 8),
                 ('variable_mod02', '42.010565 ^ 0 1 -1 0 0 0.0'),
                 ('variable_mod03', ('79.966331 STY 0 3 -1 0 0 0.0' if phospho else '0.0 X 0 3 -1 0 0 0.0'))]:
        t = set_param(t, k, v)
    if label in LABELS:
        t = set_param(t, 'add_K_lysine', LABELS[label]); t = set_param(t, 'add_Nterm_peptide', LABELS[label])
    return t


def write_mgf(spectra, path):
    with open(path, 'w', encoding='utf-8') as f:
        for i, s in enumerate(spectra):
            f.write('BEGIN IONS\nTITLE=idx%d\nSCANS=%d\nRTINSECONDS=%.2f\nPEPMASS=%.6f\n' % (i, i + 1, s['rt'], s['pmz']))
            if s['z']:
                f.write('CHARGE=%d+\n' % s['z'])
            f.write(''.join('%.5f %.2f\n' % (m, x) for m, x in zip(s['mz'], s['it'])))
            f.write('END IONS\n')


def read_results(txt_path, n):
    """Best hit per spectrum (lowest E-value), 1 % FDR by target-decoy competition."""
    lin = open(txt_path, encoding='utf-8').read().split('\n')
    k = next(i for i, l in enumerate(lin) if l.startswith('scan\t'))
    J = {c: i for i, c in enumerate(lin[k].split('\t'))}
    best = np.full(n, np.inf); target = np.zeros(n, bool)
    for l in lin[k + 1:]:
        c = l.split('\t')
        if len(c) < len(J):
            continue
        i = int(c[J['scan']]) - 1; e = float(c[J['e-value']])
        if 0 <= i < n and e < best[i]:
            best[i] = e; target[i] = not c[J['protein']].startswith('DECOY_')
    ok = np.isfinite(best); idx = np.where(ok)[0][np.argsort(best[ok])]
    t = np.cumsum(target[idx]); d = np.cumsum(~target[idx])
    q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    ident = np.zeros(n, bool); ident[idx] = target[idx] & (q <= 0.01)
    return ident


if __name__ == '__main__':
    A = [a for a in tab('identifications_1p.tsv') if float(a['id_rate_all']) >= 0.02 and a['group'] == 'contaminated'
         and int(a['n_nucleobase_spectra']) >= 20]
    E = collections.defaultdict(dict)
    for r in tab('ids_spectra_1p.tsv'):
        E[(r['project'], r['file'])][int(r['index'])] = (float(r['pmz']), float(r['nb']), r['id'] == '1')
    OUT = 'comet_ids.tsv'
    done = {tuple(l.split('\t')[:2]) for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
    if not done:
        open(OUT, 'w', encoding='utf-8').write('project\tfile\tindex\tid_comet\n')
    print('%d contaminated files' % len(A), flush=True)
    for a in A:
        p, fn = a['project'], a['file']
        if (p, fn) in done:
            continue
        name = '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000)
        spectra = pickle.load(open(os.path.join(NC.WORK, name + '.pkl'), 'rb'))
        mgf = os.path.join(DIR, name + '.mgf'); par = os.path.join(DIR, name + '.params')
        write_mgf(spectra, mgf)
        open(par, 'w', encoding='utf-8').write(params(a['phospho_enrichment'] == '1', a['labelling']))
        r = subprocess.run([COMET, '-P' + par, mgf], capture_output=True, text=True, cwd=DIR)
        txt = os.path.join(DIR, name + '.txt')
        if r.returncode != 0 or not os.path.exists(txt):
            print('  Comet failed', p, fn, r.stdout[-300:], r.stderr[-300:], flush=True); continue
        ident = read_results(txt, len(spectra))
        with open(OUT, 'a', encoding='utf-8') as f:
            f.write(''.join('%s\t%s\t%d\t%d\n' % (p, fn, i, ident[i]) for i in range(len(spectra))))
        sg = np.array([E[(p, fn)][i][2] for i in range(len(spectra))])
        print('%s %-36s Comet %5.1f %% | Sage %5.1f %% | both %5.1f %%' % (p, fn[:36], 100 * ident.mean(), 100 * sg.mean(),
              100 * (ident & sg).mean()), flush=True)

    # ---- within-file comparison, as with Sage
    C = collections.defaultdict(dict)
    for r in tab(OUT):
        C[(r['project'], r['file'])][int(r['index'])] = r['id_comet'] == '1'
    lines = []
    out = lambda s: (print(s, flush=True), lines.append(s))
    out('\nWithin-file comparison (nb >= 10 %% vs neighbours with nb < 10 %% within +-10 Th), Sage versus Comet')
    res = {'Sage': ([], []), 'Comet': ([], [])}
    for (p, fn), cm in C.items():
        v = E[(p, fn)]; n = len(v)
        pmz = np.array([v[i][0] for i in range(n)]); nb = np.array([v[i][1] for i in range(n)])
        ids = {'Sage': np.array([v[i][2] for i in range(n)]), 'Comet': np.array([cm[i] for i in range(n)])}
        m = nb >= 0.10; c = np.where(~m)[0]; o = np.argsort(pmz[c]); c, cp = c[o], pmz[c][o]
        row = []
        for engine, idn in ids.items():
            obs, exp_ = [], []
            for i in np.where(m)[0]:
                lo, hi = np.searchsorted(cp, [pmz[i] - 10, pmz[i] + 10])
                if hi - lo >= 5:
                    obs.append(idn[i]); exp_.append(idn[c[lo:hi]].mean())
            res[engine][0].append(np.mean(obs)); res[engine][1].append(np.mean(exp_))
            row.append('%s %5.1f %% vs %5.1f %%' % (engine, 100 * np.mean(obs), 100 * np.mean(exp_)))
        out('   %s %-34s %s' % (p, fn[:34], ' | '.join(row)))
    for engine, (oo, ee) in res.items():
        oo, ee = np.array(oo), np.array(ee)
        out('   %-5s %d files: observed %5.1f %% vs expected %5.1f %% (medians) | median ratio %.2f | lower in %d of %d | Wilcoxon p = %.2g'
            % (engine, len(oo), 100 * np.median(oo), 100 * np.median(ee), np.median(oo / np.maximum(ee, 1e-9)), int((oo < ee).sum()), len(oo),
               stats.wilcoxon(oo, ee).pvalue))
    tot_s = sum(len(v) for v in C.values())
    agreement = []
    for (p, fn), cm in C.items():
        v = E[(p, fn)]
        s = np.array([v[i][2] for i in range(len(v))]); c_ = np.array([cm[i] for i in range(len(v))])
        agreement.append(((s & c_).sum(), s.sum(), c_.sum()))
    a_ = np.array(agreement)
    out('   total identifications: Sage %d, Comet %d, in both %d (%d spectra)' % (a_[:, 1].sum(), a_[:, 2].sum(), a_[:, 0].sum(), tot_s))
    open('comet_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
