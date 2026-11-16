#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Table S10: extended ion panel (A, G, C, U, T, ribose; with +0.035 Da control) and RNA- or DNA-specific nucleotide
fragments for each of the 38 files searched with Sage. Reads the stored spectra from the search work folder
(nucleic_census.WORK). -> table_s10.tsv"""
import os, sys, zlib, pickle
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from nucleoqc import metrics as M
import nucleic_census as NC

H = 1.007276
BASES = {'A': 135.05450 + H, 'G': 151.04941 + H, 'C': 111.04326 + H, 'U': 112.02728 + H, 'T': 126.04293 + H, 'ribose': 97.02841}
L = [l.split('\t') for l in open('identifications_1p.tsv', encoding='utf-8').read().split('\n') if l]; header = L[0]
cols = ['project', 'file', 'group', 'n_ms2', 'nb10'] + ['%s_2pct' % k for k in BASES] + ['%s_control_2pct' % k for k in BASES] + \
       ['%s_1pct' % k for k in list(M.RNA_FRAGMENTS) + list(M.DNA_FRAGMENTS)] + ['type']
with open('table_s10.tsv', 'w', encoding='utf-8') as f:
    f.write('\t'.join(cols) + '\n')
    for r in L[1:]:
        a = dict(zip(header, r)); path = os.path.join(NC.WORK, '%s_%d.pkl' % (a['project'], zlib.crc32(a['file'].encode()) % 100000))
        if not os.path.exists(path):
            continue
        spectra = [{'mz': s['mz'], 'it': s['it'], 'pmz': s['pmz']} for s in pickle.load(open(path, 'rb'))]
        nb = np.mean([max(M.rel(s, M.ADENINE), M.rel(s, M.GUANINE)) >= 0.10 for s in spectra])
        v = [np.mean([M.rel(s, m) >= 0.02 for s in spectra]) for m in BASES.values()]
        c = [np.mean([M.rel(s, m + 0.035) >= 0.02 for s in spectra]) for m in BASES.values()]
        fr = [np.mean([M.rel(s, m) >= 0.01 for s in spectra]) for m in list(M.RNA_FRAGMENTS.values()) + list(M.DNA_FRAGMENTS.values())]
        kind = M.nucleic_acid_type(spectra)[0]
        f.write('\t'.join([a['project'], a['file'], a['group'], str(len(spectra)), '%.4f' % nb] + ['%.4f' % x for x in v + c + fr] + [kind]) + '\n')
        print(a['project'], a['file'][:40], kind, flush=True)
