#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extension of the independent validation of nucleoqc (manuscript: nucleoqc section, Figure 5).

Files from projects never used for calibration or in the first validation: contaminated files from the cross-sectional
extension (nuclease_cross_section.tsv) and from the atlas with an available proteome, plus clean files as negative
controls. Same procedure as validation_independent.py: the loss is estimated BEFORE searching and then measured with Sage
(sharded search, 1 % FDR, within-file comparison).
Usage: python validation_extended.py SAGE_EXE HUMAN_FASTA  -> validation_extended.tsv"""
import sys, os, zlib, pickle, random
SAGE, FASTA_H = sys.argv[1], sys.argv[2]
import validation_independent as V          # reads sys.argv when imported
import numpy as np
import atlas

V.PROTEOMES.update({
    'mouse_sp': 'organism_id:10090 AND reviewed:true',
    'house_fly': 'organism_id:7370 AND keyword:KW-1185',
    'toxoplasma_me49': 'organism_id:508771 AND keyword:KW-1185',
    'pseudomonas_pao1': 'organism_id:208964 AND keyword:KW-1185',
    'msmegmatis_mc2': 'organism_id:246196 AND keyword:KW-1185',
    'rabbit': 'organism_id:9986 AND keyword:KW-1185',
    'pig': 'organism_id:9823 AND keyword:KW-1185',
    'sheep': 'organism_id:9940 AND keyword:KW-1185',
})
# project, exact file name, proteome, phospho, note
CONTAMINATED = [
    ('PXD035527', '0__1hour_PP3_r1.raw.-1.mgf', 'human', True, 'phosphoproteome'),
    ('PXD058484', 'LN231114_pAVIC_MLB2_trypsin_Nterm_bRP5.mzML', 'human', False, 'infected cells (astrovirus)'),
    ('PXD068482', 'QE1_SMMG_EXP03_PP_Sample_05.raw.-1.mgf', 'human', True, 'phosphoproteome'),
    ('PXD073003', 'FemaleHousefly_Phos.mzML', 'house_fly', True, 'phosphoproteome'),
    ('PXD042938', 'LFQ_CD19_6vs6_area_ttest.mzML', 'human', False, ''),
    ('PXD035007', 'XWU_MCF7_EDX_TMTPro_pST_041420_F01-12.mzML', 'human', True, 'TMTpro phosphoproteome'),
    ('PXD044361', 'PKGAID_phosphoproteome_V2_multiconsensus.mzML', 'toxoplasma_me49', True, 'phosphoproteome'),
    ('PXD042268', 'TP_TMT6_Gel.mzML', 'mouse_sp', True, ''),
    ('PXD053240', '20211013_SARS_Cells_TMT16plex_Fosfo_Time_ttes', 'human', True, 'TMT phosphoproteome'),
    ('PXD063068', 'OE_20240808_02_070324-WT1-2v2_Ubenrichpeaklis', 'human', False, 'ubiquitin enrichment'),
    ('PXD042640', 'JCP_051422_AC4_pY_PDX11_erl_ag_18plex.mzML', 'human', True, 'pY TMTpro'),
    ('PXD053061', 'Adipocytes_Phospho_HpHF02.mzML', 'mouse_sp', True, 'phosphoproteome'),
    ('PXD044080', '20230501_AWC_IP_tgAWC028_306920-HA_Analysis_3', 'toxoplasma_me49', False, 'immunoprecipitation'),
    ('PXD054611', 'P775-R100_Hek293T_pcDNA3-1_HA_UT.mzML', 'human', False, ''),
    ('PXD047971', 'SE_515_K-3__5_0ul.mgf', 'saureus_nctc8325', False, ''),
    ('PXD014222', 'qExHF01_06004.mgf', 'human', False, ''),
    ('PXD073126', '250501_ZMM_KMI_sFtsk_4.mzML', 'msmegmatis_mc2', False, ''),
    ('PXD059573', '20230325-GJY-11-C.mgf', 'pseudomonas_pao1', False, ''),
    ('PXD055077', 'SUMOylation_PhoX_TBDSPP_FT_FAIMS_140823.mzML', 'rabbit', False, 'flow-through'),
    ('PXD070403', 'Imp70-1.mzML', 'pig', False, ''),
    ('PXD015689', '20191001Sheep_Fat_bandb.MSMS.mgf', 'sheep', False, 'TripleTOF'),
]
OUT = 'validation_extended.tsv'
tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])


def clean_files(n=10):
    """Clean files (nb10 < 0.1 %) from new projects, human or mouse, high resolution."""
    T = [t for t in tab('nuclease_cross_section.tsv') if t['valid'] == 'True' and float(t['nb10']) < 0.001
         and (t['organism'].startswith('Homo sapiens') or t['organism'].startswith('Mus musculus')) and ';' not in t['organism']]
    rng = random.Random(3); rng.shuffle(T)
    return [(t['project'], t['file'], 'human' if t['organism'].startswith('Homo') else 'mouse_sp', t['phospho_enrichment'] == '1', 'clean control')
            for t in T[:n]]


if __name__ == '__main__':
    done = {tuple(l.split('\t')[:2]) for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
    if not done:
        open(OUT, 'w', encoding='utf-8').write('group\tproject\tfile\tproteome\tnote\tn_ms2\tn_nucleobase_spectra\tnb10\tid_rate_all\t'
                                               'id_rate_nucleobase\tid_rate_expected_nucleobase\testimated_loss\tmeasured_loss\n')
    tasks = [('contaminated',) + t for t in CONTAMINATED] + [('clean',) + t for t in clean_files()]
    print('%d files (%d contaminated, %d clean)' % (len(tasks), len(CONTAMINATED), len(tasks) - len(CONTAMINATED)), flush=True)
    for g, p, pref, prot, phospho, note in tasks:
        try:
            fs = [f for f in atlas.eligible_files(p) if f[0].startswith(pref)]
        except Exception as ex:
            print('%s: error listing files %s' % (p, ex), flush=True); continue
        if not fs:
            print('%s %s: not found' % (p, pref), flush=True); continue
        fn = fs[0][0]
        if (g, p) in done:                       # columns 1-2 of the TSV: group, project
            continue
        base = os.path.join(V.I.WORK, '%s_%d' % (p, zlib.crc32(fn.encode()) % 100000))
        try:
            if os.path.exists(base + '.pkl'):
                spectra = pickle.load(open(base + '.pkl', 'rb'))
            else:
                spectra = V.I.download(p, fn); pickle.dump(spectra, open(base + '.pkl', 'wb'))
        except Exception as ex:
            print('%s %s: download error %s' % (p, fn, ex), flush=True); continue
        if len(spectra) < 200:
            print('%s %s: too few spectra (%d)' % (p, fn, len(spectra)), flush=True); continue
        nb = np.array([max(V.C.rel(s, V.C.BASES['A']), V.C.rel(s, V.C.BASES['G'])) for s in spectra]); pmz = np.array([s['pmz'] for s in spectra])
        est = V.NQ.estimated_loss(nb)                     # the estimate is fixed before searching
        lab, mass = V.I.labelling(spectra)
        if not os.path.exists(base + '.mzML'):
            V.I.write_mzml(spectra, base + '.mzML')
        try:
            ident = V.search(base, len(spectra), phospho, mass, V.proteome(prot))
        except Exception as ex:
            print('%s %s: error %s' % (p, fn, ex), flush=True); continue
        if ident is None:
            continue
        m = nb >= 0.10
        if m.sum() == 0:
            meas, o_, e_ = 0.0, np.nan, np.nan
        else:
            meas, o_, e_ = V.measured_loss(nb, pmz, ident)
        with open(OUT, 'a', encoding='utf-8') as f:
            f.write('%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (g, p, fn, prot, note, len(spectra), m.sum(), m.mean(),
                    ident.mean(), o_, e_, est, meas))
        print('%-12s %s %-36s %-16s nb10 %5.1f %% | ID %5.1f %% | estimated %5.1f %% measured %5.1f %%'
              % (g, p, fn[:36], prot, 100 * m.mean(), 100 * ident.mean(), 100 * est, 100 * meas), flush=True)
    print('DONE')
