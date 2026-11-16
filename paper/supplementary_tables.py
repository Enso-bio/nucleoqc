#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Supplementary tables as one Excel workbook (one sheet per table) plus a README sheet.
-> Supplementary_Tables.xlsx"""
import json, re, collections
import pandas as pd

rd = lambda f: pd.read_csv(f, sep='\t', dtype=str, on_bad_lines='skip')

# display names used in the workbook
PROTEOME = {'mouse_sp': 'mouse (Swiss-Prot)', 'neurospora': 'N. crassa', 'bcereus_atcc14579': 'B. cereus ATCC 14579', 'candida': 'C. albicans',
            'saureus_nctc8325': 'S. aureus NCTC 8325', 'salmonella_lt2': 'S. enterica LT2', 'house_fly': 'house fly',
            'toxoplasma_me49': 'T. gondii ME49', 'pseudomonas_pao1': 'P. aeruginosa PAO1', 'msmegmatis_mc2': 'M. smegmatis MC2 155'}
VARIANT = {'peaks_removed': 'clean'}

sheets = collections.OrderedDict()
a = rd('atlas.tsv')
sheets['S1_atlas_files'] = (a, 'Atlas: every PRIDE file measured (random sample of projects). nb10 / nb2 = fraction of MS2 spectra with adenine or guanine >=10 % / >=2 % of base peak; valid = passes coverage and resolution tests; protocol flags from text rules (1 = mentioned).')
m = json.load(open('annotation_manual.json', encoding='utf-8'))
prot = json.load(open('annotation_protocols.json', encoding='utf-8'))
SON = r'sonicat|sonif|bioruptor|ultrason|covaris|probe tip|tip sonic'
proj = sorted(set(a.loc[a['valid'] == 'True', 'project']))
son = [int(bool(re.search(SON, prot.get(p, {}).get('sample', ''), re.I)) and p not in m['sonication_false_positives']) for p in proj]
sheets['S2_manual_annotation'] = (pd.DataFrame({'project': proj, 'sonication_lysate': son,
    'nuclease': [int(p in m['nuclease_confirmed']) for p in proj], 'phospho_enrichment': [int(p in m['phospho_confirmed']) for p in proj]}),
    'Sample-processing protocol annotation (valid atlas projects): rule match on sampleProcessingProtocol followed by manual review of every matched snippet.')
sheets['S3_sage_files'] = (rd('identifications_1p.tsv'), 'Peptide identification with Sage (full human proteome, sharded search, 1 % FDR) per file; id_rate_expected_nucleobase = expected rate from neighbouring spectra without nucleobase (+-10 Th).')
sheets['S4_comet'] = (pd.read_csv('comet_ids.tsv', sep='\t').groupby(['project', 'file'])['id_comet'].mean().reset_index().rename(columns={'id_comet': 'comet_id_rate'}),
                      'Comet 2026.03.0 identification rate per file (same spectra, 1 % FDR by E-value).')
sheets['S5_rnase_natural_exp'] = (rd('nuclease_experiment.tsv'), 'Natural experiment: files from projects with and without RNase, and nuclear vs cytoplasmic fractions.')
sheets['S6_nuclease_cross_section'] = (rd('nuclease_cross_section.tsv'), 'Extension: PRIDE projects found by nuclease/sonication keywords, annotated automatically and measured (one file per project).')
r = rd('rescue.tsv'); r['variant'] = r['variant'].map(lambda v: VARIANT.get(v, v))
sheets['S7_rescue'] = (r, 'Software rescue: identification rates with nucleic-acid peaks removed, chimeric search, both, and wide window.')
v1 = rd('validation_independent.tsv'); v1.insert(0, 'group', 'contaminated'); v1.insert(0, 'set', 'validation 1')
v2 = rd('validation_extended.tsv'); v2.insert(0, 'set', 'extension')
v = pd.concat([v1, v2], ignore_index=True); v['proteome'] = v['proteome'].map(lambda x: PROTEOME.get(x, x))
sheets['S8_independent_valid'] = (v, 'Independent validation of nucleoqc: 44 files from projects never used for calibration (34 contaminated, 10 clean controls); loss estimated before searching, then measured with Sage. Files with id_rate_all < 0.02 are failed searches and were excluded from the analysis.')
sheets['S9_representativeness'] = (rd('representativeness.tsv'), 'New random PRIDE sample: projects with peak-list files vs raw-only projects, protocol flags by the same text rules.')
sheets['S10_ion_panel_RNA_DNA'] = (rd('table_s10.tsv'), 'Extended ion panel in the 38 files searched with Sage: fraction of MS2 spectra with adenine, guanine, cytosine, uracil, thymine or ribose ion >=2 % of base peak, the same at control masses shifted +0.035 Da, and nucleotide-monophosphate fragments (-H2O, >=1 %) specific to RNA (GMP, CMP, UMP) or DNA (dAMP, dCMP, dTMP), with the nucleoqc call.')
with pd.ExcelWriter('Supplementary_Tables.xlsx') as w:
    pd.DataFrame({'sheet': list(sheets), 'description': [v[1] for v in sheets.values()], 'rows': [len(v[0]) for v in sheets.values()]}).to_excel(w, sheet_name='README', index=False)
    for k, (df, _) in sheets.items():
        df.to_excel(w, sheet_name=k[:31], index=False)
print('\n'.join('%-28s %6d rows' % (k, len(v[0])) for k, v in sheets.items()))
