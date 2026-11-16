#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Cross-sectional extension: is declared nuclease treatment associated with less contamination?

PRIDE projects found by nuclease or sonication keywords (nuclease_candidates.json, 507 with MGF/mzML) that are not in the
atlas (447). The sample protocol is annotated with the same rules as the manual annotation and one file per project is
measured with the atlas metric. Studies with nucleic acid by design are excluded.
Usage: python nuclease_cross_section.py [analyse]  -> nuclease_cross_section.tsv, nuclease_cross_section_results.txt"""
import sys, os, re, json, random, math, collections, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import stats
import nucleic_census as C
import atlas

OUT = 'nuclease_cross_section.tsv'
NUC = r'benzonase|nuclease(?! inhibitor)|dnase|rnase(?! ?(inhibitor|out|in\b))|micrococcal|turbonuclease|denarase'
SON = r'sonicat|sonif|bioruptor|ultrason|covaris|tip sonic'
# broader than atlas.BY_DESIGN: these projects were found by nuclease keywords, so studies of nucleic-acid binding proteins,
# ribosomes, chromatin and RNPs are excluded as well
BY_DESIGN = (r'chip-ms|chip-mass|\brime\b|\bchip\b|nucleosome|adp-ribos|protein-rna|rna-protein|protein-dna|rna[- ]binding|'
             r'uv[- ]?(induced |inducible )?cross-?link|rna interactome|interactome capture|\b(par-|i|e)?clip\b|oligonucleotide|aptamer|'
             r'limited proteolysis|lip-ms|in vitro.{0,40}\batp\b|ribosom|polysom|chromatin|histone|nucleosom|rnp\b|ribonucleoprotein')
PHOS = r'phospho(peptide|-?enrich)|tio2|titanium dioxide|\bimac\b|fe-nta|fe\(iii\)-?nta|zr-imac|ti-imac|polymac'


def annotate(acc):
    j = C.get_json(C.API + 'projects/' + acc)
    if not isinstance(j, dict) or not j:
        return None
    m = str(j.get('sampleProcessingProtocol', '')); alltext = ' '.join(str(j.get(k, '')) for k in ('title', 'projectDescription', 'sampleProcessingProtocol', 'keywords')).lower()
    ins = atlas.as_text(j.get('instruments'))
    # nuclease applied to the material: negated mentions or mentions as a reagent for something else are discarded
    nuc = any(not re.search(r'(without|no|omit|not)\W+(\w+\W+){0,3}$', m[max(0, x.start() - 40):x.start()], re.I)
              for x in re.finditer(NUC, m, re.I))
    son = bool(re.search(SON, m, re.I)) and not re.search(r'(peptides?|gel)[^.]{0,60}sonicat', m, re.I)
    return {'project': acc, 'year': str(j.get('submissionDate', ''))[:4], 'instrument': ins[:60], 'organism': atlas.as_text(j.get('organisms'))[:50],
            'nuclease': int(nuc), 'sonication': int(son), 'phospho_enrichment': int(bool(re.search(PHOS, m, re.I))),
            'nucleic_acid_by_design': int(bool(re.search(BY_DESIGN, alltext))), 'tissue': int(bool(re.search(atlas.FLAGS['tissue'], alltext))),
            'bacteria': int(bool(re.search(atlas.FLAGS['bacteria'], alltext))), 'ip_ap_ms': int(bool(re.search(atlas.FLAGS['ip_ap_ms'], alltext))),
            'high_res': int(bool(re.search(atlas.INSTR, ins, re.I)))}


COLS = ['project', 'year', 'instrument', 'organism', 'nuclease', 'sonication', 'phospho_enrichment', 'nucleic_acid_by_design', 'tissue',
        'bacteria', 'ip_ap_ms', 'file', 'n_ms2', 'valid', 'reason', 'nb10', 'nb2']

if __name__ == '__main__' and len(sys.argv) < 2:
    cand = json.load(open('nuclease_candidates.json'))['with_peak_lists']
    in_atlas = {l.split('\t')[0] for l in open('atlas.tsv', encoding='utf-8').read().split('\n')[1:] if l}
    done = {l.split('\t')[0] for l in open(OUT, encoding='utf-8').read().split('\n')[1:] if l} if os.path.exists(OUT) else set()
    if not done:
        open(OUT, 'w', encoding='utf-8').write('\t'.join(COLS) + '\n')
    todo = [a for a in cand if a not in in_atlas and a not in done]
    print('%d projects to annotate and measure' % len(todo), flush=True)
    lock = threading.Lock()

    def one(acc):
        try:
            a = annotate(acc)
            if not a or not a['high_res'] or a['nucleic_acid_by_design'] or not (a['nuclease'] or a['sonication']):
                r = dict(a or {'project': acc}, file='', n_ms2=0, valid='skipped', reason='not eligible', nb10='', nb2='')
            else:
                fs = atlas.eligible_files(acc)
                if not fs:
                    r = dict(a, file='', n_ms2=0, valid='skipped', reason='no eligible file', nb10='', nb2='')
                else:
                    n, t, fmt, url = random.Random(acc).choice(fs)
                    m = atlas.measure(url, t, fmt)
                    r = dict(a, file=n, n_ms2=m.get('n_ms2', 0), valid=str(m.get('valid')), reason=m.get('reason', ''),
                             nb10=('%.4f' % m['nb10']) if 'nb10' in m else '', nb2=('%.4f' % m['nb2']) if 'nb2' in m else '')
        except Exception as ex:
            r = {'project': acc, 'valid': 'error', 'reason': str(ex)[:60]}
        with lock:
            open(OUT, 'a', encoding='utf-8').write('\t'.join(str(r.get(k, '')).replace('\t', ' ') for k in COLS) + '\n')
            if r.get('valid') not in ('skipped',):
                print('%s nuc=%s son=%s %s nb10=%s' % (acc, r.get('nuclease'), r.get('sonication'), r.get('reason') or 'ok', r.get('nb10')), flush=True)
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(one, todo))
    print('DONE')

if __name__ == '__main__' and len(sys.argv) >= 2:
    L = [l.split('\t') for l in open(OUT, encoding='utf-8').read().split('\n') if l]
    R = [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])]
    V = [r for r in R if r['valid'] == 'True']
    lines = []
    out = lambda s='': (print(s, flush=True), lines.append(s))
    out('Extension: %d projects measured and valid (of %d annotated)' % (len(V), len(R)))
    for name, sel in (('nuclease declared', lambda r: r['nuclease'] == '1'), ('sonication without nuclease', lambda r: r['sonication'] == '1' and r['nuclease'] != '1'),
                      ('sonication with nuclease', lambda r: r['sonication'] == '1' and r['nuclease'] == '1')):
        s = [float(r['nb10']) for r in V if sel(r)]
        if s:
            s = np.array(s)
            out('   %-27s n=%3d | nb10 > 1 %%: %2d (%.1f %%) | > 5 %%: %2d | median %.3f %%' % (name, len(s), (s > 0.01).sum(), 100 * (s > 0.01).mean(), (s > 0.05).sum(), 100 * np.median(s)))
    y = np.array([float(r['nb10']) > 0.01 for r in V], float)
    FEAT = ('nuclease', 'sonication', 'phospho_enrichment', 'tissue', 'bacteria', 'ip_ap_ms')
    X = np.array([[float(r[k] == '1') for k in FEAT] + [float(r['year'] >= '2018')] for r in V])
    a = int(((X[:, 0] == 1) & (y == 1)).sum()); b = int(((X[:, 0] == 1) & (y == 0)).sum()); c = int(((X[:, 0] == 0) & (y == 1)).sum()); d = int(((X[:, 0] == 0) & (y == 0)).sum())
    out('   Fisher nuclease vs not (among projects with sonication or nuclease): %d/%d vs %d/%d, OR %.2f, p = %.3g' % (a, a + b, c, c + d, *stats.fisher_exact([[a, b], [c, d]])))
    from scipy.optimize import minimize
    Xc = np.column_stack([np.ones(len(y)), X]); pen = np.r_[0, np.full(X.shape[1], 0.5)]
    f = lambda bb: -(y * (Xc @ bb) - np.logaddexp(0, Xc @ bb)).sum() + 0.5 * (pen * bb * bb).sum()
    bb = minimize(f, np.zeros(Xc.shape[1]), method='BFGS').x
    pr = 1 / (1 + np.exp(-Xc @ bb)); H = Xc.T @ (Xc * (pr * (1 - pr))[:, None]) + np.diag(pen); se = np.sqrt(np.diag(np.linalg.inv(H)))
    out('   multivariable: ' + ', '.join('%s OR %.2f (95%% CI %.2f-%.2f; p %.3f)' % (k, math.exp(bb[j + 1]), math.exp(bb[j + 1] - 1.96 * se[j + 1]), math.exp(bb[j + 1] + 1.96 * se[j + 1]),
                                         2 * stats.norm.sf(abs(bb[j + 1] / se[j + 1]))) for j, k in enumerate(FEAT)))
    open('nuclease_cross_section_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
