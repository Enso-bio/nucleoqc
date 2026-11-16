#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Natural experiment, step 1: PRIDE projects that compare, within the same project, samples treated and untreated with
nuclease (or different lyses). Keyword search and file classification by name -> nuclease_candidates.json (reviewed by
hand before measuring)."""
import re, json, collections
from concurrent.futures import ThreadPoolExecutor
import nucleic_census as C
import atlas

KEYWORDS = ['benzonase', 'nuclease', 'DNase', 'RNase', 'micrococcal', 'sonication', 'universal nuclease', 'turbonuclease']
# in the file name: treated with nuclease / without nuclease
WITH = re.compile(r'benzo|bzn|benz\b|nuclease|dnase|rnase|mnase|(?<![a-z])nuc(?![a-z])|turbonuc', re.I)
WITHOUT = re.compile(r'(no|wo|w-o|without|minus|ctrl|untreated|non)[-_ ]?(benzo|bzn|benz|nuclease|dnase|rnase|mnase|nuc)', re.I)


def search(kw):
    acc = set()
    for page in range(20):
        r = C.get_json(C.API + 'search/projects?keyword=%s&pageSize=100&page=%d' % (kw.replace(' ', '%20'), page))
        if not isinstance(r, list) or not r:
            break
        acc.update(x.get('accession') for x in r if x.get('accession'))
        if len(r) < 100:
            break
    return acc


def classify(acc):
    try:
        fs = atlas.eligible_files(acc)
    except Exception:
        return None
    if not fs:
        return None
    g = collections.defaultdict(list)
    for f in fs:
        n = f[0]
        g['without' if WITHOUT.search(n) else ('with' if WITH.search(n) else 'other')].append(f)
    return acc, {k: v for k, v in g.items()}


if __name__ == '__main__':
    found = set()
    for kw in KEYWORDS:
        a = search(kw); print('%-20s %d projects' % (kw, len(a)), flush=True); found |= a
    print('total %d projects' % len(found), flush=True)
    with ThreadPoolExecutor(10) as ex:
        res = [r for r in ex.map(classify, sorted(found)) if r]
    print('%d with MGF/mzML' % len(res))
    cand = {}
    for acc, g in res:
        if g.get('with') and (g.get('without') or g.get('other')):
            cand[acc] = {k: [(f[0], f[1], f[2], f[3]) for f in v] for k, v in g.items()}
            print('%s  with=%d without=%d other=%d  e.g. with: %s | e.g. other: %s' % (acc, len(g.get('with', [])), len(g.get('without', [])),
                  len(g.get('other', [])), g['with'][0][0][:45], (g.get('without') or g.get('other'))[0][0][:45]))
    json.dump({'searched': sorted(found), 'with_peak_lists': [r[0] for r in res], 'candidates': cand}, open('nuclease_candidates.json', 'w'), indent=1)
    print('%d candidates with files marked as nuclease-treated and others unmarked' % len(cand))
