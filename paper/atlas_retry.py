#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Re-measures atlas.tsv rows that failed because of the network ('too few spectra' or 'error')."""
from concurrent.futures import ThreadPoolExecutor
import nucleic_census as C
import atlas

rows = [l.split('\t') for l in open(atlas.OUTPUT, encoding='utf-8').read().split('\n') if l]
header, rows = rows[0], rows[1:]
i = {k: j for j, k in enumerate(header)}
print('%d malformed rows (left as they are)' % sum(len(f) != len(header) for f in rows))
failed = [f for f in rows if len(f) == len(header) and (f[i['reason']] == 'too few spectra' or f[i['reason']].startswith('error'))]
print('%d rows to retry' % len(failed), flush=True)


def one(f):
    fs = [x for x in atlas.eligible_files(f[i['project']]) if x[0] == f[i['file']]]
    if not fs:
        return f, None
    n, t, fmt, url = fs[0]
    try:
        return f, atlas.measure(url, t, fmt)
    except Exception as ex:
        return f, {'n_ms2': 0, 'valid': False, 'reason': 'error: %s' % str(ex)[:50]}


with ThreadPoolExecutor(4) as ex:
    for f, m in ex.map(one, failed):
        if m is None:
            continue
        for k in atlas.COLS_MEAS[3:]:
            v = m.get(k, '')
            f[i[k]] = ('%.4f' % v) if isinstance(v, float) else str(v)
        print(f[i['project']], f[i['file']][:40], m.get('reason') or 'ok', m.get('nb10', '-'), 'MB', m.get('mb', 0), flush=True)
with open(atlas.OUTPUT, 'w', encoding='utf-8') as out:
    out.write('\t'.join(header) + '\n')
    for f in rows:
        out.write('\t'.join(f) + '\n')
