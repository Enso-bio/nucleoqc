#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calibration of the loss estimated by nucleoqc on all independent files (first validation + extension), with the
project as resampling unit.
 - correlation (Pearson, Spearman) with 95 % CI by project bootstrap
 - least-squares line measured = a + b * estimated with project-bootstrap CI
 - median absolute error, mean bias; clean controls
 - sensitivity: leaving one project out, and without the most influential file
 - replication of the effect: observed/expected ratio per file and per project
-> nucleoqc_calibration_results.txt, calibration_data.tsv"""
import numpy as np
from scipy import stats

tab = lambda f: (lambda L: [dict(zip(L[0], r)) for r in L[1:] if len(r) == len(L[0])])(
    [l.split('\t') for l in open(f, encoding='utf-8').read().split('\n') if l])
F = []
for r in tab('validation_independent.tsv'):
    F.append(dict(r, group='contaminated', source='validation 1'))
for r in tab('validation_extended.tsv'):
    F.append(dict(r, source='extension'))
F = [r for r in F if float(r['id_rate_all']) >= 0.02]                 # valid searches
for r in F:
    r['est'] = 100 * float(r['estimated_loss']); r['meas'] = 100 * float(r['measured_loss'])
projects = sorted({r['project'] for r in F})
lines = []
out = lambda s='': (print(s, flush=True), lines.append(s))
cont = [r for r in F if r['group'] == 'contaminated']; clean = [r for r in F if r['group'] == 'clean']
out('Independent files with a valid search: %d (%d contaminated from %d projects, %d clean); species/proteomes: %d'
    % (len(F), len(cont), len({r['project'] for r in cont}), len(clean), len({r['proteome'] for r in F})))


def fit(R):
    e = np.array([r['est'] for r in R]); m = np.array([r['meas'] for r in R])
    b, a = np.polyfit(e, m, 1)
    return np.corrcoef(e, m)[0, 1], stats.spearmanr(e, m)[0], a, b, np.median(np.abs(e - m)), np.mean(m - e)


def boot(R, n=4000, seed=1):
    rng = np.random.default_rng(seed); P = sorted({r['project'] for r in R}); by_project = {p: [r for r in R if r['project'] == p] for p in P}
    res = []
    for _ in range(n):
        s = [x for p in rng.choice(P, len(P)) for x in by_project[p]]
        if len({r['est'] for r in s}) > 2:
            res.append(fit(s))
    return np.array(res)


for name, R in (('all (contaminated + clean)', F), ('contaminated only', cont)):
    r_, rho, a, b, mae, bias = fit(R); B = boot(R)
    ci = lambda j: (np.nanpercentile(B[:, j], 2.5), np.nanpercentile(B[:, j], 97.5))
    out('\n%s, n = %d files, %d projects:' % (name, len(R), len({r['project'] for r in R})))
    out('   Pearson r %.2f (95%% CI by projects %.2f-%.2f) | Spearman %.2f (%.2f-%.2f)' % (r_, *ci(0), rho, *ci(1)))
    out('   measured = %.2f + %.2f x estimated  (intercept 95%% CI %.2f to %.2f points; slope 95%% CI %.2f-%.2f)' % (a, b, *ci(2), *ci(3)))
    out('   median absolute error %.1f points; mean bias (measured - estimated) %+.1f points (95%% CI %+.1f to %+.1f)' % (mae, bias, *ci(5)))

# sensitivity
out('\nSensitivity (contaminated + clean):')
rs, bs = [], []
for p in projects:
    R = [r for r in F if r['project'] != p]
    x = fit(R); rs.append((x[0], p)); bs.append(x[3])
out('   leaving one project out at a time: r between %.2f (without %s) and %.2f (without %s); slope between %.2f and %.2f'
    % (min(rs)[0], min(rs)[1], max(rs)[0], max(rs)[1], min(bs), max(bs)))
worst = max(F, key=lambda r: abs(r['est'] - r['meas']))
x = fit([r for r in F if r is not worst])
out('   without the file with the largest error (%s, estimated %.1f, measured %.1f): r %.2f, slope %.2f, intercept %.2f'
    % (worst['project'], worst['est'], worst['meas'], x[0], x[3], x[2]))
pre = [r for r in F if float(r['id_rate_all']) >= 0.85]
if pre:
    x = fit([r for r in F if float(r['id_rate_all']) < 0.85])
    out('   without %d files with >=85 %% identified (peak lists probably prefiltered: %s): r %.2f, slope %.2f, intercept %.2f'
        % (len(pre), ', '.join(sorted({r['project'] for r in pre})), x[0], x[3], x[2]))
if clean:
    out('   clean controls (n = %d): maximum estimated %.2f %%, maximum measured %.2f %%' % (len(clean), max(r['est'] for r in clean), max(abs(r['meas']) for r in clean)))

# replication of the effect
E = [r for r in cont if int(r['n_nucleobase_spectra']) >= 20 and r['id_rate_nucleobase'] not in ('nan', '') and r['id_rate_expected_nucleobase'] not in ('nan', '')]
o = np.array([float(r['id_rate_nucleobase']) for r in E]); x = np.array([float(r['id_rate_expected_nucleobase']) for r in E])
out('\nReplication of the effect (contaminated files with >=20 nucleobase spectra, n = %d files, %d projects):'
    % (len(E), len({r['project'] for r in E})))
out('   median observed/expected ratio %.2f; lower in %d of %d; Wilcoxon p = %.2g' % (np.median(o / x), int((o < x).sum()), len(o), stats.wilcoxon(o, x).pvalue))
pp = {}
for r, oo, xx in zip(E, o, x):
    pp.setdefault(r['project'], []).append(oo / xx)
v = np.array([np.mean(q) for q in pp.values()])
out('   per project (mean of its files): median ratio %.2f; lower in %d of %d projects; sign test p = %.2g'
    % (np.median(v), int((v < 1).sum()), len(v), stats.binomtest(int((v < 1).sum()), len(v)).pvalue))
open('nucleoqc_calibration_results.txt', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
with open('calibration_data.tsv', 'w', encoding='utf-8') as f:
    f.write('group\tsource\tproject\tfile\tproteome\testimated\tmeasured\n')
    for r in F:
        f.write('%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\n' % (r['group'], r['source'], r['project'], r['file'], r['proteome'], r['est'], r['meas']))
