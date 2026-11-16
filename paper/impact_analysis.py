#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analysis of impact_spectra.npz: MS2 time lost and quality of peptide spectra (Figure S2)."""
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Z = np.load('impact_spectra.npz')
F = []
for k in Z.files:
    g, p, fn = k.split('|'); a = Z[k]
    a = a[~np.isnan(a[:, 1])]
    F.append({'g': g, 'p': p, 'fn': fn, 'nb': a[:, 0], 'aa': a[:, 1], 'aa0': a[:, 2], 'nuc': a[:, 3], 'nuc0': a[:, 4], 'pmz': a[:, 5]})
clean = [f for f in F if f['g'] == 'clean']; cont = [f for f in F if f['g'] == 'contaminated']

# thresholds from the NULL ladder (masses shifted +0.37 Da): a spectrum is "peptide-like" if its AA ladder exceeds the
# 99th percentile of the null ladder; "oligonucleotide" likewise with the NUC ladder (and no AA ladder)
THR_AA = float(np.percentile(np.concatenate([f['aa0'] for f in F]), 99))
THR_NUC = float(np.percentile(np.concatenate([f['nuc0'] for f in F]), 99))
print('Files: %d contaminated (%d projects), %d clean (%d projects)' % (len(cont), len({f['p'] for f in cont}), len(clean), len({f['p'] for f in clean})))
print('Thresholds (p99 of the null): peptide spectrum if AA ladder >= %.3f; oligonucleotide if NUC ladder >= %.3f and AA below threshold\n' % (THR_AA, THR_NUC))

for f in F:
    nb = f['nb'] >= 0.10; pep = f['aa'] >= THR_AA; oli = (f['nuc'] >= THR_NUC) & ~pep
    f.update(fnb=nb.mean(), fpep=pep.mean(), foli=oli.mean(), fnopep_nb=(nb & ~pep).mean(),
             pep_in_nb=pep[nb].mean() if nb.sum() >= 20 else np.nan,
             aa_pep_nb=np.median(f['aa'][nb & pep]) if (nb & pep).sum() >= 20 else np.nan,
             aa_pep_without=np.median(f['aa'][~nb & pep]))


def summary(name, xs, pct=True):
    xs = np.array([x for x in xs if not np.isnan(x)])
    return '%s median %.1f%s (p25-p75 %.1f-%.1f)' % (name, (100 if pct else 1) * np.median(xs), ' %' if pct else '',
                                                      (100 if pct else 1) * np.percentile(xs, 25), (100 if pct else 1) * np.percentile(xs, 75))


print('1) Fraction of MS2 spectra that look like peptides')
for g, fs in (('clean', clean), ('contaminated', cont)):
    print('   %-12s %s' % (g, summary('', [f['fpep'] for f in fs])))
print('   Mann-Whitney p = %.2g' % stats.mannwhitneyu([f['fpep'] for f in cont], [f['fpep'] for f in clean]).pvalue)
r = stats.spearmanr([f['fnb'] for f in cont], [f['fpep'] for f in cont])
print('   in contaminated files: Spearman(fraction with nucleobase, peptide-like fraction) = %.2f (p %.2g)' % (r.statistic, r.pvalue))

print('\n2) Spectra with nucleobase in contaminated files')
print('   ' + summary('still look like peptides:', [f['pep_in_nb'] for f in cont]))
print('   ' + summary('fraction of all MS2 with nucleobase and NO peptide ladder (time lost):', [f['fnopep_nb'] for f in cont]))
print('   ' + summary('fraction of all with an oligonucleotide ladder (contaminated):', [f['foli'] for f in cont]))
print('   ' + summary('fraction of all with an oligonucleotide ladder (clean):', [f['foli'] for f in clean]))

print('\n3) Quality of peptide spectra with co-isolated nucleobase (AA ladder, median per file)')
pairs = [(f['aa_pep_nb'], f['aa_pep_without']) for f in cont if not np.isnan(f['aa_pep_nb'])]
a = np.array(pairs)
print('   %d files: with nucleobase %.3f vs without %.3f; Wilcoxon p = %.2g; lower in %d of %d'
      % (len(a), np.median(a[:, 0]), np.median(a[:, 1]), stats.wilcoxon(a[:, 0], a[:, 1]).pvalue, int((a[:, 0] < a[:, 1]).sum()), len(a)))


# precursor m/z control: each peptide spectrum with nucleobase is compared with the peptide spectra WITHOUT nucleobase
# of the same file and precursor within +-10 Th
diff = []
for f in cont:
    nb = f['nb'] >= 0.10; pep = f['aa'] >= THR_AA
    ref = np.where(~nb & pep)[0]; o = np.argsort(f['pmz'][ref]); ref, rp = ref[o], f['pmz'][ref][o]
    d = []
    for i in np.where(nb & pep)[0]:
        a_, b_ = np.searchsorted(rp, [f['pmz'][i] - 10, f['pmz'][i] + 10])
        if b_ - a_ >= 5:
            d.append(f['aa'][i] - np.mean(f['aa'][ref[a_:b_]]))
    if len(d) >= 20:
        diff.append((np.mean(d), len(d)))
diff = np.array(diff)
print('   with precursor m/z control (+-10 Th, same file; means): median difference per file %.3f; lower in %d of %d; Wilcoxon p = %.2g'
      % (np.median(diff[:, 0]), int((diff[:, 0] < 0).sum()), len(diff), stats.wilcoxon(diff[:, 0]).pvalue))
# in clean files the same comparison with weak nucleobase (>=2 %) as a "pseudo-group" does not apply: negative control =
# split the spectra of clean files at random and repeat
rng = np.random.default_rng(1); neg = []
for f in clean:
    pep = f['aa'] >= THR_AA; idx = np.where(pep)[0]
    if len(idx) < 200:
        continue
    g = rng.random(len(idx)) < 0.2
    neg.append(np.mean(f['aa'][idx[g]]) - np.mean(f['aa'][idx[~g]]))
print('   negative control (random 20/80 split in clean files, means): median difference %.3f (p25-p75 %.3f-%.3f)'
      % (np.median(neg), np.percentile(neg, 25), np.percentile(neg, 75)))

print('\nPer file (contaminated, most to least):')
for f in sorted(cont, key=lambda f: -f['fnb']):
    print('   %s %-38s nb %5.1f %% | peptide-like %5.1f %% | nb without peptide %5.1f %% | oligo %4.1f %%'
          % (f['p'], f['fn'][:38], 100 * f['fnb'], 100 * f['fpep'], 100 * f['fnopep_nb'], 100 * f['foli']))

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
ax[0].scatter([100 * f['fnb'] for f in clean], [100 * f['fpep'] for f in clean], s=14, c='#9bbbd9', label='clean')
ax[0].scatter([100 * f['fnb'] for f in cont], [100 * f['fpep'] for f in cont], s=14, c='#c0504d', label='contaminated')
ax[0].set_xscale('symlog', linthresh=0.1); ax[0].set_xlabel('% of spectra with nucleobase'); ax[0].set_ylabel('% of spectra with a peptide ladder')
ax[0].legend(fontsize=8); ax[0].set_title('Peptide-like spectra vs contamination', fontsize=10)
ax[1].scatter(a[:, 1], a[:, 0], s=14, c='#3a6ea5'); lim = [a.min() * .9, a.max() * 1.05]; ax[1].plot(lim, lim, 'k--', lw=.7)
ax[1].set_xlabel('amino-acid ladder, peptides without nucleobase'); ax[1].set_ylabel('amino-acid ladder, peptides with nucleobase')
ax[1].set_title('Quality of peptide spectra (per file)', fontsize=10)
plt.tight_layout(); plt.savefig('impact_figure.png', dpi=150); print('\nfigure: impact_figure.png')
