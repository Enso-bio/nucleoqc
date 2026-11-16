#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Main manuscript figures as vector PDF (same data as the figures in the online document).
Usage (from this folder): python figures.py  -> figures/fig2.pdf ... fig5.pdf"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('figures', exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
BLUE, GREY, REF = '#2f6fbf', '#b9b9b0', '#6e6e6e'

# ---------------- Figure 2: protocol factors (multivariable OR, manual annotation)
f = [('Phosphopeptide enrichment', 21, 9.89, 3.22, 30.37, True), ('Lysate sonication', 61, 3.00, 1.25, 7.19, False),
     ('Bacterial sample', 48, 2.78, 1.11, 6.96, False), ('Nuclease declared', 16, 0.28, 0.04, 1.90, False),
     ('Tissue sample', 107, 0.26, 0.09, 0.75, False)]
fig, ax = plt.subplots(figsize=(6.5, 2.6))
for i, (n, k, o, lo, hi, key) in enumerate(f):
    c = BLUE if key else GREY
    ax.plot([lo, hi], [i, i], color=c, lw=2); ax.plot(o, i, 'o', color=c, ms=7)
    ax.text(45, i, '%.2f (%.2f–%.2f)' % (o, lo, hi), va='center', fontsize=7.5)
ax.axvline(1, color=REF, ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlim(0.03, 40); ax.set_xticks([0.05, 0.1, 0.3, 1, 3, 10, 30]); ax.set_xticklabels(['0.05', '0.1', '0.3', '1', '3', '10', '30'])
ax.set_yticks(range(len(f))); ax.set_yticklabels(['%s (n = %d)' % (n, k) for n, k, *_ in f]); ax.invert_yaxis()
ax.set_xlabel('Odds ratio of contamination (log scale)'); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
plt.tight_layout(); plt.savefig('figures/fig2.pdf'); plt.close()

# ---------------- Figure 3: dose-response
b = [('<0.1 %', 34765, 1.00), ('0.1–2 %', 1205, 1.34), ('2–10 %', 3551, 1.02), ('10–30 %', 4325, 0.70), ('≥30 %', 5859, 0.38)]
fig, ax = plt.subplots(figsize=(5.2, 3.0))
x = np.arange(len(b))
ax.bar(x, [r for *_, r in b], width=0.55, color=[GREY, GREY, GREY, BLUE, BLUE])
for i, (l, n, r) in enumerate(b):
    ax.text(i, r + 0.03, '%.2f' % r, ha='center', fontsize=8, fontweight='bold')
ax.axhline(1, color=REF, ls='--', lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(['%s\nn = %s' % (l, format(n, ',')) for l, n, _ in b])
ax.set_ylabel('Identified: observed / expected'); ax.set_xlabel('Adenine or guanine, % of the spectrum base peak'); ax.set_ylim(0, 1.5)
plt.tight_layout(); plt.savefig('figures/fig3.pdf'); plt.close()

# ---------------- Figure 4: robustness
r = [('Sage, neighbours ±10 Th (main)', 0.39, '17/18', False), ('+ same charge', 0.42, '18/18', False),
     ('+ charge + RT ±5 min', 0.40, '15/16', False), ('+ charge + precursor intensity*', 0.81, '10/12', True),
     ('+ charge + RT + intensity* (10 % matched)', 0.90, '11/12', True), ('Comet instead of Sage', 0.58, '15/18', False),
     ('Independent files, 21 projects, 12 species', 0.65, '21/23', False)]
ypos = [0, 1, 2, 3, 4, 5.6, 6.6]
fig, ax = plt.subplots(figsize=(6.5, 3.0))
for y, (l, v, w, intensity) in zip(ypos, r):
    c = BLUE if intensity else GREY
    ax.plot([0, v], [y, y], color=c, lw=1.5, alpha=0.6); ax.plot(v, y, 'o', color=c, ms=7)
    ax.text(v + 0.03, y, '%.2f' % v, va='center', fontsize=7.5, fontweight='bold'); ax.text(1.27, y, w, va='center', ha='right', fontsize=7.5, color='#555')
ax.text(1.27, -0.9, 'files worse', ha='right', fontsize=7, color='#555')
ax.axvline(1, color=REF, ls='--', lw=0.8)
ax.set_yticks(ypos); ax.set_yticklabels([l for l, *_ in r]); ax.invert_yaxis(); ax.set_xlim(0, 1.28); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.text(-0.02, -0.9, 'MATCHING', transform=ax.get_yaxis_transform(), ha='right', fontsize=7, fontweight='bold', color='#555')
ax.text(-0.02, 4.9, 'REPLICATION', transform=ax.get_yaxis_transform(), ha='right', fontsize=7, fontweight='bold', color='#555')
ax.set_xlabel('Identified: observed / expected (median per file)'); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
plt.tight_layout(); plt.savefig('figures/fig4.pdf'); plt.close()

# ---------------- Figure 5: nucleoqc validation
cal = [(53.5, 66.3), (55.2, 68.7), (1.9, 0.3), (13.4, 11.8), (2.7, 5.3), (12.7, 17.1), (13.7, 24.0), (9.0, 8.3), (0.6, 1.2), (46.4, 60.7),
       (27.0, 5.8), (28.2, 28.9), (41.2, 46.5), (0.5, 0.1), (1.7, -2.2), (1.7, 1.5), (0.4, 1.7), (5.1, 12.9), (1.2, 2.0), (45.6, 59.8)]
D = [l.split('\t') for l in open('calibration_data.tsv', encoding='utf-8').read().split('\n')[1:] if l]
ind = [(float(r[5]), float(r[6])) for r in D if r[0] == 'contaminated']
cln = [(float(r[5]), float(r[6])) for r in D if r[0] == 'clean']
both = ind + cln
rc = np.corrcoef(*zip(*cal))[0, 1]; ri = np.corrcoef(*zip(*both))[0, 1]
b, a = np.polyfit(*zip(*both), 1)
fig, ax = plt.subplots(figsize=(3.8, 3.6))
ax.plot([-5, 75], [-5, 75], color=REF, ls='--', lw=0.8)
ax.plot([0, 70], [a, a + 70 * b], color=BLUE, lw=0.8, alpha=0.6)
ax.scatter(*zip(*cal), s=16, color=GREY, label='Calibration files (n = 20), r = %.2f' % rc, zorder=2)
ax.scatter(*zip(*ind), s=24, color=BLUE, label='Independent contaminated (n = %d)' % len(ind), zorder=3)
ax.scatter(*zip(*cln), s=24, facecolors='white', edgecolors=BLUE, label='Independent clean (n = %d)' % len(cln), zorder=4)
ax.text(38, -2, 'independent: r = %.2f\nmeasured = %.2f + %.2f × estimated' % (ri, a, b), fontsize=6.5, color=BLUE)
for p, (e, m), pos in (('PXD012818', (11.78, 0.67), (17, 8)), ('PXD068482', (30.5, 14.1), (33, 22))):
    ax.annotate(p, (e, m), pos, fontsize=6.5, color='#555', arrowprops=dict(arrowstyle='-', lw=0.5, color='#555'))
ax.set_xlim(-5, 75); ax.set_ylim(-5, 75); ax.set_aspect('equal')
ax.set_xlabel('Loss estimated by nucleoqc (%)'); ax.set_ylabel('Loss measured with Sage (%)'); ax.legend(fontsize=6.5, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig('figures/fig5.pdf'); plt.close()
print('figures/fig2.pdf ... fig5.pdf; r calibration %.2f, r independent %.2f' % (rc, ri))
