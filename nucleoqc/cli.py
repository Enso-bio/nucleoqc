# -*- coding: utf-8 -*-
"""nucleoqc: quality control for co-isolated nucleic acid in proteomics MS2 spectra.

Examples
  nucleoqc run1.mgf run2.mzML                    local files (reads three 4 MB chunks of each)
  nucleoqc folder/ --whole                       every MGF/mzML in the folder, read entirely
  nucleoqc PXD023674 --max-files 3               public files of a PRIDE project
  nucleoqc *.mgf --tsv qc.tsv --figure qc.png    results table and figure
"""
import argparse
import glob
import json
import os
import sys
import numpy as np
from . import __version__, reading, metrics

CITATION = ('Torres Alegre EO, Mora Jimenez DE. Co-isolated nucleic acid in public proteomics data: prevalence, '
            'association with peptide identification, and a quality-control tool. 2026.')

COLS = ['source', 'file', 'verdict', 'nb10', 'nb30', 'nb2', 'estimated_loss', 'pride_percentile', 'n', 'valid',
        'reason', 'type', 'rna_fragments', 'dna_fragments', 'ppm', 'cov136', 'hires', 'A2', 'G2', 'C2', 'U2', 'ribose2']


def _targets(inputs, max_files):
    """Yield (source, name, function returning the chunks, format)."""
    for e in inputs:
        if e.upper().startswith('PXD') and not os.path.exists(e):
            from . import pride
            fs = pride.files(e.upper())
            if not fs:
                print('%s: no uncompressed MGF/mzML files in PRIDE' % e, file=sys.stderr)
            if max_files and len(fs) > max_files:      # spread over the list, not just the first ones
                fs = [fs[int(i)] for i in np.linspace(0, len(fs) - 1, max_files)]
            for n, t, fmt, url in fs:
                yield e.upper(), n, (lambda url=url, t=t: pride.read_remote(url, t)), fmt
        elif os.path.isdir(e):
            paths = sorted(glob.glob(os.path.join(e, '*.mgf')) + glob.glob(os.path.join(e, '*.mzML'))
                           + glob.glob(os.path.join(e, '*.mzml')))
            for p in dict.fromkeys(paths):
                yield 'local', p, None, reading.file_format(p)
        else:
            for p in glob.glob(e) or [e]:
                yield 'local', p, None, reading.file_format(p)


def evaluate(source, name, load, fmt, whole=False):
    try:
        chunks = load() if load else reading.read_local(name, whole=whole)[0]
        r = metrics.measure(reading.parse(chunks, fmt))
    except reading.NumpressError as ex:
        r = {'n': 0, 'valid': False, 'reason': str(ex), 'verdict': 'NOT EVALUABLE'}
    except (OSError, ValueError, ConnectionError) as ex:
        r = {'n': 0, 'valid': False, 'reason': 'error: %s' % ex, 'verdict': 'NOT EVALUABLE'}
    r['source'] = source; r['file'] = os.path.basename(name) if source == 'local' else name
    return r


def _pc(x):
    return '%5.1f %%' % (100 * x) if isinstance(x, float) else '   -   '


def line(r):
    if not r['valid']:
        return '%-14s %-45s %s' % (r['verdict'], r['file'][:45], r['reason'])
    kind = '' if r['verdict'] == 'CLEAN' else '  type (tentative): %s' % r['type']
    return ('%-14s %-45s nb10 %s  nb2 %s  estimated loss ~%s  PRIDE percentile %3.0f%s'
            % (r['verdict'], r['file'][:45], _pc(r['nb10']), _pc(r['nb2']), _pc(r['estimated_loss']).strip(),
               100 * r['pride_percentile'], kind))


def write_tsv(results, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(COLS) + '\n')
        for r in results:
            f.write('\t'.join(('%.4f' % r[k]) if isinstance(r.get(k), float) else str(r.get(k, '')).replace('\t', ' ')
                              for k in COLS) + '\n')


def figure(results, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    v = [r for r in results if r['valid']]
    if not v:
        print('figure: no valid file', file=sys.stderr); return
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.2 + 0.18 * len(v)), gridspec_kw={'width_ratios': [1.1, 1]})
    y = np.arange(len(v))
    col = {'SEVERE': '#b2302c', 'MODERATE': '#e08b2c', 'MILD': '#d9c24a', 'CLEAN': '#5b8f5b'}
    ax[0].barh(y, [100 * r['nb10'] for r in v], color=[col[r['verdict']] for r in v])
    ax[0].set_yticks(y); ax[0].set_yticklabels([r['file'][:38] for r in v], fontsize=7); ax[0].invert_yaxis()
    for u, _ in metrics.THRESHOLDS:
        ax[0].axvline(100 * u, color='0.6', lw=0.6, ls='--')
    ax[0].set_xlabel('% of MS2 with adenine/guanine >= 10 % of base peak (nb10)')
    ax[0].set_xscale('symlog', linthresh=1)
    total = np.concatenate([r['nb_per_spectrum'] for r in v])
    b = np.r_[0, np.logspace(-3, 0, 31)]
    for r in v[:12]:
        h, _ = np.histogram(r['nb_per_spectrum'], bins=b)
        ax[1].step(b[1:], h / max(h.sum(), 1), where='pre', lw=1, label=r['file'][:24])
    ax[1].set_xscale('log'); ax[1].set_xlabel('nucleobase / base peak, per spectrum'); ax[1].set_ylabel('fraction of MS2 spectra')
    ax[1].axvline(0.10, color='k', lw=0.6, ls='--')
    if len(v) <= 12:
        ax[1].legend(fontsize=6)
    fig.suptitle('nucleoqc %s - %d files (%d spectra)' % (__version__, len(v), len(total)), fontsize=9)
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='nucleoqc', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='*', help='.mgf/.mzML files, folders, wildcards or PRIDE PXD accessions')
    ap.add_argument('--cite', action='store_true', help='print how to cite nucleoqc and exit')
    ap.add_argument('--gui', action='store_true', help='open the desktop application')
    ap.add_argument('--whole', action='store_true', help='read local files entirely (default: three 4 MB chunks)')
    ap.add_argument('--max-files', type=int, default=5, help='files per PRIDE project (0 = all; default 5)')
    ap.add_argument('--tsv', help='save the results table')
    ap.add_argument('--json', help='save the results as JSON')
    ap.add_argument('--figure', help='save a PNG figure')
    ap.add_argument('--version', action='version', version='nucleoqc ' + __version__)
    a = ap.parse_args(argv)
    if a.cite:
        print('If you use nucleoqc, please cite:\n  ' + CITATION); return 0
    if a.gui:
        from . import gui
        return gui.main()
    if not a.inputs:
        ap.error('give at least one file, folder or PRIDE accession (or --cite)')

    results = []
    for source, name, load, fmt in _targets(a.inputs, a.max_files):
        r = evaluate(source, name, load, fmt, a.whole)
        results.append(r); print(line(r), flush=True)
    if not results:
        print('no file to evaluate', file=sys.stderr); return 2

    v = [r for r in results if r['valid']]
    worst = max(v, key=lambda r: r['nb10']) if v else None
    print('\n%d files, %d evaluable.' % (len(results), len(v)))
    if worst:
        print('Worst: %s (%s). %s' % (worst['file'], worst['verdict'], metrics.ADVICE[worst['verdict']]))
        print('The estimated loss applies the dose-response measured with Sage in PRIDE files; it is approximate.')
    if a.tsv:
        write_tsv(results, a.tsv); print('table: %s' % a.tsv)
    if a.json:
        clean = [{k: v_ for k, v_ in r.items() if not isinstance(v_, np.ndarray)} for r in results]
        json.dump(clean, open(a.json, 'w', encoding='utf-8'), indent=1); print('json: %s' % a.json)
    if a.figure:
        figure(results, a.figure); print('figure: %s' % a.figure)
    print('\nIf you use nucleoqc, please cite: ' + CITATION)
    return 1 if any(r.get('verdict') == 'SEVERE' for r in results) else 0


if __name__ == '__main__':
    sys.exit(main())
