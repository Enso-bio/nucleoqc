# -*- coding: utf-8 -*-
"""Metrics of co-isolated nucleic acid in MS2 spectra.

nb10 : fraction of MS2 spectra in which protonated adenine (136.06177) or guanine (152.05669) reaches at least
       10 % of the base peak (+-0.005 Th after recalibration).
nb2  : the same with a 2 % threshold.
Validity: the file must record low masses (>=50 % of MS2 spectra start below m/z 136) and have high-resolution
MS2 (>=70 % of reference ions within +-0.005 Th of their exact mass).

The estimated identification loss uses the dose-response measured with Sage in 18 PRIDE files (observed /
expected identification against control spectra, nucleobase < 0.1 % of base peak, in the same file and 20 Th
precursor window): nucleobase 10-30 % of
base peak -> 0.70; >=30 % -> 0.38; below 10 % no appreciable loss. It assumes that, without contamination, all
intensity bands would be identified at the same rate as the rest of the file, and it uses unadjusted ratios.
"""
import numpy as np

PROTON = 1.007276466
ADENINE, GUANINE, CYTOSINE, URACIL = 136.06177, 152.05669, 112.05054, 113.03456
RIBOSE = 97.02841            # dehydrated ribose ion (C5H5O2+)
TOL = 0.005
REFS_CAL = np.array([70.06513, 72.08078, 86.09643, 110.07127, 120.08078, 129.10224, 136.07569, 147.11280, 175.11895])
REFS_RES = np.array([86.09643, 110.07127, 120.08078, 147.11280, 175.11895])   # immonium L/I, H, F; y1 K, R
MIN_SPECTRA = 200

# Nucleotide fragments (monophosphate - H2O, protonated) specific to RNA or DNA. AMP-H2O (330.060) is not used:
# it has the same formula as dGMP-H2O. Free pyrimidine bases are rarely seen as protonated ions, but their
# nucleotides are.
RNA_FRAGMENTS = {'GMP-H2O': 346.05472, 'CMP-H2O': 306.04857, 'UMP-H2O': 307.03259}
DNA_FRAGMENTS = {'dAMP-H2O': 314.06488, 'dCMP-H2O': 290.05366, 'dTMP-H2O': 305.05332}
FRAGMENT_THRESHOLD = 0.01    # fragment >= 1 % of base peak
MIN_TYPE = 0.005             # at least 0.5 % of MS2 spectra with a fragment to call a type

# dose-response: lower bound of relative intensity -> observed/expected ratio
DOSE = ((0.30, 0.38), (0.10, 0.70))

# verdict from nb10. The atlas calls a project contaminated above 1 %; above ~20 % the deficits measured with
# Sage were 12-69 % of identifications (unadjusted).
THRESHOLDS = ((0.20, 'SEVERE'), (0.05, 'MODERATE'), (0.01, 'MILD'))


def pride_percentile(nb10):
    """Fraction of PRIDE atlas projects (random sample, n = 406) whose maximum nb10 is below nb10."""
    from .reference import NB10_PRIDE
    return float(np.mean(np.asarray(NB10_PRIDE) < nb10))


def recalibrate(spectra):
    """Global ppm offset from ubiquitous immonium and y1 ions; corrects the spectra in place."""
    err = []
    for s in spectra:
        m = s['mz'][s['mz'] < 380]
        if len(m) > 1:
            j = np.clip(np.searchsorted(m, REFS_CAL), 1, len(m) - 1)
            for jj in (j - 1, j):
                e = (m[jj] - REFS_CAL) / REFS_CAL * 1e6
                err.extend(e[np.abs(e) < 20])
    if len(err) < 200:
        return 0.0
    ppm = float(np.median(err))
    f = 1 - ppm * 1e-6
    for s in spectra:
        s['mz'] = s['mz'] * f; s['pmz'] *= f
    return ppm


def rel(s, m, tol=TOL):
    """Intensity of the most intense peak within +-tol of m, relative to the spectrum base peak."""
    a = int(np.searchsorted(s['mz'], m - tol, side='left'))
    b = int(np.searchsorted(s['mz'], m + tol, side='right'))
    return float(s['it'][a:b].max()) / s['it'].max() if b > a else 0.0


def resolution(spectra):
    """Fraction of reference ions (most intense peak within +-0.2 Th) that fall within +-0.005 Th of their mass."""
    ok = tot = 0
    for s in spectra:
        for m in REFS_RES:
            a, b = np.searchsorted(s['mz'], [m - 0.2, m + 0.2])
            if b > a:
                k = a + int(np.argmax(s['it'][a:b]))
                tot += 1; ok += abs(s['mz'][k] - m) <= TOL
    return ok / tot if tot >= 100 else 0.0


def nucleobase_intensity(spectra):
    """Per spectrum: max(adenine, guanine) relative to the base peak."""
    return np.array([max(rel(s, ADENINE), rel(s, GUANINE)) for s in spectra])


def estimated_loss(nb):
    """Approximate fraction of identifications lost, from the dose-response (see module docstring)."""
    nb = np.asarray(nb); f = 0.0; upper = np.inf
    for low, ratio in DOSE:
        f += np.mean((nb >= low) & (nb < upper)) * (1 - ratio)
        upper = low
    return float(f)


def nucleic_acid_type(spectra):
    """'RNA', 'DNA', 'mixed' or 'undetermined' from the fraction of MS2 spectra with type-specific fragments."""
    rna = max(np.mean([rel(s, m) >= FRAGMENT_THRESHOLD for s in spectra]) for m in RNA_FRAGMENTS.values())
    dna = max(np.mean([rel(s, m) >= FRAGMENT_THRESHOLD for s in spectra]) for m in DNA_FRAGMENTS.values())
    if max(rna, dna) < MIN_TYPE:
        t = 'undetermined'
    elif rna >= 2 * dna:
        t = 'RNA'
    elif dna >= 2 * rna:
        t = 'DNA'
    else:
        t = 'mixed'
    return t, float(rna), float(dna)


def verdict(nb10):
    for u, name in THRESHOLDS:
        if nb10 >= u:
            return name
    return 'CLEAN'


def measure(spectra):
    """Return a dict with all metrics for a set of MS2 spectra."""
    r = {'n': len(spectra)}
    if len(spectra) < MIN_SPECTRA:
        r.update(valid=False, reason='too few MS2 spectra (%d < %d)' % (len(spectra), MIN_SPECTRA), verdict='NOT EVALUABLE')
        return r
    r['ppm'] = recalibrate(spectra)
    r['cov136'] = float(np.mean([s['mz'][0] < 136.0 for s in spectra]))
    r['hires'] = resolution(spectra)
    nb = nucleobase_intensity(spectra)
    A = np.array([rel(s, ADENINE) for s in spectra]); G = np.array([rel(s, GUANINE) for s in spectra])
    C = np.array([rel(s, CYTOSINE) for s in spectra]); U = np.array([rel(s, URACIL) for s in spectra])
    R = np.array([rel(s, RIBOSE) for s in spectra])
    r.update(nb10=float(np.mean(nb >= 0.10)), nb30=float(np.mean(nb >= 0.30)), nb2=float(np.mean(nb >= 0.02)),
             A2=float(np.mean(A >= 0.02)), G2=float(np.mean(G >= 0.02)), C2=float(np.mean(C >= 0.02)),
             U2=float(np.mean(U >= 0.02)), ribose2=float(np.mean(R >= 0.02)))
    r['estimated_loss'] = estimated_loss(nb)
    r['valid'] = bool(r['cov136'] >= 0.5 and r['hires'] >= 0.7)
    if not r['valid']:
        r['reason'] = 'low-resolution MS2 (ion trap?)' if r['hires'] < 0.7 else 'MS2 spectra do not record m/z < 136'
        r['verdict'] = 'NOT EVALUABLE'
    else:
        r['reason'] = ''
        r['verdict'] = verdict(r['nb10'])
        r['pride_percentile'] = pride_percentile(r['nb10'])
        r['type'], r['rna_fragments'], r['dna_fragments'] = nucleic_acid_type(spectra)
    r['nb_per_spectrum'] = nb
    r['pmz'] = np.array([s['pmz'] for s in spectra])
    return r


ADVICE = {
    'SEVERE': ('A large share of MS2 spectra carry co-isolated nucleic acid. Spectra with intense nucleobase ions are '
               'identified at roughly 40-65 % of the expected rate (unadjusted; much of this coincides with weaker '
               'precursors), and the search strategies we tested (removing '
               'nucleic-acid peaks, chimeric, wide-window and open search) did not recover them. Review the sample '
               'preparation: nuclease digestion in a compatible buffer plus protein precipitation before enrichment. In '
               'the PRIDE atlas, contamination is associated with phosphopeptide enrichment and lysate sonication.'),
    'MODERATE': ('Appreciable contamination (above a level reached by only ~6 % of PRIDE projects). Consider adding a '
                 'nuclease to the lysate and re-checking the first runs.'),
    'MILD': 'Low contamination; the expected effect on identifications is small.',
    'CLEAN': 'No appreciable nucleic-acid contamination.',
    'NOT EVALUABLE': 'This file cannot be measured.',
}
