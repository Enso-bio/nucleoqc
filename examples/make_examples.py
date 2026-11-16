# -*- coding: utf-8 -*-
"""Writes two small synthetic MGF files for trying nucleoqc offline (no download needed):
  example_clean.mgf         400 HCD-like spectra without nucleobase ions          -> CLEAN
  example_contaminated.mgf  400 spectra, 30 % with adenine at 50 % of base peak   -> SEVERE
The spectra are synthetic (exact immonium/y1 reference ions plus random background), not real data.
Usage: python examples/make_examples.py"""
import os
import numpy as np

ADENINE = 136.06177
REF = [70.06513, 72.08078, 86.09643, 110.07127, 120.08078, 129.10224, 136.07569, 147.11280, 175.11895]


def spectra(n, frac, rng):
    out = []
    for i in range(n):
        mz = REF + list(rng.uniform(200, 1500, 60)) + [60.0, 1000.0]
        it = list(rng.uniform(0.05, 0.4, len(REF))) + list(rng.uniform(0.01, 0.3, 60)) + [0.01, 1.0]
        if i < int(n * frac):
            mz.append(ADENINE); it.append(0.5)
        o = np.argsort(mz)
        out.append((float(rng.uniform(400, 1200)), np.array(mz)[o], 1e6 * np.array(it)[o]))
    return out


def write(path, sp):
    with open(path, 'w') as f:
        for k, (pmz, mz, it) in enumerate(sp):
            f.write('BEGIN IONS\nTITLE=synthetic_%d\nPEPMASS=%.5f\nCHARGE=2+\n' % (k, pmz))
            f.writelines('%.5f %.1f\n' % (m, t) for m, t in zip(mz, it))
            f.write('END IONS\n')


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(1)
    write(os.path.join(here, 'example_clean.mgf'), spectra(400, 0.0, rng))
    write(os.path.join(here, 'example_contaminated.mgf'), spectra(400, 0.3, rng))
    print('wrote example_clean.mgf and example_contaminated.mgf in', here)
