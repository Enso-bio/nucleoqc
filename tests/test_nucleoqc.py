# -*- coding: utf-8 -*-
"""Tests with synthetic spectra: metrics, validity, MGF/mzML readers, nucleic-acid type and CLI."""
import base64
import zlib
import numpy as np
import pytest
from nucleoqc import reading, metrics, cli

REF = list(metrics.REFS_CAL)


def spectrum(rng, nb=0.0, base=None, shift=0.0, tyr=False, mz_min=60.0):
    """HCD-like spectrum: exact reference ions + background peaks; optionally a nucleobase ion."""
    mz = [m + shift for m in REF] + list(rng.uniform(200, 1500, 60))
    it = list(rng.uniform(0.05, 0.4, len(REF))) + list(rng.uniform(0.01, 0.3, 60))
    mz.append(mz_min); it.append(0.01)
    mz.append(1000.0); it.append(1.0)                      # base peak
    if nb > 0:
        mz.append((base or metrics.ADENINE) + shift); it.append(nb)
    if tyr:
        mz.append(136.07569 + shift); it.append(0.8)       # Tyr immonium ion, NOT adenine
    return {'mz': np.array(mz), 'it': np.array(it), 'pmz': float(rng.uniform(400, 1200))}


def batch(n=400, frac=0.0, level=0.5, **kw):
    rng = np.random.default_rng(1)
    k = int(round(n * frac))
    return [reading._spectrum(**spectrum(rng, level if i < k else 0.0, **kw)) for i in range(n)]


def test_clean():
    r = metrics.measure(batch())
    assert r['valid'] and r['nb10'] == 0 and r['verdict'] == 'CLEAN' and r['estimated_loss'] == 0


def test_contaminated_and_loss():
    r = metrics.measure(batch(frac=0.5, level=0.5))
    assert r['valid'] and abs(r['nb10'] - 0.5) < 1e-9 and r['verdict'] == 'SEVERE'
    assert abs(r['estimated_loss'] - 0.5 * (1 - 0.38)) < 1e-9
    r = metrics.measure(batch(frac=0.5, level=0.15))
    assert abs(r['estimated_loss'] - 0.5 * (1 - 0.70)) < 1e-9 and r['nb30'] == 0


def test_guanine_counts():
    r = metrics.measure(batch(frac=0.1, level=0.5, base=metrics.GUANINE))
    assert abs(r['nb10'] - 0.1) < 1e-9 and r['verdict'] == 'MODERATE'


def test_thresholds_2_and_10():
    r = metrics.measure(batch(frac=0.3, level=0.05))
    assert r['nb10'] == 0 and abs(r['nb2'] - 0.3) < 1e-9


def test_tyr_is_not_adenine():
    r = metrics.measure(batch(tyr=True))
    assert r['valid'] and r['nb10'] == 0


def test_recalibration():
    # whole spectrum shifted +18 ppm: the offset is estimated and adenine is still counted
    rng = np.random.default_rng(2)
    spectra = []
    for i in range(400):
        s = spectrum(rng, 0.5 if i < 200 else 0.0)
        s['mz'] = s['mz'] * (1 + 18e-6)
        spectra.append(reading._spectrum(s['mz'], s['it'], s['pmz']))
    r = metrics.measure(spectra)
    assert abs(r['ppm'] - 18) < 1 and abs(r['nb10'] - 0.5) < 1e-9
    assert abs(spectra[0]['mz'][np.argmin(np.abs(spectra[0]['mz'] - metrics.ADENINE))] - metrics.ADENINE) < 1e-4


def test_low_resolution_not_evaluable():
    r = metrics.measure(batch(frac=0.5, shift=0.05))
    assert not r['valid'] and r['verdict'] == 'NOT EVALUABLE'


def test_no_low_masses_not_evaluable():
    spectra = [dict(s, mz=s['mz'][s['mz'] > 140], it=s['it'][s['mz'] > 140]) for s in batch(frac=0.5)]
    r = metrics.measure(spectra)
    assert not r['valid'] and '136' in r['reason']


def test_too_few_spectra():
    r = metrics.measure(batch(n=50))
    assert not r['valid'] and r['verdict'] == 'NOT EVALUABLE'


def test_rna_dna_type():
    rng = np.random.default_rng(5)

    def with_fragment(frag):
        out = []
        for i in range(400):
            s = spectrum(rng, 0.5)
            if i < 100:
                s['mz'] = np.append(s['mz'], frag); s['it'] = np.append(s['it'], 0.2)
            out.append(reading._spectrum(s['mz'], s['it'], s['pmz']))
        return out
    assert metrics.measure(with_fragment(metrics.RNA_FRAGMENTS['UMP-H2O']))['type'] == 'RNA'
    assert metrics.measure(with_fragment(metrics.DNA_FRAGMENTS['dCMP-H2O']))['type'] == 'DNA'
    assert metrics.measure(batch(frac=0.5))['type'] == 'undetermined'


def test_percentile():
    assert metrics.pride_percentile(0.0) == 0.0
    assert 0.85 < metrics.pride_percentile(0.02) < 0.95
    assert metrics.pride_percentile(1.0) == 1.0


# ---------------------------------------------------------------- readers
def write_mgf(spectra, path):
    with open(path, 'w') as f:
        for i, s in enumerate(spectra):
            f.write('BEGIN IONS\nTITLE=s%d\nPEPMASS=%.5f\nCHARGE=2+\n' % (i, s['pmz']))
            for m, t in zip(s['mz'], s['it']):
                f.write('%.6f %.4f\n' % (m, t))
            f.write('END IONS\n')


def write_mzml(spectra, path, numpress=False):
    def arr(a, acc):
        b = base64.b64encode(zlib.compress(np.asarray(a, '<f8').tobytes())).decode()
        cv = '<cvParam accession="MS:1002312" name="numpress"/>' if numpress else ''
        return ('<binaryDataArray><cvParam accession="MS:1000523" name="64-bit float"/>'
                '<cvParam accession="MS:1000574" name="zlib compression"/>%s<cvParam accession="%s"/>'
                '<binary>%s</binary></binaryDataArray>' % (cv, acc, b))
    with open(path, 'w') as f:
        f.write('<?xml version="1.0"?><mzML><run><spectrumList>\n')
        f.write('<spectrum index="0"><cvParam accession="MS:1000511" name="ms level" value="1"/>'
                '<binaryDataArrayList></binaryDataArrayList></spectrum>\n')
        for i, s in enumerate(spectra):
            f.write('<spectrum index="%d"><cvParam accession="MS:1000511" name="ms level" value="2"/>'
                    '<cvParam accession="MS:1000127" name="centroid spectrum"/>'
                    '<precursorList><precursor><selectedIonList><selectedIon>'
                    '<cvParam accession="MS:1000744" name="selected ion m/z" value="%.5f"/>'
                    '</selectedIon></selectedIonList></precursor></precursorList><binaryDataArrayList>%s%s'
                    '</binaryDataArrayList></spectrum>\n' % (i + 1, s['pmz'], arr(s['mz'], 'MS:1000514'), arr(s['it'], 'MS:1000515')))
        f.write('</spectrumList></run></mzML>\n')


@pytest.mark.parametrize('fmt', ['mgf', 'mzml'])
def test_readers(tmp_path, fmt):
    spectra = batch(frac=0.25)
    path = str(tmp_path / ('a.' + ('mgf' if fmt == 'mgf' else 'mzML')))
    (write_mgf if fmt == 'mgf' else write_mzml)(spectra, path)
    read = reading.parse(reading.read_local(path)[0], fmt)
    assert len(read) == len(spectra)                            # the mzML MS1 spectrum is ignored
    assert abs(metrics.measure(read)['nb10'] - 0.25) < 1e-9


def test_ms3_skipped(tmp_path):
    path = str(tmp_path / 'sps.mzML')
    write_mzml(batch(n=50), path)
    txt = open(path).read()
    ms3 = txt.split('\n')[2].replace('value="2"', 'value="3"', 1).replace('index="1"', 'index="999"', 1)
    open(path, 'w').write(txt.replace('</spectrumList>', ms3 + '\n</spectrumList>'))
    assert len(reading.parse(reading.read_local(path)[0], 'mzml')) == 50     # the MS3 scan is not read as MS2


def test_rel_takes_maximum_within_tolerance():
    # three peaks within +-0.005 Th of adenine: the most intense is the middle one, not one of the two neighbours
    mz = np.array([100.0, metrics.ADENINE - 0.004, metrics.ADENINE - 0.001, metrics.ADENINE + 0.001, metrics.ADENINE + 0.004, 300.0])
    it = np.array([1000.0, 10.0, 20.0, 500.0, 30.0, 1.0])
    o = np.argsort(mz)
    s = {'mz': mz[o], 'it': it[o]}
    assert metrics.rel(s, metrics.ADENINE) == 0.5
    mz2 = np.array([100.0, metrics.ADENINE - 0.003, metrics.ADENINE - 0.002, metrics.ADENINE - 0.001, metrics.ADENINE + 0.0005, 300.0])
    it2 = np.array([1000.0, 800.0, 10.0, 10.0, 10.0, 1.0])
    assert metrics.rel({'mz': mz2, 'it': it2}, metrics.ADENINE) == 0.8   # outside the two nearest peaks
    assert metrics.rel({'mz': np.array([100.0, 136.07]), 'it': np.array([1.0, 1.0])}, metrics.ADENINE) == 0.0


def test_numpress(tmp_path):
    path = str(tmp_path / 'np.mzML')
    write_mzml(batch(n=5), path, numpress=True)
    with pytest.raises(reading.NumpressError):
        reading.parse(reading.read_local(path)[0], 'mzml')


def test_chunks_cut_spectra(tmp_path):
    path = str(tmp_path / 'large.mgf')
    write_mgf(batch(n=3000, frac=0.5), path)
    chunks, size = reading.read_local(path, chunk=200_000)
    assert len(chunks) == 3
    read = reading.parse(chunks, 'mgf')
    assert 0 < len(read) < 3000                                 # spectra cut at chunk edges are discarded
    assert all(len(s['mz']) >= reading.MIN_PEAKS for s in read)


def test_cli(tmp_path, capsys):
    write_mgf(batch(frac=0.4), str(tmp_path / 'dirty.mgf'))
    write_mgf(batch(), str(tmp_path / 'clean.mgf'))
    tsv = str(tmp_path / 'qc.tsv')
    code = cli.main([str(tmp_path), '--tsv', tsv, '--json', str(tmp_path / 'qc.json')])
    out = capsys.readouterr().out
    assert code == 1 and 'SEVERE' in out and 'CLEAN' in out
    rows = open(tsv, encoding='utf-8').read().strip().split('\n')
    assert len(rows) == 3 and rows[0].startswith('source\tfile\tverdict')


def test_unknown_format():
    with pytest.raises(ValueError):
        reading.file_format('x.raw')
