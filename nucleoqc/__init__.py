# -*- coding: utf-8 -*-
"""nucleoqc: detection of co-isolated nucleic acid in proteomics MS2 spectra."""
__version__ = '0.3.0'
__authors__ = (('Enso Onill Torres Alegre', 'University of Notre Dame'),
               ('Diana E. Mora Jimenez', 'University of Notre Dame'))
__url__ = 'https://github.com/Enso-bio/nucleoqc'

from .metrics import measure, estimated_loss, verdict, nucleic_acid_type   # noqa: E402,F401
from .reading import mgf_spectra, mzml_spectra, read_local, parse          # noqa: E402,F401
