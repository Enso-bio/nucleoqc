# Changelog

## 0.3.0 (2026-11-16)
- Desktop application (`nucleoqc-gui`, `nucleoqc --gui` or `python -m nucleoqc.gui`): add files, folders or PRIDE
  projects, run the QC, colour-coded results with per-file advice, TSV and PNG export, quick guide, About and
  citation dialogs. It uses the same measurement as the command line; Tkinter only, no new dependency.
- Thermo .raw files are recognised in the application and the user is told how to convert them.
- `paper/`: analysis code and tables in English (file names, column names, comments and messages); the search work
  folder is set with `NUCLEOQC_SEARCH_WORK`. Re-running the offline analyses gives the same numbers as before.

## 0.2.1 (2026-11-14)
- mzML: only MS2 scans are read; SPS-MS3 reporter scans were previously counted as MS2. Re-reading the 145 valid
  mzML files of the PRIDE atlas, together with the next change this altered nb10 in 4 files and no file changed
  status at the 1 % cut-off.
- Peak search: the most intense peak within ±0.005 Th is used, not only the two peaks nearest the target mass.
- RNA/DNA type is reported as a tentative call.
- PRIDE reference distribution: 406 projects (one project previously lost to a malformed table row restored).
- `--cite` option; the citation is printed at the end of each run.
- Two synthetic example files (`examples/`) for trying the tool offline.
- 19 unit tests.

## 0.2.0 (2026-11-14)
- English interface, field names and documentation; modules renamed to `reading`, `metrics`, `reference`.
- Tentative RNA/DNA call from nucleotide fragments.

## 0.1.0
- First version used for the analyses of the manuscript (atlas measurement, loss estimate, CLI).
