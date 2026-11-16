# -*- coding: utf-8 -*-
"""nucleoqc desktop application (Tkinter, no extra dependencies).

Start it with `nucleoqc-gui`, `nucleoqc --gui` or `python -m nucleoqc.gui`.
It runs exactly the same measurement as the command line (cli.evaluate).
"""
import os
import queue
import sys
import threading
import webbrowser
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from . import __version__, __authors__, __url__, metrics, cli

# ---------------------------------------------------------------- palette
GREEN_DARK = '#1b5e20'
GREEN = '#2e7d32'
GREEN_HOVER = '#388e3c'
GREEN_LIGHT = '#e8f5e9'
GREEN_PALE = '#f4faf4'
LINE = '#c8e6c9'
TEXT = '#1f2d1f'
MUTED = '#5f6f5f'
VERDICT_BG = {'SEVERE': '#f6d5d3', 'MODERATE': '#fbe3c6', 'MILD': '#f7f0c4', 'CLEAN': '#d7ecd7',
              'NOT EVALUABLE': '#e6e6e6', 'ERROR': '#e6e6e6'}
VERDICT_FG = {'SEVERE': '#8c1d18', 'MODERATE': '#8a4b00', 'MILD': '#6b5b00', 'CLEAN': '#1b5e20',
              'NOT EVALUABLE': '#555555', 'ERROR': '#555555'}

ACCEPTED = ('.mgf', '.mzml')
RAW_HELP = ('Thermo .raw files cannot be read directly. Convert them to mzML first, for example with\n'
            'ThermoRawFileParser:  ThermoRawFileParser -i file.raw -f 2\n'
            'or with msconvert (ProteoWizard), without numpress compression.')

GUIDE = """How to use nucleoqc

1. Add your data
   - "Add files" for .mgf or .mzML files (uncompressed; centroid or profile).
   - "Add folder" to add every .mgf/.mzML file in a folder.
   - Or type a PRIDE accession (for example PXD023674) and press "Add PRIDE project".
     Public files are read in part over the internet, without downloading them entirely.
   Thermo .raw files must be converted to mzML first (ThermoRawFileParser or msconvert).

2. Press "Run QC". Each file takes a few seconds (three 4 MB chunks are read).
   Tick "Read whole files" when the exact figure matters.

3. Read the verdict
   CLEAN     nb10 < 1 %      no appreciable nucleic-acid contamination
   MILD      1-5 %           small expected effect on identifications
   MODERATE  5-20 %          consider adding a nuclease to the lysate
   SEVERE    >= 20 %         review the sample preparation
   NOT EVALUABLE             the MS2 spectra are low-resolution or do not record m/z < 136

   nb10 = % of MS2 spectra in which adenine (136.0618) or guanine (152.0567) reaches
   at least 10 % of the base peak. Click a row to see all metrics and the advice.

4. Export the table (TSV) or a figure (PNG) for your records.

Note: ADP-ribosylation, protein-RNA cross-linking, ChIP-MS and ATP assays give signal
by design; in those studies a high nb10 is not contamination.
"""


def citation():
    return cli.CITATION


# ---------------------------------------------------------------- helpers
def _resource(name):
    base = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    p = os.path.join(base, name)
    if not os.path.exists(p):
        p = os.path.join(base, 'nucleoqc', name)
    return p


def _pc(x):
    return '%.1f %%' % (100 * x) if isinstance(x, float) else '-'


def _open_path(path):
    try:
        if sys.platform.startswith('win'):
            os.startfile(path)                      # noqa: S606
        elif sys.platform == 'darwin':
            import subprocess; subprocess.Popen(['open', path])
        else:
            import subprocess; subprocess.Popen(['xdg-open', path])
    except OSError:
        pass


class GreenButton(tk.Label):
    """Flat green button that looks the same on Windows, macOS and Linux."""

    def __init__(self, master, text, command, big=False, secondary=False, **kw):
        self._bg = GREEN_LIGHT if secondary else GREEN
        self._hover = LINE if secondary else GREEN_HOVER
        self._fg = GREEN_DARK if secondary else 'white'
        super().__init__(master, text=text, bg=self._bg, fg=self._fg, cursor='hand2',
                         font=('Segoe UI', 13 if big else 10, 'bold' if big or not secondary else 'normal'),
                         padx=22 if big else 12, pady=9 if big else 5, **kw)
        self._command = command
        self._enabled = True
        self.bind('<Enter>', lambda e: self._enabled and self.config(bg=self._hover))
        self.bind('<Leave>', lambda e: self.config(bg=self._bg if self._enabled else '#bdbdbd'))
        self.bind('<Button-1>', lambda e: self._enabled and self._command())

    def enable(self, on=True):
        self._enabled = on
        self.config(bg=self._bg if on else '#bdbdbd', fg=self._fg if on else 'white', cursor='hand2' if on else 'arrow')


# ---------------------------------------------------------------- application
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('nucleoqc %s - nucleic acid QC for proteomics' % __version__)
        self.geometry('1120x780')
        self.minsize(900, 600)
        self.configure(bg=GREEN_PALE)
        try:
            self._icon = tk.PhotoImage(file=_resource('icon64.png'))
            self.iconphoto(True, self._icon)
        except tk.TclError:
            self._icon = None
        self.inputs = []          # (label, value) ; value is a path or a PXD accession
        self.results = []
        self.q = queue.Queue()
        self.running = False
        self._style()
        self._menu()
        self._build()
        self.protocol('WM_DELETE_WINDOW', self.destroy)

    # ------------------------------------------------------------ look
    def _style(self):
        s = ttk.Style(self)
        try:
            s.theme_use('clam')
        except tk.TclError:
            pass
        s.configure('.', background=GREEN_PALE, foreground=TEXT, font=('Segoe UI', 10))
        s.configure('Treeview', background='white', fieldbackground='white', rowheight=26, bordercolor=LINE)
        s.configure('Treeview.Heading', background=GREEN, foreground='white', font=('Segoe UI', 10, 'bold'),
                    relief='flat', padding=5)
        s.map('Treeview.Heading', background=[('active', GREEN_HOVER)])
        s.map('Treeview', background=[('selected', GREEN_HOVER)], foreground=[('selected', 'white')])
        s.configure('green.Horizontal.TProgressbar', troughcolor=GREEN_LIGHT, background=GREEN, bordercolor=LINE,
                    lightcolor=GREEN, darkcolor=GREEN)
        s.configure('TCheckbutton', background='white')
        s.configure('TSpinbox', arrowsize=12)
        s.configure('Vertical.TScrollbar', background=GREEN_LIGHT, troughcolor='white', bordercolor=LINE,
                    arrowcolor=GREEN_DARK)

    def _menu(self):
        m = tk.Menu(self)
        f = tk.Menu(m, tearoff=0)
        f.add_command(label='Add files...', command=self.add_files, accelerator='Ctrl+O')
        f.add_command(label='Add folder...', command=self.add_folder)
        f.add_separator()
        f.add_command(label='Export table (TSV)...', command=self.export_tsv)
        f.add_command(label='Save figure (PNG)...', command=self.save_figure)
        f.add_separator()
        f.add_command(label='Exit', command=self.destroy)
        m.add_cascade(label='File', menu=f)
        h = tk.Menu(m, tearoff=0)
        h.add_command(label='Quick guide', command=self.show_guide, accelerator='F1')
        h.add_command(label='Try the example files', command=self.load_examples)
        h.add_command(label='How to cite', command=self.show_cite)
        if 'YOUR-GITHUB-USER' not in __url__:
            h.add_command(label='Project website', command=lambda: webbrowser.open(__url__))
        h.add_separator()
        h.add_command(label='About nucleoqc', command=self.show_about)
        m.add_cascade(label='Help', menu=h)
        self.config(menu=m)
        self.bind('<Control-o>', lambda e: self.add_files())
        self.bind('<F1>', lambda e: self.show_guide())

    def _card(self, parent, title):
        outer = tk.Frame(parent, bg=LINE, padx=1, pady=1)
        inner = tk.Frame(outer, bg='white', padx=14, pady=10)
        inner.pack(fill='both', expand=True)
        tk.Label(inner, text=title, bg='white', fg=GREEN_DARK, font=('Segoe UI', 11, 'bold')).pack(anchor='w', pady=(0, 6))
        return outer, inner

    def _build(self):
        # header
        hd = tk.Frame(self, bg=GREEN_DARK, padx=18, pady=10)
        hd.pack(fill='x')
        if self._icon:
            self._icon_small = self._icon.subsample(1)
            tk.Label(hd, image=self._icon_small, bg=GREEN_DARK).pack(side='left', padx=(0, 12))
        t = tk.Frame(hd, bg=GREEN_DARK); t.pack(side='left')
        tk.Label(t, text='nucleoqc', bg=GREEN_DARK, fg='white', font=('Segoe UI', 20, 'bold')).pack(anchor='w')
        tk.Label(t, text='Quality control for co-isolated nucleic acid in proteomics MS2 spectra',
                 bg=GREEN_DARK, fg='#c8e6c9', font=('Segoe UI', 10)).pack(anchor='w')
        r = tk.Frame(hd, bg=GREEN_DARK); r.pack(side='right')
        tk.Label(r, text='version %s' % __version__, bg=GREEN_DARK, fg='#c8e6c9', font=('Segoe UI', 9)).pack(anchor='e')
        a = tk.Label(r, text='About  |  Guide  |  Cite', bg=GREEN_DARK, fg='white', cursor='hand2',
                     font=('Segoe UI', 9, 'underline'))
        a.pack(anchor='e'); a.bind('<Button-1>', lambda e: self.show_about())

        # footer (authors always visible)
        ft = tk.Frame(self, bg=GREEN_LIGHT, padx=18, pady=5)
        ft.pack(side='bottom', fill='x')
        names = ' & '.join(n for n, _ in __authors__)
        tk.Label(ft, text='Developed by %s  -  %s  -  MIT licence' % (names, __authors__[0][1]),
                 bg=GREEN_LIGHT, fg=MUTED, font=('Segoe UI', 9)).pack(side='left')
        c = tk.Label(ft, text='How to cite', bg=GREEN_LIGHT, fg=GREEN, cursor='hand2', font=('Segoe UI', 9, 'underline'))
        c.pack(side='right'); c.bind('<Button-1>', lambda e: self.show_cite())

        body = tk.Frame(self, bg=GREEN_PALE, padx=14, pady=12)
        body.pack(fill='both', expand=True)
        body.columnconfigure(0, weight=0, minsize=330); body.columnconfigure(1, weight=1); body.rowconfigure(0, weight=1)

        # ---- left: inputs
        left = tk.Frame(body, bg=GREEN_PALE); left.grid(row=0, column=0, sticky='nsew', padx=(0, 12))
        card, c1 = self._card(left, '1   Add your data')
        card.pack(fill='both', expand=True)
        bb = tk.Frame(c1, bg='white'); bb.pack(fill='x')
        GreenButton(bb, 'Add files', self.add_files).pack(side='left')
        GreenButton(bb, 'Add folder', self.add_folder).pack(side='left', padx=6)
        tk.Label(c1, text='.mgf or .mzML  (convert .raw first)', bg='white', fg=MUTED,
                 font=('Segoe UI', 8)).pack(anchor='w', pady=(4, 6))
        pr = tk.Frame(c1, bg='white'); pr.pack(fill='x', pady=(2, 2))
        tk.Label(pr, text='PRIDE', bg='white', fg=TEXT).pack(side='left')
        self.pxd = tk.StringVar()
        e = tk.Entry(pr, textvariable=self.pxd, width=12, relief='solid', bd=1, highlightthickness=1,
                     highlightcolor=GREEN, highlightbackground=LINE, font=('Segoe UI', 10))
        e.pack(side='left', padx=6, ipady=3); e.insert(0, ''); e.bind('<Return>', lambda ev: self.add_pride())
        tk.Label(pr, text='files', bg='white', fg=MUTED).pack(side='left')
        self.maxf = tk.IntVar(value=5)
        ttk.Spinbox(pr, from_=1, to=200, textvariable=self.maxf, width=4).pack(side='left', padx=4)
        GreenButton(c1, 'Add PRIDE project', self.add_pride, secondary=True).pack(anchor='w', pady=(4, 8))

        lf = tk.Frame(c1, bg='white'); lf.pack(fill='both', expand=True)
        self.listbox = tk.Listbox(lf, relief='flat', bd=0, highlightthickness=1, highlightbackground=LINE,
                                  selectbackground=GREEN_HOVER, selectmode='extended', font=('Segoe UI', 9),
                                  activestyle='none')
        sb = ttk.Scrollbar(lf, orient='vertical', command=self.listbox.yview)
        self.listbox.configure(yscrollcommand=sb.set)
        self.listbox.pack(side='left', fill='both', expand=True); sb.pack(side='right', fill='y')
        self.listbox.bind('<Delete>', lambda e: self.remove_selected())
        rb = tk.Frame(c1, bg='white'); rb.pack(fill='x', pady=(6, 0))
        GreenButton(rb, 'Remove', self.remove_selected, secondary=True).pack(side='left')
        GreenButton(rb, 'Clear', self.clear_inputs, secondary=True).pack(side='left', padx=6)
        self.count = tk.Label(rb, text='0 inputs', bg='white', fg=MUTED, font=('Segoe UI', 9))
        self.count.pack(side='right')

        card2, c2 = self._card(left, '2   Run')
        card2.pack(fill='x', pady=(12, 0))
        self.whole = tk.BooleanVar(value=False)
        ttk.Checkbutton(c2, text='Read whole files (slower; exact figure)', variable=self.whole).pack(anchor='w')
        self.run_btn = GreenButton(c2, 'Run QC', self.run, big=True)
        self.run_btn.pack(fill='x', pady=(8, 6))
        self.progress = ttk.Progressbar(c2, style='green.Horizontal.TProgressbar', mode='determinate')
        self.progress.pack(fill='x')
        self.status = tk.Label(c2, text='Add files and press Run QC.', bg='white', fg=MUTED, font=('Segoe UI', 9),
                               anchor='w', justify='left', wraplength=290)
        self.status.pack(fill='x', pady=(4, 0))

        # ---- right: results
        right = tk.Frame(body, bg=GREEN_PALE); right.grid(row=0, column=1, sticky='nsew')
        card3, c3 = self._card(right, '3   Results')
        card3.pack(fill='both', expand=True)
        self.summary = tk.Label(c3, text='', bg='white', fg=TEXT, font=('Segoe UI', 10), anchor='w', justify='left')
        self.summary.pack(fill='x')
        eb = tk.Frame(c3, bg='white'); eb.pack(side='bottom', fill='x', pady=(8, 0))
        self.tsv_btn = GreenButton(eb, 'Export table (TSV)', self.export_tsv, secondary=True)
        self.tsv_btn.pack(side='left')
        self.fig_btn = GreenButton(eb, 'Save figure (PNG)', self.save_figure, secondary=True)
        self.fig_btn.pack(side='left', padx=6)
        self.tsv_btn.enable(False); self.fig_btn.enable(False)
        df = tk.Frame(c3, bg='white'); df.pack(side='bottom', fill='x', pady=(8, 0))
        self.detail = tk.Text(df, height=10, relief='flat', bg=GREEN_PALE, fg=TEXT, wrap='word', font=('Segoe UI', 9),
                              padx=10, pady=8, highlightthickness=1, highlightbackground=LINE)
        ds = ttk.Scrollbar(df, orient='vertical', command=self.detail.yview)
        self.detail.configure(yscrollcommand=ds.set)
        self.detail.pack(side='left', fill='x', expand=True); ds.pack(side='right', fill='y')
        tf = tk.Frame(c3, bg='white'); tf.pack(fill='both', expand=True, pady=(4, 0))
        cols = ('file', 'verdict', 'nb10', 'loss', 'pct', 'type', 'source')
        heads = ('File', 'Verdict', 'nb10', 'Est. loss', 'PRIDE pct.', 'RNA/DNA', 'Source')
        widths = (200, 105, 62, 72, 74, 108, 70)
        self.tree = ttk.Treeview(tf, columns=cols, show='headings', selectmode='browse')
        for c, h, w in zip(cols, heads, widths):
            self.tree.heading(c, text=h, command=lambda c=c: self.sort(c))
            self.tree.column(c, width=w, anchor='w' if c in ('file', 'source', 'type') else 'center', stretch=c == 'file')
        for v in VERDICT_BG:
            self.tree.tag_configure(v, background=VERDICT_BG[v], foreground=VERDICT_FG[v])
        ts = ttk.Scrollbar(tf, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscrollcommand=ts.set)
        self.tree.pack(side='left', fill='both', expand=True); ts.pack(side='right', fill='y')
        self.tree.bind('<<TreeviewSelect>>', lambda e: self.show_detail())

        self.detail.tag_configure('h', font=('Segoe UI', 11, 'bold'), foreground=GREEN_DARK)
        self.detail.tag_configure('b', font=('Segoe UI', 9, 'bold'))
        self._set_detail([('Click a result to see all metrics and what to do.\n', ''),
                          ('New here? Help > Try the example files runs two small demo files.', '')])

    # ------------------------------------------------------------ inputs
    def _add(self, label, value):
        if any(v == value for _, v in self.inputs):
            return
        self.inputs.append((label, value)); self.listbox.insert('end', label)
        self.count.config(text='%d input%s' % (len(self.inputs), '' if len(self.inputs) == 1 else 's'))

    def add_files(self):
        paths = filedialog.askopenfilenames(
            title='Choose MGF or mzML files',
            filetypes=[('MS/MS spectra', '*.mgf *.mzML *.mzml'), ('Thermo raw (needs conversion)', '*.raw'),
                       ('All files', '*.*')])
        raw = []
        for p in paths:
            if p.lower().endswith('.raw'):
                raw.append(os.path.basename(p)); continue
            if not p.lower().endswith(ACCEPTED):
                messagebox.showwarning('nucleoqc', 'Not an .mgf or .mzML file:\n%s' % p); continue
            self._add(os.path.basename(p), p)
        if raw:
            messagebox.showinfo('Thermo .raw files', '%d .raw file(s) were not added.\n\n%s' % (len(raw), RAW_HELP))

    def add_folder(self):
        d = filedialog.askdirectory(title='Choose a folder with MGF or mzML files')
        if not d:
            return
        n = 0
        for f in sorted(os.listdir(d)):
            if f.lower().endswith(ACCEPTED):
                self._add(f, os.path.join(d, f)); n += 1
        raw = sum(f.lower().endswith('.raw') for f in os.listdir(d))
        if not n:
            messagebox.showinfo('nucleoqc', 'No .mgf or .mzML files in this folder.' + ('\n\n' + RAW_HELP if raw else ''))
        elif raw:
            self.status.config(text='%d files added; %d .raw files skipped (convert them to mzML first).' % (n, raw))

    def add_pride(self):
        acc = self.pxd.get().strip().upper()
        if not (acc.startswith('PXD') and acc[3:].isdigit() and len(acc) == 9):
            messagebox.showwarning('nucleoqc', 'Type a PRIDE accession such as PXD023674.'); return
        self._add('%s  (PRIDE, up to %d files)' % (acc, self.maxf.get()), acc)
        self.pxd.set('')

    def remove_selected(self):
        for i in reversed(self.listbox.curselection()):
            self.listbox.delete(i); del self.inputs[i]
        self.count.config(text='%d inputs' % len(self.inputs))

    def clear_inputs(self):
        self.listbox.delete(0, 'end'); self.inputs = []; self.count.config(text='0 inputs')

    def load_examples(self):
        d = _resource('examples')
        if not os.path.isdir(d):
            d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'examples')
        files = [os.path.join(d, f) for f in ('example_clean.mgf', 'example_contaminated.mgf')]
        if not all(os.path.exists(f) for f in files):
            messagebox.showinfo('nucleoqc', 'Example files are not included in this installation.'); return
        self.clear_inputs()
        for f in files:
            self._add(os.path.basename(f), f)
        self.run()

    # ------------------------------------------------------------ run
    def run(self):
        if self.running:
            return
        if not self.inputs:
            messagebox.showinfo('nucleoqc', 'Add at least one file, folder or PRIDE project first.'); return
        self.running = True; self.run_btn.enable(False); self.tsv_btn.enable(False); self.fig_btn.enable(False)
        self.results = []; self.tree.delete(*self.tree.get_children()); self.summary.config(text='')
        self.progress.config(value=0, maximum=1)
        self.status.config(text='Preparing...')
        args = ([v for _, v in self.inputs], int(self.maxf.get() or 5), bool(self.whole.get()))
        threading.Thread(target=self._worker, args=args, daemon=True).start()
        self.after(100, self._poll)

    def _worker(self, inputs, max_files, whole):
        try:
            targets = []
            for e in inputs:
                if e.upper().startswith('PXD') and not os.path.exists(e):
                    self.q.put(('status', 'Listing public files of %s in PRIDE...' % e))
                try:
                    found = list(cli._targets([e], max_files))
                except Exception as ex:                       # noqa: BLE001  network or format problems
                    self.q.put(('result', {'source': e if e.upper().startswith('PXD') else 'local',
                                           'file': os.path.basename(e), 'verdict': 'NOT EVALUABLE', 'valid': False,
                                           'n': 0, 'reason': str(ex)}))
                    continue
                if not found and e.upper().startswith('PXD'):
                    self.q.put(('result', {'source': e, 'file': e, 'verdict': 'NOT EVALUABLE', 'valid': False, 'n': 0,
                                           'reason': 'no uncompressed MGF/mzML files in PRIDE'}))
                targets += found
            self.q.put(('total', len(targets)))
            for i, (source, name, load, fmt) in enumerate(targets):
                self.q.put(('status', 'Measuring %d of %d: %s' % (i + 1, len(targets), os.path.basename(name))))
                try:
                    r = cli.evaluate(source, name, load, fmt, whole)
                except Exception as ex:                       # noqa: BLE001  keep the window alive
                    r = {'source': source, 'file': os.path.basename(name), 'verdict': 'NOT EVALUABLE', 'valid': False,
                         'n': 0, 'reason': 'error: %s' % ex}
                self.q.put(('result', r))
                self.q.put(('step', i + 1))
        finally:
            self.q.put(('done', None))

    def _poll(self):
        try:
            while True:
                kind, x = self.q.get_nowait()
                if kind == 'status':
                    self.status.config(text=x)
                elif kind == 'total':
                    self.progress.config(maximum=max(x, 1))
                elif kind == 'result':
                    self._add_result(x)
                elif kind == 'step':
                    self.progress.config(value=min(x, float(self.progress.cget('maximum'))))
                elif kind == 'done':
                    self.progress.config(value=float(self.progress.cget('maximum')))
                    self._finish(); return
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def _add_result(self, r):
        self.results.append(r)
        v = r['verdict']
        valid = r.get('valid')
        self.tree.insert('', 'end', iid=str(len(self.results) - 1), tags=(v,), values=(
            r['file'], v, _pc(r.get('nb10')) if valid else '-', ('~' + _pc(r.get('estimated_loss'))) if valid else '-',
            '%.0f' % (100 * r['pride_percentile']) if valid else '-',
            (r.get('type', '') if v != 'CLEAN' else '') if valid else r.get('reason', '')[:40], r['source']))

    def _finish(self):
        self.running = False; self.run_btn.enable(True)
        v = [r for r in self.results if r.get('valid')]
        self.tsv_btn.enable(bool(self.results)); self.fig_btn.enable(bool(v))
        counts = {}
        for r in self.results:
            counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
        order = ['SEVERE', 'MODERATE', 'MILD', 'CLEAN', 'NOT EVALUABLE']
        self.summary.config(text='%d files, %d evaluable:   %s' % (
            len(self.results), len(v), '   '.join('%s %d' % (k, counts[k]) for k in order if k in counts)))
        self.status.config(text='Done.' if self.results else 'Nothing to evaluate.')
        if self.results:
            worst = max(range(len(self.results)),
                        key=lambda i: self.results[i].get('nb10', -1) if self.results[i].get('valid') else -1)
            self.tree.selection_set(str(worst)); self.tree.see(str(worst))

    # ------------------------------------------------------------ detail
    def _set_detail(self, parts):
        self.detail.config(state='normal'); self.detail.delete('1.0', 'end')
        for text, tag in parts:
            self.detail.insert('end', text, tag)
        self.detail.config(state='disabled')

    def show_detail(self):
        sel = self.tree.selection()
        if not sel:
            return
        r = self.results[int(sel[0])]
        p = [('%s   %s\n' % (r['verdict'], r['file']), 'h')]
        if not r.get('valid'):
            p.append(('Not evaluable: %s.\n' % (r.get('reason') or 'unknown reason'), ''))
            if r.get('n'):
                p.append(('%d MS2 spectra read.' % r['n'], ''))
            self._set_detail(p); return
        p += [(metrics.ADVICE[r['verdict']] + '\n\n', ''),
              ('nb10 ', 'b'), ('%s   ' % _pc(r['nb10']), ''), ('nb30 ', 'b'), ('%s   ' % _pc(r['nb30']), ''),
              ('nb2 ', 'b'), ('%s   ' % _pc(r['nb2']), ''),
              ('estimated loss ', 'b'), ('~%s   ' % _pc(r['estimated_loss']), ''),
              ('PRIDE percentile ', 'b'), ('%.0f\n' % (100 * r['pride_percentile']), ''),
              ('type (tentative) ', 'b'), ('%s   (RNA fragments %s, DNA fragments %s)\n' % (
                  r['type'], _pc(r['rna_fragments']), _pc(r['dna_fragments'])), ''),
              ('ions >= 2 %% of base peak:  A %s   G %s   C %s   U %s   ribose %s\n' % (
                  _pc(r['A2']), _pc(r['G2']), _pc(r['C2']), _pc(r['U2']), _pc(r['ribose2'])), ''),
              ('%d MS2 spectra read;  mass offset corrected %.1f ppm.  ' % (r['n'], r['ppm']), ''),
              ('The estimated loss is approximate (dose-response measured with Sage in PRIDE files).', '')]
        self._set_detail(p)

    def sort(self, col):
        key = {'file': lambda r: r['file'].lower(), 'source': lambda r: r['source'],
               'type': lambda r: r.get('type', ''), 'verdict': lambda r: r['verdict'],
               'nb10': lambda r: r.get('nb10', -1) if r.get('valid') else -1,
               'loss': lambda r: r.get('estimated_loss', -1) if r.get('valid') else -1,
               'pct': lambda r: r.get('pride_percentile', -1) if r.get('valid') else -1}[col]
        rev = getattr(self, '_sort', None) == col
        self._sort = None if rev else col
        items = sorted(range(len(self.results)), key=lambda i: key(self.results[i]), reverse=rev or col in ('nb10', 'loss', 'pct'))
        for pos, i in enumerate(items):
            self.tree.move(str(i), '', pos)

    # ------------------------------------------------------------ export
    def export_tsv(self):
        if not self.results:
            return
        p = filedialog.asksaveasfilename(title='Save results table', defaultextension='.tsv',
                                         initialfile='nucleoqc_results.tsv',
                                         filetypes=[('Tab-separated (Excel, R, Python)', '*.tsv'), ('All files', '*.*')])
        if p:
            cli.write_tsv(self.results, p)
            self.status.config(text='Table saved: %s' % p)

    def save_figure(self):
        if not any(r.get('valid') for r in self.results):
            return
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            messagebox.showwarning('nucleoqc', 'The figure needs matplotlib:\n\npip install matplotlib'); return
        p = filedialog.asksaveasfilename(title='Save figure', defaultextension='.png', initialfile='nucleoqc_figure.png',
                                         filetypes=[('PNG image', '*.png')])
        if p:
            cli.figure(self.results, p)
            self.status.config(text='Figure saved: %s' % p)
            _open_path(p)

    # ------------------------------------------------------------ dialogs
    def _dialog(self, title, width=560):
        w = tk.Toplevel(self); w.title(title); w.configure(bg='white'); w.transient(self); w.resizable(False, False)
        hd = tk.Frame(w, bg=GREEN_DARK, padx=16, pady=10); hd.pack(fill='x')
        tk.Label(hd, text=title, bg=GREEN_DARK, fg='white', font=('Segoe UI', 13, 'bold')).pack(anchor='w')
        body = tk.Frame(w, bg='white', padx=18, pady=14, width=width); body.pack(fill='both', expand=True)
        return w, body

    def _center(self, w):
        w.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() - w.winfo_width()) // 2
        y = self.winfo_rooty() + (self.winfo_height() - w.winfo_height()) // 3
        w.geometry('+%d+%d' % (max(x, 0), max(y, 0))); w.grab_set(); w.focus_set()
        w.bind('<Escape>', lambda e: w.destroy())

    def show_about(self):
        w, b = self._dialog('About nucleoqc')
        top = tk.Frame(b, bg='white'); top.pack(fill='x')
        if self._icon:
            tk.Label(top, image=self._icon, bg='white').pack(side='left', padx=(0, 14))
        t = tk.Frame(top, bg='white'); t.pack(side='left')
        tk.Label(t, text='nucleoqc', bg='white', fg=GREEN_DARK, font=('Segoe UI', 18, 'bold')).pack(anchor='w')
        tk.Label(t, text='version %s' % __version__, bg='white', fg=MUTED).pack(anchor='w')
        tk.Label(b, text='Quality control for co-isolated nucleic acid in proteomics MS2 spectra.', bg='white',
                 fg=TEXT, wraplength=520, justify='left').pack(anchor='w', pady=(12, 8))
        tk.Label(b, text='Developed by', bg='white', fg=GREEN_DARK, font=('Segoe UI', 10, 'bold')).pack(anchor='w')
        for name, aff in __authors__:
            tk.Label(b, text='%s  -  %s' % (name, aff), bg='white', fg=TEXT).pack(anchor='w')
        tk.Label(b, text='Citation', bg='white', fg=GREEN_DARK, font=('Segoe UI', 10, 'bold')).pack(anchor='w', pady=(10, 0))
        tk.Label(b, text=citation(), bg='white', fg=TEXT, wraplength=520, justify='left').pack(anchor='w')
        tk.Label(b, text='Licence: MIT.  Detection by nucleic-acid marker ions follows Potel et al. 2018 '
                         '(Mol Cell Proteomics).', bg='white', fg=MUTED, wraplength=520, justify='left',
                 font=('Segoe UI', 9)).pack(anchor='w', pady=(10, 0))
        if 'YOUR-GITHUB-USER' not in __url__:
            link = tk.Label(b, text=__url__, bg='white', fg=GREEN, cursor='hand2', font=('Segoe UI', 9, 'underline'))
            link.pack(anchor='w'); link.bind('<Button-1>', lambda e: webbrowser.open(__url__))
        bb = tk.Frame(b, bg='white'); bb.pack(fill='x', pady=(14, 0))
        GreenButton(bb, 'Close', w.destroy).pack(side='right')
        GreenButton(bb, 'Copy citation', lambda: self._copy(citation()), secondary=True).pack(side='right', padx=6)
        self._center(w)

    def show_cite(self):
        w, b = self._dialog('How to cite nucleoqc')
        tk.Label(b, text='If you use nucleoqc, please cite the article:', bg='white', fg=TEXT).pack(anchor='w')
        tk.Label(b, text=citation(), bg=GREEN_PALE, fg=TEXT, wraplength=520, justify='left', padx=10, pady=8).pack(
            fill='x', pady=8)
        bb = tk.Frame(b, bg='white'); bb.pack(fill='x')
        GreenButton(bb, 'Close', w.destroy).pack(side='right')
        GreenButton(bb, 'Copy citation', lambda: self._copy(citation()), secondary=True).pack(side='right', padx=6)
        self._center(w)

    def show_guide(self):
        w, b = self._dialog('Quick guide')
        t = tk.Text(b, width=86, height=31, relief='flat', bg='white', fg=TEXT, font=('Consolas', 9), wrap='word')
        t.insert('1.0', GUIDE); t.config(state='disabled'); t.pack()
        GreenButton(b, 'Close', w.destroy).pack(anchor='e', pady=(8, 0))
        self._center(w)

    def _copy(self, text):
        self.clipboard_clear(); self.clipboard_append(text)
        self.status.config(text='Citation copied to the clipboard.')


def main():
    if sys.platform.startswith('win'):
        try:                                           # sharp text on high-DPI screens
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except (AttributeError, OSError):
            pass
    App().mainloop()
    return 0


if __name__ == '__main__':
    sys.exit(main())
