#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_gf_spectrum.py -- spectral functions A(w) of the retarded G (DYSON route), from gf_spectrum.py files.

Make the files first, e.g.:   python gf_spectrum.py v_gbg/ndb.G --out bg
Everything you may want to change is in the SETTINGS block below.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ======================= SETTINGS: edit here =======================

TITLE  = r'WS2 9$\times$9, $\Gamma$: spectral functions of the PPA $G_0W_0$ Green function'
XCOL   = 0           # x axis: 0 = w (eV), 4 = w - Eo (eV)
XMIN   = -40.0       # x axis from (eV)
XMAX   = 45.0        # x axis to (eV)
LOGY   = True        # True: log scale (shows the satellites); False: linear (shows the QP peaks)
YMIN   = 1e-4        # lower y limit with LOGY
OUT    = 'gf_spectrum.png'

# Vertical dotted lines (energy, label): e.g. the Newton QP energies (o-v_nbg.qp: Eo + E-Eo)
LINES = [
    (-0.5992, 'VB QP (Newton)'),
    ( 3.4226, 'CB QP (Newton)'),
]

# One line per curve:  (file, column, label, colour, line style, line width)
#   column: 1 = A(w) (1/eV), 2 = Re Sigma_c (eV), 3 = Im Sigma_c (eV)   (plot one kind per figure)
#   Put a # in front of a line to leave that curve out.
CURVES = [
    ('bg_k1_b13.dat',    1, 'band 13 (VB), BG',            'tab:blue',   '-',  1.6),
    ('bg_k1_b14.dat',    1, 'band 14 (CB), BG',            'tab:red',    '-',  1.6),
#   ('bg_k1_b15.dat',    1, 'band 15, BG',                 'tab:orange', '-',  1.2),
    ('none_k1_b13.dat',  1, 'band 13, no terminator',      'navy',       ':',  1.4),
    ('none_k1_b14.dat',  1, 'band 14, no terminator',      'darkred',    ':',  1.4),
]

# ===================================================================

fig, ax = plt.subplots(figsize=(10.5, 4.6))
top = 0.0
for file, col, label, colour, style, width in CURVES:
    try:
        data = np.loadtxt(file, comments='#')
    except OSError:
        print(f'  not found, skipped: {file}')
        continue
    w, y = data[:, XCOL], data[:, col]
    keep = (w >= XMIN) & (w <= XMAX)
    if LOGY:
        y = np.where(y > 0, y, np.nan)
    ax.plot(w[keep], y[keep], label=label, color=colour, ls=style, lw=width)
    top = max(top, np.nanmax(y[keep]))

for i, (x, name) in enumerate(LINES):
    ax.axvline(x, color='k', lw=0.7, ls=':')
    y = top * 2 * 0.2 ** i if LOGY else top * (1.04 - 0.12 * i)     # staggered labels
    ax.text(x + 0.3, y, f'{name}\n{x:.3f} eV', fontsize=8, va='top')

ax.set_xlim(XMIN, XMAX)
if LOGY:
    ax.set_yscale('log')
    ax.set_ylim(YMIN, top * 3)
else:
    ax.set_ylim(0, top * 1.1)
ax.set_xlabel(r'$\omega$ (eV)' if XCOL == 0 else r'$\omega - \varepsilon_{KS}$ (eV)')
ax.set_ylabel(r'$A(\omega)$ (1/eV)')
ax.set_title(TITLE, fontsize=11)
ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=False)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
