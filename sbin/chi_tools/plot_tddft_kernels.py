#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_tddft_kernels.py -- one absorption panel: BSE, TDDFT with the exact kernel, ALDA and RPA.

Run it in prod/:   python plot_tddft_kernels.py
Everything you may want to change is in the SETTINGS block below.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ======================= SETTINGS: edit here =======================

TITLE   = 'WS2: absorption, BSE vs TDDFT kernels'
EMIN    = 2.2            # x axis from (eV)
EMAX    = 3.3            # x axis to (eV)
EXCITON = 2.3774         # dotted line "BSE A" (set to None to hide it)
GAP     = 2.9238         # dotted line "QP gap" (set to None to hide it)
OUT     = 'tddft_kernels.png'

# One line per curve:  (file, column, label, colour, line style, line width)
#
#   column: in Yambo o-* files, 1 = Im alpha (the absorption)
#           in casida_spectra.dat: 1 = BSE, 2 = chi0, 3 = chi (RPA), 4 = exact fxc at w=0,
#                                  5 = exact fxc at E_A, 6 = exact fxc(w) at every w
#   colour: any matplotlib name ('black', 'tab:blue', 'purple', ...) or a hex code ('#2a78d6')
#   line style: '-' solid, '--' dashed, '-.' dash-dot, ':' dotted
#
# Put a # in front of a line to leave that curve out; remove the # to put it back.

CURVES = [
    ('casida_spectra.dat',               1, 'BSE',                              'black',      '-',  3.0),
    ('casida_spectra.dat',               5, 'TDDFT, exact $f_{xc}$ at $E_A$',   'tab:blue',   '-',  1.8),
    ('o-tal.alpha_q1_diago_alda_bse',    1, 'TDDFT, ALDA',                      'tab:purple', '-.', 1.6),
    ('o-thr.alpha_q1_diago_bse',         1, 'RPA (no kernel)',                  'gray',       '--', 1.4),
#   ('casida_spectra.dat',               6, r'TDDFT, exact $f_{xc}(\omega)$',   'navy',       '--', 1.2),
#   ('o-tal_ks.alpha_q1_diago_alda_bse', 1, 'TDDFT, ALDA (DFT energies)',       'plum',       '-.', 1.3),
#   ('o-thr_ks.alpha_q1_diago_bse',      1, 'RPA (DFT energies)',               'silver',     '--', 1.3),
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
    w, y = data[:, 0], data[:, col]
    keep = (w >= EMIN) & (w <= EMAX)
    ax.plot(w[keep], y[keep], label=label, color=colour, ls=style, lw=width)
    top = max(top, y[keep].max())

for x, name in ((EXCITON, 'BSE A'), (GAP, 'QP gap')):
    if x is not None:
        ax.axvline(x, color='k', lw=0.7, ls=':')
        ax.text(x + 0.008, top * 1.04, f'{name}\n{x:.3f} eV', fontsize=8, va='top')

ax.set_xlim(EMIN, EMAX)
ax.set_ylim(0, top * 1.1)
ax.set_xlabel(r'$\omega$ (eV)')
ax.set_ylabel(r'Im $\alpha$ (bohr)')
ax.set_title(TITLE, fontsize=11)
ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=False)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
