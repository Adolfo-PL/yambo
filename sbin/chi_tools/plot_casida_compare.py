#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_casida_compare.py -- absorption and emission of the two Casida solvers on the same exact kernel
(ours, from Yambo's matrices, and Volodymyr's k-space program), with Yambo's own spectra for reference.

Run it in prod/:   python plot_casida_compare.py
Everything you may want to change is in the SETTINGS block below.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ======================= SETTINGS: edit here =======================

TITLE = 'WS2: Casida with the exact kernel, coupled (solid) vs TDA (dashed)'
EMIN  = 2.2              # x axis from (eV), all panels
EMAX  = 3.2              # x axis to (eV), all panels
BSE_A = 2.3774           # grey line "BSE A" in every panel (None to hide it)
OUT   = 'casida_compare.png'

# Each curve is divided by its own maximum inside EMIN..EMAX (their units differ). To show the real
# height of a curve relative to another one, put that other file in the last field ("scale by"):
# e.g. 'abs_coupled.dat' on the TDA line shows the TDA peak relative to the coupled one.

# The four panels: name -> title. LAYOUT places them (rows, columns).
PANELS = {
    'abs_fix':  'absorption: Yambo vs Volo with the fixes (dipoles)',
    'em_fix':   'emission: Yambo vs Volo with the fixes',
    'abs_orig': 'absorption: Yambo vs Volo, original program',
    'em_orig':  'emission: Yambo vs Volo, original program',
}
LAYOUT = [['abs_fix',  'em_fix'],
          ['abs_orig', 'em_orig']]

# One line per curve: (panel, file, column, label, colour, line style, line width, scale by)
#   column 1 for all these files (ours: abs_*.dat/em_*.dat; Volo: absorption_0/emission_0;
#   Yambo o-* files: 1 = Im alpha)
#   Put a # in front of a line to leave that curve out.
CURVES = [
    # --- top row: ours and Volo's program with the Yambo dipoles ---
    ('abs_fix', 'abs_coupled.dat',                 1, 'TDDFT: Yambo, coupled',          '#2a78d6', '-',  1.8, None),
    ('abs_fix', 'abs_tda.dat',                     1, 'TDDFT: Yambo, TDA',              '#eb6834', '--', 1.8, None),
    ('abs_fix', 'kq_dip/absorption_0',             1, 'TDDFT: Volo + Yambo dipoles (TDA)', '#1baf7a', ':', 2.6, None),
    ('abs_fix', 'o-bsp.alpha_q1_diago_bse',        1, 'BSE (Yambo)',                    'black',   '-',  1.2, None),
#   ('abs_fix', 'o-tal.alpha_q1_diago_alda_bse',   1, 'TDDFT: ALDA (Yambo)',            '#8a3ffc', '-.', 1.6, None),
    ('em_fix',  'em_coupled.dat',                  1, 'TDDFT: Yambo, coupled',          '#2a78d6', '-',  1.8, None),
    ('em_fix',  'em_tda.dat',                      1, 'TDDFT: Yambo, TDA',              '#eb6834', '--', 1.8, None),
    ('em_fix',  'kq_dip/emission_0',               1, 'TDDFT: Volo + Yambo dipoles (TDA)', '#1baf7a', ':', 2.6, None),
    # --- bottom row: ours and Volo's original program ---
    ('abs_orig', 'abs_coupled.dat',                1, 'TDDFT: Yambo, coupled',          '#2a78d6', '-',  1.8, None),
    ('abs_orig', 'abs_tda.dat',                    1, 'TDDFT: Yambo, TDA',              '#eb6834', '--', 1.8, None),
    ('abs_orig', '../VOLO/absorption_0',           1, 'TDDFT: Volo (TDA)',              '#1baf7a', ':',  2.6, None),
    ('abs_orig', 'o-bsp.alpha_q1_diago_bse',       1, 'BSE (Yambo)',                    'black',   '-',  1.2, None),
    ('em_orig',  'em_coupled.dat',                 1, 'TDDFT: Yambo, coupled',          '#2a78d6', '-',  1.8, None),
    ('em_orig',  'em_tda.dat',                     1, 'TDDFT: Yambo, TDA',              '#eb6834', '--', 1.8, None),
    ('em_orig',  '../VOLO/emission_0',             1, 'TDDFT: Volo (TDA)',              '#1baf7a', ':',  2.6, None),
]

# ===================================================================


def load(file, col):
    """x and the chosen column inside EMIN..EMAX, or None if the file is missing."""
    try:
        data = np.loadtxt(file, comments='#')
    except OSError:
        return None
    keep = (data[:, 0] >= EMIN) & (data[:, 0] <= EMAX)
    return data[keep, 0], data[keep, col]


nrow, ncol = len(LAYOUT), len(LAYOUT[0])
fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 3.6 * nrow), sharex=True, squeeze=False)
where = {name: axes[r][c] for r, row in enumerate(LAYOUT) for c, name in enumerate(row)}

for panel, file, col, label, colour, style, width, scale_by in CURVES:
    d = load(file, col)
    if d is None:
        print(f'  not found, skipped: {file}')
        continue
    x, y = d
    ref = load(scale_by, 1) if scale_by else d
    height = ref[1].max()
    if scale_by:
        print(f'  {file}: maximum {y.max() / height:.3f} of {scale_by}')
    where[panel].plot(x, y / height, label=label, color=colour, ls=style, lw=width)

for name, ax in where.items():
    if BSE_A is not None:
        ax.axvline(BSE_A, color='grey', lw=0.8)
        ax.text(BSE_A + 0.005, 0.95, 'BSE A', color='grey', fontsize=8, transform=ax.get_xaxis_transform())
    ax.set_title(PANELS[name], fontsize=10)
    ax.set_xlim(EMIN, EMAX)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('normalized')
    if ax.lines:
        ax.legend(fontsize=7, loc='upper right')
for ax in axes[-1]:
    ax.set_xlabel(r'$\omega$ (eV)')
fig.suptitle(TITLE, fontsize=12)
fig.tight_layout()
fig.savefig(OUT, dpi=170)
print('wrote', OUT)
