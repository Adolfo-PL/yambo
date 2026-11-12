#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_bands.py -- band structures on one plot: ypp's interpolated bands (o-<job>.bands_interpolated) and
exciting's BAND*.OUT files (blocks of "k E[Ha]" separated by blank lines).

Run it in the folder with the ypp outputs:   python plot_bands.py
Everything you may want to change is in the SETTINGS block below.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ======================= SETTINGS: edit here =======================

TITLE = 'WS2 monolayer: GW vs DFT, Yambo vs exciting (no SOC)'
EMIN  = -4.0             # energy axis from (eV, relative to each set's mid-gap)
EMAX  = 5.0              # energy axis to
OUT   = 'bands.png'

# Positions of the high-symmetry points along the path (the |k| column of the file, a.u.) and their labels.
# For the path Gamma - M - K - Gamma of this cell (a = 3.1842 A): read from the file, where the k_x, k_y
# columns are (0,0), (0.5,0), (1/3,1/3), (0,0).
TICKS = [(0.0000, r'$\Gamma$'), (0.6029, 'M'), (0.9509, 'K'), (1.6470, r'$\Gamma$')]

# One line per band set: (file, format, first column, last column, shift, label, colour, line style, width)
#   format:  'ypp' (o-*.bands_interpolated, eV) or 'exciting' (BAND*.OUT, Hartree; columns are ignored)
#   columns: ypp files: column 0 is |k|; columns 1..14 are bands 7..20 (b7 ... b20), so 7 = band 13 (top
#            valence band) and 8 = band 14 (lowest conduction band)
#   shift:   subtracted from the energies (eV): the mid-gap of that set, (VBM + CBM)/2, so every set has
#            its gap centred on 0 and you see how GW opens it on both sides.
#            Yambo: VBM and CBM from ypp's report ("Valence Band Max", "Conduction Band Min"; GW: "[G+QP] ...")
#            exciting GW: VBM -1.3690 eV (K), CBM 1.3592 eV (Q, on K-Gamma), from its BAND file
#   Put a # in front of a line to leave that set out.
BANDS = [
    ('o-bands_dft.bands_interpolated', 'ypp',      1, 14, -0.067823, 'Yambo DFT (PBE)',      'gray',     '--', 1.1),
    ('o-bands_qp.bands_interpolated',  'ypp',      1, 14,  0.000240, 'Yambo GW (G0W0, PPA)', 'tab:blue', '-',  1.6),
    ('BAND_WANNIER_GW.OUT',            'exciting', 0,  0, -0.004892, 'exciting GW (Wannier)', 'tab:red', ':',  1.8),
]

# ===================================================================

HARTREE = 27.211386     # eV

fig, ax = plt.subplots(figsize=(6.5, 6.0))
for file, fmt, first, last, shift, label, colour, style, width in BANDS:
    try:
        lines = open(file).read().split('\n')
    except OSError:
        print(f'  not found, skipped: {file}')
        continue
    curves = []                                   # one (k, E) array per band, E in eV
    if fmt == 'ypp':
        data = np.loadtxt(file, comments='#', usecols=range(0, last + 1))
        curves = [(data[:, 0], data[:, col]) for col in range(first, last + 1)]
    else:                                         # exciting: blocks "k E[Ha]" separated by blank lines
        block = []
        for line in lines + ['']:
            if line.split():
                block.append([float(v) for v in line.split()[:2]])
            elif block:
                b = np.array(block)
                curves.append((b[:, 0], b[:, 1] * HARTREE))
                block = []
    for i, (k, e) in enumerate(curves):
        ax.plot(k, e - shift, color=colour, ls=style, lw=width, label=label if i == 0 else None)

for x, name in TICKS:
    ax.axvline(x, color='k', lw=0.6)
ax.axhline(0.0, color='k', lw=0.6, ls=':')
ax.set_xticks([x for x, name in TICKS])
ax.set_xticklabels([name for x, name in TICKS])
ax.set_xlim(TICKS[0][0], TICKS[-1][0])
ax.set_ylim(EMIN, EMAX)
ax.set_ylabel('E - E$_{mid-gap}$ (eV)')
ax.set_title(TITLE, fontsize=11)
ax.legend(fontsize=9, loc='center right')    # sits in the empty gap region
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
