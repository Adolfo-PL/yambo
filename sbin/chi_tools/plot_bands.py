#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_bands.py -- band structure from ypp's interpolated bands (o-<job>.bands_interpolated), e.g. GW and DFT.

Run it in the folder with the ypp outputs:   python plot_bands.py
Everything you may want to change is in the SETTINGS block below.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ======================= SETTINGS: edit here =======================

TITLE = 'WS2 monolayer, 18x18: GW vs DFT (no SOC)'
EMIN  = -4.0             # energy axis from (eV, relative to each set's VBM)
EMAX  = 5.0              # energy axis to
OUT   = 'bands.png'

# Positions of the high-symmetry points along the path (the |k| column of the file, a.u.) and their labels.
# For the path Gamma - M - K - Gamma of this cell (a = 3.1842 A): read from the file, where the k_x, k_y
# columns are (0,0), (0.5,0), (1/3,1/3), (0,0).
TICKS = [(0.0000, r'$\Gamma$'), (0.6029, 'M'), (0.9509, 'K'), (1.6470, r'$\Gamma$')]

# One line per band set: (file, first column, last column, shift, label, colour, line style, width)
#   columns: column 0 is |k|; columns 1..14 are bands 7..20 (b7 ... b20), so 7 = band 13 (top valence
#            band) and 8 = band 14 (lowest conduction band)
#   shift:   subtracted from the energies, eV: the VBM of that set, so both start from 0 at the VBM
#            (ypp's report: DFT "Valence Band Max", GW "[G+QP] Valence Band Max")
#   Put a # in front of a line to leave that set out.
BANDS = [
    ('o-bands_dft.bands_interpolated', 1, 14, -0.977324, 'DFT (PBE)', 'gray',     '--', 1.2),
    ('o-bands_qp.bands_interpolated',  1, 14, -1.408245, 'GW (G0W0, PPA)', 'tab:blue', '-',  1.6),
]

# ===================================================================

fig, ax = plt.subplots(figsize=(6.5, 6.0))
for file, first, last, shift, label, colour, style, width in BANDS:
    try:
        data = np.loadtxt(file, comments='#', usecols=range(0, last + 1))
    except OSError:
        print(f'  not found, skipped: {file}')
        continue
    k = data[:, 0]
    for col in range(first, last + 1):
        ax.plot(k, data[:, col] - shift, color=colour, ls=style, lw=width,
                label=label if col == first else None)

for x, name in TICKS:
    ax.axvline(x, color='k', lw=0.6)
ax.axhline(0.0, color='k', lw=0.6, ls=':')
ax.set_xticks([x for x, name in TICKS])
ax.set_xticklabels([name for x, name in TICKS])
ax.set_xlim(TICKS[0][0], TICKS[-1][0])
ax.set_ylim(EMIN, EMAX)
ax.set_ylabel('E - E$_{VBM}$ (eV)')
ax.set_title(TITLE, fontsize=11)
ax.legend(fontsize=9, loc='upper right')
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
