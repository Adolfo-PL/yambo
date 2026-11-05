#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
plot_casida_compare.py -- absorption and emission of the two Casida solvers on the same exact kernel.

Run in the folder with the files of (all optional; what exists is drawn):

  abs_coupled.dat, em_coupled.dat   ktt_casida.py ...           --spectrum --emission   (ours, coupled)
  abs_tda.dat, em_tda.dat           ktt_casida.py ... --tda     --spectrum --emission   (ours, TDA)
  volo_tda/                         ktt_casida.py ... --tda --volo volo_tda   (his formulas, TDA)
  volo_coupled/                     ktt_casida.py ...       --volo volo_coupled  (his formulas, coupled)
  --his DIR                         his program's own output (absorption_0, emission_0), plain input
  --his-fold DIR                    his program's own output on the ktt_to_kq.py --fold input
  --his-dip DIR                     output of the copy of his program with dipole weights
                                    (eigenvalues_SPN_dipoles.f90, 0.03 eV width), plain input: top row
  --his-dip-fold DIR                the same on the --fold input: top row

  python plot_casida_compare.py [--his ../TEST_CAS] [--his-fold kq_fold] [--his-dip kq_dip]
                                [--his-dip-fold kq_dip_fold] [--bse 2.3774] [--label WS2]
                                [--emin 2.2 --emax 3.2] [--out casida_compare.png]

Every curve is normalized to its maximum inside the plotted range. Upper panels: our formulas
(oscillator strengths, eta = 0.03 eV, Boltzmann emission); lower panels: his formulas (every state
weight 1/N, 0.1 eV width, Bose emission). Solid = coupled, dashed = TDA.
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(path):
    try:
        x = np.loadtxt(path, comments='#', ndmin=2)
        return x[:, 0], x[:, 1]
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--his', help="folder with his program's absorption_0 and emission_0 (plain, TDA input)")
    ap.add_argument('--his-fold', help="folder with his program's output on the --fold input")
    ap.add_argument('--his-dip', help="folder with the dipole-weighted copy of his program's output (plain input)")
    ap.add_argument('--his-dip-fold', help="folder with the dipole-weighted copy's output on the --fold input")
    ap.add_argument('--bse', type=float, default=None, help='BSE A exciton (eV), marked')
    ap.add_argument('--label', default='')
    ap.add_argument('--emin', type=float, default=2.2)
    ap.add_argument('--emax', type=float, default=3.2)
    ap.add_argument('--out', default='casida_compare.png')
    a = ap.parse_args()

    C_COUP, C_TDA, C_HIS = '#2a78d6', '#eb6834', '#1baf7a'
    rows = [  # (row, column, file, colour, style, width, label)
        (0, 0, 'abs_coupled.dat', C_COUP, '-', 1.8, 'ours, coupled'),
        (0, 0, 'abs_tda.dat', C_TDA, '--', 1.8, 'ours, TDA'),
        (0, 1, 'em_coupled.dat', C_COUP, '-', 1.8, 'ours, coupled'),
        (0, 1, 'em_tda.dat', C_TDA, '--', 1.8, 'ours, TDA'),
        (1, 0, 'volo_coupled/absorption_0', C_COUP, '-', 1.8, 'his formulas, coupled (our run)'),
        (1, 0, 'volo_tda/absorption_0', C_TDA, '--', 1.8, 'his formulas, TDA (our run)'),
        (1, 1, 'volo_coupled/emission_0', C_COUP, '-', 1.8, 'his formulas, coupled (our run)'),
        (1, 1, 'volo_tda/emission_0', C_TDA, '--', 1.8, 'his formulas, TDA (our run)'),
    ]
    if a.his:
        rows += [(1, 0, os.path.join(a.his, 'absorption_0'), C_HIS, ':', 2.6, 'his program, TDA input'),
                 (1, 1, os.path.join(a.his, 'emission_0'), C_HIS, ':', 2.6, 'his program, TDA input')]
    if a.his_fold:
        rows += [(1, 0, os.path.join(a.his_fold, 'absorption_0'), '#0d366b', ':', 2.6, 'his program, --fold input'),
                 (1, 1, os.path.join(a.his_fold, 'emission_0'), '#0d366b', ':', 2.6, 'his program, --fold input')]

    # the dipole-weighted copy of his program has our absorption formula: top row (its emission keeps
    # his Bose factor)
    if a.his_dip:
        rows += [(0, 0, os.path.join(a.his_dip, 'absorption_0'), C_HIS, ':', 2.6, 'his program + dipoles, TDA input'),
                 (0, 1, os.path.join(a.his_dip, 'emission_0'), C_HIS, ':', 2.6, 'his program + dipoles (Bose), TDA input')]
    if a.his_dip_fold:
        rows += [(0, 0, os.path.join(a.his_dip_fold, 'absorption_0'), '#0d366b', ':', 2.6,
                  'his program + dipoles, --fold input'),
                 (0, 1, os.path.join(a.his_dip_fold, 'emission_0'), '#0d366b', ':', 2.6,
                  'his program + dipoles (Bose), --fold input')]

    fig, ax = plt.subplots(2, 2, figsize=(12, 7.2), sharex=True)
    drawn = 0
    for r, c, f, col, ls, lw, lab in rows:
        d = load(f)
        if d is None:
            print(f'  (no {f}: skipped)')
            continue
        w, y = d
        m = (w >= a.emin) & (w <= a.emax)
        if not m.any() or y[m].max() <= 0:
            continue
        ax[r, c].plot(w[m], y[m] / y[m].max(), ls=ls, color=col, lw=lw, label=lab)
        drawn += 1
    titles = [['absorption, our formulas (oscillator strengths)', 'emission, our formulas (Boltzmann)'],
              ['absorption, his formulas (every state weight 1/N)', 'emission, his formulas (Bose)']]
    for r in range(2):
        for c in range(2):
            x = ax[r, c]
            if a.bse:
                x.axvline(a.bse, color='grey', lw=0.8)
                x.text(a.bse + 0.005, 0.92, 'BSE A', color='grey', fontsize=8)
            x.set_title(titles[r][c], fontsize=10)
            x.set_xlim(a.emin, a.emax)
            x.set_ylim(0, 1.08)
            x.set_ylabel('normalized')
            if x.lines:
                x.legend(fontsize=7, loc='upper right')
        ax[1, c].set_xlabel('ω (eV)')
    if a.label:
        fig.suptitle(f'{a.label}: Casida with the exact kernel, coupled (solid) vs TDA (dashed)', fontsize=12)
    fig.tight_layout()
    fig.savefig(a.out, dpi=170)
    print(f'wrote {a.out} ({drawn} curves)')


if __name__ == '__main__':
    main()
