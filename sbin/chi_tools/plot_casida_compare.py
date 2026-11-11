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
  --his DIR                         his program's own output (absorption_0, emission_0), plain input
  --his-fold DIR                    his program's own output on the ktt_to_kq.py --fold input
  --his-dip DIR                     output of the copy of his program with dipole weights
                                    (eigenvalues_SPN_dipoles.f90, 0.03 eV width), plain input: top row
  --his-dip-fold DIR                the same on the --fold input: top row

  --bse-spec FILE                    the full BSE spectrum (Yambo o-*.alpha/eps file, e.g. prod/o-bsp.alpha_q1_*
                                    of exact_kernel.sh STAGE=spectra): absorption panels, both rows
  --alda FILE                       TDDFT-ALDA (prod/o-tal.* of STAGE=tddft): absorption panels, both rows
  --lrc FILE[:LABEL]                TDDFT-LRC (prod/o-tlr<alpha>.* of STAGE=tddft), repeatable; the label
                                    defaults to the alpha in the file name

  python plot_casida_compare.py [--his ../TEST_CAS] [--his-fold kq_fold] [--his-dip kq_dip]
                                [--his-dip-fold kq_dip_fold] [--bse 2.3774] [--label WS2]
                                [--emin 2.2 --emax 3.2] [--out casida_compare.png] [--common-scale]
                                [--bse-spec o-bsp.alpha_q1_...] [--alda o-tal.alpha_q1_...]
                                [--lrc o-tlr2.alpha_q1_... --lrc o-tlr0.6860.alpha_q1_...:"LRC, exact head"]

Every curve is normalized to its maximum inside the plotted range; with --common-scale the Yambo TDA
curves are divided by the Yambo coupled maximum instead (same formula and constant: real relative
heights), his curves still by their own maximum (his units differ). Both rows show the two Yambo-kernel
results (ktt_casida.py: oscillator strengths, eta = 0.03 eV, Boltzmann emission; solid = coupled,
dashed = TDA); the upper row adds his program with the fixes (--his-dip, --his-dip-fold), the lower row
his original program (--his, --his-fold). The Yambo spectra (--bse-spec, --alda, --lrc: column 2, Im alpha
or Im eps, same broadening as STAGE=spectra) go in both absorption panels, each normalized to its own maximum
(their units differ from the Casida curves); they have no emission counterpart.
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
    ap.add_argument('--common-scale', action='store_true',
                    help='scale the Yambo TDA curves by the Yambo coupled maximum (same formula and constant), '
                         'so their real relative heights show; his curves stay normalized to their own maximum')
    ap.add_argument('--bse-spec', help='full BSE spectrum (Yambo o-*.alpha/eps file): absorption, both rows')
    ap.add_argument('--alda', help='TDDFT-ALDA spectrum (Yambo o-tal.* file): absorption, both rows')
    ap.add_argument('--lrc', action='append', default=[], metavar='FILE[:LABEL]',
                    help='TDDFT-LRC spectrum (Yambo o-tlr<alpha>.* file), repeatable')
    a = ap.parse_args()

    C_COUP, C_TDA, C_HIS = '#2a78d6', '#eb6834', '#1baf7a'
    rows = []  # (row, column, file, colour, style, width, label)
    # both rows: the two Yambo-kernel results with the oscillator strengths (ktt_casida.py)
    for r in (0, 1):
        rows += [(r, 0, 'abs_coupled.dat', C_COUP, '-', 1.8, 'TDDFT:Yambo, coupled'),
                 (r, 0, 'abs_tda.dat', C_TDA, '--', 1.8, 'TDDFT:Yambo, TDA'),
                 (r, 1, 'em_coupled.dat', C_COUP, '-', 1.8, 'TDDFT:Yambo, coupled'),
                 (r, 1, 'em_tda.dat', C_TDA, '--', 1.8, 'TDDFT:Yambo, TDA')]
    if a.his:
        rows += [(1, 0, os.path.join(a.his, 'absorption_0'), C_HIS, ':', 2.6, 'TDDFT:Volo, (TDA)'),
                 (1, 1, os.path.join(a.his, 'emission_0'), C_HIS, ':', 2.6, 'TDDFT:Volo, (TDA)')]
    if a.his_fold:
        rows += [(1, 0, os.path.join(a.his_fold, 'absorption_0'), '#0d366b', ':', 2.6, 'TDDFT:Volo, (--fold input)'),
                 (1, 1, os.path.join(a.his_fold, 'emission_0'), '#0d366b', ':', 2.6, 'TDDFT:Volo, (--fold input)')]

    # the dipole-weighted copy of his program has our absorption formula: top row (its emission keeps
    # his Bose factor)
    if a.his_dip:
        rows += [(0, 0, os.path.join(a.his_dip, 'absorption_0'), C_HIS, ':', 2.6, 'TDDFT:Volo + Yambo dipoles, (TDA)'),
                 (0, 1, os.path.join(a.his_dip, 'emission_0'), C_HIS, ':', 2.6, 'TDDFT:Volo + Yambo dipoles, (TDA)')]
    if a.his_dip_fold:
        rows += [(0, 0, os.path.join(a.his_dip_fold, 'absorption_0'), '#0d366b', ':', 2.6,
                  'TDDFT:Volo + Yambo dipoles, (--fold input)'),
                 (0, 1, os.path.join(a.his_dip_fold, 'emission_0'), '#0d366b', ':', 2.6,
                  'TDDFT:Volo + Yambo dipoles, (--fold input)')]

    # Yambo's own spectra (full BSE, ALDA, LRC): absorption panels of both rows, own normalization
    import re
    ref_curves = []
    if a.bse_spec:
        ref_curves.append((a.bse_spec, 'black', '-', 1.2, 'BSE (Yambo, full spectrum)'))
    if a.alda:
        ref_curves.append((a.alda, '#8a3ffc', '-.', 1.6, 'TDDFT:Yambo ALDA'))
    lrc_cols = ['#a0522d', '#6b7280', '#c71e5d', '#4a7c59']
    for i, item in enumerate(a.lrc):
        f, _, lab = item.partition(':') if not os.path.exists(item) else (item, '', '')
        if not lab:
            m = re.search(r'tlr([-\d.]+?)\.(alpha|eps)', os.path.basename(f))
            lab = f'TDDFT:Yambo LRC, α = {m.group(1)}' if m else 'TDDFT:Yambo LRC'
        ref_curves.append((f, lrc_cols[i % len(lrc_cols)], (0, (1, 1.5)), 1.8, lab))
    for r in (0, 1):
        rows += [(r, 0, f, col, ls, lw, lab) for f, col, ls, lw, lab in ref_curves]

    # --common-scale: the Yambo curves of one quantity share the coupled curve's maximum as their scale
    ref = {'abs_coupled.dat': 'abs_coupled.dat', 'abs_tda.dat': 'abs_coupled.dat',
           'em_coupled.dat': 'em_coupled.dat', 'em_tda.dat': 'em_coupled.dat'}

    def peak(f):
        d = load(f)
        if d is None:
            return None
        m = (d[0] >= a.emin) & (d[0] <= a.emax)
        return d[1][m].max() if m.any() and d[1][m].max() > 0 else None

    fig, ax = plt.subplots(2, 2, figsize=(12, 7.2), sharex=True)
    top = np.ones((2, 2))
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
        scale = y[m].max()
        if a.common_scale and f in ref and peak(ref[f]):
            scale = peak(ref[f])
        ax[r, c].plot(w[m], y[m] / scale, ls=ls, color=col, lw=lw, label=lab)
        top[r, c] = max(top[r, c], (y[m] / scale).max())
        drawn += 1
        if a.common_scale and f in ref and f != ref[f]:
            print(f'  {f}: maximum {(y[m] / scale).max():.3f} of the coupled one ({ref[f]})')
    titles = [['absorption: Yambo vs Volo with the fixes (dipoles)', 'emission: Yambo vs Volo with the fixes'],
              ['absorption: Yambo vs Volo, original program', 'emission: Yambo vs Volo, original program']]
    for r in range(2):
        for c in range(2):
            x = ax[r, c]
            if a.bse:
                x.axvline(a.bse, color='grey', lw=0.8)
                x.text(a.bse + 0.005, 0.92, 'BSE A', color='grey', fontsize=8)
            x.set_title(titles[r][c], fontsize=10)
            x.set_xlim(a.emin, a.emax)
            x.set_ylim(0, 1.08 * top[r, c])
            x.set_ylabel('normalized' if not a.common_scale else 'relative height')
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
