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

Run in prod/ (or give --dir) after exact_kernel.sh STAGE=spectra and STAGE=tddft. All curves are Im alpha
(2D polarizability, bohr) on one common axis, no per-curve normalization:

  BSE                         casida_spectra.dat, column BSE (Yambo, o-bsp)
  TDDFT, exact fxc at E_A     casida_spectra.dat, column fxc(E_A): Casida, coupled, kernel frozen at the A
                              exciton (scaled to Yambo's units by casida_spectra.py)
  TDDFT, exact fxc(w)         casida_spectra.dat, column fxc(w) at every w (--dyn)
  ALDA (QP energies)          o-tal.alpha_q1_*
  RPA  (QP energies)          o-thr.alpha_q1_*
  ALDA, RPA (DFT energies)    o-tal_ks.alpha_q1_*, o-thr_ks.alpha_q1_* (--ks)

Missing files are skipped. The QP gap and the A exciton are marked when given (--gap, --exciton).

  python plot_tddft_kernels.py [--dir prod] [--emin 2.2 --emax 3.3] [--ks] [--dyn] [--label WS2]
                               [--exciton 2.3774] [--gap 2.9238] [--out tddft_kernels.png]
"""
import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def yambo_file(d, job):
    f = sorted(glob.glob(os.path.join(d, f'o-{job}.alpha_q1_*')) or glob.glob(os.path.join(d, f'o-{job}.eps_q1_*')))
    if not f:
        return None
    x = np.loadtxt(f[0], comments='#', ndmin=2)
    return x[:, 0], x[:, 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--dir', default='.', help='the prod/ folder')
    ap.add_argument('--emin', type=float, default=2.2)
    ap.add_argument('--emax', type=float, default=3.3)
    ap.add_argument('--ks', action='store_true', help='also ALDA and RPA on DFT energies')
    ap.add_argument('--dyn', action='store_true', help='also the exact kernel fxc(w) at every w')
    ap.add_argument('--exciton', type=float, help='A exciton of the BSE [eV] (marked)')
    ap.add_argument('--gap', type=float, help='QP gap [eV] (marked)')
    ap.add_argument('--label', default='')
    ap.add_argument('--out', default='tddft_kernels.png')
    a = ap.parse_args()

    curves = []     # (w, y, style, label)
    cs = os.path.join(a.dir, 'casida_spectra.dat')
    if os.path.exists(cs):
        x = np.loadtxt(cs, comments='#', ndmin=2)
        curves.append((x[:, 0], x[:, 1], dict(color='#1f1f1f', lw=4, alpha=0.30), 'BSE'))
        curves.append((x[:, 0], x[:, 5], dict(color='#2a78d6', lw=1.8), r'TDDFT, exact $f_{xc}$ at $E_A$ (Casida)'))
        if a.dyn and x.shape[1] > 6:
            curves.append((x[:, 0], x[:, 6], dict(color='#0d366b', lw=1.2, ls='--'),
                           r'TDDFT, exact $f_{xc}(\omega)$ at every $\omega$'))
    else:
        print(f'  (no {cs}: BSE and exact-kernel curves skipped; run STAGE=spectra)')
    for job, st, lab in (('tal', dict(color='#8a3ffc', lw=1.6, ls='-.'), 'TDDFT, ALDA'),
                         ('thr', dict(color='#6b7280', lw=1.4, ls='--'), 'RPA (no kernel)')):
        d = yambo_file(a.dir, job)
        if d:
            curves.append((d[0], d[1], st, lab))
        else:
            print(f'  (no o-{job}.alpha_q1_*: skipped)')
    if a.ks:
        for job, st, lab in (('tal_ks', dict(color='#c4a6ff', lw=1.3, ls='-.'), 'TDDFT, ALDA (DFT energies)'),
                             ('thr_ks', dict(color='#b8bcc4', lw=1.3, ls='--'), 'RPA (DFT energies)')):
            d = yambo_file(a.dir, job)
            if d:
                curves.append((d[0], d[1], st, lab))
            else:
                print(f'  (no o-{job}.alpha_q1_*: skipped)')

    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    top = 0.0
    for w, y, st, lab in curves:
        m = (w >= a.emin) & (w <= a.emax)
        if not m.any():
            continue
        ax.plot(w[m], y[m], label=lab, **st)
        top = max(top, y[m].max())
    for xv, lab in ((a.exciton, 'BSE A'), (a.gap, 'QP gap')):
        if xv and a.emin <= xv <= a.emax:
            ax.axvline(xv, color='k', lw=0.7, ls=':')
            ax.text(xv + 0.008, top * 1.04, f'{lab}\n{xv:.3f} eV', fontsize=8, va='top')
    ax.set_xlim(a.emin, a.emax)
    ax.set_ylim(0, top * 1.1 if top > 0 else 1)
    ax.set_xlabel(r'$\omega$ (eV)')
    ax.set_ylabel(r'Im $\alpha$ (bohr)')
    ax.set_title(f'{a.label + ": " if a.label else ""}absorption, BSE vs TDDFT kernels' + ('' if a.ks else ' (QP energies)'), fontsize=11)
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0, frameon=False)
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)
    print(f'wrote {a.out} ({len(curves)} curves)')


if __name__ == '__main__':
    main()
