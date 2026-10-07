#!/usr/bin/env python3
"""
xq_summary.py -- exciton energies and binding energies at every q: BSE, Casida with the
exported kernel solved self-consistently, and Casida with the static kernel fxc(q, w=0).

Inputs (written by xq_dispersion.sh):
  r-<bse job>*        report of the exporting BSE run: per momentum index the gap of the
                      window, the lowest excitations, the lowest bright one
  sc_q<iq>.dat        sc_scan.py table of that q (lambda(w) on undamped real frequencies)
  static_q<iq>.log    ktt_casida.py output for the z = 0 matrix of that q
  <npz>               fxc_export.py output of the same ndb.Chi (only for |q|)

For each q the gap is the lowest transition energy of the window, E_c(k) - E_v(k-q)
minimised over k, and the binding energy is gap - exciton energy.

Writes <out>.dat (table), <out>_omega.dat ("iq w": the self-consistent energies, input of
fxc_export.py --at) and <out>.png.

Usage
  python xq_summary.py --bse r-xq_bse_* --npz fxc_xq_static.npz [--dir .] [--out xq_dispersion]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sc_scan import crossings   # noqa: E402


def read_bse(paths):
    """Per momentum index: gap, lowest energies, their strengths, lowest bright."""
    out, iq = {}, None
    for p in paths:
        for line in open(p, errors='ignore'):
            if 'BSE momentum index' in line:
                iq = int(line.split()[-1])
                out[iq] = {}
            elif iq is None:
                continue
            elif 'lowest excitation energies' in line:
                out[iq]['E'] = [float(x) for x in line.split(':')[1].split()]
            elif 'their relative optical strength' in line:
                out[iq]['s'] = [float(x) for x in line.split(':')[1].split()]
            elif 'lowest bright excitation (' in line:
                out[iq]['bright'] = float(line.split()[-1])
            elif 'lowest transition energy of the window' in line:
                out[iq]['gap'] = float(line.split()[-1])
    return {k: v for k, v in out.items() if 'gap' in v and 'E' in v}


def read_static(path):
    """Lowest eigenvalue and lowest bright eigenvalue of M(0) from a ktt_casida.py log."""
    low = bright = None
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        if 'lowest Re E' in line and i + 1 < len(lines):
            low = float(lines[i + 1].split()[0])
        m = re.search(r'lowest bright:\s+([-\d.]+)', line)
        if m:
            bright = float(m.group(1))
    return low, bright


def nearest(xs, x0):
    return min(xs, key=lambda x: abs(x - x0)) if xs else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--bse', nargs='+', required=True, help='r- report(s) of the exporting BSE')
    ap.add_argument('--npz', help='fxc_export.py npz of the same ndb.Chi (for |q|)')
    ap.add_argument('--dir', default='.', help='where sc_q<iq>.dat and static_q<iq>.log are')
    ap.add_argument('--out', default='xq_dispersion')
    a = ap.parse_args()

    bse = read_bse(sorted(set(sum([glob.glob(p) for p in a.bse], []))))
    if not bse:
        sys.exit('no "[BSE]" exciton lines in ' + ' '.join(a.bse))
    qabs = {}
    if a.npz:
        z = np.load(a.npz)
        qabs = {int(i): float(p[0]) for i, p in zip(z['iq'], z['qpg'])}

    rows, omega = [], []
    print(' iq   |q|      gap    BSE low  BSE brt   sc low   sc brt  static brt   Eb BSE   Eb sc  Eb static  sc-BSE')
    for iq in sorted(bse):
        b = bse[iq]
        gap, blow, bbr = b['gap'], b['E'][0], b.get('bright', np.nan)
        sc_low = sc_br = st_br = np.nan
        f = os.path.join(a.dir, f'sc_q{iq}.dat')
        if os.path.exists(f):
            t = np.loadtxt(f, ndmin=2)
            if len(t) >= 2:
                sc_br = nearest(crossings(t[:, 0], t[:, 1])[0], bbr)
                sc_low = nearest(crossings(t[:, 0], t[:, 3])[0], blow)
        f = os.path.join(a.dir, f'static_q{iq}.log')
        if os.path.exists(f):
            st_br = read_static(f)[1] or np.nan
        qa = qabs.get(iq, np.nan)
        rows.append((iq, qa, gap, blow, bbr, sc_low, sc_br, st_br, gap - bbr, gap - sc_br, gap - st_br))
        if np.isfinite(sc_br):
            omega.append((iq, sc_br))
        print(f'{iq:3d} {qa:7.4f} {gap:8.4f} {blow:8.4f} {bbr:8.4f} {sc_low:8.4f} {sc_br:8.4f} {st_br:10.4f}'
              f' {gap - bbr:8.4f} {gap - sc_br:7.4f} {gap - st_br:9.4f} {1000 * (sc_br - bbr):+7.2f} meV')

    rows = np.array(rows)
    np.savetxt(a.out + '.dat', rows, fmt=['%4d'] + ['%10.6f'] * 10,
               header='energies in eV, |q| in bohr^-1; gap = lowest transition E_c(k)-E_v(k-q); '
                      'brt = lowest bright (> 1% strength); sc = Casida with fxc(q,w) solved '
                      'self-consistently; static = Casida with fxc(q,0); Eb = gap - bright\n'
                      'iq |q| gap BSE_low BSE_bright sc_low sc_bright static_bright '
                      'Eb_BSE Eb_sc Eb_static')
    if omega:
        np.savetxt(a.out + '_omega.dat', np.array(omega), fmt=['%4d', '%12.6f'],
                   header='iq  w [eV]: self-consistent bright exciton (input of fxc_export.py --at)')
    print(f'\nwrote {a.out}.dat' + (f', {a.out}_omega.dat' if omega else ''))
    d = rows[:, 6] - rows[:, 4]
    d = d[np.isfinite(d)]
    if len(d):
        print(f'self-consistent Casida vs BSE (bright): max |difference| {1000 * np.abs(d).max():.3f} meV '
              f'over {len(d)} q')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib missing: no plot')
        return
    x = rows[:, 1] if np.all(np.isfinite(rows[:, 1])) else rows[:, 0]
    xl = '|q| (bohr$^{-1}$)' if x is rows[:, 1] else 'q index'
    o = np.argsort(x)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].plot(x[o], rows[o, 2], 'k^', ms=4, label='gap (lowest transition)')
    ax[0].plot(x[o], rows[o, 4], 'o', mfc='none', ms=7, label='BSE')
    ax[0].plot(x[o], rows[o, 6], 'x', ms=6, label='Casida, fxc(q, $\\omega$) self-consistent')
    ax[0].plot(x[o], rows[o, 7], 's', mfc='none', ms=4, label='Casida, static fxc(q, 0)')
    ax[0].set_xlabel(xl)
    ax[0].set_ylabel('energy (eV)')
    ax[0].legend(fontsize=7)
    ax[1].plot(x[o], rows[o, 8], 'o', mfc='none', ms=7, label='BSE')
    ax[1].plot(x[o], rows[o, 9], 'x', ms=6, label='Casida, fxc(q, $\\omega$)')
    ax[1].plot(x[o], rows[o, 10], 's', mfc='none', ms=4, label='Casida, fxc(q, 0)')
    ax[1].set_xlabel(xl)
    ax[1].set_ylabel('binding energy (eV)')
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(a.out + '.png', dpi=150)
    print(f'wrote {a.out}.png')


if __name__ == '__main__':
    main()
