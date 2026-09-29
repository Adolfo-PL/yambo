#!/usr/bin/env python3
"""
sc_scan.py -- self-consistent exciton energy from a scan of Casida matrices M(w).

Input: the BSEChiKoutW files of a BSKmod="CHI" + BSEChiDyn run on UNDAMPED real
frequencies (BDmRange 0 | 0, kernel exported by a BSEChiOut run on the same grid):

  o-<job>.Ktt_q1_transitions  and  o-<job>.Ktt_q1_w<iw>, one per frequency w

At every w the script diagonalizes M(w) and records the lowest bright eigenvalue
lambda(w). The exact kernel reproduces the BSE exciton E_x only where

    lambda(E_x) = E_x                      (as the quasiparticle equation E = e + Sigma(E))

so the script looks for the crossing of Re lambda(w) - w and interpolates it linearly.
The whole transition space must be in the files (BSEChiKoutB 0 | 0): a sub-matrix does
not have the BSE exciton as an eigenvalue.

Output: a table  w  Re/Im lambda_bright  Re lambda_lowest  (eV) and the crossing(s).

Usage
  python sc_scan.py o-sc_cas.Ktt_q1_transitions o-sc_cas.Ktt_q1_w* --out sc.dat
         [--exciton 2.2342] [--delete] [--bright 0.01] [--tda]

--delete removes each w-file once it is read (they are ~1 GB each for 2000 transitions);
--append merges the rows with those already in --out (a scan split into batches).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ktt_casida as kc   # noqa: E402


def crossings(w, lam, jump=10.0):
    """Zeros of lam(w) - w by linear interpolation, and the sign changes that are jumps.

    Near a pole of the kernel lambda(w) changes branch: a sign change with
    |d lambda / d w| > jump between two grid points is such a jump, not a solution."""
    g = lam - w
    out, poles = [], []
    for j in range(len(w) - 1):
        if g[j] == 0.0:
            out.append(w[j])
        elif g[j] * g[j + 1] < 0.0:
            x = w[j] - g[j] * (w[j + 1] - w[j]) / (g[j + 1] - g[j])
            steep = abs(lam[j + 1] - lam[j]) > jump * abs(w[j + 1] - w[j])
            (poles if steep else out).append(x)
    return out, poles


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('transitions')
    ap.add_argument('wfiles', nargs='+')
    ap.add_argument('--out', default='sc_scan.dat')
    ap.add_argument('--exciton', type=float, help='BSE exciton energy [eV] to compare with')
    ap.add_argument('--bright', type=float, default=0.01, help='relative strength of a bright state')
    ap.add_argument('--tda', action='store_true')
    ap.add_argument('--delete', action='store_true', help='remove each w-file after reading it')
    ap.add_argument('--append', action='store_true',
                    help='merge with the rows already in --out (a scan run in several batches)')
    a = ap.parse_args()

    _, tr = kc.read_table(a.transitions)
    n = len(tr)
    E = tr[:, 9]
    res = tr[:, 11] + 1j * tr[:, 12]
    pairs = {(int(v), int(c)) for v, c in tr[:, 7:9]}
    print(f'{n} transitions, {len(pairs)} band pair(s), lowest transition {E.min():.4f} eV')
    if len(pairs) == 1:
        print('  !!! one band pair only: the crossing is not the BSE exciton (use BSEChiKoutB 0 | 0)')

    files = sorted(a.wfiles, key=lambda f: int(f.rsplit('_w', 1)[1]))
    rows = []
    print('     w [eV]   Re lam_bright  Im lam_bright   Re lam_lowest   lam_bright - w')
    for f in files:
        z, M, coupled = kc.read_matrix(f, n, E, a.tda, False)
        if a.delete:
            os.remove(f)
        if abs(z.imag) > 1e-6:
            print(f'  !!! {f}: z has damping {z.imag:.3g} eV; the crossing needs undamped real w')
        if coupled:
            A = np.concatenate([res, 1j * np.conj(res)])
            B = np.concatenate([np.conj(res), 1j * res])
        else:
            A, B = res, np.conj(res)
        lam, R = np.linalg.eig(M)
        del M
        st = np.abs((B @ R) * np.linalg.solve(R, A))
        del R
        pos = np.where(lam.real > 1e-6)[0]
        pos = pos[np.argsort(lam[pos].real)]
        srel = st[pos] / max(st[pos].max(), 1e-300)
        b = pos[np.argmax(srel > a.bright)]
        rows.append((z.real, lam[b].real, lam[b].imag, lam[pos[0]].real))
        print(f'  {z.real:9.4f}  {lam[b].real:13.5f}  {lam[b].imag:13.2e}  {lam[pos[0]].real:14.5f}'
              f'  {lam[b].real - z.real:+14.5f}', flush=True)

    rows = np.array(rows).reshape(-1, 4)
    if a.append and os.path.exists(a.out):
        old = np.loadtxt(a.out, ndmin=2)
        keep = ~np.isin(np.round(old[:, 0], 6), np.round(rows[:, 0], 6))
        rows = np.vstack([old[keep], rows])
    rows = rows[np.argsort(rows[:, 0])]
    np.savetxt(a.out, rows, fmt='%.6f',
               header='omega_eV lambda_bright_re lambda_bright_im lambda_lowest_re')
    print(f'\ntable written to {a.out}')
    xb, pb = crossings(rows[:, 0], rows[:, 1])
    xl, _ = crossings(rows[:, 0], rows[:, 3])
    ref = f'   (BSE {a.exciton:.4f} eV)' if a.exciton else ''
    if xb:
        print(f'lambda_bright(w) = w at  ' + ', '.join(f'{x:.4f}' for x in xb) + ' eV' + ref)
        if a.exciton:
            print(f'  self-consistent - BSE = {min(xb, key=lambda x: abs(x - a.exciton)) - a.exciton:+.4f} eV')
    else:
        print('no crossing of lambda_bright(w) = w in the scanned window: widen or move it')
    if pb:
        print('lambda_bright jumps across w near ' + ', '.join(f'{x:.3f}' for x in pb) +
              ' eV: a pole of fxc(w), not a solution')
    if xl:
        print(f'lambda_lowest(w) = w at  ' + ', '.join(f'{x:.4f}' for x in xl) + ' eV (dark or bright)')


if __name__ == '__main__':
    main()
