#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""Check a retarded ndb.G (GRetarded) against the conditions ChiGMode="DYSON" relies on.

Per state:
  causality   A(w) = -Im G(w)/pi >= 0 on the whole grid, as the loader: dips down to 1e-2 of the peak,
              adding up to 1e-3 of the norm,
              (PPA Im Sigma > 0 at isolated frequencies) are set to zero and reported
  norm        trapezoid sum of A, as in Chi_spectral_weights (the loader stops if |norm-1| > --norm-tol)
  Sc(Eo)      Re Sigma_c at the KS energy, Sigma_c = SE_Operator - G_SIGMA_X + G_VXC_REFERENCE,
              interpolated on the grid (whose points carry +i*GDmRnge)
With --newton o-<job>.qp (a Newton run of the same states with ExtendOut and GDamping = GDmRnge):
  Re Sc(Eo) must equal its Sc|Eo column. Both are the same PPA sum at w = Eo + i*eta, and the real
  part is the same for retarded and time-ordered denominators, so they differ only by round-off.
With --ref <ndb.G of a second run> and --ref-newton <its o-*.qp> (e.g. the run without terminator):
  the change of Re Sc(Eo) between the two ndb.G must equal the change of Sc|Eo between the Newton runs
  (this isolates the terminator contribution on the grid).
Prints one line per state; exits 1 if any check fails.

  gf_check.py gw/ndb.G [--newton o-gw.qp] [--ref ref/ndb.G --ref-newton o-ref.qp] [--norm-tol 0.05] [--sc-tol 0.002]
"""
import argparse
import re
import sys

import numpy as np
from netCDF4 import Dataset

HA2EV = 27.211386245988
CLIP_REL = 1e-3    # Chi_clip_rel in src/modules/mod_Chi_Green.F: largest clipped weight (of the norm)
CLIP_PEAK = 1e-2   # Chi_clip_peak: deepest clipped dip (of the peak)


def cplx(ds, name):
    """A Yambo complex array. Yambo defines netCDF variables with the Fortran VAR_SZ order
    (2, n[, steps]), so Python reads it reversed, ([steps,] n, 2). Returns complex (n[, steps])."""
    a = np.array(ds.variables[name][:], dtype=float)
    return (a[..., 0] + 1j * a[..., 1]).T


def load_g(path):
    with Dataset(path) as ds:
        if 'G_RETARDED' not in ds.variables:
            sys.exit(f'!!! {path}: no G_RETARDED (not a GRetarded ndb.G)')
        eo = np.array(ds.variables['G_REFERENCE_ENERGIES'][:], dtype=float).ravel()
        n = eo.size
        t = np.array(ds.variables['QP_table'][:]).T     # Fortran (n, 3+n_sp_pol-1), read reversed
        if t.shape[0] != n:
            sys.exit(f'!!! {path}: QP_table has {t.shape[0]} states, G_REFERENCE_ENERGIES {n}')
        g = {'eo': eo * HA2EV, 'band': t[:, 0].astype(int), 'k': t[:, 2].astype(int),
             'w': cplx(ds, 'Green_Functions_Energies') * HA2EV,
             'G': cplx(ds, 'Green_Functions') / HA2EV,
             'se': cplx(ds, 'SE_Operator') * HA2EV,
             'sx': cplx(ds, 'G_SIGMA_X') * HA2EV,
             'vxc': cplx(ds, 'G_VXC_REFERENCE') * HA2EV}
    g['sc'] = g['se'] - g['sx'][:, None] + g['vxc'][:, None]
    return g


def sc_at_eo(g, i):
    e = g['w'][i].real
    return np.interp(g['eo'][i], e, g['sc'][i].real), np.min(np.abs(e - g['eo'][i]))


def norm_and_minA(g, i):
    e = g['w'][i].real
    a = -g['G'][i].imag / np.pi
    de = np.empty_like(e)
    de[0] = (e[1] - e[0]) / 2
    de[-1] = (e[-1] - e[-2]) / 2
    de[1:-1] = (e[2:] - e[:-2]) / 2
    # as the loader (Chi_spectral_weights): negative A tolerated down to round-off or CLIP_PEAK of the
    # peak, adding up to at most CLIP_REL of the weight; tolerated negatives are set to zero (A in 1/Ha)
    tol = max(100 * np.finfo(float).eps * max(1.0, np.abs(a).max() * HA2EV) / HA2EV, CLIP_PEAK * max(a.max(), 0.0))
    clipped = float(-np.sum(np.minimum(a, 0.0) * de))
    neg = bool(a.min() < -tol or clipped > CLIP_REL * max(float(np.sum(np.maximum(a, 0.0) * de)), 0.0))
    return float(np.sum(np.maximum(a, 0.0) * de)), float(a.min() / a.max()), float(e[np.argmax(a)]), neg, clipped


def load_newton(path):
    """(k, band) -> Sc|Eo [eV] from an o-*.qp written with ExtendOut."""
    titles, rows = None, []
    for line in open(path):
        if line.startswith('#'):
            if 'K-point' in line and 'Band' in line:
                # one title = a name and its optional [unit]; names can be followed by several blanks
                titles = re.findall(r'[A-Za-z][^\s\[]*(?:\s*\[[^\]]*\])?', line.lstrip('#'))
            continue
        if line.strip():
            rows.append([float(x) for x in line.split()])
    if titles is None or not any(t.startswith('Sc|Eo') for t in titles):
        sys.exit(f'!!! {path}: no Sc|Eo column (run the Newton GW with ExtendOut)')
    ic = [j for j, t in enumerate(titles) if t.startswith('Sc|Eo')][0]
    return {(int(r[0]), int(r[1])): r[ic] for r in rows}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('ndb')
    p.add_argument('--newton')
    p.add_argument('--ref')
    p.add_argument('--ref-newton')
    p.add_argument('--norm-tol', type=float, default=0.05)
    p.add_argument('--sc-tol', type=float, default=0.002, help='eV')
    a = p.parse_args()
    g = load_g(a.ndb)
    nw = load_newton(a.newton) if a.newton else None
    r = load_g(a.ref) if a.ref else None
    rn = load_newton(a.ref_newton) if a.ref_newton else None
    if r is not None and (r['eo'].size != g['eo'].size or np.any(r['band'] != g['band']) or np.any(r['k'] != g['k'])):
        sys.exit('!!! --ref has other states')
    bad = []
    head = '#    k  band     Eo[eV]   peak[eV]      norm    minA/maxA   ReSc(Eo)[eV]'
    head += '   -Newton[meV]' if nw else ''
    head += '   dSc_grid[meV] dSc_Newton[meV] diff[meV]' if r is not None and nw and rn else ''
    print(head)
    worst = {'norm': 0.0, 'neg': np.inf, 'sc': 0.0, 'shift': 0.0, 'clip': 0.0}
    n_clipped = 0
    for i in range(g['eo'].size):
        k, b = g['k'][i], g['band'][i]
        nrm, mina, peak, neg, clipped = norm_and_minA(g, i)
        sc, off = sc_at_eo(g, i)
        line = f'{k:6d}{b:6d}{g["eo"][i]:11.4f}{peak:11.4f}{nrm:10.5f}{mina:13.2e}{sc:15.5f}'
        worst['norm'] = max(worst['norm'], abs(nrm - 1))
        worst['neg'] = min(worst['neg'], mina)
        if abs(nrm - 1) > a.norm_tol:
            bad.append(f'k {k} band {b}: norm {nrm:.4f}')
        if neg:
            bad.append(f'k {k} band {b}: negative spectral weight ({mina:.2e} of the peak, {clipped:.1e} in all)')
        elif clipped > 0:
            n_clipped += 1
            worst['clip'] = max(worst['clip'], clipped)
        if nw:
            if (k, b) not in nw:
                sys.exit(f'!!! k {k} band {b} not in {a.newton}')
            d = sc - nw[(k, b)]
            line += f'{1000 * d:15.3f}'
            worst['sc'] = max(worst['sc'], abs(d))
            if abs(d) > a.sc_tol:
                bad.append(f'k {k} band {b}: Re Sc(Eo) differs from Newton by {1000 * d:.2f} meV')
            if r is not None and rn:
                dg = sc - sc_at_eo(r, i)[0]
                dn = nw[(k, b)] - rn[(k, b)]
                line += f'{1000 * dg:16.3f}{1000 * dn:16.3f}{1000 * (dg - dn):10.3f}'
                worst['shift'] = max(worst['shift'], abs(dg - dn))
                if abs(dg - dn) > a.sc_tol:
                    bad.append(f'k {k} band {b}: grid and Newton terminator shifts differ by {1000 * (dg - dn):.2f} meV')
        if off > 1e-3:
            line += f'   (Eo {1000 * off:.1f} meV from a grid point: interpolated)'
        print(line)
    print(f'# worst: |norm-1| {worst["norm"]:.4f} (tol {a.norm_tol}), min A/max A {worst["neg"]:.2e}'
          + (f', |Sc - Newton| {1000 * worst["sc"]:.3f} meV' if nw else '')
          + (f', |shift difference| {1000 * worst["shift"]:.3f} meV' if r is not None and nw and rn else '')
          + f' (tol {1000 * a.sc_tol:.1f} meV)')
    if n_clipped:
        print(f'# negative weight set to zero (as the loader) in {n_clipped} states, largest {worst["clip"]:.1e} of the norm')
    if bad:
        print('!!! FAILED:\n  ' + '\n  '.join(bad[:20]))
        sys.exit(1)
    print('# all checks passed')


if __name__ == '__main__':
    main()
