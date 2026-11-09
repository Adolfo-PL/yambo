#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
qp_gaps.py -- KS and QP band gaps from a Yambo o-*.qp file (or KS gaps alone from a pw.x output).

The o-*.qp columns are k, band, Eo (KS), E-Eo (QP correction), ...; the QP energy is Eo + (E-Eo).
Yambo's report prints the gaps of the KS and Hartree-Fock energies, not of the QP ones: this does.

  python qp_gaps.py o-gw.qp [--vb 13]
  python qp_gaps.py nscf.out [--vb 13]     the DFT (KS) gaps only, before any GW run: the band energies
                                           pw.x prints at each k (needs verbosity='high' when the run
                                           has more than 100 k-points, or pw.x does not print them)

--vb is the highest occupied band (13 for WS2/WSe2 without SOC); the next band is the conduction band.
Prints the VBM and CBM, the direct gap at K and the indirect gap K -> Q (VB at K to CB at Q), and the
gap at each k (Yambo's k indices: match them to the "QP [eV] @ K [n] (iku)" lines of the r-file).
K is the k of the smallest KS direct gap; Q the k of the lowest KS conduction band away from the K
valley (its KS valence band at least 0.3 eV below the one at K, which leaves out K and its neighbours).
--k and --q set them by hand. Q lies between grid points on coarse grids: the nearest one is used.
"""
import argparse
import re
import sys

import numpy as np


NUM = re.compile(r'-?\d+\.\d+')


def read_pw(path):
    """k-points (cartesian, 2pi/alat) and band energies (eV) of a pw.x output: the last listing."""
    txt = open(path, errors='replace').read()
    m = re.search(r'number of Kohn-Sham states\s*=\s*(\d+)', txt)
    if not m:
        sys.exit(f'{path}: not a pw.x output (no "number of Kohn-Sham states")')
    nb = int(m.group(1))
    last = txt.rfind('End of band structure calculation')
    last = last if last >= 0 else txt.rfind('End of self-consistent calculation')
    if last < 0:
        sys.exit(f'{path}: no band listing (unfinished run?)')
    kpts, bands = [], []
    for blk in re.split(r'\n\s*k =', txt[last:])[1:]:
        head, _, rest = blk.partition('\n')
        if 'bands (ev)' not in head:
            continue
        kpts.append([float(x) for x in NUM.findall(head.split('(')[0])][:3])
        e = []
        for line in rest.split('\n'):
            if 'occupation numbers' in line or 'highest' in line or 'Fermi' in line or 'k =' in line:
                break
            e += [float(x) for x in NUM.findall(line)]
            if len(e) >= nb:
                break
        bands.append(e[:nb])
    if not bands:
        sys.exit(f'{path}: no band energies printed (more than 100 k-points: rerun with verbosity=\'high\')')
    return np.array(kpts), np.array(bands)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('qpfile')
    ap.add_argument('--vb', type=int, default=13, help='highest occupied band')
    ap.add_argument('--k', type=int, help='k index of K (default: the smallest KS direct gap)')
    ap.add_argument('--q', type=int, help='k index of Q (default: the lowest KS CB away from the K valley)')
    a = ap.parse_args()
    kc_pw = None
    with open(a.qpfile, errors='replace') as f:
        is_pw = 'Program PWSCF' in f.read(4000)
    if is_pw:                       # KS only: E-Eo = 0, k numbered as pw.x lists them
        kc_pw, B = read_pw(a.qpfile)
        nk, nb = B.shape
        k = np.repeat(np.arange(1, nk + 1), nb)
        b = np.tile(np.arange(1, nb + 1), nk)
        e0, de = B.ravel(), np.zeros(B.size)
    else:
        d = np.loadtxt(a.qpfile, comments='#', ndmin=2)
        k, b, e0, de = d[:, 0].astype(int), d[:, 1].astype(int), d[:, 2], d[:, 3]
    ks = sorted(set(k[b == a.vb]) & set(k[b == a.vb + 1]))
    if not ks:
        sys.exit(f'bands {a.vb} and {a.vb + 1} are not both in {a.qpfile} (QPkrange)')
    v0 = {q: e0[(k == q) & (b == a.vb)][0] for q in ks}
    c0 = {q: e0[(k == q) & (b == a.vb + 1)][0] for q in ks}
    K = a.k if a.k else min(ks, key=lambda q: c0[q] - v0[q])
    away = [q for q in ks if v0[q] < v0[K] - 0.3]
    Q = a.q if a.q else (min(away, key=c0.get) if away else None)
    if K not in ks or (Q is not None and Q not in ks):
        sys.exit(f'--k/--q must be k indices in {a.qpfile}: {ks[0]}..{ks[-1]}')
    print(f'K = k {K}, Q = k {Q}' + ('' if a.k or a.q else '  (automatic; set --k/--q to override)'))
    if kc_pw is not None:
        for lab, q in (('K', K), ('Q', Q)):
            if q is not None:
                x = kc_pw[q - 1]
                print(f'  {lab}: k = {x[0]:.4f} {x[1]:.4f} {x[2]:.4f} (2pi/alat), |k| = {np.linalg.norm(x):.4f}'
                      + ('  (K points have |k| = 2/3 = 0.6667)' if lab == 'K' else ''))
    for lab, E in ((('KS',  e0),) if is_pw else (('KS', e0), ('QP', e0 + de))):
        v = {q: E[(k == q) & (b == a.vb)][0] for q in ks}
        c = {q: E[(k == q) & (b == a.vb + 1)][0] for q in ks}
        kv, kc = max(v, key=v.get), min(c, key=c.get)
        line = (f'{lab}: VBM {v[kv]:9.4f} eV at k {kv}, CBM {c[kc]:9.4f} eV at k {kc}; '
                f'direct gap (K) {c[K] - v[K]:.4f} eV')
        if Q is not None:
            line += f', indirect gap (K -> Q) {c[Q] - v[K]:.4f} eV'
        print(line)
    print('\n   k    VB(KS)    CB(KS)  gap(KS)' + ('' if is_pw else '     VB(QP)    CB(QP)  gap(QP)'))
    for q in ks:
        m = (k == q)
        vk, ck = e0[m & (b == a.vb)][0], e0[m & (b == a.vb + 1)][0]
        vq, cq = vk + de[m & (b == a.vb)][0], ck + de[m & (b == a.vb + 1)][0]
        print(f'{q:4d} {vk:9.4f} {ck:9.4f} {ck - vk:8.4f}'
              + ('' if is_pw else f'  {vq:9.4f} {cq:9.4f} {cq - vq:8.4f}'))


if __name__ == '__main__':
    main()
