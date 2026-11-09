#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
qp_gaps.py -- KS and QP band gaps from a Yambo o-*.qp file.

The o-*.qp columns are k, band, Eo (KS), E-Eo (QP correction), ...; the QP energy is Eo + (E-Eo).
Yambo's report prints the gaps of the KS and Hartree-Fock energies, not of the QP ones: this does.

  python qp_gaps.py o-gw.qp [--vb 13]

--vb is the highest occupied band (13 for WS2/WSe2 without SOC); the next band is the conduction band.
Prints the VBM and CBM, the direct gap at K and the indirect gap K -> Q (VB at K to CB at Q), and the
gap at each k (Yambo's k indices: match them to the "QP [eV] @ K [n] (iku)" lines of the r-file).
K is the k of the smallest KS direct gap; Q the k of the lowest KS conduction band away from the K
valley (its KS valence band at least 0.3 eV below the one at K, which leaves out K and its neighbours).
--k and --q set them by hand. Q lies between grid points on coarse grids: the nearest one is used.
"""
import argparse
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('qpfile')
    ap.add_argument('--vb', type=int, default=13, help='highest occupied band')
    ap.add_argument('--k', type=int, help='k index of K (default: the smallest KS direct gap)')
    ap.add_argument('--q', type=int, help='k index of Q (default: the lowest KS CB away from the K valley)')
    a = ap.parse_args()
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
    for lab, E in (('KS', e0), ('QP', e0 + de)):
        v = {q: E[(k == q) & (b == a.vb)][0] for q in ks}
        c = {q: E[(k == q) & (b == a.vb + 1)][0] for q in ks}
        kv, kc = max(v, key=v.get), min(c, key=c.get)
        line = (f'{lab}: VBM {v[kv]:9.4f} eV at k {kv}, CBM {c[kc]:9.4f} eV at k {kc}; '
                f'direct gap (K) {c[K] - v[K]:.4f} eV')
        if Q is not None:
            line += f', indirect gap (K -> Q) {c[Q] - v[K]:.4f} eV'
        print(line)
    print('\n   k    VB(KS)    CB(KS)  gap(KS)     VB(QP)    CB(QP)  gap(QP)')
    for q in ks:
        m = (k == q)
        vk, ck = e0[m & (b == a.vb)][0], e0[m & (b == a.vb + 1)][0]
        vq, cq = vk + de[m & (b == a.vb)][0], ck + de[m & (b == a.vb + 1)][0]
        print(f'{q:4d} {vk:9.4f} {ck:9.4f} {ck - vk:8.4f}  {vq:9.4f} {cq:9.4f} {cq - vq:8.4f}')


if __name__ == '__main__':
    main()
