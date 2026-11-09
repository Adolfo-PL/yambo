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
Prints the VBM and CBM k index and energy, the indirect and the smallest direct gap, and the gap at
each k (Yambo's k indices: match them to the "QP [eV] @ K [n] (iku)" lines of the r-file).
"""
import argparse
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('qpfile')
    ap.add_argument('--vb', type=int, default=13, help='highest occupied band')
    a = ap.parse_args()
    d = np.loadtxt(a.qpfile, comments='#', ndmin=2)
    k, b, e0, de = d[:, 0].astype(int), d[:, 1].astype(int), d[:, 2], d[:, 3]
    ks = sorted(set(k[b == a.vb]) & set(k[b == a.vb + 1]))
    if not ks:
        sys.exit(f'bands {a.vb} and {a.vb + 1} are not both in {a.qpfile} (QPkrange)')
    for lab, E in (('KS', e0), ('QP', e0 + de)):
        v = {q: E[(k == q) & (b == a.vb)][0] for q in ks}
        c = {q: E[(k == q) & (b == a.vb + 1)][0] for q in ks}
        kv, kc = max(v, key=v.get), min(c, key=c.get)
        kd = min(ks, key=lambda q: c[q] - v[q])
        print(f'{lab}: VBM {v[kv]:9.4f} eV at k {kv}, CBM {c[kc]:9.4f} eV at k {kc}; '
              f'indirect gap {c[kc] - v[kv]:.4f} eV, direct gap {c[kd] - v[kd]:.4f} eV at k {kd}')
    print('\n   k    VB(KS)    CB(KS)  gap(KS)     VB(QP)    CB(QP)  gap(QP)')
    for q in ks:
        m = (k == q)
        vk, ck = e0[m & (b == a.vb)][0], e0[m & (b == a.vb + 1)][0]
        vq, cq = vk + de[m & (b == a.vb)][0], ck + de[m & (b == a.vb + 1)][0]
        print(f'{q:4d} {vk:9.4f} {ck:9.4f} {ck - vk:8.4f}  {vq:9.4f} {cq:9.4f} {cq - vq:8.4f}')


if __name__ == '__main__':
    main()
