#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
gf_spectrum.py -- spectral functions of a retarded ndb.G (GRetarded) as plain text, one file per state.

  python gf_spectrum.py <ndb.G> [--out gf]

Writes <out>_k<k>_b<band>.dat with the columns
  0  w (eV)              real part of the state's grid (Yambo's energy zero)
  1  A(w) (1/eV)         -Im G/pi; negative values are left as they are (the loader clips small ones)
  2  Re Sigma_c (eV)     Sigma_c = SE_Operator - G_SIGMA_X + G_VXC_REFERENCE, at w + i*GDmRnge
  3  Im Sigma_c (eV)
  4  w - Eo (eV)         energy from the state's KS energy
and prints, per state, Eo, the peak of A, the integral of A on the grid and the weight within
+-2 eV of the peak. Plot with plot_gf_spectrum.py.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_check import load_g  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('ndb')
    p.add_argument('--out', default='gf')
    a = p.parse_args()
    g = load_g(a.ndb)
    print('#    k  band     Eo[eV]   peak[eV]   integral   weight(peak+-2eV)   file')
    for i in range(g['eo'].size):
        w = g['w'][i].real
        A = -g['G'][i].imag / np.pi
        sc = g['sc'][i]
        de = np.gradient(w)
        ip = int(np.argmax(A))
        near = np.abs(w - w[ip]) <= 2.0
        name = f"{a.out}_k{g['k'][i]}_b{g['band'][i]}.dat"
        with open(name, 'w') as o:
            o.write(f"# {os.path.abspath(a.ndb)}: k {g['k'][i]} band {g['band'][i]}, Eo = {g['eo'][i]:.5f} eV,"
                    f" peak of A at {w[ip]:.5f} eV, broadening {g['w'][i].imag[0] * 1000:.1f} meV\n")
            o.write('# w(eV)  A(1/eV)  ReSigma_c(eV)  ImSigma_c(eV)  w-Eo(eV)\n')
            for j in range(w.size):
                o.write(f'{w[j]:12.6f} {A[j]:14.6e} {sc[j].real:14.6e} {sc[j].imag:14.6e} {w[j] - g["eo"][i]:12.6f}\n')
        print(f"{g['k'][i]:6d}{g['band'][i]:6d}{g['eo'][i]:11.4f}{w[ip]:11.4f}{np.sum(A * de):11.4f}"
              f"{np.sum((A * de)[near]):20.4f}   {name}")


if __name__ == '__main__':
    main()
