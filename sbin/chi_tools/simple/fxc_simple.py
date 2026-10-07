#!/usr/bin/env python3
# The kernel f_xc(G,G',q,w) of an ndb.Chi as plain text: one file per q, at one frequency.
#
#   f_xc = chi0^-1 - P^-1    (computed by Yambo in the BSE run, stored in ndb.Chi)
#
# usage:   python3 fxc_simple.py kx_fq [iw]
#          iw = 1 (default) is w = 0, the static kernel
#
# output:  G_list.dat     iG  G_1 G_2 G_3  (Yambo units: component i in 2 pi / alat_i)
#          fxc_q<iq>.dat  iG  iG'  Re f_xc  Im f_xc   (Hartree atomic units), one file per q
#
# Only the G where the kernel is not zero are written (the 7 in-plane G with BSEGinplane).
import os
import sys

import numpy as np
from netCDF4 import Dataset

folder = sys.argv[1]
iw = int(sys.argv[2]) if len(sys.argv) > 2 else 1

main = Dataset(os.path.join(folder, 'ndb.Chi'))
G = np.array(main.variables['CHI_RL_vecs'][:]).T            # (nG, 3)
nq = np.array(main.variables['HEAD_QPT'][:]).shape[1]

keep = None
for iq in range(1, nq + 1):
    # the kernel of q is in ndb.Chi_fragment_<iq> (or in ndb.Chi itself)
    f = os.path.join(folder, f'ndb.Chi_fragment_{iq}')
    ds = Dataset(f) if os.path.exists(f) else main
    if f'FXC_Q_{iq}' not in ds.variables:
        continue
    x = np.array(ds.variables[f'FXC_Q_{iq}'][:])               # (nw, nG, nG, 2), G' first
    fxc = (x[..., 0] + 1j * x[..., 1]).transpose(0, 2, 1)[iw - 1]   # fxc[G, G']
    w = np.array(ds.variables[f'CHI_FREQ_Q_{iq}'][:])[iw - 1] * 27.211386245988
    if keep is None:
        keep = [i for i in range(len(G)) if np.abs(fxc[i]).max() > 0]
        with open('G_list.dat', 'w') as o:
            for n, i in enumerate(keep, 1):
                o.write(f'{n:4d} {G[i, 0]:12.6f} {G[i, 1]:12.6f} {G[i, 2]:12.6f}\n')
    with open(f'fxc_q{iq}.dat', 'w') as o:
        o.write(f'# f_xc(G,G\') at q index {iq}, w = {w[0]:.6f} {w[1]:+.6f}i eV, Hartree a.u.\n')
        for n, i in enumerate(keep, 1):
            for m, j in enumerate(keep, 1):
                o.write(f'{n:4d} {m:4d} {fxc[i, j].real:18.9e} {fxc[i, j].imag:18.9e}\n')
    print(f'fxc_q{iq}.dat: {len(keep)} G, w = {w[0]:.4f} eV')
