#!/usr/bin/env python3
"""
fxc_export.py -- write the kernel fxc_GG'(q, w) of an ndb.Chi as plain text (and .npz).

  python fxc_export.py <dir with ndb.Chi> [--out fxc_GGq] [--static] [--sym]

For every q fragment present (a BSEChiOut export holds only the optical q -> 0) the text
file lists q, the G vectors with |q+G|, and the full G x G' matrix at each frequency
(--static: w = 0 only). --sym also writes fxc |q+G||q+G'| / 4 pi, the kernel in units of
the bare Coulomb interaction; its head is -alpha/(4 pi), alpha = -fxc_00 q^2.

Units: fxc as stored (Hartree atomic units); |q+G| in bohr^-1; frequencies in eV.
G as stored in CHI_RL_vecs (Yambo's Cartesian components, units of 2 pi/alat per axis).
The .npz holds the same arrays: G, qpg[iq], q[iq], freqs_eV[iq], fxc[iq] (nw, nG, nG).
"""
import argparse
import os

import numpy as np
from netCDF4 import Dataset

HA2EV = 27.211386245988


def text(v):
    a = np.asarray(v[:]).ravel()
    return b''.join(x if isinstance(x, bytes) else bytes(str(x), 'ascii') for x in a).decode(errors='ignore').strip()


def cplx(v):
    x = np.asarray(v[:])
    return x[..., 0] + 1j * x[..., 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('dir')
    ap.add_argument('--out', default='fxc_GGq')
    ap.add_argument('--static', action='store_true', help='w = 0 only')
    ap.add_argument('--sym', action='store_true', help='also fxc |q+G||q+G\'|/4pi')
    a = ap.parse_args()

    with Dataset(os.path.join(a.dir, 'ndb.Chi')) as ds:
        G = np.asarray(ds.variables['CHI_RL_vecs'][:]).T
        qpg_all = cplx(ds.variables['CHI_QPG'])              # (ng, nq)
        qpts = np.asarray(ds.variables['HEAD_QPT'][:]).T     # (nq, 3)
        mode = text(ds.variables['CHI_G_MODE'])
        conv = text(ds.variables['CHI_CONVENTION'])
        action = text(ds.variables['CHI_QP_ACTION'])
        q0 = np.asarray(ds.variables['CHI_OPTICAL_Q0'][:])
    ng = G.shape[0]

    data = {}
    for iq in range(1, qpts.shape[0] + 1):
        for f in (f'ndb.Chi_fragment_{iq}', 'ndb.Chi'):
            p = os.path.join(a.dir, f)
            if not os.path.exists(p):
                continue
            with Dataset(p) as ds:
                if f'FXC_Q_{iq}' not in ds.variables:
                    continue
                w = cplx(ds.variables[f'CHI_FREQ_Q_{iq}']) * HA2EV
                fx = np.transpose(cplx(ds.variables[f'FXC_Q_{iq}']), (0, 2, 1))
                data[iq] = (w, fx)
            break
    if not data:
        raise SystemExit('no FXC_Q_* in ' + a.dir)

    with open(a.out + '.dat', 'w') as o:
        o.write(f'# fxc_GG\'(q, w) from {os.path.abspath(a.dir)}/ndb.Chi\n')
        o.write(f'# mode: {mode}\n# {action}\n# convention: {conv}\n')
        o.write(f'# units: fxc Hartree atomic units; |q+G| bohr^-1; w eV. G: CHI_RL_vecs as stored\n')
        o.write(f'# q fragments present: {sorted(data)} of {qpts.shape[0]} q points in the header\n')
        o.write(f'# nG = {ng}\n')
        for iq, (w, fx) in data.items():
            qpg = np.abs(qpg_all[:, iq - 1])
            qv = q0 if (iq == 1 and np.any(q0)) else qpts[iq - 1]
            o.write(f'\n# ===== iq = {iq}   q = {qv[0]:.8e} {qv[1]:.8e} {qv[2]:.8e} (q0 direction for iq=1)\n')
            o.write('#   iG      G_1          G_2          G_3        |q+G|\n')
            for i in range(ng):
                o.write(f'  {i + 1:4d} {G[i, 0]:12.6f} {G[i, 1]:12.6f} {G[i, 2]:12.6f} {qpg[i]:14.6e}\n')
            iws = [0] if a.static else range(len(w))
            for iw in iws:
                o.write(f'# --- iw = {iw + 1}   w = {w[iw].real:.6f} {w[iw].imag:+.6f}i eV\n')
                o.write('#   iG  iG\'      Re fxc            Im fxc' +
                        ('          Re fxc*|q+G||q+G\'|/4pi  Im' if a.sym else '') + '\n')
                for i in range(ng):
                    for j in range(ng):
                        v = fx[iw, i, j]
                        line = f'  {i + 1:4d} {j + 1:4d} {v.real:17.9e} {v.imag:17.9e}'
                        if a.sym:
                            s = v * qpg[i] * qpg[j] / (4 * np.pi)
                            line += f' {s.real:17.9e} {s.imag:17.9e}'
                        o.write(line + '\n')
    np.savez(a.out + '.npz', G=G, q=np.array([q0 if iq == 1 else qpts[iq - 1] for iq in data]),
             iq=np.array(sorted(data)), qpg=np.array([np.abs(qpg_all[:, iq - 1]) for iq in data]),
             freqs_eV=np.array([data[iq][0] for iq in data]),
             fxc=np.array([data[iq][1][:1] if a.static else data[iq][1] for iq in data]))
    for iq, (w, fx) in data.items():
        qpg = np.abs(qpg_all[:, iq - 1])
        print(f'iq={iq}: nG={ng}, {len(w)} frequencies; alpha(w=0) = {-fx[0, 0, 0].real * qpg[0] ** 2:.5f}')
    print(f'wrote {a.out}.dat and {a.out}.npz')


if __name__ == '__main__':
    main()
