#!/usr/bin/env python3
"""
ktt_casida.py -- diagonalize the Casida matrix written by BSEChiKoutW.

Input: the files of a BSKmod="CHI" + BSEChiDyn run with BSEChiKoutW (and BSEChiKoutB):

  o-<job>.Ktt_q1_transitions   t, T, ik_bz, ik_ibz, k (rlu), v, c, E_t [eV], f_t, Re/Im a_t
  o-<job>.Ktt_q1_w<iw>         t t' and the X (exchange) and F (fxc) parts of the rr, rc, cr, cc
                               blocks in eV (only rr in TDA)

The matrix at the frequency z of each w-file is

  M(z) = [[ diag(E) + X_rr + F_rr ,           X_rc + F_rc ],
          [           X_cr + F_cr , -diag(E) + X_cc + F_cc ]]

and the response Resp(w) = -Co B^T [w - M]^-1 A with A = (a, i b), B = (b, i a), b = conj(a).
The script reports its eigenvalues with Re > 0 (lowest first), their oscillator strength
|(B^T R_l)(L_l^T A)| relative to the strongest, the lowest bright one (> 1% by default) and the
binding energy against the lowest transition of the file, and optionally writes eps2 up to the
constant Co (a Lorentzian sum over the eigenvalues).

Notes
 * The kernel is exact only on the whole transition space it was built for. A file restricted
   to one band pair (BSEChiKoutB) is a sub-matrix: the binding it gives is not the BSE one.
 * A fixed-z eigenvalue equals the exciton only at z = E_exciton with zero damping. The
   ndb.Chi grid has damping eta and a finite step, so expect deviations of order eta.
 * --tda keeps the rr block only; --no-fxc drops F (exchange only, the RPA+LF matrix).

Usage
  python ktt_casida.py o-cas_dyn_exc.Ktt_q1_transitions o-cas_dyn_exc.Ktt_q1_w225 [more w-files]
         [--tda] [--no-fxc] [--nlow 8] [--bright 0.01] [--spectrum eps.dat --eta 0.05]
"""
import argparse
import sys

import numpy as np


def read_table(path):
    """Numeric rows of a Yambo text file ('#' comments skipped), fast for large files."""
    with open(path) as fh:
        header = []
        for line in fh:
            if line.startswith('#'):
                header.append(line)
                continue
            first = line
            break
        else:
            return header, np.zeros((0, 0))
        ncol = len(first.split())
        body = first + fh.read()
    data = np.array(body.split(), dtype=float)
    if data.size % ncol:
        sys.exit(f'{path}: ragged table ({data.size} numbers, {ncol} columns)')
    return header, data.reshape(-1, ncol)


def read_frequency(header, path):
    for line in header:
        if 'z =' in line:
            re_z, im_z = line.split('z =')[1].split()[:2]
            return complex(float(re_z), float(im_z))
    sys.exit(f'{path}: no "z =" line in the header')


def build_matrix(rows, n, E, tda, no_fxc):
    """M(z) from the w-file rows; TDA files (6 columns) or --tda give the rr block only."""
    i = rows[:, 0].astype(int) - 1
    j = rows[:, 1].astype(int) - 1
    nblk = (rows.shape[1] - 2) // 4
    if nblk not in (1, 4):
        sys.exit(f'unexpected number of columns: {rows.shape[1]}')
    blocks = []
    for b in range(nblk):
        blk = np.zeros((n, n), complex)
        blk[i, j] = rows[:, 2 + 4 * b] + 1j * rows[:, 3 + 4 * b]              # X
        if not no_fxc:
            blk[i, j] += rows[:, 4 + 4 * b] + 1j * rows[:, 5 + 4 * b]         # F
        blocks.append(blk)
    if nblk == 1 or tda:
        return np.diag(E).astype(complex) + blocks[0], False
    M = np.block([[np.diag(E) + blocks[0], blocks[1]],
                  [blocks[2], -np.diag(E) + blocks[3]]])
    return M, True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('transitions')
    ap.add_argument('wfiles', nargs='+')
    ap.add_argument('--tda', action='store_true', help='resonant block only')
    ap.add_argument('--no-fxc', action='store_true', help='drop the fxc part F (exchange only)')
    ap.add_argument('--nlow', type=int, default=8, help='eigenvalues to print')
    ap.add_argument('--bright', type=float, default=0.01, help='relative strength of a bright state')
    ap.add_argument('--spectrum', help='write eps2 up to a constant (Lorentzians) to this file')
    ap.add_argument('--eta', type=float, default=0.05, help='Lorentzian width for --spectrum [eV]')
    ap.add_argument('--wmax', type=float, default=4.0, help='upper energy of --spectrum [eV]')
    a = ap.parse_args()

    _, tr = read_table(a.transitions)
    n = len(tr)
    E = tr[:, 9]
    res = tr[:, 11] + 1j * tr[:, 12]              # a_t
    pairs = sorted({(int(v), int(c)) for v, c in tr[:, 7:9]})
    gap = E.min()
    print(f'{n} transitions, band pairs {pairs[:6]}{" ..." if len(pairs) > 6 else ""}; '
          f'lowest transition {gap:.4f} eV')
    if len(pairs) == 1:
        print('  !!! one band pair only: a sub-matrix of the Casida problem, not the BSE exciton')

    for wf in a.wfiles:
        header, rows = read_table(wf)
        z = read_frequency(header, wf)
        M, coupled = build_matrix(rows, n, E, a.tda, a.no_fxc)
        if coupled:
            A = np.concatenate([res, 1j * np.conj(res)])
            B = np.concatenate([np.conj(res), 1j * res])
        else:
            A, B = res, np.conj(res)
        lam, R = np.linalg.eig(M)
        L = np.linalg.inv(R)                       # rows: left eigenvectors, L R = 1
        strength = (B @ R) * (L @ A)               # residue of each pole in B^T (w-M)^-1 A
        pos = np.where(lam.real > 1e-6)[0]
        pos = pos[np.argsort(lam[pos].real)]
        s = np.abs(strength[pos])
        s_rel = s / max(s.max(), 1e-300)
        kind = ('coupled' if coupled else 'TDA') + (', exchange only' if a.no_fxc else '')
        print(f'\n{wf}: z = {z.real:.4f} {z.imag:+.4f}i eV, {kind}, matrix {M.shape[0]}')
        print('    lowest Re E [eV]   Im E [eV]   rel. strength')
        for k in range(min(a.nlow, len(pos))):
            print(f'    {lam[pos[k]].real:12.5f} {lam[pos[k]].imag:11.5f} {s_rel[k]:14.4e}')
        bright = np.where(s_rel > a.bright)[0]
        if len(bright):
            eb = lam[pos[bright[0]]].real
            print(f'    lowest bright: {eb:.4f} eV   binding (lowest transition - bright) = {gap - eb:.4f} eV')
        if a.spectrum:
            w = np.linspace(0.0, a.wmax, 801)
            eps2 = np.zeros_like(w)
            for k in pos:
                eps2 += (-(strength[k] / (w + 1j * a.eta - lam[k]))).imag
            out = a.spectrum if len(a.wfiles) == 1 else f'{a.spectrum}.{wf.rsplit("_w", 1)[-1]}'
            np.savetxt(out, np.column_stack([w, eps2]), fmt='%14.6e',
                       header=f'eps2 up to a constant, Lorentzians eta={a.eta} eV, M at z={z}')
            print(f'    spectrum written to {out}')


if __name__ == '__main__':
    main()
