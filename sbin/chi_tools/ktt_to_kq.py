#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
ktt_to_kq.py -- the Casida matrix of the exact kernel as input files of the k-space Casida program
(bindingenergies: H = diag(E_c - E_v) + Hartree * F_kq * sqrt(w_k w_q), diagonalized with ZHEEV).

From a transitions file and one w-file of a BSKmod="CHI" + BSEChiDyn run with BSEChiKoutW (the
kernel frozen at the frequency z of that file, e.g. the A exciton: prod/o-xspA.Ktt_q1_w<i> of
exact_kernel.sh STAGE=spectra), it writes the resonant block (Tamm-Dancoff, the program solves no
coupling) of the whole transition space, every transition t = (v, c, k) one "k" of the program:

  energy_all    E(t, spin, band) in the order the program reads it (t, then band, then spin):
                E_t on (spin 1, band NC), 0 on (spin 1, band NV) and everywhere else
  kpoints.dat   kx ky kz weight: k of the transition (crystal units), weight 1/N
  F_kq_LR.dat   t t' (re,im): (X + F)_tt' / (Hartree * sqrt(w_t w_t')), row by row, so that the
                program's Hartree * F_kq * sqrt(w w) is the Casida matrix element in eV again
  params.txt    the parameter values to set in the program (n_k_max = N_m = N, ...)

and diagonalizes the same matrix itself (numpy, Hermitian): its eigenvalues minus the lowest
transition are what the program writes to Eq0, for a direct check of the hand-off.

  python ktt_to_kq.py o-xspA.Ktt_q1_transitions o-xspA.Ktt_q1_w37 --out kq_EA [--no-x]

--no-x leaves the exchange (Hartree) part X out of F_kq: the kernel part alone.
--fold writes the coupled problem instead: H_eff = sum_n lambda_n |X~_n><X~_n| over its N positive
eigenvalues, with X~_n the orthonormalized excitation parts of the eigenvectors. It is Hermitian and
N x N, so the program diagonalizes it unchanged and returns the coupled energies (its absorption
weights stay 1/N per state).
The exact kernel reproduces the BSE only in the coupled problem (ktt_casida.py without --tda);
ktt_casida.py --tda on the same files gives the TDA numbers this hand-off should reproduce.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ktt_casida as kc   # noqa: E402

HARTREE = 2.0 * 13.605698   # eV, the program's own value (Ry = 13.605698 eV)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('transitions')
    ap.add_argument('wfile', help='w-file of the frequency to freeze the kernel at')
    ap.add_argument('--out', default='kq', help='folder for the program input files')
    ap.add_argument('--no-x', action='store_true', help='kernel part only, without the exchange X')
    ap.add_argument('--fold', action='store_true', help='write the coupled problem folded into an N x N Hermitian '
                    'matrix with the same N positive eigenvalues (the program then gives the coupled energies)')
    ap.add_argument('--n-bands', type=int, default=17, help='bands in energy_all (program n_bands)')
    ap.add_argument('--nv', type=int, default=13, help='valence band index of the program (n_v)')
    ap.add_argument('--nc', type=int, default=14, help='conduction band index of the program (n_c)')
    ap.add_argument('--n-spin', type=int, default=2, help='spin entries per band in energy_all (n_spin)')
    ap.add_argument('--nlow', type=int, default=8)
    a = ap.parse_args()
    if not (1 <= a.nv <= a.n_bands and 1 <= a.nc <= a.n_bands and a.nv != a.nc):
        sys.exit('--nv and --nc must be different bands within --n-bands')

    tr, E, Ea, res, Aa, Ba = kc.read_rows(a.transitions)
    n = len(tr)
    z, M, coupled = kc.read_matrix(a.wfile, n, E, not a.fold, False, E_ares=Ea)   # TDA: resonant block, X + F
    if a.no_x:
        _, MX, _ = kc.read_matrix(a.wfile, n, E, not a.fold, True, E_ares=Ea)    # X only
        M = M - MX + np.diag(np.concatenate([E, Ea]) if M.shape[0] == 2 * n else E)
    if a.fold:
        if M.shape[0] != 2 * n:
            sys.exit('--fold needs a w-file with the four blocks (a coupled run)')
        # the N eigenvalues with Re > 0 of the coupled problem, and the excitation parts X_n of their
        # eigenvectors, orthonormalized (Loewdin): H_eff = Xt diag(lambda) Xt^dagger has exactly
        # those eigenvalues and Hermitian form, as the program's ZHEEV needs
        lam, R = np.linalg.eig(M)
        pos = np.where(lam.real > 1e-6)[0]
        pos = pos[np.argsort(lam[pos].real)]
        if len(pos) != n:
            sys.exit(f'--fold: {len(pos)} eigenvalues with Re > 0, expected {n}')
        print(f'  coupled problem: largest |Im lambda| of the excitations {np.abs(lam[pos].imag).max():.1e} eV')
        Xn = R[:n, pos]
        s, V = np.linalg.eigh(Xn.conj().T @ Xn)
        if s.min() <= 1e-12 * s.max():
            sys.exit('--fold: the excitation parts are linearly dependent; cannot fold')
        Xt = Xn @ (V * s ** -0.5) @ V.conj().T
        H = (Xt * lam[pos].real) @ Xt.conj().T
        H = 0.5 * (H + H.conj().T)
        K = H - np.diag(E)
    else:
        K = M - np.diag(E)                                                   # X + F
    herm = np.abs(K - K.conj().T).max() / max(np.abs(K).max(), 1e-300)
    pairs = sorted({(int(v), int(c)) for v, c in tr[:, 7:9]})
    print(f'{n} transitions (band pairs {pairs}), kernel frozen at z = {z.real:.4f} {z.imag:+.4f}i eV, '
          f'{"F only" if a.no_x else "X + F"}, '
          + ('coupled problem folded into N x N' if a.fold else 'resonant block (TDA)'))
    print(f'  non-Hermitian part of the block: {herm:.1e} of its largest element'
          + ('' if herm < 1e-6 else '  !!! ZHEEV reads the upper triangle only: the hand-off is not exact'))
    if abs(z.imag) > 1e-6:
        print('  !!! damped frequency: use a w-file of an undamped run (the kernel at a real frequency)')

    os.makedirs(a.out, exist_ok=True)
    w = np.full(n, 1.0 / n)
    with open(os.path.join(a.out, 'energy_all'), 'w') as f:      # t, then band, then spin
        for t in range(n):
            for b in range(1, a.n_bands + 1):
                for s in range(1, a.n_spin + 1):
                    f.write(f'{E[t] if (b == a.nc and s == 1) else 0.0:.10f}\n')
    with open(os.path.join(a.out, 'kpoints.dat'), 'w') as f:
        for t in range(n):
            f.write(f'{tr[t, 4]:.8f} {tr[t, 5]:.8f} {tr[t, 6]:.8f} {w[t]:.12e}\n')
    F = K / (HARTREE * np.sqrt(np.outer(w, w)))
    with open(os.path.join(a.out, 'F_kq_LR.dat'), 'w') as f:     # row by row: l1, then m1
        for i in range(n):
            f.write(''.join(f'{i + 1} {j + 1} ({F[i, j].real:.12e},{F[i, j].imag:.12e})\n' for j in range(n)))
    with open(os.path.join(a.out, 'params.txt'), 'w') as f:
        f.write(f'n_k_max = {n}   (N_m = n_k_max)\nn_bands = {a.n_bands}\nn_spin = {a.n_spin}\n'
                f'n_v = {a.nv}\nn_c = {a.nc}\nn_spin_v = 1\nn_spin_c = 1\nepsilon = 1.0d0\n'
                f'kernel frozen at z = {z.real:.6f} eV; {"F only" if a.no_x else "X + F"}; '
                f'{"coupled, folded" if a.fold else "TDA"}\n')

    lam = np.linalg.eigvalsh(np.diag(E) + K)          # what the program diagonalizes
    gap = E.min()
    print(f'  wrote {a.out}/energy_all, kpoints.dat, F_kq_LR.dat, params.txt')
    print(f'  set in the program: n_k_max = {n} (N_m), n_bands = {a.n_bands}, n_v = {a.nv}, n_c = {a.nc}, '
          f'n_spin = {a.n_spin}, epsilon = 1')
    print(f'  expected program output: gap = {gap:.6f}; Eq0 (eigenvalue - gap), lowest {a.nlow}:')
    for i, l in enumerate(lam[:a.nlow], 1):
        print(f'    {i:4d} {l - gap:12.6f}   (eigenvalue {l:.6f} eV)')


if __name__ == '__main__':
    main()
