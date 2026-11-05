#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
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
At finite q (explicit anti-resonant transitions) o-<job>.Ktt_q<iq>_antiresonant lists rows N+t
with their own energies E_ares and A, B; it is read when present, and the lower-right block is
diag(E_ares) + X_cc + F_cc.
The script reports its eigenvalues with Re > 0 (lowest first), their oscillator strength
|(B^T R_l)(L_l^T A)| relative to the strongest, the lowest bright one (> 1% by default) and the
binding energy against the lowest transition of the file, and optionally writes eps2 up to the
constant Co (a Lorentzian sum over the eigenvalues).

Notes
 * The kernel is exact only on the whole transition space it was built for. A file restricted
   to one band pair (BSEChiKoutB) is a sub-matrix: the binding it gives is not the BSE one.
 * A fixed-z eigenvalue equals the exciton only at z = E_exciton with zero damping. The
   ndb.Chi grid has damping eta and a finite step, so expect deviations of order eta.
 * Only bright eigenvalues are excitons. With a frequency-dependent kernel, M(w) also has
   eigenvalues that are dark in G space and are not BSE eigenvalues (doc/greens_function_fxc.md,
   BSEChiKoutW): do not take dark states off these eigenvectors without checking them.
 * --tda keeps the rr block only; --no-fxc drops F (exchange only, the RPA+LF matrix).
 * Memory: the files are read in chunks; the peak is the dense eigenproblem, about
   0.4 GB for 1024 transitions (coupled 2048) and ~2 GB for 2025 (4050). Run it on a
   compute node rather than a login node.

Usage
  python ktt_casida.py o-cas_dyn_exc.Ktt_q1_transitions o-cas_dyn_exc.Ktt_q1_w225 [more w-files]
         [--tda] [--no-fxc] [--nlow 8] [--bright 0.01] [--spectrum eps.dat --eta 0.05]
         [--eq0 Eq0.dat] [--emission em.dat --temp 0.01 --occupation boltzmann|bose] [--volo DIR]
"""
import argparse
import os
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


def read_rows(path):
    """Transitions file (and its _antiresonant companion at finite q): table, E, E_ares, a, A_ares, B_ares."""
    _, tr = read_table(path)
    E = tr[:, 9]
    res = tr[:, 11] + 1j * tr[:, 12]
    ares = path.replace('_transitions', '_antiresonant')
    if ares != path and os.path.exists(ares):
        _, ta = read_table(ares)
        if len(ta) != len(tr):
            sys.exit(f'{ares}: {len(ta)} rows, {path}: {len(tr)}')
        return tr, E, ta[:, 9], res, ta[:, 11] + 1j * ta[:, 12], ta[:, 13] + 1j * ta[:, 14]
    return tr, E, -E, res, 1j * np.conj(res), 1j * res


def full_weight(tr, R, nbf):
    """Weight of each eigenvector (columns of R) on transitions out of a full band (lower band <= nbf):
    with partly filled bands, the interband states are those with more than half of it."""
    m = tr[:, 7:9].min(axis=1) <= nbf
    if R.shape[0] == 2 * len(m):
        m = np.concatenate([m, m])
    w = np.abs(R) ** 2
    return w[m].sum(axis=0) / np.maximum(w.sum(axis=0), 1e-300)


def read_matrix(path, n, E, tda, no_fxc, chunk_lines=500000, E_ares=None):
    """Stream a w-file into M(z) chunk by chunk (the text is never held in memory whole).

    TDA files (6 columns) or --tda give the rr block only. Returns z, M, coupled."""
    z = None
    blocks = None
    nblk = None
    with open(path) as fh:
        buf = []
        for line in fh:
            if line.startswith('#'):
                if 'z =' in line:
                    re_z, im_z = line.split('z =')[1].split()[:2]
                    z = complex(float(re_z), float(im_z))
                continue
            if nblk is None:
                ncol = len(line.split())
                nblk = (ncol - 2) // 4
                if nblk not in (1, 4) or ncol != 2 + 4 * nblk:
                    sys.exit(f'{path}: unexpected number of columns: {ncol}')
                nkeep = 1 if (tda or nblk == 1) else 4
                blocks = [np.zeros((n, n), complex) for _ in range(nkeep)]
            buf.append(line)
            if len(buf) >= chunk_lines:
                _scatter(buf, ncol, blocks, no_fxc, path)
                buf = []
        if buf:
            _scatter(buf, ncol, blocks, no_fxc, path)
    if z is None:
        sys.exit(f'{path}: no "z =" line in the header')
    if blocks is None:
        sys.exit(f'{path}: no matrix rows')
    if len(blocks) == 1:
        return z, np.diag(E).astype(complex) + blocks[0], False
    M = np.empty((2 * n, 2 * n), complex)
    M[:n, :n] = np.diag(E) + blocks[0]
    M[:n, n:] = blocks[1]
    M[n:, :n] = blocks[2]
    M[n:, n:] = np.diag(-E if E_ares is None else E_ares) + blocks[3]
    return z, M, True


def volo_files(W, gap, out):
    """Eq0, absorption_0 and emission_0 of the k-space Casida program (eigenvalues_SPN.f90) from the
    eigenvalues W, with its formulas: every eigenvalue weighted d = 1/N (its sum_t |U_tn|^2 w_t with
    w_t = 1/N), A(w) = sum_n d / ((w - W_n)^2 + 0.01) on w = 0, 0.02, ..., 20 eV, and
    B(w) = 4 pi |w^4 A / (exp((w - W_1)/T) - 1)| for w > W_1, T = 0.01 eV, written as 0.8 B / 1600."""
    W = np.sort(np.asarray(W, float))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'Eq0'), 'w') as f:
        for i, x in enumerate(W, 1):
            f.write(f'{i:21d} {x - gap:25.17g}\n')
    w = 0.02 * np.arange(1001)
    A = (1.0 / len(W) / ((w[:, None] - W[None, :]) ** 2 + 0.01)).sum(1)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        B = np.where(w > W[0], 4.0 * np.pi * np.abs(w ** 4 * A / (np.exp((w - W[0]) / 0.01) - 1.0)), 0.0)
    B = np.nan_to_num(B, nan=0.0, posinf=0.0)
    np.savetxt(os.path.join(out, 'absorption_0'), np.column_stack([w, A]), fmt='%25.17g')
    np.savetxt(os.path.join(out, 'emission_0'), np.column_stack([w, 0.8 * B / 1600.0]), fmt='%25.17g')


def _scatter(lines, ncol, blocks, no_fxc, path):
    rows = np.fromstring(''.join(lines), sep=' ')
    if rows.size % ncol:
        sys.exit(f'{path}: ragged table')
    rows = rows.reshape(-1, ncol)
    i = rows[:, 0].astype(int) - 1
    j = rows[:, 1].astype(int) - 1
    for b, blk in enumerate(blocks):
        blk[i, j] = rows[:, 2 + 4 * b] + 1j * rows[:, 3 + 4 * b]              # X
        if not no_fxc:
            blk[i, j] += rows[:, 4 + 4 * b] + 1j * rows[:, 5 + 4 * b]         # F


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('transitions')
    ap.add_argument('wfiles', nargs='+')
    ap.add_argument('--tda', action='store_true', help='resonant block only')
    ap.add_argument('--no-fxc', action='store_true', help='drop the fxc part F (exchange only)')
    ap.add_argument('--nlow', type=int, default=8, help='eigenvalues to print')
    ap.add_argument('--bright', type=float, default=0.01, help='relative strength of a bright state')
    ap.add_argument('--full-bands', type=int, help='partly filled bands: also report the lowest bright interband '
                    'state (> 50%% weight on transitions out of bands <= this one) and its binding against the '
                    'lowest such transition')
    ap.add_argument('--spectrum', help='write eps2 up to a constant (Lorentzians) to this file')
    ap.add_argument('--eta', type=float, default=0.05, help='Lorentzian width for --spectrum [eV]')
    ap.add_argument('--wmax', type=float, default=4.0, help='upper energy of --spectrum [eV]')
    ap.add_argument('--eq0', help='write every eigenvalue (Re > 0) to this file as in the k-space program\'s Eq0: '
                    'index, eigenvalue - lowest transition, eigenvalue, relative strength')
    ap.add_argument('--emission', help='write an emission spectrum to this file: w^4 * absorption * occupation of '
                    'the states above the lowest eigenvalue E1 at temperature --temp')
    ap.add_argument('--temp', type=float, default=0.01, help='temperature of --emission [eV] (0.01 eV ~ 116 K)')
    ap.add_argument('--volo', help='also write Eq0, absorption_0 and emission_0 in this folder with the formulas '
                    'of the k-space Casida program (equal weights, 0.1 eV width, its grid and emission)')
    ap.add_argument('--occupation', choices=['boltzmann', 'bose'], default='boltzmann',
                    help='--emission factor: exp(-(w-E1)/T), or 1/(exp((w-E1)/T)-1) as the k-space program '
                    '(the latter diverges at w = E1)')
    a = ap.parse_args()

    tr, E, Ea, res, Aa, Ba = read_rows(a.transitions)
    n = len(tr)
    pairs = sorted({(int(v), int(c)) for v, c in tr[:, 7:9]})
    gap = E.min()
    print(f'{n} transitions, band pairs {pairs[:6]}{" ..." if len(pairs) > 6 else ""}; '
          f'lowest transition {gap:.4f} eV')
    if len(pairs) == 1:
        print('  !!! one band pair only: a sub-matrix of the Casida problem, not the BSE exciton')

    for wf in a.wfiles:
        z, M, coupled = read_matrix(wf, n, E, a.tda, a.no_fxc, E_ares=Ea)
        if coupled:
            A = np.concatenate([res, Aa])
            B = np.concatenate([np.conj(res), Ba])
        else:
            A, B = res, np.conj(res)
        dim = M.shape[0]
        lam, R = np.linalg.eig(M)
        del M
        # residue of each pole in B^T (w-M)^-1 A: (B^T R)_l (R^-1 A)_l
        strength = (B @ R) * np.linalg.solve(R, A)
        wfull = full_weight(tr, R, a.full_bands) if a.full_bands else None
        del R
        pos = np.where(lam.real > 1e-6)[0]
        pos = pos[np.argsort(lam[pos].real)]
        s = np.abs(strength[pos])
        s_rel = s / max(s.max(), 1e-300)
        kind = ('coupled' if coupled else 'TDA') + (', exchange only' if a.no_fxc else '')
        print(f'\n{wf}: z = {z.real:.4f} {z.imag:+.4f}i eV, {kind}, matrix {dim}')
        print('    lowest Re E [eV]   Im E [eV]   rel. strength')
        for k in range(min(a.nlow, len(pos))):
            print(f'    {lam[pos[k]].real:12.5f} {lam[pos[k]].imag:11.5f} {s_rel[k]:14.4e}')
        bright = np.where(s_rel > a.bright)[0]
        if len(bright):
            eb = lam[pos[bright[0]]].real
            print(f'    lowest bright: {eb:.4f} eV   binding (lowest transition - bright) = {gap - eb:.4f} eV')
        if a.full_bands:
            gi = E[tr[:, 7:9].min(axis=1) <= a.full_bands].min()
            bi = [k for k in range(len(pos)) if s_rel[k] > a.bright and wfull[pos[k]] > 0.5]
            if bi:
                ei = lam[pos[bi[0]]].real
                print(f'    lowest bright interband: {ei:.4f} eV   binding (interband gap {gi:.4f} - it) = {gi - ei:.4f} eV')
            else:
                print(f'    no bright interband state (interband gap {gi:.4f} eV)')
        tag = '' if len(a.wfiles) == 1 else '.' + wf.rsplit('_w', 1)[-1]
        if a.eq0:
            with open(a.eq0 + tag, 'w') as f:
                f.write(f'# {kind}, M at z = {z.real:.6f} eV; lowest transition {gap:.6f} eV\n'
                        f'#   i   E_i - gap [eV]    E_i [eV]   rel. strength\n')
                for i, k in enumerate(pos, 1):
                    f.write(f'{i:6d} {lam[k].real - gap:14.6f} {lam[k].real:12.6f} {s_rel[i - 1]:14.4e}\n')
            print(f'    eigenvalues written to {a.eq0 + tag} ({len(pos)})')
        if a.volo:
            volo_files(lam[pos].real, gap, a.volo + tag)
            print(f'    Eq0, absorption_0, emission_0 of the k-space program written to {a.volo + tag}/')
        if a.spectrum or a.emission:
            w = np.linspace(0.0, a.wmax, 801)
            eps2 = np.zeros_like(w)
            for k in pos:
                eps2 += (-(strength[k] / (w + 1j * a.eta - lam[k]))).imag
        if a.spectrum:
            np.savetxt(a.spectrum + tag, np.column_stack([w, eps2]), fmt='%14.6e',
                       header=f'eps2 up to a constant, Lorentzians eta={a.eta} eV, M at z={z}')
            print(f'    spectrum written to {a.spectrum + tag}')
        if a.emission:
            # van Roosbroeck-Shockley-type emission from the absorption: w^4 x absorption x the thermal
            # occupation of the excitations above the lowest one, E1 (the k-space program's W(1))
            e1 = lam[pos[0]].real
            x = (w - e1) / a.temp
            if a.occupation == 'boltzmann':
                occ = np.exp(-np.clip(x, 0.0, 700.0))
            else:
                occ = np.where(x > 0, 1.0 / np.expm1(np.clip(x, 1e-12, 700.0)), 0.0)
            em = np.where(w >= e1, w ** 4 * np.abs(eps2) * occ, 0.0)
            np.savetxt(a.emission + tag, np.column_stack([w, em]), fmt='%14.6e',
                       header=f'emission up to a constant: w^4 * absorption * {a.occupation} occupation, '
                              f'T = {a.temp} eV, E1 = {e1:.6f} eV, eta = {a.eta} eV')
            print(f'    emission written to {a.emission + tag} (E1 = {e1:.4f} eV, T = {a.temp} eV)')


if __name__ == '__main__':
    main()
