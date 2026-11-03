#!/usr/bin/env python3
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
"""
casida_spectra.py -- absorption at q1 from the Casida matrices of the exact kernel, against the BSE.

From one transitions file and the w-files of two frequencies (BSEChiKoutW, BSKmod="CHI" + BSEChiDyn):

  chi0      independent QP transitions                         M = diag(E)
  chi       + exchange (Hartree local fields, in-plane G)      M = diag(E) + X
  fxc(0)    + the exact kernel frozen at w = 0                 M = diag(E) + X + F(0)
  fxc(E_A)  + the exact kernel frozen at the A exciton         M = diag(E) + X + F(E_A)

each a sum of Lorentzians over the poles of (w - M)^-1 with the residues of ktt_casida.py, and the
BSE spectrum Yambo wrote for the same model (o-<job>.alpha_q1_* or eps). The Casida curves are scaled
to Yambo's units by fitting their chi0 to the independent-particle column of the Yambo file (the same
transitions and vertices), so all five share one axis.

  python casida_spectra.py o-xsp0.Ktt_q1_transitions o-xsp0.Ktt_q1_w1 o-xspA.Ktt_q1_w<i> \
         --bse o-bsp.alpha_q1_diago_bse [--eta 0.05] [--exciton 2.03] [--label WSe2] [--out casida_spectra]

Writes <out>.dat (w, Im of BSE, chi0, chi, fxc(0), fxc(E_A), ...) and three figures: <out>.png (all curves),
<out>_chi.png (chi0 on DFT and QP energies, chi) and <out>_tddft.png (BSE and the exact-kernel TDDFT curves).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ktt_casida as kc   # noqa: E402


def lorentz(w, eta, lam, strength):
    """Im of -sum_l s_l / (w + i eta - lam_l) over the poles with Re lam > 0."""
    out = np.zeros_like(w)
    for l, s in zip(lam, strength):
        if l.real > 1e-6:
            out += (-(s / (w + 1j * eta - l))).imag
    return out


def casida(M, A, B, w, eta):
    lam, R = np.linalg.eig(M)
    s = (B @ R) * np.linalg.solve(R, A)
    return lorentz(w, eta, lam, s)


def yambo_spectrum(path):
    """E, Im, Im of the independent-particle column, from a Yambo alpha/eps file (titles in the header)."""
    titles = None
    for line in open(path):
        if line.startswith('#') and ('E/ev' in line or 'E/eV' in line):
            titles = line.lstrip('#').split()
    x = np.loadtxt(path, comments='#', ndmin=2)
    io = 3
    if titles:
        for i, t in enumerate(titles):
            if '_o' in t and 'Im' in t:
                io = i
                break
    return x[:, 0], x[:, 1], x[:, io]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('transitions')
    ap.add_argument('w0', help='w-file at z = 0')
    ap.add_argument('wA', help='w-file at the frequency nearest the A exciton')
    ap.add_argument('--bse', required=True, help="Yambo's BSE spectrum of the same model (broadened)")
    ap.add_argument('--dyn', help="Yambo's spectrum of Casida with fxc(w) at every w (BSEChiDyn, same grid)")
    ap.add_argument('--ks', help="a Yambo BSE spectrum on KS energies: its independent-particle column is chi0 with DFT energies")
    ap.add_argument('--eta', type=float, default=0.05, help='Lorentzian width [eV], as the Yambo run')
    ap.add_argument('--exciton', type=float, help='A exciton of the BSE [eV] (marked)')
    ap.add_argument('--gap', type=float, help='gap [eV] (marked)')
    ap.add_argument('--emin', type=float, help='plot from (default: exciton - 0.35 eV)')
    ap.add_argument('--emax', type=float, help='plot to (default: gap + 0.75 eV)')
    ap.add_argument('--label', default='')
    ap.add_argument('--out', default='casida_spectra')
    a = ap.parse_args()

    tr, E, Ea, res, Aa, Ba = kc.read_rows(a.transitions)
    n = len(tr)
    w, bse, ip_y = yambo_spectrum(a.bse)
    dyn = ks = None
    if a.ks:
        wk, _, ksip = yambo_spectrum(a.ks)
        ks = np.interp(w, wk, ksip)
    if a.dyn:
        wd, dyn, _ = yambo_spectrum(a.dyn)
        dyn = np.interp(w, wd, dyn)
        dev = np.abs(dyn - bse).max() / max(np.abs(bse).max(), 1e-300)
        print(f'Casida with fxc(w) at every w vs BSE: max deviation {dev:.2%} of the BSE maximum'
              + ('' if dev < 0.01 else '  !!! should be ~0: same model, same grid?'))

    curves = {}
    zA = None
    for name, f in (('fxc0', a.w0), ('fxcA', a.wA)):
        z, M, coupled = kc.read_matrix(f, n, E, False, False, E_ares=Ea)
        if name == 'fxcA':
            zA = z.real
        A = np.concatenate([res, Aa]) if coupled else res
        B = np.concatenate([np.conj(res), Ba]) if coupled else np.conj(res)
        curves[name] = casida(M, A, B, w, a.eta)
        if name == 'fxc0':          # the same file without the kernel: chi; without exchange too: chi0
            _, MX, _ = kc.read_matrix(f, n, E, False, True, E_ares=Ea)
            curves['chi'] = casida(MX, A, B, w, a.eta)
            M0 = np.diag(np.concatenate([E, Ea]) if coupled else E).astype(complex)
            curves['chi0'] = casida(M0, A, B, w, a.eta)
        del M
    sc = float(np.dot(ip_y, curves['chi0']) / max(np.dot(curves['chi0'], curves['chi0']), 1e-300))
    for k in curves:
        curves[k] *= sc
    resid = np.abs(curves['chi0'] - ip_y).max() / max(np.abs(ip_y).max(), 1e-300)
    print(f'Casida chi0 scaled by {sc:.4e} onto the Yambo independent-particle column: max deviation {resid:.1%}'
          + ('' if resid < 0.05 else '  !!! the Casida chi0 does not match Yambo\'s: check the files are of one model'))
    print(f'kernel frozen at w = {zA:.4f} eV' + (f' (A exciton {a.exciton:.4f} eV)' if a.exciton else ''))
    for k, lab in (('chi0', 'chi0'), ('chi', 'chi'), ('fxc0', 'fxc(0)'), ('fxcA', 'fxc(E_A)')):
        i = np.argmax(curves[k])
        print(f'  strongest peak {lab:9s} {w[i]:.4f} eV')
    print(f'  strongest peak {"BSE":9s} {w[np.argmax(bse)]:.4f} eV')

    cols = [w, bse, curves['chi0'], curves['chi'], curves['fxc0'], curves['fxcA']] + ([dyn] if dyn is not None else []) \
        + ([ks] if ks is not None else [])
    np.savetxt(f'{a.out}.dat', np.column_stack(cols), fmt='%.6e',
               header=f'w[eV] Im BSE, chi0, chi, fxc(0), fxc(w={zA:.4f})' + (', fxc(w) at every w' if dyn is not None else '')
               + (', chi0 on DFT energies' if ks is not None else '')
               + ' (Casida curves scaled to the Yambo units)')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print('no plot:', e)
        return
    # every curve once: key -> (y, style, label); three figures draw subsets of them
    S = {}
    if ks is not None:
        S['ks'] = (ks, dict(color='#a8a39a', lw=1.0), r'$\chi_0$ (DFT energies)')
    S['chi0'] = (curves['chi0'], dict(color='#898781', lw=1.2, ls=':'), r'$\chi_0$ (independent QP transitions)')
    S['chi'] = (curves['chi'], dict(color='#52514e', lw=1.2, ls='--'), r'$\chi$ (+ Hartree local fields)')
    S['bse'] = (bse, dict(color='#2a78d6', lw=4, alpha=0.45), 'BSE')
    if dyn is not None:
        S['dyn'] = (dyn, dict(color='#0d366b', lw=1.3, ls='--'), r'TDDFT, exact $f_{xc}(\omega)$ at every $\omega$')
    S['fxc0'] = (curves['fxc0'], dict(color='#1baf7a', lw=1.5), r'TDDFT, exact $f_{xc}$ at $\omega=0$')
    S['fxcA'] = (curves['fxcA'], dict(color='#eb6834', lw=1.5, ls='-.'),
                 rf'TDDFT, exact $f_{{xc}}$ at $\omega=E_A$ ({zA:.3f} eV)')
    figs = (('', list(S), 'in-plane absorption'),
            ('_chi', [k for k in ('ks', 'chi0', 'chi') if k in S], 'independent particles and RPA'),
            ('_tddft', [k for k in ('bse', 'dyn', 'fxc0', 'fxcA') if k in S], 'BSE and TDDFT with the exact kernel'))
    for sfx, keys, what in figs:
        fig, ax = plt.subplots(figsize=(9.5, 4.6))
        for k in keys:
            y, st, lab = S[k]
            ax.plot(w, y, label=lab, **st)
        lo = a.emin if a.emin is not None else (a.exciton - 0.35 if a.exciton else w.min())
        if a.emin is None and 'ks' in keys:     # the DFT chi0 starts at the KS gap, far below the QP one
            on = np.where(ks > 0.05 * ks.max())[0]
            if len(on):
                lo = min(lo, max(w[on[0]] - 0.15, w.min()))
        hi = a.emax if a.emax is not None else (a.gap + 0.75 if a.gap else w.max())
        m = (w >= lo) & (w <= hi)
        top = max(S[k][0][m].max() for k in keys)
        for x, lab in ((a.exciton, 'exciton'), (a.gap, 'QP gap')):
            if x:
                ax.axvline(x, color='k', lw=0.8, ls=':')
                ax.annotate(f'{lab}\n{x:.3f} eV', (x, top * 1.02), xytext=(4, 0), textcoords='offset points',
                            va='top', fontsize=9)
        ax.set_xlabel(r'$\omega$ (eV)')
        ax.set_ylabel(r'Im $\alpha(\omega)$ ($q \to 0$)')
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, top * 1.12)
        ax.set_title(f'{a.label + ": " if a.label else ""}{what}, $q \\to 0$')
        ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
        fig.tight_layout()
        fig.savefig(f'{a.out}{sfx}.png', dpi=200)
        plt.close(fig)
        print('wrote', f'{a.out}{sfx}.png')
    print('wrote', f'{a.out}.dat')

if __name__ == '__main__':
    main()
