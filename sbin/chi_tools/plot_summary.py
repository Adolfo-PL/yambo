#!/usr/bin/env python3
"""
plot_summary.py -- figures for the Sigma -> G -> chi -> fxc -> Casida chain.

Each figure is made only when its inputs are given; all are PNG files in --outdir.

  energies.png         KS and QP energies (o-<gw>.qp): the QP correction of every state,
                       and the vb->cb gap at every k with the exciton and the static-fxc
                       Casida level drawn on top                     --qp [--vb 13]
  spectra.png          eps2: independent QP transitions, Casida with the static fxc(0),
                       Casida with fxc(w) (= BSE)                     --eps LABEL=FILE ...
  fxc_omega.png        the exported kernel fxc(G,G';w) of ndb.Chi, in units of the Coulomb
                       interaction, ft(G,G') = fxc(G,G') |q+G||q+G'| / 4pi:
                       the head as alpha(w) = -fxc(0,0) q^2 (the LRC convention fxc = -alpha/q^2,
                       as read_ndb_Chi.py reports it), one diagonal body element per shell of
                       |q+G|, and ft(G,G') at w = 0                           --chi DIR
  fkq.png              F_kq = the fxc part of the Casida matrix K_tt' (BSEChiKoutW, rr block):
                       the matrix ordered by transition energy, its diagonal next to the
                       exchange X_tt, and the diagonal over the Brillouin zone
                                                       --ktt TRANSITIONS WFILE [--pair 13 14]
  selfconsistency.png  lambda(w) of sc_scan.py with the line lambda = w      --scan sc_scan.dat

Reference levels (eV) come from the [BSE] lines of r- files, or are given by hand:
  --bse-report r-bse_sex_exc_*   (exciton, QP gap)   --static-report r-cas_st_exc_*
  --exciton 2.2342 --gap 2.5433 --static 2.528

Example (06_TEST_VAL, after binding.sh with KOUTW="2.23 | 2.23"):
  python3 plot_summary.py --qp g0w0/o-g0w0.qp --vb 13 \\
     --eps "QP transitions=o-cas_st_exc.eps_q1_diago_bse:ip" \\
     --eps "Casida, static fxc(0)=o-cas_st_exc.eps_q1_diago_bse" \\
     --eps "Casida, fxc(w) = BSE=o-cas_dyn_exc.eps_q1_chidyn_bse" \\
     --chi kx_exc --ktt o-cas_dyn_exc.Ktt_q1_transitions o-cas_dyn_exc.Ktt_q1_w225 \\
     --bse-report r-bse_sex_exc_* --static-report r-cas_st_exc_* --scan sc_scan.dat

A file given as FILE:ip plots its independent-particle columns (4, 5 of a yambo eps file).
--hex draws the zone of a hexagonal lattice (b1, b2 at 60 degrees) instead of reduced units.
"""
import argparse
import glob
import os
import sys

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
except ImportError:
    sys.exit('need matplotlib:  pip install --user matplotlib')

# ---- one palette for every figure -------------------------------------------------------
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']      # categorical, fixed order
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
SURFACE = '#fcfcfb'
DIVERGING = LinearSegmentedColormap.from_list(
    'blue_red', ['#184f95', '#5598e7', '#f0efec', '#ec8a89', '#b52f2e'])
SEQUENTIAL = LinearSegmentedColormap.from_list('blue', ['#cde2fb', '#6da7ec', '#256abf', '#0d366b'])

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'axes.titlecolor': INK,
    'axes.titlesize': 11, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
    'axes.labelsize': 9.5, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2, 'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6,
    'legend.frameon': False, 'legend.fontsize': 8.5, 'legend.labelcolor': INK2,
    'font.family': 'sans-serif', 'text.color': INK, 'axes.axisbelow': True,
})


def refline(ax, x, label, vertical=True, ls='--', color=MUTED, y=0.97, left=False):
    """A labelled reference level in muted ink (left=True puts a vertical label on the left)."""
    if x is None:
        return
    if vertical:
        ax.axvline(x, color=color, lw=1.0, ls=ls, zorder=1)
        ax.text(x, y, f'{label} ' if left else f' {label}', transform=ax.get_xaxis_transform(),
                color=INK2, fontsize=8, va='top', ha='right' if left else 'left', rotation=90)
    else:
        ax.axhline(x, color=color, lw=1.0, ls=ls, zorder=1)
        ax.text(0.99, x, f'{label} ', transform=ax.get_yaxis_transform(), color=INK2,
                fontsize=8, va='bottom', ha='right')


def save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {path}')


# ---- inputs --------------------------------------------------------------------------------
def report_value(patterns, key):
    for pat in patterns or []:
        for f in sorted(glob.glob(pat)):
            for line in open(f, errors='replace'):
                if key in line:
                    return float(line.split(':')[-1].split()[0])
    return None


def load_cols(path):
    return np.loadtxt(path, comments='#', ndmin=2)


def read_rr(path, n, chunk_lines=500000):
    """X_rr and F_rr of a BSEChiKoutW w-file (the other blocks are skipped), and z."""
    X = np.zeros((n, n), complex)
    F = np.zeros((n, n), complex)
    z = None
    ncol = None

    def scatter(lines):
        rows = np.fromstring(''.join(lines), sep=' ').reshape(-1, ncol)
        i = rows[:, 0].astype(int) - 1
        j = rows[:, 1].astype(int) - 1
        X[i, j] = rows[:, 2] + 1j * rows[:, 3]
        F[i, j] = rows[:, 4] + 1j * rows[:, 5]

    with open(path) as fh:
        buf = []
        for line in fh:
            if line.startswith('#'):
                if 'z =' in line:
                    re_z, im_z = line.split('z =')[1].split()[:2]
                    z = complex(float(re_z), float(im_z))
                continue
            if ncol is None:
                ncol = len(line.split())
            buf.append(line)
            if len(buf) >= chunk_lines:
                scatter(buf)
                buf = []
        if buf:
            scatter(buf)
    return z, X, F


def read_chi(dirname, iq):
    """Header and q fragment of an ndb.Chi (layout of src/io/io_Chi.F).

    netCDF reverses Fortran's (2,ng,ng,nw): arrays arrive as (nw,ng,ng,2) with (Re,Im)
    last and each frequency slice transposed, which is undone here ([G,G'] order)."""
    from netCDF4 import Dataset

    def cplx(v):
        x = np.asarray(v[:])
        return x[..., 0] + 1j * x[..., 1]

    with Dataset(os.path.join(dirname, 'ndb.Chi')) as ds:
        pars = np.asarray(ds.variables['CHI_PARS_1'][:]).ravel()
        h = {'ng': int(round(pars[0])), 'qpg': cplx(ds.variables['CHI_QPG']).T}
    for path in (os.path.join(dirname, f'ndb.Chi_fragment_{iq}'), os.path.join(dirname, 'ndb.Chi')):
        if not os.path.exists(path):
            continue
        with Dataset(path) as ds:
            if f'FXC_Q_{iq}' not in ds.variables:
                continue
            q = {'freqs': cplx(ds.variables[f'CHI_FREQ_Q_{iq}'])}
            for key, name in (('chi0', 'CHI0'), ('P', 'P'), ('fxc', 'FXC')):
                q[key] = np.transpose(cplx(ds.variables[f'{name}_Q_{iq}']), (0, 2, 1))
            if f'CHI_RCOND_Q_{iq}' in ds.variables:
                q['rcond'] = np.asarray(ds.variables[f'CHI_RCOND_Q_{iq}'][:])
            return h, q
    raise FileNotFoundError(f'no FXC_Q_{iq} in {dirname}')


# ---- figures -------------------------------------------------------------------------------
def fig_energies(a, lev, outdir):
    d = load_cols(a.qp)                       # k, band, Eo, E-Eo, (Sc)
    k, b, e0, de = d[:, 0].astype(int), d[:, 1].astype(int), d[:, 2], d[:, 3]
    vb = a.vb or int(b[e0 <= max(e0[e0 <= 0.0], default=0.0) + 1e-9].max())
    cb = vb + 1
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={'width_ratios': [1, 1.25]})

    # states within 8 eV of the gap (the semicore levels would squeeze the rest)
    vbm = e0[b == vb].max()
    cbm = e0[b == cb].min()
    near = (e0 > vbm - 8.0) & (e0 < cbm + 8.0)
    val, con = (b <= vb) & near, (b >= cb) & near
    ax1.scatter(e0[val], de[val], s=10, color=SERIES[0], label='valence', lw=0)
    ax1.scatter(e0[con], de[con], s=10, color=SERIES[1], label='conduction', lw=0)
    ax1.axhline(0, color=AXIS, lw=0.8)
    ax1.set_xlabel('Kohn-Sham energy $E_0$ [eV]')
    ax1.set_ylabel('QP correction $E-E_0$ [eV]')
    ax1.set_title('GW correction of the states near the gap')
    ax1.legend(loc='best')

    ks = np.array(sorted(set(k)))
    g0, g = [], []
    for kk in ks:
        v = (k == kk) & (b == vb)
        c = (k == kk) & (b == cb)
        if v.any() and c.any():
            g0.append(e0[c][0] - e0[v][0])
            g.append(e0[c][0] + de[c][0] - e0[v][0] - de[v][0])
    g0, g = np.array(g0), np.array(g)
    order = np.argsort(g0)
    x = np.arange(len(order))
    ax2.plot(x, g0[order], 'o-', color=SERIES[0], ms=4, label=f'KS {vb}$\\to${cb}')
    ax2.plot(x, g[order], 'o-', color=SERIES[1], ms=4, label=f'QP {vb}$\\to${cb}')
    ax2.set_xticks(x)
    ax2.set_xticklabels([str(ks[i]) for i in order], fontsize=7)
    ax2.set_xlabel('IBZ k point (sorted by the KS gap)')
    ax2.set_ylabel('transition energy [eV]')
    ax2.set_title(f'Band gap {vb}$\\to${cb} at each k, and the excitations')
    refline(ax2, lev.get('exciton'), f'BSE exciton {lev["exciton"]:.3f}' if lev.get('exciton') else '',
            vertical=False, ls='-', color=SERIES[2])
    refline(ax2, lev.get('static'), f'Casida static fxc(0) {lev["static"]:.3f}' if lev.get('static') else '',
            vertical=False)
    ax2.legend(loc='upper left')
    print(f'  energies: gap {vb}->{cb}: KS min {g0.min():.4f} eV, QP min {g.min():.4f} eV '
          f'(k {ks[np.argmin(g)]}), QP shift at that k {g[np.argmin(g)] - g0[np.argmin(g)]:+.4f} eV')
    save(fig, outdir, 'energies.png')


def fig_spectra(a, lev, outdir):
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for i, item in enumerate(a.eps):
        label, path = item.rsplit('=', 1)
        ip = path.endswith(':ip')
        path = path[:-3] if ip else path
        d = load_cols(path)
        y = d[:, 3] if ip else d[:, 1]
        ls = ':' if ip else '-'
        ax.plot(d[:, 0], y, ls, color=SERIES[i % len(SERIES)], label=label)
        pk = [j for j in range(1, len(y) - 1) if y[j] > y[j - 1] and y[j] >= y[j + 1] and y[j] > 0.1 * y.max()]
        print(f'  spectra: {label:40s} first peaks ' + ', '.join(f'{d[j, 0]:.3f}' for j in pk[:3]) + ' eV')
    refline(ax, lev.get('gap'), f'QP gap {lev["gap"]:.3f}' if lev.get('gap') else '')
    refline(ax, lev.get('exciton'), f'exciton {lev["exciton"]:.3f}' if lev.get('exciton') else '', left=True)
    if a.emin is not None or a.emax is not None:
        ax.set_xlim(a.emin, a.emax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('photon energy [eV]')
    ax.set_ylabel(r'$\varepsilon_2$')
    ax.set_title('Absorption')
    ax.legend(loc='upper right')
    save(fig, outdir, 'spectra.png')


def fig_fxc(a, lev, outdir):
    try:
        h, q = read_chi(a.chi, a.iq)
    except ImportError:
        print('  !!! fxc_omega: needs netCDF4 (pip install --user netCDF4)')
        return
    HA = 27.211386
    w = q['freqs'] * HA
    qpg = h['qpg'][a.iq - 1][:h['ng']]
    ft = q['fxc'] * (qpg[None, :, None] * qpg[None, None, :]) / (4 * np.pi)
    o = np.argsort(w.real)
    w, ft = w[o], ft[o]
    eta = np.median(np.abs(w.imag[w.real > 1e-6])) if np.any(w.real > 1e-6) else 0.0
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw={'width_ratios': [1.2, 1.2, 1]})

    alpha = -q['fxc'][o][:, 0, 0] * qpg[0] ** 2      # LRC convention: fxc(0,0) = -alpha/q^2
    ax1.plot(w.real, alpha.real, color=SERIES[0], label=r'Re $\alpha$')
    ax1.plot(w.real, alpha.imag, color=SERIES[1], label=r'Im $\alpha$')
    ax1.axhline(0, color=AXIS, lw=0.8)
    refline(ax1, lev.get('gap'), 'QP gap')
    refline(ax1, lev.get('exciton'), 'exciton', left=True)
    ax1.set_xlabel(r'$\omega$ [eV]' + (f'  (z = $\\omega$ + {eta:.2g}i)' if eta > 0 else ''))
    ax1.set_ylabel(r'$\alpha(\omega) = -f_{xc}(0,0)\,q^2$')
    ax1.set_title(r'Head: long-range $-\alpha/q^2$ part')
    ax1.legend(loc='upper left')
    lim = np.percentile(np.abs(alpha), 97) * 1.4
    ax1.set_ylim(-lim, lim)

    ng = ft.shape[1]
    # the first G of each of the first three shells of |q+G| (a shell shares one curve)
    shells, seen = [], []
    for G in range(1, ng):
        m = abs(qpg[G])
        if all(abs(m - s0) > 1e-3 * m for s0 in seen):
            seen.append(m)
            shells.append(G)
    shells = shells[:3]
    for i, G in enumerate(shells):
        ax2.plot(w.real, ft[:, G, G].real, color=SERIES[i], label=f'G = {G + 1} (|q+G| = {abs(qpg[G]):.2f})')
    ax2.axhline(0, color=AXIS, lw=0.8)
    refline(ax2, lev.get('gap'), 'QP gap')
    refline(ax2, lev.get('exciton'), 'exciton', left=True)
    if shells:
        body = np.array([ft[:, G, G].real for G in shells])
        lo, hi = np.percentile(body, 3), np.percentile(body, 97)
        pad = 0.4 * (hi - lo)
        ax2.set_ylim(lo - pad, hi + pad)
        ax2.legend(loc='lower left')
    else:
        ax2.text(0.5, 0.5, 'one plane wave: no body', transform=ax2.transAxes, ha='center', color=INK2)
    ax2.set_xlabel(r'$\omega$ [eV]')
    ax2.set_ylabel(r'Re $f_{xc}(G,G;\omega)\,|q+G|^2/4\pi$')
    ax2.set_title('Body: diagonal, one G per shell, in units of v(G)')

    i0 = int(np.argmin(np.abs(w)))
    m = ft[i0].real
    v = np.abs(m).max()
    im = ax3.imshow(m, cmap=DIVERGING, norm=TwoSlopeNorm(0, -v, v), origin='upper')
    ax3.grid(False)
    ax3.set_xticks(range(ng))
    ax3.set_yticks(range(ng))
    ax3.set_xticklabels(range(1, ng + 1))
    ax3.set_yticklabels(range(1, ng + 1))
    ax3.set_xlabel("G'")
    ax3.set_ylabel('G')
    ax3.set_title(f'Re ft(G,G\') at $\\omega$ = {w[i0].real:.2f} eV')
    cb = fig.colorbar(im, ax=ax3, shrink=0.85)
    cb.outline.set_visible(False)
    print(f'  fxc: alpha(0) = {alpha[i0].real:.4f}; ng = {ng}; {len(w)} frequencies '
          f'{w.real.min():.2f}..{w.real.max():.2f} eV')
    # conditioning: f_xc ~ delta/lambda^2 along small eigenvalues lambda of D chi0 D
    d = np.sqrt(4 * np.pi) / qpg
    lam = np.sort(np.abs(np.linalg.eigvals(d[:, None] * q['chi0'][o][i0] * d[None, :])))
    msg = f'  fxc: |eigenvalues| of D chi0 D at w = 0 span {lam[0]:.1e} .. {lam[-1]:.1e}'
    if 'rcond' in q:
        msg += f'; stored min rcond (chi0, P) {q["rcond"][:, 0].min():.1e}, {q["rcond"][:, 1].min():.1e}'
    print(msg)
    if lam[0] < 1e-6 * lam[-1]:
        print('  !!! fxc: nearly dependent plane waves: the body of fxc (up to '
              f'{np.abs(ft[i0][1:, 1:]).max():.1e} v) is set by near-null directions of chi0.')
        print('      The head and the optical response are unaffected; for the kernel itself use fewer G.')
    save(fig, outdir, 'fxc_omega.png')


def fig_fkq(a, lev, outdir):
    tr = load_cols(a.ktt[0])
    n = len(tr)
    z, X, F = read_rr(a.ktt[1], n)
    E, v, c = tr[:, 9], tr[:, 7].astype(int), tr[:, 8].astype(int)
    kr = tr[:, 4:7]
    MEV = 1000.0
    fig = plt.figure(figsize=(15, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 1])
    ax1, ax2, ax3 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

    o = np.argsort(E)
    m = F.real[np.ix_(o, o)] * MEV
    vmax = np.percentile(np.abs(m), 99.5)
    im = ax1.imshow(m, cmap=DIVERGING, norm=TwoSlopeNorm(0, -vmax, vmax), interpolation='nearest')
    ax1.grid(False)
    ax1.set_xlabel("t' (sorted by $E_{t'}$)")
    ax1.set_ylabel('t (sorted by $E_t$)')
    ax1.set_title("Re F$_{kq}$(t,t') [meV]")
    fig.suptitle(f'Casida matrix at z = {z.real:.3f} {z.imag:+.3f}i eV: fxc part F$_{{kq}}$ = K$_{{tt\'}}$ - X$_{{tt\'}}$',
                 x=0.01, ha='left', fontsize=11, color=INK2, y=1.02)
    cb = fig.colorbar(im, ax=ax1, shrink=0.85)
    cb.outline.set_visible(False)

    dF, dX = np.diag(F).real * MEV, np.diag(X).real * MEV
    ax2.scatter(E, dF, s=12, color=SERIES[0], lw=0, label=r'Re F$_{tt}$ (fxc, attraction)')
    ax2.scatter(E, dX, s=12, color=SERIES[1], lw=0, label=r'Re X$_{tt}$ (exchange)')
    ax2.axhline(0, color=AXIS, lw=0.8)
    refline(ax2, lev.get('exciton'), 'exciton')
    ax2.set_xlabel('transition energy $E_t$ [eV]')
    ax2.set_ylabel('diagonal element [meV]')
    ax2.set_title('Diagonal elements')
    ax2.legend(loc='best')

    pair = tuple(a.pair) if a.pair else (v.min() if len(set(v)) else 0, c.min() if len(set(c)) else 0)
    if a.pair is None and len(set(zip(v, c))) > 1:
        # the pair of the lowest transition
        j = np.argmin(E)
        pair = (v[j], c[j])
    sel = (v == pair[0]) & (c == pair[1])
    kk = (kr[sel] + 0.5) % 1.0 - 0.5                 # fold to [-1/2, 1/2)
    if np.abs(kr).max() < 1e-8:
        # dump written before the k coordinates were added (ec600a8): per IBZ point instead
        ik = tr[sel, 3].astype(int)
        ax3.scatter(ik, dF[sel], s=14, color=SERIES[0], lw=0)
        ax3.axhline(0, color=AXIS, lw=0.8)
        ax3.set_xlabel('IBZ k point of the transition')
        ax3.set_ylabel(r'Re F$_{tt}$ [meV]')
        ax3.set_title(f'Diagonal per IBZ k, {pair[0]}$\\to${pair[1]}')
        print('  fkq: no k coordinates in the transitions file (dump older than ec600a8): '
              'zone map replaced by F_tt per IBZ k')
        save_fkq(fig, a, lev, outdir, n, z, E, dF, dX, F)
        return
    if a.hex:
        b1, b2 = np.array([1.0, 0.0]), np.array([0.5, np.sqrt(3) / 2])
        xy = kk[:, :1] * b1 + kk[:, 1:2] * b2
        xl, yl = r'$k_x$ [$b$]', r'$k_y$ [$b$]'
    else:
        xy = kk[:, :2]
        xl, yl = '$k_1$ [rlu]', '$k_2$ [rlu]'
    val = dF[sel]
    vm = np.abs(val).max()
    sc = ax3.scatter(xy[:, 0], xy[:, 1], c=val, cmap=DIVERGING, norm=TwoSlopeNorm(0, -vm, vm), s=38,
                     marker='h' if a.hex else 's', lw=0)
    ax3.set_aspect('equal')
    ax3.set_xlabel(xl)
    ax3.set_ylabel(yl)
    ax3.set_title(f'Re F$_{{tt}}$ over the zone, {pair[0]}$\\to${pair[1]} [meV]', fontsize=10)
    cb = fig.colorbar(sc, ax=ax3, shrink=0.85)
    cb.outline.set_visible(False)
    save_fkq(fig, a, lev, outdir, n, z, E, dF, dX, F)


def save_fkq(fig, a, lev, outdir, n, z, E, dF, dX, F):
    ex = lev.get('exciton')
    print(f'  fkq: {n} transitions, z = {z}; diagonal F: mean {dF.mean():.2f} meV, min {dF.min():.2f} meV '
          f'at E_t = {E[np.argmin(dF)]:.3f} eV; diagonal X: mean {dX.mean():.2f} meV'
          + (f'; exciton line {ex:.3f} eV' if ex else ''))
    ev = np.sort(np.linalg.eigvals(F).real) * 1000.0
    nz = int(np.sum(np.abs(ev) > 1e-3 * np.abs(ev).max()))
    print(f'  fkq: trace of Re F {np.trace(F).real * 1000:.1f} meV carried by {nz} non-zero modes; '
          f'most attractive {", ".join(f"{e:.1f}" for e in ev[:min(4, nz)])} meV (rank <= BSENGfxc)')
    herm = np.abs(F - F.conj().T).max() / max(np.abs(F).max(), 1e-300)
    print(f'  fkq: |F - F^H| / max|F| = {herm:.2e} (non-zero: fxc(z) at complex z is not Hermitian)')
    save(fig, outdir, 'fkq.png')


def fig_scan(a, lev, outdir):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sc_scan import crossings
    r = np.loadtxt(a.scan, ndmin=2)
    w, lam = r[:, 0], r[:, 1]
    xs, poles = crossings(w, lam)
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(w, w, '--', color=MUTED, lw=1.0)
    ax.text(w[0], w[0], r'  $\lambda=\omega$', color=INK2, fontsize=8.5, va='bottom')
    # do not draw a line across a jump of lambda (a pole of the kernel): a steep rise, or a
    # steep step right after one (the hop back to the branch); a steep fall alone is the crossing
    s = np.diff(lam) / np.diff(w)
    up = s > 10.0
    cut = [j + 1 for j in range(len(s)) if up[j] or (j > 0 and up[j - 1] and abs(s[j]) > 10.0)]
    seg = np.split(np.arange(len(w)), cut)
    for s in seg:
        ax.plot(w[s], lam[s], 'o-', color=SERIES[0], ms=4.5)
    ax.plot([], [], 'o-', color=SERIES[0], ms=4.5, label=r'lowest bright eigenvalue $\lambda(\omega)$ of M($\omega$)')
    refline(ax, lev.get('exciton'), f'BSE exciton {lev["exciton"]:.4f}' if lev.get('exciton') else '',
            vertical=False, ls='-', color=SERIES[2])
    for x in xs:
        ax.plot([x], [x], 'o', ms=9, mfc='none', mec=INK, mew=1.4)
        ax.annotate(f'self-consistent\n{x:.4f} eV', (x, x), xytext=(30, 45), textcoords='offset points',
                    fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
    lo, hi = np.percentile(lam, 2), np.percentile(lam, 98)
    lo, hi = min(lo, w.min()), max(hi, w.max())
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    ax.set_xlabel(r'frequency $\omega$ of the kernel $f_{xc}(\omega)$ [eV]')
    ax.set_ylabel(r'$\lambda(\omega)$ [eV]')
    ax.set_title(r'Self-consistency: the exciton is where $\lambda(\omega)=\omega$')
    ax.legend(loc='best')
    print(f'  scan: crossing(s) {", ".join(f"{x:.4f}" for x in xs) or "none"} eV'
          + (f'; jumps (kernel poles) near {", ".join(f"{x:.3f}" for x in poles)} eV' if poles else ''))
    save(fig, outdir, 'selfconsistency.png')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--outdir', default='summary_plots')
    ap.add_argument('--qp', help='o-<gw>.qp')
    ap.add_argument('--vb', type=int, help='highest valence band (default: from E0 <= 0)')
    ap.add_argument('--eps', action='append', help='LABEL=FILE[:ip], repeatable (split at the last =)')
    ap.add_argument('--emin', type=float)
    ap.add_argument('--emax', type=float)
    ap.add_argument('--chi', help='directory with ndb.Chi')
    ap.add_argument('--iq', type=int, default=1)
    ap.add_argument('--ktt', nargs=2, metavar=('TRANSITIONS', 'WFILE'))
    ap.add_argument('--pair', type=int, nargs=2, metavar=('V', 'C'), help='band pair of the zone map')
    ap.add_argument('--hex', action='store_true', help='hexagonal zone (b1, b2 at 60 degrees)')
    ap.add_argument('--scan', help='sc_scan.dat of sc_scan.py')
    ap.add_argument('--bse-report', nargs='+', help='r- file(s) of the BSE run')
    ap.add_argument('--static-report', nargs='+', help='r- file(s) of the static-fxc Casida run')
    ap.add_argument('--exciton', type=float)
    ap.add_argument('--gap', type=float)
    ap.add_argument('--static', type=float)
    a = ap.parse_args()

    lev = {
        'exciton': a.exciton or report_value(a.bse_report, 'lowest bright excitation'),
        'gap': a.gap or report_value(a.bse_report, 'lowest transition energy of the window'),
        'static': a.static or report_value(a.static_report, 'lowest bright excitation'),
    }
    print('levels [eV]: ' + ', '.join(f'{k} {v:.4f}' for k, v in lev.items() if v is not None))
    os.makedirs(a.outdir, exist_ok=True)
    done = False
    for flag, fn in (('qp', fig_energies), ('eps', fig_spectra), ('chi', fig_fxc),
                     ('ktt', fig_fkq), ('scan', fig_scan)):
        if getattr(a, flag):
            fn(a, lev, a.outdir)
            done = True
    if not done:
        ap.print_help()


if __name__ == '__main__':
    main()
