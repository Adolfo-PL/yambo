#!/usr/bin/env python3
# Static kernel f_xc(q, w=0) next to the kernel at the exciton frequency f_xc(q, w=E_X(q)).
#
# usage: python3 fxc_q_compare.py fxc_xq_static.npz fxc_GGq_scw.npz fxc_q_compare.png
#        (both from xq_dispersion.sh: the first is f_xc at w=0, the second at w=E_X(q))
#
# panels: (1) head: -f_xc^{00}(q), log scale (q1, at |q| = 1e-5, is printed, not plotted)
#         (2) body: -f_xc^{GG} averaged over the six in-plane G
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


def load(path):
    d = np.load(path)
    out = {}
    for iq, qpg, fx in zip(d['iq'], d['qpg'], d['fxc']):
        out[int(iq)] = (qpg[0], fx[0].real)          # |q| (G = 0), f_xc(G, G') at the one frequency
    return out


s, w = load(sys.argv[1]), load(sys.argv[2])
iqs = sorted(set(s) & set(w))
q = np.array([s[i][0] for i in iqs])
fs = np.array([s[i][1] for i in iqs])
fw = np.array([w[i][1] for i in iqs])
fin = q > 1e-4                                       # q1 is q -> 0 (|q| = 1e-5)

alpha_s = -fs[:, 0, 0] * q**2
alpha_w = -fw[:, 0, 0] * q**2
a0_s, a0_w = alpha_s[np.argmin(q)], alpha_w[np.argmin(q)]
body_s = -np.array([np.mean(np.diag(x)[1:]) for x in fs])
body_w = -np.array([np.mean(np.diag(x)[1:]) for x in fw])

fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
# q1 sits at |q| = 1e-5 where the head is ~1e10: left out of the plot, printed below
ax[0].plot(q[fin], -fs[fin, 0, 0], 's', mfc='none', label=r'$\omega=0$')
ax[0].plot(q[fin], -fw[fin, 0, 0], 'o', label=r'$\omega=E_X(q)$')
ax[0].set_yscale('log')
ax[0].set_xlabel(r'$|q|$ (bohr$^{-1}$)')
ax[0].set_ylabel(r'$-f_{xc}^{00}(q)$ (Ha bohr$^3$)')
ax[0].set_title('head (G = G\' = 0)')
ax[0].legend()
ax[1].plot(q, body_s, 's', mfc='none', label=r'$\omega=0$')
ax[1].plot(q, body_w, 'o', label=r'$\omega=E_X(q)$')
ax[1].set_xlabel(r'$|q|$ (bohr$^{-1}$)')
ax[1].set_ylabel(r'$-f_{xc}^{GG}$, mean of 6 in-plane $G$ (Ha bohr$^3$)')
ax[1].set_title('body diagonal')
ax[1].legend()
fig.tight_layout()
fig.savefig(sys.argv[3], dpi=150)

print(f'alpha_0: w=0 {a0_s:.4f}   w=E_X {a0_w:.4f}   (ratio {a0_w / a0_s:.2f})')
print(' iq    |q|     head(w=0)    head(E_X)   body(w=0)  body(E_X)')
for i, iq in enumerate(iqs):
    print(f'{iq:3d} {q[i]:7.4f} {-fs[i, 0, 0]:12.4e} {-fw[i, 0, 0]:12.4e} {body_s[i]:11.2f} {body_w[i]:10.2f}')
