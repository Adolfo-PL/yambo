# Finite-q `ndb.Chi` to transition-space kernel

`BSKmod="CHI"` reads the zero-frequency, unsymmetrized `FXC_Q_iq` from an
existing `ndb.Chi` and adds its projection to the Hartree part of Yambo's
transition-space matrix. The number of G components resolved from `BSENGfxc`
must equal the stored `CHI_PARS_1` G-basis size. The stored G vectors and
q-point database must match the current run.
Use Kohn–Sham transition energies: `KfnQPdb` is rejected because the current
`fxc = chi0^-1 - P_G^-1` already represents the quasiparticle bubble shift.

For a transition Fourier column `d_t(G)` built from Yambo's bare wavefunction
overlaps, the added matrix element is

    spin_occ / (cell_volume * Nk) * d_i^† fxc(q, 0) d_j .

This is the double real-space contraction with the inverse Fourier kernel

    f_q(r,r') = (1/cell_volume) sum_G,G'
                 exp(i(q+G)r) fxc_GG'(q,0) exp(-i(q+G')r') .

The production path uses Yambo's Fourier transition densities and contracts
the G coefficients without allocating an `Nr × Nr` array. The test below also
applies the inverse transform explicitly, one transition at a time, and checks
it against both reciprocal projection and direct double-grid quadrature.

Current limits: finite q only; scalar, unpolarized, serial CPU calculations;
first frequency exactly zero; Hermitian static kernel; matrix BSE solver.
The `ChiGMode="LEVELS"` response has a separate optical dipole-limit route,
but `G0`, `COHSEX`, and `DYSON` response bubbles do not yet export the q=0
head and wings needed for an optical Casida calculation. A static W can still
produce a frequency-dependent `fxc`; this mode deliberately uses only its
zero-frequency slice. It does not turn the current bubble-derived kernel into
an excitonic GW vertex.

Run the isolated numerical check with `python run_tests.py --fc <gfortran>`.
It checks complex off-diagonal finite-q kernels at two frequencies.
