<!-- License-Identifier: GPL. Copyright (C) 2026 Adolfo Partida. Authors: Adolfo Partida, with AI assistance (Claude Code). -->

# Two Casida solvers on the same exact kernel

This note compares the two ways the exact kernel of the fork has been used to compute excitons:

- **Solver Y**: `sbin/chi_tools/ktt_casida.py`, on the Casida matrices Yambo writes (`BSEChiKoutW`);
- **Solver V**: Volodymyr Turkowski's Fortran program (`eigenvalues_SPN.f90`, "bindingenergies"), fed by
  `sbin/chi_tools/ktt_to_kq.py`.

It writes down what each solver computes, shows that both diagonalize the same numbers, and derives
each difference between their results. Numbers: WS2 monolayer, no SOC, 24×24 k-grid, BSE bands 13–15
(N = 1152 transitions), kernel frozen at $z$ = 2.3795 eV (`o-xspA.Ktt_q1_w98`).

---

## 1. Notation

A transition is $t = (v \to c, \mathbf k)$, with

$$
E_t = \varepsilon_{c\mathbf k} - \varepsilon_{v\mathbf k}, \qquad
\rho_t(\mathbf G) = \langle c\mathbf k \,|\, e^{i(\mathbf q+\mathbf G)\cdot\mathbf r} \,|\, v\mathbf k\rangle ,
$$

the QP transition energy and the pair density. $a_t$ is the optical vertex of transition $t$
(the $q \to 0$ dipole, column "Re/Im $a_t$" of the transitions file). $N$ is the number of
transitions.

## 2. The exact kernel

The fork computes, at each frequency $z$ of a grid (`K_Chi_export.F`),

$$
f_{xc}(z) = \chi_0^{-1}(z) - P^{-1}(z), \qquad P^{-1} = \bar\chi_{\rm BSE}^{-1} + \bar v ,
$$

where both response functions are the **full** ones, resonant plus anti-resonant:

$$
\chi_0(z) = \sum_t \left[ \frac{\rho_t \rho_t^\dagger}{z - E_t + i\eta} - \frac{\rho_t^* \rho_t^{T}}{z + E_t + i\eta} \right]
\equiv \chi_0^{\rm res}(z) + \chi_0^{\rm ares}(z).
$$

With this kernel the TDDFT Dyson equation $\bar\chi^{-1} = \chi_0^{-1} - \bar v - f_{xc}$ gives
$\bar\chi^{-1} = P^{-1} - \bar v = \bar\chi_{\rm BSE}^{-1}$ exactly, at every $z$.

## 3. The kernel in the transition basis

The Casida run of the fork projects $\bar v$ and $f_{xc}(z)$ on the transitions and writes four blocks:

$$
X^{rr}_{tt'} + F^{rr}_{tt'} = \sum_{\mathbf G\mathbf G'} \rho_t^*(\mathbf G)\,\big[\bar v(\mathbf G)\delta_{\mathbf G\mathbf G'} + f_{xc}(\mathbf G,\mathbf G';z)\big]\,\rho_{t'}(\mathbf G') ,
$$

$$
X^{rc}_{tt'} + F^{rc}_{tt'} = \sum_{\mathbf G\mathbf G'} \rho_t^*(\mathbf G)\,\big[\bar v(\mathbf G)\delta_{\mathbf G\mathbf G'} + f_{xc}(\mathbf G,\mathbf G';z)\big]\,\rho_{t'}^*(\mathbf G') ,
$$

and the mirror blocks $cr$, $cc$. $X$ is the exchange (Hartree local fields), $F$ the kernel. The
$rr$ block couples excitations, the $rc$ block couples excitations to de-excitations; it comes from
$\chi_0^{\rm ares}$.

---

## 4. Solver Y (ours): coupled Casida

**Matrix** ($2N \times 2N$, non-Hermitian):

$$
M(z) = \begin{pmatrix} A & B \\ C & D \end{pmatrix}, \qquad
\begin{aligned}
A &= {\rm diag}(E) + X^{rr} + F^{rr}, & B &= X^{rc} + F^{rc},\\
C &= X^{cr} + F^{cr}, & D &= -{\rm diag}(E) + X^{cc} + F^{cc}
\end{aligned}
$$

(at $q \to 0$ with time reversal $C = -B^*$, $D = -A^*$).

**Eigenproblem** (`np.linalg.eig`): $M R = R\,\Lambda$, eigenvalues $\lambda_n$; the excitations
are the $N$ eigenvalues with ${\rm Re}\,\lambda_n > 0$.

**Response** (`ktt_casida.py`, the same form as Yambo's):

$$
\chi(\omega) = -C_0\, \mathcal B^{T} \left[\omega - M\right]^{-1} \mathcal A,
\qquad \mathcal A = \begin{pmatrix} a \\ i\,a^* \end{pmatrix},\;
\mathcal B = \begin{pmatrix} a^* \\ i\,a \end{pmatrix},
$$

so each eigenvalue carries the oscillator strength

$$
s_n = \big(\mathcal B^{T} R\big)_n \big(R^{-1}\mathcal A\big)_n .
$$

**Absorption** (`--spectrum`), Lorentzians of half-width $\eta$ (0.03 eV):

$$
\alpha_2(\omega) \propto \sum_{{\rm Re}\,\lambda_n>0} {\rm Im}\left[ \frac{-s_n}{\omega + i\eta - \lambda_n} \right]
= \sum_n \frac{s_n\,\eta}{(\omega - \lambda_n)^2 + \eta^2} \quad (\text{real } s_n, \lambda_n).
$$

**Emission** (`--emission`), from the lowest eigenvalue $E_1$ at temperature $T$:

$$
I(\omega) \propto \omega^4\, \alpha_2(\omega)\, e^{-(\omega - E_1)/T}\, \theta(\omega - E_1).
$$

**Eigenvalue list** (`--eq0`): $\lambda_n - \min_t E_t$.

---

## 5. Solver V (Volodymyr's program): Tamm–Dancoff

**Program equations** (k-space form, two bands $v$, $c$):

$$
H_{kq} = \big(E_c(k) - E_v(k)\big)\,\delta_{kq} + \frac{\rm Ha}{\epsilon}\, F_{kq}\,\sqrt{w_k w_q},
\qquad H U = U\,{\rm diag}(W), \quad (\texttt{ZHEEV})
$$

with ${\rm Ha} = 2 \times 13.605698$ eV and $\epsilon = 1$.

**What the converter gives it** (`ktt_to_kq.py`): every transition $t$ is one "$k$", with

$$
E_c(t) - E_v(t) = E_t, \qquad w_t = \frac1N, \qquad
F_{tt'} = \frac{X^{rr}_{tt'} + F^{rr}_{tt'}}{{\rm Ha}\,\sqrt{w_t w_{t'}}} ,
$$

so the program builds

$$
H = {\rm diag}(E) + X^{rr} + F^{rr} = A ,
$$

**the $rr$ block of Solver Y, and nothing else**: Tamm–Dancoff.

**Eigenvalue list** (`Eq0`): $W_n - \min_t E_t$.

**"Dipole" and absorption** (`absorption_0`), on $\omega = 0, 0.02, \dots, 20$ eV:

$$
d_n = \sum_t |U_{tn}|^2\, w_t , \qquad
A(\omega) = \sum_n \frac{d_n}{(\omega - W_n)^2 + \gamma^2}, \quad \gamma^2 = 0.01\ {\rm eV}^2 .
$$

Because the columns of $U$ are normalized and $w_t = 1/N$,

$$
d_n = \frac1N \sum_t |U_{tn}|^2 = \frac1N \quad \text{for every } n:
$$

**every eigenvalue gets the same weight.** $A(\omega)$ is the density of eigenvalues broadened by
$\gamma = 0.1$ eV, not an optical spectrum. (Checked against the program's output: equal weights
reproduce `absorption_0` to $5\times10^{-14}$.)

**Emission** (`emission_0`), with $W_1$ the lowest eigenvalue and $T = 0.01$ eV:

$$
B(\omega) = 4\pi \left| \frac{\omega^4 A(\omega)}{e^{(\omega - W_1)/T} - 1} \right| \theta(\omega - W_1)
\times \frac{0.8}{1600}.
$$

(Checked against `emission_0` to $5\times10^{-16}$.)

---

## 6. The differences, one at a time

### 6.1 Same matrix, same eigenvalues

Solver V diagonalizes $A$; Solver Y with `--tda` diagonalizes $A$ too. Their eigenvalues agree
exactly (rebuilding $H$ from the program's input files: $|W_n - \lambda_n^{\rm TDA}| < 5\times10^{-14}$
eV for all 1152). **The two codes are the same solver; every other difference below is a choice of
equations.**

### 6.2 Tamm–Dancoff vs coupled: $-0.21$ eV

Dropping $B$, $C$, $D$ is the same as solving the Dyson equation with $\chi_0^{\rm res}$ only:

$$
\bar\chi_{\rm TDA}^{-1} = \big(\chi_0^{\rm res}\big)^{-1} - \bar v - f_{xc}.
$$

Insert the exact kernel of §2:

$$
\bar\chi_{\rm TDA}^{-1}
= \big(\chi_0^{\rm res}\big)^{-1} - \bar v - \chi_0^{-1} + P^{-1}
= \bar\chi_{\rm BSE}^{-1} + \Big[\big(\chi_0^{\rm res}\big)^{-1} - \chi_0^{-1}\Big].
$$

The bracket vanishes only if $\chi_0^{\rm ares} = 0$. The kernel was built to cancel the full
$\chi_0^{-1}$; in TDA it over-corrects by the bracket. For the BSE itself TDA is a good approximation
because the attraction $W$ sits in the $rr$ block ($W^{rc}$ is small); a density kernel enters $rr$
and $rc$ through the same pair densities, so $F^{rc}$ is not small.

| | lowest bright | $\lambda - \min E_t$ |
|---|---|---|
| TDA (Solver V; Solver Y `--tda`) | 2.5785 eV | −0.345 eV |
| coupled (Solver Y) | 2.3689 eV | −0.555 eV |

### 6.3 Frozen frequency: +8.5 meV

The exact solution is the self-consistent $\omega$ with $\lambda(\omega) = \omega$: $E_A$ = 2.3774 eV.
The kernel here is frozen at $z$ = 2.3795 eV. Near the crossing,

$$
\lambda(z) \approx E_A + \lambda'\,(z - E_A)
\quad\Rightarrow\quad
\lambda' \approx \frac{2.3689 - 2.3774}{2.3795 - 2.3774} \approx -4 .
$$

A 2.1 meV offset in $z$ moves the eigenvalue by $|\lambda'| \times 2.1 \approx 8.5$ meV, in the
opposite direction. The same slope sets the weight of the A peak in the exact (frequency-dependent)
response, $s_A / (1 - \lambda')$: a kernel frozen at one $z$ puts A at the right energy only for
$z = E_A$, and always with the weight $s_A$, about $1 - \lambda' \approx 5$ times too large.

### 6.4 Oscillator strengths

In TDA ($R = U$ unitary, $\mathcal A = a$, $\mathcal B = a^*$), Solver Y weights eigenvalue $n$ with

$$
s_n = \Big| \sum_t a_t^*\, U_{tn} \Big|^2 ,
$$

the coherent sum of the transition dipoles: bright excitons are those whose transitions add in phase.
Solver V weights it with $d_n = 1/N$. So Solver V's peaks sit at the right energies (the same $W_n$)
but its heights carry no optical information: dark states (here the partner of A,
$s/s_A \sim 10^{-8}$) are as tall as A.

### 6.5 Line shape and grid

$$
{\rm Im}\frac{-s}{\omega + i\eta - \lambda} = \frac{s\,\eta}{(\omega-\lambda)^2 + \eta^2}
\qquad\text{vs}\qquad
\frac{d}{(\omega-W)^2 + \gamma^2}
$$

are the same Lorentzian shape. Solver Y uses $\eta$ = 0.03 eV on a 0.005 eV grid, Solver V
$\gamma$ = 0.1 eV on a 0.02 eV grid: Solver V's peaks are about three times wider.

### 6.6 Emission occupation

$$
\text{Y: } e^{-(\omega - E_1)/T}
\qquad\text{vs}\qquad
\text{V: } \frac{1}{e^{(\omega - W_1)/T} - 1} \;\xrightarrow{\;\omega \to W_1^+\;}\; \frac{T}{\omega - W_1}.
$$

The Bose form diverges at the lowest eigenvalue, so Solver V's emission spikes at its first grid
point above $W_1$; far above $W_1$ both decay as $e^{-(\omega-W_1)/T}$. Both use $E_1 = W_1$ = the
lowest eigenvalue, bright or dark. The factor $0.8/1600$ and the overall scale are arbitrary.

---

## 7. Bringing the results together without changing the program

1. **Replicate Solver V with Solver Y:** `--tda` gives its eigenvalues (§6.1); evaluating §5's
   $A(\omega)$ and $B(\omega)$ with $d_n = 1/N$ on its grid gives its two spectra.
2. **Give Solver V the coupled energies:** solve the coupled problem, orthonormalize the excitation
   parts $X_n$ of its $N$ positive eigenvectors ($\tilde X_n$), and hand the program

   $$
   H_{\rm eff} = \sum_{n=1}^{N} \lambda_n\, |\tilde X_n\rangle\langle \tilde X_n| ,
   $$

   Hermitian and $N\times N$: ZHEEV returns exactly the coupled $\lambda_n$, so `Eq0` gives A at
   −0.555 eV. The heights stay $1/N$ (§6.4): that is in the program's formula, not its input.

## 8. Summary of the WS2 numbers

| step | A | difference | cause |
|---|---|---|---|
| Solver V (TDA) | 2.5785 eV | | |
| Solver Y, `--tda` | 2.5785 eV | 0 | same matrix $A$ (§6.1) |
| Solver Y, coupled, kernel at 2.3795 eV | 2.3689 eV | −0.210 eV | anti-resonant coupling (§6.2) |
| exact kernel, $\lambda(\omega) = \omega$ | 2.3774 eV | +0.0085 eV | frozen vs self-consistent $z$ (§6.3) |
| BSE (export model) | 2.3774 eV | 0 | the kernel is exact (§2) |
