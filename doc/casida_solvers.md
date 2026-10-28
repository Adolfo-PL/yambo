<!-- License-Identifier: GPL. Copyright (C) 2026 Adolfo Partida. Authors: Adolfo Partida, with AI assistance (Claude Code). -->

# How the fork solves the Casida equations

This note covers the two Casida solvers of the fork. Each equation is tied to the file and the
lines that compute it, so the code can be checked against the literature one step at a time.

- **Static Casida** (`BSKmod= "CHI"`) projects f_xc(ω=0) and diagonalizes the matrix.
- **Dynamic Casida** (`BSKmod= "CHI"` with `BSEChiDyn`) uses f_xc(ω) and solves a linear system at each frequency.

The full reference is `doc/greens_function_fxc.md`, which has the derivations and the validation tables.
Line numbers are from commit 96ee648.

**Sources.** Each citation in §6 says whether it was checked from the PDF in this project ("checked")
or is a standard reference whose equation numbers were not checked here ("not checked"). A checked
equation number can be compared with the code directly. Look up a not-checked one before quoting it.

---

## 1. Notation

| Symbol | Meaning | Code |
| --- | --- | --- |
| $t=(v,c,\mathbf k)$ | a transition $v\,\mathbf k-\mathbf q \to c\,\mathbf k$, $N$ of them | `BS_T_grp(...)%table` |
| $E_t$ | transition energy $\varepsilon_{c\mathbf k}-\varepsilon_{v\mathbf k-\mathbf q}$ (QP or KS) | `BSS_eh_E`, `BS_T_grp%E` |
| $f_t$ | occupation factor $(f_v-f_c)/\text{spin\_occ}$, negative for anti-resonant rows | `BS_T_grp%f`, `BSS_eh_f` |
| $s_t=\sqrt{f_t}$ | $i\sqrt{\lvert f_t\rvert}$ when $f_t<0$ | `K_kernel.F:1000-1001`, `Chi_dyn_s` |
| $\rho_t(\mathbf G)$ | plane-wave transition density (vertex) of $t$ | `O_x` (exchange), `O_tddft_L`, `chi_R_res`/`chi_R_cpl`, `Chi_dyn_V` |
| $v_{\mathbf G}=4\pi/\lvert\mathbf q+\mathbf G\rvert^2$ | Coulomb; $\bar v$ drops $\mathbf G=0$ (`Lkind="bar"`) | `bare_qpg` |
| $D_{\mathbf G}=\sqrt{4\pi}/\lvert\mathbf q+\mathbf G\rvert$ | so that $v_{\mathbf G}=D_{\mathbf G}^2$ | `K_Chi_export.F:170` |
| $f(\mathbf G,\mathbf G';z)$ | the kernel f_xc read from `ndb.Chi` | `Chi_fxc(:,:,iw)` |
| $C_o=4\pi\,\text{spin\_occ}/(\Omega N_k)$ | Yambo's kernel prefactor | `K_kernel.F:134` |
| $c=C_o/4\pi=\text{spin\_occ}/(\Omega N_k)$ | prefactor of the f_xc projection | `K_kernel.F:460, 506` (`Chi_dyn_c`) |
| $d_t$ | optical dipole of $t$ | `BSS_dipoles_opt` |
| $z=\omega+i\eta$ | complex frequency | `W%p(iw)` |

$\Omega$ is the cell volume and $N_k$ the number of k points in the full Brillouin zone.

---

## 2. The matrix both methods start from

### 2.1 Casida form

In the coupled (non-Tamm–Dancoff) form, the linear-response TDDFT eigenproblem in transition space
is the Casida equation:

```math
\begin{pmatrix} A & B \\ -B^{*} & -A^{*} \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix} = \Omega_\lambda \begin{pmatrix} X \\ Y \end{pmatrix},
\qquad
A_{tt'} = E_t\,\delta_{tt'} + K^{\rm Hxc}_{tt'},\quad
B_{tt'} = K^{\rm Hxc}_{tt'}\ \text{(coupling pairs)} .
```

The BSE has the same shape, with $K^{\rm Hxc}$ replaced by $\bar v - W$. Yambo builds both with the
same code. The only difference in `BSKmod= "CHI"` is that the f_xc projection replaces $-W$.

### 2.2 One matrix element (`K_kernel.F`)

For each pair of transitions $(t,t')$ the kernel loop computes:

```math
K_{tt'} = s_t\,\Big[\,C_o\,\big(K^{v}_{tt'} + K^{f}_{tt'}/4\pi\big) - K^{W}_{tt'}\Big]\,s_{t'}
```

The terms come from these lines:

- **Prefactor and occupations.** `K_kernel.F:1003`: `BS_mat_tmp = sqrt_f_itk*(H_x*Co - H_c)*sqrt_f_itp`.
- **Exchange.** `K_exchange_kernel.F:10-65`:
  ```math
  K^{v}_{tt'}=\sum_{\mathbf G}\rho^{*}_{t'}(\mathbf G)\,\rho_{t}(\mathbf G)/\lvert\mathbf q+\mathbf G\rvert^2 .
  ```
  Since $C_o K^{v} = c\sum_{\mathbf G}\rho^{*}\,v_{\mathbf G}\,\rho$, the exchange is the bare Coulomb $\bar v$ (or $v$) between transition densities.
- **f_xc projection, static path.** `K_kernel.F:930-938` (resonant) and `:970-978` (coupling) call `TDDFT_Chi_pair_G` (`src/tddft/TDDFT_Chi_Fourier.F:20-43`). That routine returns `dot_product(d_left, matmul(fxc, d_right))`:
  ```math
  K^{f}_{tt'}=\sum_{\mathbf G\mathbf G'}\big[\rho^{L}_{t}(\mathbf G)\big]^{*}\,f(\mathbf G,\mathbf G')\,\rho^{R}_{t'}(\mathbf G') .
  ```
  The result is divided by $4\pi$ (`:938`, `:978`) because $C_o$ already contains it. The net f_xc term is therefore $c\,s_t\,\rho^{L\dagger}_t f \rho^R_{t'}\,s_{t'}$.
- **Screened interaction.** $K^{W}$ (`H_c`) is zero for `BSKmod= "CHI"`: no correlation pass runs.
- **Factor $i$ on the coupling block.** `K_kernel.F:1004`: `if (mode_now=="C") BS_mat_tmp = cI*BS_mat_tmp`.

### 2.3 Assembly into the full matrix (`K_stored_in_a_big_matrix.F`)

- **Diagonal energies.** `K_stored_in_a_big_matrix.F:71`: `Mij = real(Mij) + BSS_eh_E(...)`, which adds $E_t$.
- **Block placement.** `K_stored_in_a_big_matrix.F:79-109`. The resonant block is $R$, with $R^\dagger=R$ (`:79-81`). The anti-resonant block is $-R^{*}$ (`:84-86`). The coupling block is stored with its factor $i$, and its partner is $-(\cdot)^{*}$ (`:89-95`).

The stored matrix is therefore

```math
H_{\rm Yambo}=\begin{pmatrix} R & iC \\ iC^{*} & -R^{*} \end{pmatrix}
= T\begin{pmatrix} R & C \\ -C^{*} & -R^{*} \end{pmatrix}T^{-1},\qquad T={\rm diag}(1,-i),
```

where the second equality assumes $C^{T}=C$. This is a similarity transform of the Casida form in
§2.1, so the eigenvalues are the same. The anti-resonant rows have $f_t<0$, so their $s_t$ is
imaginary (`K_kernel.F:1000-1001`).

---

## 3. Method 1: static Casida (`BSKmod= "CHI"`)

### 3.1 Kernel

The first frequency sample of `ndb.Chi` is used:

```math
f_{\rm static}(\mathbf G,\mathbf G')=f(\mathbf G,\mathbf G';z=0)
```

- `K_kernel.F:486-487`: `Chi_static = Chi_fxc(:,:,1)`.
- `K_kernel.F:485` stops the run if this matrix is not Hermitian, because a static Casida matrix must be.
- For a `BSEChiOut` export, sample 1 is exactly $z=0$ (doc, `BSEChiOut` section).

### 3.2 Eigenproblem

The matrix of §2 is assembled once and diagonalized:

- `K_solvers.F` → `K_diago_driver.F:159`: `K_is_not_hermitian = BS_K_coupling .or. ...`.
- `K_diago_driver.F:227`: `LINEAR_ALGEBRA_driver(DIAGO, M=BS_mat, E_cmpl=BS_E, V_left=..., V_right=...)`.

This is a general complex eigensolver. It returns the
eigenvalues $E_\lambda$ and the right and left eigenvectors $V^{R}_\lambda$, $V^{L}_\lambda$:

```math
H\,V^{R}_\lambda = E_\lambda V^{R}_\lambda .
```

### 3.3 Spectrum

**Residuals.** `K_diago_right_residuals.F:11, 61, 68` compute

```math
R^{R}_\lambda=\sum_t d_t^{*}\,s_t\,V^{R}_{t\lambda},
```

and `K_diago_left_residuals.F` computes $R^{L}_\lambda$ with the overlap of the left and right
eigenvectors.

**Response.** `K_diago_response_functions.F:11, 129` give

```math
\chi(\omega) \propto -\sum_\lambda \frac{R^{L}_\lambda R^{R}_\lambda}{z-E_\lambda},
```

with $\text{Resp\_ii}=-C_o\times(\ldots)$. The output file is `o-<job>.eps_q1_diago_bse`. The
excitation energies are the eigenvalues $E_\lambda$, and the `[BSE]` lines of the `r-` file report
the lowest bright one and the binding energy.

**Cost.** One diagonalization of a $2N\times 2N$ matrix.

**Physics.** The static kernel reproduces $\varepsilon(\omega=0)$ of the reference, but it can barely
move the poles away from the transition energies. On WS2 it binds only a small part of the BSE
exciton (doc, "Frequency-dependent Casida" and `BSEChiOut`). This is the known need for a
frequency-dependent kernel (§6, rows 3 and 4).

---

## 4. Method 2: frequency-dependent Casida (`BSEChiDyn`)

With $f(z)$ the matrix depends on the frequency, so the eigenproblem becomes nonlinear in
$\Omega_\lambda$. The code does not search for those roots. Instead it computes the response at
every frequency of the grid, which needs one linear solve per frequency.

### 4.1 What is prepared first

| Step | Equation | Code |
| --- | --- | --- |
| $H_0$: the matrix of §2 without f_xc ($E_t$ plus exchange) | $H_0=H\big\rvert_{f=0}$ | `K_kernel.F` skips the f_xc term when `BS_K_Chi_dynamic` is on (`:930`, `:970`); `K_Chi_dynamic.F:128-129, 137` copy `BS_mat` |
| vertices of every transition | $\rho_A(\mathbf G)$, $A=1..2N$ | `K_kernel.F:724-737` (`Chi_dyn_V`, `Chi_dyn_s`); `K_Chi_dynamic.F:153-159` (`Vf`, `sf`, `df`) |
| routing to this solver | n/a | `K_solvers.F:59-61` |

### 4.2 Kernel block at one frequency

`TDDFT_Chi_full_dynamic_response`, `src/tddft/TDDFT_Chi_dynamic.F:295-300` (two `zgemm` calls,
then the $s$ scaling):

```math
K^{f}_{AB}(z)=c\,s_A\,s_B\sum_{\mathbf G\mathbf G'}\rho^{*}_A(\mathbf G)\,f(\mathbf G,\mathbf G';z)\,\rho_B(\mathbf G') .
```

All $2N$ rows (resonant and anti-resonant) see the same retarded $f(z)$ (header comment,
`TDDFT_Chi_dynamic.F:34-36`). For a Hermitian static $f$ this is exactly the term §3 adds to the
matrix (`TDDFT_Chi_dynamic.F:16-20`).

### 4.3 The solve

Same routine, `TDDFT_Chi_dynamic.F:275-312`:

```math
\mathrm{Resp}(\omega) = -C_o \sum_{AB} s_A\,d^{*}_A\,\Big[\,z - H_0 - K^{f}(z)\Big]^{-1}_{AB}\, s_B\,d_B .
```

- `:295-304` form the matrix $z-H_0-K^{f}(z)$ (`:301` subtracts $H_0$, `:302-304` add $z$).
- `:305-306` solve $[z-H_0-K^f]\,x = s\,d$ with LAPACK `zgesv` (an LU factorization; no inverse is formed).
- `:310` takes $\sum_A s_A d_A^* x_A$.
- The $-C_o$ is applied at `K_Chi_dynamic.F:228`.

The frequency loop is `K_Chi_dynamic.F:201-217`: the frequencies are shared across OpenMP threads,
with one solve each.

The Tamm–Dancoff variant (`TDDFT_Chi_dynamic_response`, `:79-155`) solves the resonant and
anti-resonant blocks separately. The code warns against it, because the exported kernel inverts
the full bubble.

**Output.** `o-<job>.eps_q1_chidyn_bse`. The excitons are the **peaks** of this spectrum, not
eigenvalues.

**Cost.** One $2N\times 2N$ LU factorization per frequency. The grid must be the `ndb.Chi` grid
exactly, because $f$ is not interpolated.

### 4.4 Static limit

If $f$ does not depend on $z$ and is Hermitian, use the spectral resolution of the resolvent:

```math
[z-H]^{-1}=\sum_\lambda \frac{V^R_\lambda (V^L_\lambda)^\dagger}{z-E_\lambda}
```

The formula of §4.3 then reduces to the sum over poles of §3.3. `tests/chi_dynamic` check 2
verifies this numerically, for both TDA and coupling (`tests/chi_dynamic/README.md`).

---

## 5. Why the dynamic solve reproduces the BSE (Woodbury identity)

### 5.1 The export (`K_Chi_export.F:252-265`)

- **BSE response in G space.** `TDDFT_Chi_full_response`, `TDDFT_Chi_dynamic.F:227-257`:
  ```math
  \bar\chi(\mathbf G,\mathbf G';z)=c\sum_{AB} s_A\,\rho_A(\mathbf G)\,\big[z-M\big]^{-1}_{AB}\,s_B\,\rho^{*}_B(\mathbf G') .
  ```
  Here $M$ is the BSE matrix: energies, $\bar v$ and $-W$.
- **Reference bubble.** `TDDFT_Chi_full_bubble`, `TDDFT_Chi_dynamic.F:259-273`:
  ```math
  \chi_0(\mathbf G,\mathbf G';z)=c\sum_A s_A^2\,\frac{\rho_A(\mathbf G)\rho_A^{*}(\mathbf G')}{z-E_A}
  ```
  For `BSEChiOut= "exc"` the energies are those of the BSE (QP). For `"full"` they are KS.
- **Kernel.** `TDDFT_Chi_kernel_from_response`, `TDDFT_Chi_dynamic.F:532-600`:
  ```math
  P^{-1}=\bar\chi^{-1}+\bar v,\qquad f=\chi_0^{-1}-P^{-1}
  ```
  The inversions are done on the symmetrized matrices $D\chi D$ (`:561-570`). $\bar v$ is added as $+1$ on the diagonal for $\mathbf G\neq 0$, and also for $\mathbf G=0$ with `Lkind="full"` (`:571-578`). The result is $f=D(\tilde\chi_0^{-1}-\tilde P^{-1})D$ (`:579-583`), and $P$ is inverted back (`:584-591`).

### 5.2 The identity

In G space, Casida with exchange $\bar v$ and kernel $f$ gives:

```math
\chi_{\rm Casida}=\big[\chi_0^{-1}-\bar v-f\big]^{-1}
=\big[\chi_0^{-1}-\bar v-\chi_0^{-1}+\bar\chi^{-1}+\bar v\big]^{-1}=\bar\chi .
```

The pair-space solve of §4.3 is the same quantity. By the Woodbury identity, the transition-space
resolvent $[z-E-c\,S\rho^\dagger(\bar v+f)\rho S]^{-1}$, projected onto the vertices, equals the
G-space Dyson equation above. This holds whenever the transitions rebuild the same $\chi_0$:

- same bands, k grid and energies;
- same G basis and frequencies.

`K_Chi_dynamic.F:167-189` checks the bubble condition and stops above $10^{-3}$ for BSE-exported
kernels. `tests/chi_dynamic` check 1 verifies the identity for an arbitrary complex, non-Hermitian
$f(z)$.

**What this does and does not show.** The Casida spectrum then reproduces the BSE spectrum
**by construction**. That checks the export, the file and the solver; it does not validate the
physics of the BSE. An error in $W$, the QP energies or the occupations is carried over just as
exactly (doc, `BSEChiOut` section).

---

## 6. Checking the code against the literature

| # | Statement in the code | Literature | Source status |
| --- | --- | --- | --- |
| 1 | Casida eigenproblem with $A$, $B$ blocks (§2.1) | M. E. Casida, in *Recent Advances in Density Functional Methods*, Part I (World Scientific, 1995) | not checked (standard) |
| 2 | TDDFT and BSE share the transition-space structure; $f$ replaces $\bar v-W$ in the kernel (§2.2) | G. Onida, L. Reining, A. Rubio, Rev. Mod. Phys. 74, 601 (2002) | not checked (standard) |
| 3 | $f=\chi_0^{-1}-P^{-1}$ from a BSE $P$; QP bubble ("exc") or KS bubble ("full") (§5.1) | R. Del Sole, G. Adragna, V. Olevano, L. Reining, PRB 67, 045207 (2003), Eq. 4 ($f^{\rm FQP}$, $f^{\rm FLDA}$); they evaluate the head only | checked |
| 4 | the kernel must depend on frequency to bind excitons (§3.3) | Del Sole et al. (row 3); F. Sottile, V. Olevano, L. Reining, PRL 91, 056402 (2003); A. Marini, R. Del Sole, A. Rubio, PRL 91, 256402 (2003) | checked |
| 5 | static kernel in a Tamm–Dancoff Casida equation in transition space | Shi, Turkowski, Rahman, PRB 107, 155431 (2023), Eq. 14 (Casida), Eq. 21 (kernel) | checked |
| 6 | the exact kernel has poles between excitons (relevant to the per-frequency solve, §4) | Gulevich et al., PRB 107, 165118 (2023); found on WS2 at 1.78 and 2.23 eV | checked |
| 7 | Woodbury (matrix-inversion) identity (§5.2) | W. W. Hager, SIAM Review 31, 221 (1989) | not checked (standard) |

### Suggested checks, in order

1. **§2.1 against Casida (row 1) and Onida et al. (row 2).** Match the $A$ and $B$ blocks of §2.1 to
   their equations. Then confirm that Yambo's stored form (§2.3) is the similarity transform shown
   there. The two papers may use different sign conventions or a different normalization of $X$
   and $Y$; the eigenvalues must agree.
2. **§2.2 normalization.** Check that $c=\text{spin\_occ}/(\Omega N_k)$ matches the factor in the
   $\chi_0$ of these papers. This is the factor that makes the $\chi_0$ of §5.1 equal to the RPA
   bubble.
3. **§5.1 against Del Sole et al. Eq. 4 (row 3).** Their formula is ours restricted to
   $\mathbf G=\mathbf G'=0$. `doc/greens_function_fxc.md`, "Literature checks", records that comparison.
4. **§3 against Shi et al. Eq. 14 (row 5).** This compares the static TDA Casida structure. Their
   kernel (Eq. 21) is a different, static approximation (`"nq"` mode, doc), so compare the
   equation's structure, not the numbers.
5. **§4 has no single literature equation.** It is the pair-space Dyson equation solved per
   frequency, and its correctness rests on §4.4 and §5.2. The tests that check it are listed in
   §7. A literature comparison is possible for the result: frequency-dependent-kernel spectra
   (rows 4 and 6).

---

## 7. Tests that check these equations

| Test | Checks | Where |
| --- | --- | --- |
| pair-space identity, coupling and TDA, arbitrary complex $f(z)$ | §4.3 = G-space Dyson | `tests/chi_dynamic`, check 1 |
| static limit equals the diagonalization sums | §4.4 | `tests/chi_dynamic`, check 2 |
| `"mdr"` and `"nq"` kernels against term-by-term references | §5.1 variants | `tests/chi_dynamic`, checks 5 and 6 |
| bulk Si: `BSEChiDyn` against the QP RPA+LF response of `ndb.Chi`, $3\times10^{-4}$ | §4 + §5 on a material | doc, "Frequency-dependent Casida" table |
| WS2 SOC: exported kernel, Casida against its BSE, Im ε $7\times10^{-5}$, Re ε $1\times10^{-4}$ | §5.2 on a 2D material | doc, `BSEChiOut` test table |

Run them with:

```
python tests/chi_dynamic/run_tests.py --fc gfortran [--openmp]
```
