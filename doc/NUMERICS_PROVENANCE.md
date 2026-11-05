<!-- License-Identifier: GPL. Copyright (C) 2026 Adolfo Partida. Authors: Adolfo Partida, with AI assistance (Claude Code). -->

# Provenance of the numerics of the exact kernel

This note records, for every numerical piece the fork added to Yambo for the exact kernel and its
Casida solvers, **where it comes from**: Yambo, the literature, or neither. For anything not taken from
Yambo it gives the justification. It is meant for checking the work item by item, and as the credit
section of any write-up.

## 0. How to read this note

### 0.1 Scope

**Covered, in depth:** the exact kernel exported from a BSE (`BSEChiOut = "exc"` or `"full"`, any q), and
the solvers that use it: static Casida (`BSKmod = "CHI"`), frequency-dependent Casida (`BSEChiDyn`), the
text dump of the Casida matrix (`BSEChiKoutW`).

**Not covered in depth** (classified file by file in section 7, nothing more): the `"mdr"` and `"nq"` kernels, the
Green-function route (`ChiFxc`, `ChiGMode`, `GRetarded`), the screened-potential database `ndb.W`, and the
scripts, documentation and tests.

### 0.2 Pinned versions and reference notation

| | |
| --- | --- |
| upstream | Yambo 5.3.0, revision 23959, commit `19c1241` |
| fork | branch `claude/nbd-file-generation-path-yidcr2`; no source file cited here changed after commit `260ecff` |
| `path:line` | a line of the fork |
| `U:path:line` | a line of upstream Yambo 5.3.0 |

Every line reference was checked against the source when this note was written (section 9).

### 0.3 The tiers

| tier | meaning |
| --- | --- |
| **A** | Yambo code called unchanged |
| **B** | adapted from a Yambo routine: the same construction, generalized or reorganized |
| **C** | not in Yambo; taken from the literature (equation numbers checked against the papers where marked) |
| **D** | no source in Yambo or the literature: numerical design, consistency checks, bookkeeping |

"Not in Yambo" means that the upstream routines named were read for it. It is not a search of all of Yambo.

### 0.4 Status words

Each claim carries one of these, where it matters:

- **verified in source**: read in the upstream or fork file at the line cited;
- **stated in a code comment**: the fork's author wrote the reason; not independently verified here;
- **tested**: a run or test in the project supports it (section 9);
- **not verified**: said plainly, and listed again in section 8.

### 0.5 No claim of novelty

This note says what is taken from where. It does not claim that anything is new. The literature checked is the
set of papers the project holds (Del Sole et al.; Adragna, Del Sole, Marini; Sottile, Olevano, Reining;
Marini, Del Sole, Rubio; Shi, Turkowski, Rahman; Gulevich et al.; and the doping papers listed in
`doc/greens_function_fxc.md`). A paper outside that set may contain any of Tier C or D.

---

## 1. The path in one picture

```
 Yambo (A)                                fork                                           Yambo (A)
 ---------                                ----                                           ---------
 em1s, W, exchange --> BS_mat  ---------> K_Chi_export.F ----> ndb.Chi (chi0, P, fxc)
 K_Transition_check    (H0, 2N x 2N)        |  vertices rho_A(G)   [B1-B3]                    |
 scatter_Bamp, dipoles                      |  bubble chi0         [B5]                       |
                                            |  response chib       [B6]                       v
                                            |  P^-1 = chib^-1+vbar [B7, C2]        K_kernel.F  (BSKmod="CHI")
                                            |  fxc = chi0^-1-P^-1  [C1]               |  reads ndb.Chi
                                            |  LU inversion        [B8]               |  static:  H = H0 + K^f(0)  [B4]
                                            v                                         |  -> Yambo diagonalization (A)
                                       checks, guards  [D]                            |
                                                                                      |  dynamic: K_Chi_dynamic.F
                                                                                      |  one solve per z       [B9]
                                                                                      v
                                                                          o-*.eps_q1_diago_bse / o-*.eps_q1_chidyn_bse
```

The mathematics of each arrow is derived in the LaTeX notes (`bse_to_casida`) and, in short form, in
`doc/casida_solvers.md`.

---

## 2. Tier A: Yambo code called unchanged

| id | what it provides | upstream routine | used in the fork at |
| --- | --- | --- | --- |
| A1 | transitions, `E_t`, `f_t=(f_v-f_c)/spin_occ`, the window, dropped transitions | `U:src/bse/K_Transition_check.F:64` (`f_eh_REF`) | the BSE run itself |
| A2 | plane-wave matrix elements `<in\|e^{i(q+G)r}\|out>` (FFT, spinors, USPP augmentation) | `U:src/wf_and_fft/scatter_Bamp.F:10` | `src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:156`, `_L.F:137` |
| A3 | position dipoles `<c\|r\|v>` and their rotation to the full zone | `DIPOLE_rotate`, `DIP_iR` (`U:src/bse/K_dipoles.F:60, 158`) | `src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:160` |
| A4 | exchange `Vbar`, with the `G=0` term set by `Lkind` | `U:src/bse/K_exchange_kernel.F`, `U:src/bse/K_exchange_collisions.F:102-103` | the BSE run |
| A5 | screened interaction `W` from `ndb.em1s` | `U:src/bse/K_screened_interaction.F`, `U:src/bse/K_correlation_kernel_std.F` | the BSE run |
| A6 | the matrix element `sqrt(f) (H_x Co - H_c) sqrt(f)`, `Co`, the factor `i` on coupling blocks | `U:src/bse/K_kernel.F:114` (`Co`), `U:src/bse/K_kernel.F:755` (`BS_mat_tmp = sqrt_f (H_x Co - H_c) sqrt_f`) | `src/bse/K_kernel.F:134, 1003-1004` |
| A7 | assembly of the `2N x 2N` matrix: energies on the diagonal, blocks `R`, `iB`, `-R*`, `-(iB)*` | `U:src/bse/K_components_folded_in_serial_arrays`, `U:src/bse/K_stored_in_a_big_matrix.F:71, 79-109` | `src/bse/K_Chi_export.F:130, 133`; `src/bse/K_Chi_dynamic.F:128-129` |
| A8 | diagonalization, residues, overlap, sum over poles, `eps_M` | `U:src/bse/K_diago_driver.F`, `U:src/bse/K_diago_left_residuals.F`, `K_diago_right_residuals.F`, `K_diago_response_functions.F`, `U:src/output/K_OUTPUT.F` | static Casida (unchanged); `K_OUTPUT` in `src/bse/K_Chi_dynamic.F:236` |
| A9 | the output-file machinery of a solver (`X_OUTPUT_driver`, `Resp_ii`, `para_term_w0`) | `U:src/pol_function/X_OUTPUT_driver.F` | `src/bse/K_Chi_dynamic.F:68, 228, 232` |
| A10 | the `ndb` database mechanism (`io_control`, `io_elemental`, fragments) | `U:src/io/`, `U:src/Yio/` | `src/io/io_Chi.F` follows the pattern |
| A11 | the BSE kernel loop over transition blocks, symmetry rotation, OpenMP, `BS_blk`, `BS_T_grp` | `U:src/bse/K_kernel.F` | `src/bse/K_kernel.F` (edited, B12) |
| A12 | LAPACK / BLAS (`zgesv`, `zgetrf`, `zgetri`, `zgemm`) | the libraries | `src/tddft/TDDFT_Chi_dynamic.F` |


---

## 3. Tier B: adapted from Yambo

### B1. Transition-density vertices `rho_A(G)` for `G != 0`

- **Fork:** `src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:148-157` (right vertex), `src/tddft/TDDFT_ALDA_eh_space_G_collisions_L.F:130-138` (left vertex).
- **Upstream origin:** the same two routines. `U:src/tddft/TDDFT_ALDA_eh_space_G_collisions_L.F:122-123` forms the real-space product
  `psi_1(r) conj(psi_2(r))` and takes its FFT; `U:src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:133-134` does the same
  after multiplying by the local kernel `F_xc(r)`. Both are verified in source.
- **Change:** in "bare vertex" mode both vertices are the transition density from `scatter_Bamp`, conjugated:
  `O_tddft_R(:,t,1,1) = conjg(chi_scatt%rhotw)` (line 157). No `F_xc` is multiplied in, and the wavefunction-symmetry step
  `WF_apply_symm` is skipped (`_R.F:143-146`, `_L.F:125-128`).
- **Why (justification):**
  1. a nonlocal kernel `f(G,G')` cannot be multiplied in real space beforehand, so the vertex must be bare;
  2. the vertex must be exactly the one the exchange and the bubble use, because `vbar` and `f` act on the same `rho`.
     The code comment gives the reason for `scatter_Bamp` rather than the FFT product: it matches the native `X_irredux`
     density "at every G (summed over the spinor components, USPP augmentation included)" (`src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:150-152`).
     **Status: stated in a code comment.** Whether the upstream FFT product would differ for the pseudopotentials used was not tested.
- **Evidence:** tested. The head fit (D4) compares the vertices with Yambo's own response; the exported response
  reproduces the BSE to `1e-4` (WS2) and the Casida runs reproduce it (section 9).

### B2. The optical head of the vertex at `q = 0` (the `G = 0` element)

- **Fork:** `src/tddft/TDDFT_ALDA_eh_space_G_collisions_R.F:158-164` and `src/tddft/TDDFT_ALDA_eh_space_G_collisions_L.F:139-145` (the same block).
- **Upstream origin:** `U:src/bse/K_dipoles.F:157-171`. At `q = 0` upstream builds the head vertex of the BSE for optics as
  `DIP_expanded = -i DIP_bare` (line 159), `DIP_projected = field_direction . DIP_expanded` (line 161), times `q0_def_norm` (line 164), then `dipoles_opt = conjg(DIP_projected)` (line 171),
  with `DIP_bare` from `DIPOLE_rotate(..., DIP_iR, ...)` (line 158). At finite `q` the same routine takes the vertex from `scatter_Bamp`: `DIP_projected(1) = -conjg(rhotw(i_G_shift))` (lines 166-168),
  and a comment says the `q = 0` branch is multiplied by `q0_def_norm` "to have something consistent with the definition from scatter_Bamp at finite-q" (line 163). Verified in source.
  So Yambo itself uses both constructions in one place; the fork uses the `scatter_Bamp` one for every `G` and the dipole one for the `G = 0` element at `q = 0`.
- **Change:** the fork recomputes the head in the vertex routine from the same `DIPOLE_rotate(..., 'DIP_iR', ...)`, as
  `chi_head = -dot_product(BSS_q0, chi_dipole)` (`_R.F:161`), conjugated and sign-flipped for the anti-resonant and coupling modes
  (`'A'` or `'C'`, `_R.F:162`; `'A'` or `'Q'`, `_L.F:143`). The scalar conventions differ from upstream's (`-i`, `q0_def_norm`, an
  explicit `conjg`); they are absorbed by `BSS_q0`.
- **Why:** at `q = 0` the density `<v|e^{iqr}|c>` vanishes for `v != c`, so the `G = 0` element has to be the dipole limit. The vertex must carry exactly
  the phase and normalization that Yambo's own spectrum uses.
- **Evidence:** tested. This is what the head-consistency fit (D4) checks: the `G = 0` element of `chib` must equal Yambo's response built from
  `dipoles_opt`, up to the known scale `c`, with residual below `1e-4`. Also verified: silicon and WS2 runs pass it.
- **Caveat:** a convention error here would not be seen by the identity of section 6 (C2), which holds for any vertices. D4 is what catches it.

### B3. Vertex bookkeeping for the export and the dynamic solver

- **Fork:** `src/bse/K_kernel.F:709-737` (collection), `K_kernel.F:455-460, 502-506` (allocation), `src/modules/mod_Chi.F` (`Chi_dyn_V`, `Chi_dyn_s`, `Chi_dyn_filled`, `Chi_dyn_c`).
- **Upstream origin:** the per-group cache of the ALDA right vertices, `BS_T_grp(i)%O_tddft_R`, filled in the loop over transition groups
  (`U:src/bse/K_kernel.F`, the "TDDFT oscillators" block).
- **Change:** the resonant (`"R"`) and de-excitation (`"C"`) right vertices of each group are copied, per transition, into `Chi_dyn_V(:, t, 1:2)`, together with
  `s_t = sqrt(f_t)` and a `Chi_dyn_filled` mask. The cache is invalidated per column group with `chi_grp_p` (lines 709, 742).
  At finite `q` the explicit anti-resonant transitions are groups of their own: their `"R"` vertex goes to `Chi_dyn_V(:,t,i_ra)` with
  `i_ra = BS_T_grp%i_res_ares` (lines 733-737).
- **Why:** the ALDA cache is overwritten group by group and holds one vertex type; the export and the dynamic solver need every transition's vertices
  at once, for both vertex types, after the kernel loop. The mask lets the export stop if a transition was missed (guard at `K_Chi_export.F:128`).
- **Evidence:** tested at `q = 0` (derived anti-resonant rows) and finite `q` (explicit rows) in `tests/chi_dynamic` and the hBN and WS2 runs
  (`doc/greens_function_fxc.md`, "Finite momentum").

### B4. The static projection `K^f_{tt'} = rho_L^H(t) f_xc(0) rho_R(t') / 4 pi`

- **Fork:** `src/tddft/TDDFT_Chi_Fourier.F:20-43` (`TDDFT_Chi_pair_G`), called at `src/bse/K_kernel.F:931` (resonant) and `971` (coupling),
  added to `H_x` at `938` and `978`.
- **Upstream origin:** `U:src/tddft/TDDFT_ALDA_eh_space_G_kernel.F:35-52`, which computes `sum_G conj(L(G)) R(G)` with `Vstar_dot_V_omp` (lines 35, 41) and divides by `4 pi` (line 52);
  `R` contains `F_xc` already (B1). The result is added to `H_x`, which `K_kernel` multiplies by `Co`.
- **Change:** `R` is bare and the kernel matrix sits between the vertices: `value = dot_product(d_left, matmul(fxc, d_right))`. The `1/4 pi` and the use of `Co` are unchanged.
- **Why:** the ALDA routine can only apply a kernel that is diagonal in real space. `f(G,G')` is a full matrix.
- **Evidence:** tested (`tests/chi_fxc_fourier`; the static Casida reproduces `eps(0)` of the BSE).

### B5. The reference bubble `chi0`

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:259-273` (`TDDFT_Chi_full_bubble`), called at `src/bse/K_Chi_export.F:255`.
- **Upstream origin:** the free electron-hole propagator `Lo` of the inversion solver (`U:src/bse/K_inversion_Lo.F:150-151`, via `X_GreenF_analytical`, retarded `"r"` for resonant and `"Ra"` for anti-resonant rows),
  whose residue is `dipole * conjg(dipole)` (`U:src/bse/K_inversion_Lo.F:177`); and the Adler-Wiser bubble of `X_irredux`.
- **Change:** `chi0(G,G') = c sum_A s_A^2 rho_A(G) conj(rho_A(G')) / (z - E_A)` on the BSE's own rows (resonant and anti-resonant), with the BSE's vertices and the reference energies of B11.
  Upstream's `Lo` is shifted by the diagonal of the kernel (`U:src/bse/K_inversion_Lo.F:142-143`); the fork's bubble is not.
- **Why not call `X_irredux`:** the identity needs exactly this bubble (section 6). `X_irredux` differs in the transitions (its own band window and energy filters),
  the vertices (its own optical `q0`, `LongDrXs`), the frequency ordering (time-ordered by default, `U:src/modules/mod_X.F:276`; retarded only for finite temperature or a causal Green function, `U:src/pol_function/X_pre_setup.F:90-93`),
  the precision (`complex(SP)`, `U:src/modules/mod_matrix.F:61`) and the parallel layout. **Verified in source.**
- **Evidence:** tested (D6 compares it with the stored `chi0`; the Casida runs depend on it).

### B6. The `G`-resolved BSE response `chib = U (z - H)^-1 U^H`

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:227-257` (`TDDFT_Chi_full_response`), called at `src/bse/K_Chi_export.F:252`.
- **Upstream origin:** the BSE inversion solver. `U:src/bse/K_inversion_do_it_full.F:31-38` builds `Lo^-1 - M` and inverts it;
  `U:src/bse/K_inversion_engine.F:139-141` contracts the inverse with `BSS_dipoles_opt(1,:)` on both sides and multiplies by `-Co_factor`.
- **Change:** the contraction is with the `N_G` vertex vectors instead of the dipole vector, giving the whole `G G'` matrix. The linear algebra is an LU solve with
  `N_G` right-hand sides (`zgesv`, line 251) followed by one `zgemm` (line 255), not an explicit inverse.
- **Why:** the exact kernel is a `G G'` matrix and needs the `G G'` response. Upstream's solvers return only the macroscopic response:
  the diagonalization (`U:src/bse/K_diago_left_residuals.F:69`), the inversion (`U:src/bse/K_inversion_engine.F:139-141`) and Haydock
  (`U:src/bse/K_Haydock.F:230`) all start from `dipoles_opt`. **Verified in source for the three solver families; SLEPc uses the diagonalization residuals.**
- **Evidence:** tested (the head of `chib` is compared with Yambo's response, D4).

### B7. Removing the exchange: `P^-1 = chib^-1 + vbar`, in symmetrized form

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:532-591` (`TDDFT_Chi_kernel_from_response`).
- **Upstream origin:** `U:src/pol_function/X_redux.F`. The Dyson step builds `KERNEL = 1 - chi0 v` (`U:src/pol_function/X_redux.F:418-419`: `KERNEL - Xo*4pi/bare_qpg(iq,ig2)**2`, `+ cONE` on the diagonal),
  solves `KERNEL X = Xo` (`U:...X_redux.F:253`), and symmetrizes with `4 pi / (bare_qpg(ig1) bare_qpg(ig2))` (`U:...X_redux.F:280, 295`), i.e. `D X D` with `D = sqrt(4 pi)/|q+G|`. Verified in source.
- **Change:** the same `D` (`src/bse/K_Chi_export.F:170`, complex with a Coulomb cutoff, like the exchange) and the same `G = 0` convention (`Lbar`:
  `O_x(1) = 0`, `U:src/bse/K_exchange_collisions.F:103`), but the Dyson step is **reversed**. Yambo adds `v` to go from `chi0` to `chi`; here `v` is subtracted to go from the BSE response to `P`:
  `P~^-1 = chib~^-1 + J`, `J = diag(0,1,...,1)`, with `J(1,1) = 1` for `Lkind = "full"` (lines 571-578).
- **Why symmetrized:** in `D chi D` every block is `O(1)`, while the raw head of `chib` is `O(q^2)`. Inverting the raw matrix would mix magnitudes about `1e6` apart.
- **Evidence:** tested (`tests/chi_dynamic` checks the pair-space identity; the WS2 Casida reproduces the BSE).

### B8. The LU inversion

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:599-620` (`invert`, internal to `TDDFT_Chi_kernel_from_response`): `zgetrf`, `zgecon`, `zgetri`.
- **Upstream origin:** `U:src/modules/mod_linear_algebra.F:144-172` (`LU_factorization`, `LU_inversion`): the same pair, `getrf` then `getri`.
- **Change:** double precision called directly, plus `zgecon` for the reciprocal condition number (D2).
  Yambo's wrapper works on `complex(SP)` and selects `zgetrf` only if built with `_DOUBLE` (`U:src/modules/mod_linear_algebra.F:152-156, 167-171`).
- **Why not the wrapper:** the kernel is the difference of two nearly cancelling inverses, so single precision is not acceptable (D9); and the export calls these routines from OpenMP threads
  over frequencies (`!$omp` at `src/bse/K_Chi_export.F:247`), so they are `recursive` and use only local arrays. **Whether Yambo's wrapper is thread-safe was not tested.**

### B9. The per-frequency Casida solve

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:275-312` (`TDDFT_Chi_full_dynamic_response`); driver `src/bse/K_Chi_dynamic.F:201-228`.
- **Upstream origin:** the inversion solver again: one linear-algebra step per frequency (`U:src/bse/K_inversion_engine.F:78-80, 139-141`); the residual vectors of `K_diago_*_residuals`
  (`a = d sqrt(f)`, `b = conj(d) sqrt(f)`); the scaling `Resp_ii = -Co_factor * resp` and `eps_M` from `K_OUTPUT`.
- **Change:** a frequency-dependent block `Kf(z) = c s_A s_B rho^H f_xc(z) rho` (two `zgemm` calls and the `s` scaling, lines 295-300) is added to `H0` at every `z`,
  and the system `[z - H0 - Kf(z)] x = s d` is solved for one right-hand side (`zgesv`, line 306), giving `resp = sum_A s_A conj(d_A) x_A` (line 310).
  No explicit inverse is formed.
- **Why:** with `f(omega)` there is no fixed matrix to diagonalize. Yambo's inversion solver has the right shape; it never had a frequency-dependent kernel.
- **Evidence:** tested (`tests/chi_dynamic` check 1 for an arbitrary complex `f(z)`; check 2 against the Lorentzian sums; WS2 dynamic Casida reproduces the BSE to `1e-4`).

### B10. Reading the BSE matrix, and the shape of the two drivers

- **Fork:** `src/bse/K_Chi_export.F:130-133, 212` and `src/bse/K_Chi_dynamic.F:128-129, 137`: `K_components_folded_in_serial_arrays(iq)`, `K_stored_in_a_big_matrix(1,iq,"DIAGO")`, `H0 = cmplx(BS_mat, kind=DP)`.
- **Upstream origin:** the same two calls in `U:src/bse/K_diago_driver.F:162, 171` and `U:src/bse/K_inversion_driver.F:121, 144`.
- **Change:** the matrix is read in the way the diagonalization solver does, and copied to double precision.
- **Why serial:** `BS_mat` is assembled serially by this path; the drivers refuse an MPI run (guards `K_Chi_export.F:121`, `K_Chi_dynamic.F:72`).

### B11. The reference energies `E_ref`

- **Fork:** `src/bse/K_Chi_export.F:151-165`.
- **Upstream origin:** the Yambo energy structure `levels`: `Ken%E` (the energies the BSE uses) and `Ken%Eo` (the uncorrected ones, kept when QP corrections are applied).
- **Change:** `"exc"` takes the BSE transition energies `BSS_eh_E` (line 154). `"full"` with QP corrections rebuilds `E_c(k) - E_v(k-q)` from `Ken%Eo` (lines 164-165);
  anti-resonant rows derived from the resonant ones are `-E` (lines 160-163).
- **Why:** the kernel is defined relative to the reference bubble (C1, C3); the Casida run must start from the same energies (D7).

### B12. Control-flow edits in `K_kernel.F`

- **Fork:** `src/bse/K_kernel.F:181-182` (`l_chi_vtx`, `l_chi_vtx_only`), `203-207`, `342`, `357-508` (the Chi block), `647`, `695-737`, `766`, `926-938`, `966-978`, `1104-1111`, `1155`.
- **Upstream origin:** `U:src/bse/K_kernel.F`: the `BS_K_is_ALDA` branches.
- **Change:** the ALDA machinery is reused as a vehicle: the TDDFT pass runs when `BS_K_is_ALDA` **or** a Chi vertex is needed (`l_chi_vtx`); the real-space `F_xc` setup is wrapped in
  `if (BS_K_is_ALDA)` (this re-indents 24 existing lines, so the diff looks larger than the change); the kernel is not written to disk for a CHI run (`l_write_kernel`, line 647);
  the Chi block of the loop adds the static projection (B4). **Verified in source.**
- **Why:** this is the least intrusive way to get the transition-pair loop, the symmetry rotation and the OpenMP structure of the BSE kernel, which are Yambo's (A11), without duplicating them.

### B13. Keeping the position dipoles

- **Fork:** `src/bse/K_dipoles.F:260` (`if (.not.(BS_K_is_Chi.or.trim(BS_K_Chi_out)/='none')) call DIPOLE_ALLOC_global()`), `src/bse/K_kernel.F:277, 1155` (`call DIPOLE_ALLOC_global()` after the kernel).
- **Upstream origin:** `U:src/bse/K_dipoles.F` frees the dipole arrays after building `dipoles_opt`.
- **Change:** they are kept until `K_kernel` has used them for the `q -> 0` head (B2), then released.
- **Why:** B2 needs `DIP_iR` inside the vertex routine, which runs after `K_dipoles`.

### B14. The text dump of the Casida matrix

- **Fork:** `src/bse/K_Chi_Kdump.F:166-191` (`X` from `H0`, `F` from `f_xc(z)`).
- **Upstream origin:** the storage of `U:src/bse/K_kernel.F` (factor `i` on coupling) and `U:src/bse/K_stored_in_a_big_matrix.F:79-109`.
- **Change:** `X = H0 - diag(E)` (lines 169-176); `F = c s V^H f_xc(z) V s` over all rows (lines 183-189), split into the blocks `F_rr = Kf_rr`, `F_rc = i Kf_rc`, `F_cr = i Kf_cr`, `F_cc = -Kf_cc` (line 190), written in eV.
- **Why:** the dumped matrix must be the one Yambo would build.
- **Evidence:** tested. Silicon: the eigenvalues of `M(0)` rebuilt from the files match the static run to `4e-5` eV and the response at `z = 0.1i` eV matches the dynamic solver to `1e-6`; hBN at finite q to `2e-7`
  (`doc/greens_function_fxc.md`, "The Casida matrix as text").

### B15. Routing

- **Fork:** `src/bse/K_solvers.F:56` (call the export), `59-67` (send `BSEChiDyn` to `K_Chi_dynamic` and return).
- **Upstream origin:** `U:src/bse/K_solvers.F:38-46`, where a `BSfxc` branch is commented out (B-note below).
- **Why:** a place after the BSE matrix exists and before the requested solver.

*Note on upstream `BSfxc`.* Yambo's own BSE-derived kernel (`U:src/tddft/TDDFT_BSK_2_FXC.F`, Marini 2015) is reachable only through
`tddft_BSK_disk_2_FXC`, whose single call is commented out at `U:src/bse/K_solvers.F:43`; the matching branch of `X_redux` is empty (`U:src/pol_function/X_redux.F:196-197`). **Verified in source.**
It is therefore disabled in 5.3, and it is a different construction from the exact kernel.

---

## 4. Tier C: not in Yambo, from the literature

### C1. The kernel `f_xc = chi0^-1 - P^-1`

- **Fork:** `src/tddft/TDDFT_Chi_dynamic.F:579-583`, `src/bse/K_Chi_export.F:265`.
- **Source:** Del Sole, Adragna, Olevano, Reining, PRB 67, 045207 (2003), Eq. 4; Adragna, Del Sole, Marini, PRB 68, 165108 (2003), Eq. 7. **Equation numbers checked against the PDFs.**
- **Difference:** they evaluated the head only, approximately (their Eq. 6, neglecting 1 against `eps`). The fork evaluates the full `G G'` matrix exactly, by inversion.
- **Justification:** by definition this is the kernel for which TDDFT, starting from `chi0`, gives `P`. With it, `[chi0^-1 - vbar - f]^-1 = chib` is an exact algebraic identity (section 6).

### C2. `P` and the exchange: `P^-1 = chib^-1 + vbar`

- **Source:** Del Sole et al., Eq. 8 (their `P-bar` is `chib`): `P-bar^-1 = chi0^-1 - vbar - f`, with Eq. 4 giving `P^-1 = P-bar^-1 + vbar`. **Checked.**
- **Identity used for the many-body meaning:** with `U = sqrt(c) rho s` and `K_W = -W`, `P_MB = U [z - E - K_W]^-1 U^H` and `chib = U [z - E - X - K_W]^-1 U^H`
  give `chib^-1 + vbar = P_MB^-1` by the push-through (Woodbury) identity. `P_MB` is the irreducible polarizability `-i G G Gamma` with the vertex built from `W` alone
  (Adragna et al., Eqs. 3-4). The derivation is in the LaTeX notes, section "Our P is the irreducible polarizability of many-body theory".

### C3. The two reference choices

- **Source:** Del Sole et al.: `f^FQP` (quasiparticle bubble) and `f^FLDA` (Kohn-Sham bubble). In the fork: `BSEChiOut = "exc"` and `"full"`.
- **Justification:** `"exc"` measures the kernel from the BSE's own energies, so it carries only the electron-hole attraction; `"full"` measures it from the KS bubble, so it also carries the quasiparticle shift.

### C4. TDDFT Dyson equation, mapped to transitions

- **Source:** the TDDFT response equation `chi = [chi0^-1 - vbar - f]^-1` (standard; Casida, in *Recent Advances in Density Functional Methods*, 1995; Onida, Reining, Rubio, Rev. Mod. Phys. 74, 601, 2002) and the matrix-inversion
  (Woodbury) identity (Hager, SIAM Review 31, 221, 1989). **Equation numbers of these three references were not checked here: not verified.**
- **Fork:** `chi = U [z - H0 - K^f(z)]^-1 U^H`, with `H0 = E + U^H vbar U` and `K^f = U^H f U`.
- **Justification:** derived in the LaTeX notes (section "The Casida equation in transition space"); tested for an arbitrary complex, non-Hermitian `f(z)` in `tests/chi_dynamic`, check 1.

### C5. A frequency-dependent kernel in a Casida-type equation

- **Source:** the frequency dependence of the exact kernel is Del Sole et al. and Adragna et al. (the long-range part is frequency dependent, unlike the static kernel of Reining et al., PRL 88, 066404); the poles of the exact kernel between excitons are Gulevich et al., PRB 107, 165118 (2023). **Checked.**
- **Fork:** the per-frequency solve (B9) and the self-consistent root `lambda(omega) = omega` in `sbin/chi_tools/sc_scan.py`.
- **Difference:** the literature computes the head and uses the kernel in `G`-space Dyson equations. Solving a transition-space Casida equation with the full `G G'` kernel at each frequency is the fork's use of it.
  Shi, Turkowski, Rahman (PRB 107, 155431) solve a transition-space Casida equation with a **static** kernel built from `W`; the comparison is in `doc/greens_function_fxc.md`.

---

## 5. Tier D: no source in Yambo or the literature

These are numerical choices and checks. For each, the reason and the evidence.

### D1. The frequency set: `z = 0` exactly, then the BSE grid

- **Where:** `src/bse/K_Chi_export.F:206-207`; loop `244-289`; early exit `285-288`.
- **Why:** the static Casida needs `f(0)` at `z = 0`, not at `z = i eta`; the dynamic solver needs the BSE grid exactly, because `f` is not interpolated (`src/bse/K_Chi_dynamic.F:86-123`).
  Doing `z = 0` first lets a bad basis fail in seconds (on WS2 the frequency loop costs more than the diagonalization).
- **Evidence:** tested (an ill-conditioned run skipped the loop: export step 1 s instead of 77 s, `doc/greens_function_fxc.md`).

### D2. Condition-number monitoring

- **Where:** `src/tddft/TDDFT_Chi_dynamic.F:599-620` (`zgecon`); floor `ChiRcondMin` (`src/modules/mod_Chi.F:65`, default `1e-12`; `src/interface/INIT_load.F:301`); report and stop `src/bse/K_Chi_export.F:285-318`.
- **Why:** `f` is the difference of two nearly cancelling inverses, so an ill-conditioned inversion turns rounding into a spurious kernel. The reciprocal condition numbers of `chi0`, `chib` and `P^-1` are all checked. Below the floor
  the export warns, writes no `ndb.Chi`, and the BSE run still finishes.
- **Evidence:** tested. WS2 (bands 11-16, 9 G, no cutoff) reached `2e-9`: the eight pure-`G_z` plane waves are nearly dependent, and the body of `f` reaches `3e6` in units of `v`
  (`doc/greens_function_fxc.md`, `BSEChiOut`).
- **Caveat:** the default floor `1e-12` is a choice, not derived.

### D3. The in-plane `G` restriction (`BSEGinplane`)

- **Where:** exchange zeroed on `G_z != 0` at `src/bse/K_exchange_collisions.F:172-176`; kernel set `S = {G_z = 0}` at `src/bse/K_Chi_export.F:175-181`; the kernel is embedded into the full matrix at `266-271`.
- **Exactness:** if `vbar` vanishes outside `S`, then `chib = P + P vbar chib` restricted to `S` is a closed Dyson equation,
  `chib_SS = P_SS + P_SS vbar_SS chib_SS`, so `f` on `S` is exact for that BSE.
- **Why:** it removes the nearly dependent out-of-plane vectors that ruin the conditioning of a slab (D2).
- **Caveat:** it changes the BSE itself (no out-of-plane exchange). The kernel is exact for that BSE, not for the original one.
- **Guards:** the Casida run must use the same setting (`src/bse/K_kernel.F:390-392`).

### D4. The head-consistency fit

- **Where:** `src/bse/K_Chi_export.F:323-343`.
- **What:** the `G = 0` element of `chib` is compared, at nine frequencies, with Yambo's own response built from `dipoles_opt` (the `f = 0` limit of the dynamic solver). One complex scale
  `lambda = sum(conjg(rr) head) / sum(abs(rr)^2)` is fitted (line 337); it must equal `c`, with residual (line 338) below `1e-4`.
- **Why:** an independent check that the vertices, phases, anti-resonant convention and normalization match Yambo's. A mismatch here is invisible to the identity of section 6, which holds for any vertices (B2).

### D5. The static Hermiticity test

- **Where:** `src/bse/K_kernel.F:466-485` (stop); report `src/bse/K_Chi_export.F:392-403`.
- **What:** tested on `f` and on `|q+G| f |q+G'|`.
- **Why:** the static Casida matrix must be Hermitian in its resonant block. At `q = 0` the raw head is `O(1/q0^2)` and would hide an anti-Hermitian body.

### D6. The bubble check in the dynamic solver

- **Where:** `src/bse/K_Chi_dynamic.F:167-189` (stop at 185 for a BSE-exported kernel, warning at 187 otherwise).
- **What:** the Casida run's own bubble is compared with the stored `chi0` in the symmetrized metric.
- **Why:** the identity holds only with the same `chi0`; this tests it directly. The threshold `1e-3` is a choice.

### D7. Export-to-Casida guards

Every precondition of the identity is checked by the Casida run before any matrix element is computed (`src/bse/K_kernel.F`):

| precondition | line | threshold |
| --- | --- | --- |
| optical `q0` direction and normalization | 376-379 | `1e-3` relative |
| `BSENGfxc` equals the `ndb.Chi` basis | 381 | exact |
| `BSENGexx` equals the basis, exchange on | 385 | exact |
| `Lkind` equals the export's (read from the convention string) | 389 | exact |
| `BSEGinplane` equals the export's | 392 | exact |
| QP corrections present exactly when the export had them | 401, 403, 407 | exact |
| occupations (bands, k points, spins) | 415 | exact |
| occupations and temperature | 421 | `1e-6` |
| band energies of the reference bubble | 429 | `1e-5` Ha |
| the fragment belongs to this header (export serial) | 443 | exact |
| number of transitions at this `q` | 447 | exact |
| first frequency exactly zero (static) | 465 | `1e-7` |
| `f(0)` Hermitian (D5) | 485 | `1e-4` |

- **Why:** each is a precondition of the identity. Without them a mismatched Casida run gives a plausible but wrong spectrum. Several were added after a case in which the run went through silently: for example, a Casida run with a smaller transition window
  (8 of 576 transitions) ran 76% off before the transition-count check existed, and a 0.1 K difference in temperature changed the occupations (`doc/greens_function_fxc.md`, "Kernel from the BSE" and "Partly filled bands").
- **Caveat:** all thresholds are choices.

### D8. The head structure report

- **Where:** `src/bse/K_Chi_export.F:405-433`.
- **What:** the local minima of `|P_00|` (response zeros) and the local maxima of `-Im q^2 f_00` (kernel poles), up to six each, on the BSE grid.
- **Why:** the poles are physical (Gulevich et al.) and explain the jumps of `lambda(omega)` in the root search; reporting them keeps them from being taken for numerical failures. Diagnostics only. The 20% and 50% thresholds are conventions.

### D9. Precision

- **Where:** all new linear algebra is `complex(DP)`; `ndb.Chi` stores `chi0`, `P` and `f_xc` as `complex(SP)` (`src/modules/mod_Chi.F:78-80`); the input `BS_mat` is single precision.
- **Why:** to keep the two inversions accurate (D2).
- **Limit:** the result cannot be more accurate than about `1e-7` relative, the precision of `BS_mat`. The stored `f_xc` has a head up to about `5e10` on WS2 (`1/q0^2` scaling).

### D10. The `[BSE]` binding report

- **Where:** `src/bse/K_diago_driver.F:343, 377-484`.
- **What:** lowest excitations, strengths `|R_L R_R|` relative to the strongest, the lowest bright one (above 1% of the maximum), and the binding energy against the lowest transition of the window.
  For partly filled bands, the interband weight `sum_{t in full} |V_{t,lambda}|^2 / sum_t |V_{t,lambda}|^2 > 0.5` decides which excitations count.
- **Why:** with partly filled bands the lowest transition of the window is not a gap.
- **Caveat:** the 1% and 50% thresholds are conventions, not derived.

### D11. The `ndb.Chi` metadata

- **Where:** `src/bse/K_Chi_export.F:469-495`; `src/io/io_Chi.F`.
- **What is stored besides `chi0`, `P`, `f_xc`:**
  - an export serial (line 473): tenths of a second in the month, modulo `2^24 - 1 = 16777215`, so that it is exact in the single-precision record; written in the header and in every `q` fragment;
  - the occupations of the BSE levels and the temperature (lines 475-482), and the band energies of the reference bubble (lines 484-485);
  - the number of resonant transitions per `q`; the optical `q0`; the convention string (Lbar or Lfull, QP action).
- **Why:** each lets the Casida run check a precondition (D7). Older databases without them remain readable and are not checked.

### D12. The kernel embedded on the set `S`

- **Where:** `src/bse/K_Chi_export.F:266-271`.
- **What:** `fs` and `Ps` live on `S`; `f` and `Pm` are zero outside it.
- **Why:** `ndb.Chi` has a fixed `ng x ng` shape; with `BSEGinplane` the kernel is exact on `S` only (D3).

### D13. Threading and parallelism

- **Where:** `!$omp` frequency loops at `src/bse/K_Chi_export.F:247, 250`, `src/bse/K_Chi_dynamic.F:201, 203`; MPI guards `K_Chi_export.F:121`, `K_Chi_dynamic.F:72`.
- **Why:** each frequency is independent, so OpenMP over frequencies is exact; one MPI rank because the BSE matrix is assembled serially (B10).

### D14. The Lkind and compiler fixes ported from Yambo 5.4

- **Where:** `src/bse/K_driver_init.F:88-100, 224-227` (the default `Lkind` at `q /= 0`, `Lkind = "tilde"` accepted); `src/bse/K_correlation_collisions_std.F`, `_dir.F` (`intent(out)` to `intent(in)` for `j_k`, `j_p`).
- **Source:** Yambo 5.4 commits `ecd4a9667`, `1fceec807`, `d71faa3d3` (Paleari) and `81d6a687d` (Varsano). They are the Yambo authors' changes, ported unchanged; they are not fork work.
- **Fork additions:** a stop with `Lkind = "tilde"` in the export (`src/bse/K_driver_init.F:74-75`) and a check that the Casida `Lkind` equals the export's (`src/bse/K_kernel.F:389`).

---

## 6. Why the numerics hang together (the identity)

With the same transitions, vertices and energies in the export and in the Casida run (so the same `chi0`, D7):

```
chi_Casida = [chi0^-1 - vbar - f_xc]^-1 = [chi0^-1 - vbar - chi0^-1 + chib^-1 + vbar]^-1 = chib.
```

- It holds **by construction**, at every `z`, exciton poles included. A Casida spectrum that matches its BSE checks that the export, the file and the solver are consistent.
- It does **not** test the physics of the BSE: an error in `W`, the QP energies or the occupations is reproduced just as exactly.
- It holds for the dynamic solver and for the self-consistent root, not for a kernel frozen at one frequency (the static Casida: silicon 0.043 eV binding against 0.173 eV).

---

## 7. Coverage: every source file of the fork

Classified at file level. "Detailed" means the sections above; "summary" means classified here only.

### 7.1 New Fortran files (14)

| file | tier | what it is | depth |
| --- | --- | --- | --- |
| `src/tddft/TDDFT_Chi_dynamic.F` | B5-B9, C1-C2, D2 | linear algebra in transition space | detailed. Also holds `"mdr"`/`"nq"` routines, the Tamm-Dancoff solver and unused pair-space routines (section 7.4) |
| `src/tddft/TDDFT_Chi_Fourier.F` | B4 | the pair projection; also `TDDFT_Chi_project_G/R` (unused outside tests) | detailed |
| `src/bse/K_Chi_export.F` | B5-B11, C1-C3, D1-D4, D8, D11-D13 | the export driver | detailed |
| `src/bse/K_Chi_dynamic.F` | B9-B10, D6, D13 | the dynamic Casida driver | detailed |
| `src/bse/K_Chi_Kdump.F` | B14 | the text dump | detailed |
| `src/modules/mod_Chi.F` | D | data module for `ndb.Chi`; defaults (`Chi_rcond_floor`) | detailed |
| `src/io/io_Chi.F` | A10, D11 | database I/O in Yambo's pattern | summary: no numerics |
| `src/pol_function/Chi_G_bubble.F` | B, C, D | Green-function route: band-pair bubble from spectral functions; follows `X_irredux_residuals` conventions (stated in its header) | summary |
| `src/pol_function/Chi_fxc.F` | A, B | Green-function route: re-runs `X_irredux` with the KS levels for `chi0`, builds `P` and `f_xc` | summary |
| `src/modules/mod_Chi_Green.F` | C, D | numerical kernels for a diagonal dressed bubble (spectral quadrature, causality checks) | summary |
| `src/pol_function/W_potential.F` | B | `ndb.W`: the same construction as `QP_ppa_cohsex`, `K_screened_interaction`, `LSEX_potential` (stated in its header) | summary |
| `src/io/io_W.F` | A10 | `ndb.W` I/O | summary: no numerics |
| `src/modules/mod_Wpot.F` | D | data module for `ndb.W` | summary |
| `src/qp/W_Sigma_check.F` | B, D | recomputes `Sigma_c` from `ndb.W` and compares it with `QP_ppa_cohsex`: a consistency check | summary |

### 7.2 Edited upstream files (36 source files)

| file | tier | change | depth |
| --- | --- | --- | --- |
| `src/bse/K_kernel.F` | B3, B4, B12, D5, D7 | see B12 | detailed |
| `src/bse/K_solvers.F` | B15 | routing | detailed |
| `src/bse/K_dipoles.F` | B13 | keeps the dipoles | detailed |
| `src/bse/K_exchange_collisions.F` | D3 | `BSEGinplane` | detailed |
| `src/bse/K_driver_init.F` | D7, D14 | input checks; ported 5.4 `Lkind` | detailed |
| `src/bse/K_diago_driver.F` | D10 | binding report | detailed |
| `src/bse/K_correlation_collisions_std.F`, `_dir.F` | D14 | ported 5.4 `intent` fix | detailed |
| `src/tddft/TDDFT_ALDA_eh_space_G_collisions_L.F`, `_R.F` | B1, B2 | bare vertex | detailed |
| `src/tddft/TDDFT_do_X_W_typs.F` | D | skips the local TDDFT setup for the CHI kernel (`l_chi_bs_fxc`) | summary |
| `src/modules/mod_BS.F` | D | flags `BS_K_is_Chi`, `BS_K_Chi_dynamic`, `BS_K_Chi_out`, dump window, `BS_G_inplane` | detailed |
| `src/modules/mod_drivers.F` | D | `l_Wpot`, `l_chi_bs_fxc` | summary |
| `src/modules/SET_defaults.F` | D | default `BS_K_is_Chi` | summary |
| `src/output/K_OUTPUT.F`, `src/output/X_OUTPUT_driver.F` | A9 | output name `chidyn`; observables off | summary |
| `src/pol_function/X_dielectric_matrix.F` | D | calls `Chi_fxc_eval` in the q loop (`ChiFxc`) | summary |
| `src/pol_function/X_pre_setup.F` | D | `ChiFxc` input checks; forces the recalculation of the response | summary |
| `src/coulomb/cutoff_driver.F` | D | two MPI barriers around writing `ndb.cutoff` (a race) | summary |
| `src/interface/INIT_load.F`, `INIT_activate.F`, `INIT.F`, `INIT_barriers.F`, `INIT_read_command_line.F` | D | input variables, activation, runlevel `Wpot`, `-k chi` | summary: plumbing |
| `src/driver/options_yambo.c`, `driver/yambo.F` | D | the `-W` option; saves the KS levels before the QP correction | summary |
| `src/Yio/io_fragment.F`, `io_fragment_disconnect.F` | A10 | fragments for `ndb.Chi` and `ndb.W` | summary |
| `src/qp/QP_driver.F`, `QP_Green_Function.F`, `QP_ppa_cohsex.F`, `QP_newton.F`, `QP_descriptions.F`, `src/qp_control/QP_load_GFs.F`, `src/modules/mod_QP.F`, `src/io/io_QP_and_GF.F` | B, C, D | Green-function route (`GRetarded`, conventional retarded Dyson `G = 1/(z - eps - Delta Sigma)`); `Wpot` hook in `QP_newton` | summary |

The six `.objects` files (`src/{bse,io,modules,pol_function,qp,tddft}/.objects`) only register the new sources; `.gitignore` un-ignores `tests/chi_green/*.F90`.

### 7.3 Not source: scripts, documentation, tests

`sbin/chi_tools/**` (26 scripts), `doc/*.md` (4), `tests/**` (14): outside this note. They contain no Yambo numerics.

### 7.4 Written but not used by the production path

| routine | file | status |
| --- | --- | --- |
| `TDDFT_Chi_pair_response`, `TDDFT_Chi_pair_bubble` | `src/tddft/TDDFT_Chi_dynamic.F` | no production caller (the comment at `src/bse/K_Chi_export.F:15` still names `TDDFT_Chi_pair_response`); used by tests |
| `TDDFT_Chi_project_G`, `TDDFT_Chi_project_R` | `src/tddft/TDDFT_Chi_Fourier.F` | no production caller; used by tests |
| `TDDFT_Chi_dynamic_response`, `TDDFT_Chi_dynamic_kernel` | `src/tddft/TDDFT_Chi_dynamic.F` | Tamm-Dancoff dynamic solve: reachable with `BSEmod = "resonant"` and `BSEChiDyn`, warned against (`src/bse/K_Chi_dynamic.F:192`); not used in any reported run |
| `TDDFT_Chi_mdr_response`, `TDDFT_Chi_kernel_mdr`, `TDDFT_Chi_exchange_free_split`, `TDDFT_Chi_first_order_kernel` | `src/tddft/TDDFT_Chi_dynamic.F` | `"mdr"` and `"nq"` kernels (outside the scope of this note) |

---

## 8. Open points: what is not verified, and what to double-check

1. **B1, the `scatter_Bamp` vertex.** The reason given in the code (USPP augmentation, spinor sum) is a code comment. Not tested: whether upstream's FFT product would give a different vertex for the pseudopotentials used.
2. **B8, thread safety of Yambo's LU wrapper.** The fork calls LAPACK directly; it was not tested whether the wrapper would also work inside the OpenMP frequency loop.
3. **C4, textbook references.** Casida (1995), Onida-Reining-Rubio (2002) and Hager (1989): equation numbers not checked against the sources.
4. **B6, "only the macroscopic response".** Verified for the diagonalization, inversion and Haydock solvers; not for the SLEPc path beyond its use of the diagonalization residuals.
5. **Novelty.** Not claimed (0.5). A paper outside the project's set may contain Tier C or D items.
6. **D-tier thresholds** (`1e-3`, `1e-4`, `1e-5`, `1e-6`, `1e-12`, 1%, 20%, 50%) are choices, not derivations. D2's default is the only one tied to a measured case (WS2).
7. **Not tested:** `"full"` at finite q; partly filled bands with `"full"` (a guard stops it, `src/bse/K_driver_init.F:81-83`); the Tamm-Dancoff dynamic path; more than one MPI rank (refused).
8. **Tier B items rest on structural resemblance.** "Adapted from" means the same construction with identified differences. It does not mean the fork's code was derived line by line from Yambo's.

---

## 9. Verification record

### 9.1 Line references

Every `path:line` and `U:path:line` in this note was checked against the source (fork at `df32790`, upstream at `19c1241`) by a script that prints the cited line for review.

### 9.2 Tests of the numerics

| test | what it shows |
| --- | --- |
| `tests/chi_dynamic` (check 1) | the transition-space response equals the `G`-space Dyson solution for an arbitrary complex, non-Hermitian `f(z)` |
| `tests/chi_dynamic` (check 2) | the static limit equals the Lorentzian sums of `K_diago_response_functions` (TDA and coupling) |
| `tests/chi_fxc_fourier` | the pair projection against a small grid |
| silicon, `BSEChiDyn` | against the QP RPA+LF response of `ndb.Chi`: `3e-4`; "exc" and "full" dynamic Casida against the BSE: `1.4e-4` |
| WS2 with SOC | dynamic Casida against its BSE: `Im eps` `7e-5`, `Re eps` `1e-4` |
| silicon, text dump | eigenvalues of `M(0)` against the static run `4e-5` eV; response at `z = 0.1i` eV `1e-6` |
| silicon, root search | `lambda(omega) = omega` at 3.4276 eV against the BSE 3.4278 eV; WS2 2 meV grid: exact to `< 1 ueV` |

### 9.3 The reading copy of the numerics

The branch `numerics-reading` holds the numerics of this path without comments, guards or messages. It was checked three ways:

1. its `TDDFT_Chi_dynamic.F` against the fork's, on random data: bit for bit;
2. all its changed files type-checked against the fork's module files;
3. a Yambo binary built with nine of its files in place of the fork's reproduced the fork's WS2 `ndb.Chi` (`chi0`, `P`, `f_xc`, 322 x 21 x 21), the BSE spectrum, the static Casida spectrum and the dynamic Casida spectrum, with a maximum difference of 0.

### 9.4 How to check an item yourself

```
git diff 19c1241 HEAD -- <file>                      # what the fork changed in an upstream file
git show 19c1241:<file> | sed -n '<a>,<b>p'          # an upstream line range
sed -n '<a>,<b>p' <file>                              # a fork line range
```
