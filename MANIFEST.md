# numerics-reading: the numerics of the exact kernel, without the fluff

This branch is Yambo 5.3.0 (base commit `19c1241`) plus the **numerics** of the exact kernel,
written to be read, not to be built or run as a product. It is generated from the fork branch
`claude/nbd-file-generation-path-yidcr2` (commit `611894a`).

`git diff 19c1241 numerics-reading` shows everything the fork adds on this path.

## What the path covers

The exact kernel exported from a BSE, `BSEChiOut = "exc"` or `"full"`, at any q (explicit
anti-resonant rows included), and the Casida solvers that use it: static (`BSKmod = "CHI"`) and
frequency-dependent (`BSEChiDyn`), with coupling, length gauge, scalar or spinor states.

## Size

| | fork | this branch |
| --- | --- | --- |
| new Fortran files on this path | 5 files, 1,667 lines | 5 files, 480 lines (322 code) |
| edits to upstream files | 43 files | 7 files, +215 / -52 lines |

(The fork has 14 new Fortran files in all; the others belong to the routes left out below.)

## What each file does

| file | role |
| --- | --- |
| `src/tddft/TDDFT_Chi_dynamic.F` | linear algebra in transition space: bare bubble, G-resolved BSE response, kernel from the response (`fxc = chi0^-1 - (chib^-1 + vbar)`), per-frequency Casida solve |
| `src/tddft/TDDFT_Chi_Fourier.F` | projection of `fxc(G,G')` onto two transition vertices (static Casida) |
| `src/bse/K_Chi_export.F` | `BSEChiOut`: builds H0, vertices, reference energies, runs the frequency loop, writes `ndb.Chi` |
| `src/bse/K_Chi_dynamic.F` | `BSEChiDyn`: one linear solve per frequency with `fxc(z)` read from `ndb.Chi` |
| `src/modules/mod_Chi.F` | the arrays of `ndb.Chi` and the transition vertices |
| `src/bse/K_kernel.F` (edit) | reads `ndb.Chi`; collects the vertices; adds the static projection to the matrix elements |
| `src/tddft/TDDFT_ALDA_eh_space_G_collisions_L/R.F` (edit) | "bare vertex" mode: transition density from `scatter_Bamp`, optical head from the dipoles |
| `src/bse/K_solvers.F` (edit) | calls the export; sends `BSEChiDyn` to `K_Chi_dynamic` |
| `src/bse/K_dipoles.F` (edit) | keeps the position dipoles for the q -> 0 head |
| `src/bse/K_exchange_collisions.F` (edit) | `BSEGinplane`: exchange vertex zero on G_z /= 0 |
| `src/modules/mod_BS.F` (edit) | four flags: `BS_K_is_Chi`, `BS_K_Chi_dynamic`, `BS_K_Chi_out`, `BS_G_inplane` |

Variable names are those of the fork. Each file starts with a legend; each routine with one equation.

## What was removed, file by file

New files:

- `TDDFT_Chi_dynamic.F` (622 -> 149 lines)
  - routines: `TDDFT_Chi_dynamic_response` and `_kernel` (Tamm-Dancoff dynamic solve), `TDDFT_Chi_pair_response`,
    `TDDFT_Chi_pair_bubble`, `TDDFT_Chi_mdr_response`, `TDDFT_Chi_kernel_mdr`, `TDDFT_Chi_exchange_free_split`,
    `TDDFT_Chi_first_order_kernel` (the last four belong to `"mdr"` and `"nq"`);
  - inside the kept routines: array-size checks and the `ierr` codes, the reciprocal condition numbers
    (`zgecon`, `rcond`, `rcond_P`), and the `exch` argument (always `.TRUE.` on this path).
- `TDDFT_Chi_Fourier.F` (122 -> 31): `TDDFT_Chi_project_G` and `_project_R` (used only by tests), size checks.
- `K_Chi_export.F` (507 -> 154):
  - `"mdr"` and `"nq"` (exchange-free split, shifted bubble, first-order kernel and their checks);
  - the two-pass loop (z = 0 first) with the condition-number floor and early exit;
  - the head-consistency fit against Yambo's own response (nine frequencies);
  - the static Hermiticity report, the head zero and pole report;
  - the `ndb.Chi` metadata: `Chi_G_mode`, convention string, QP action, optical q0, occupations, band energies,
    number of transitions, export serial;
  - MPI, gauge, widths, G=0 and dimension guards; section, messages, timing.
- `K_Chi_dynamic.F` (245 -> 86): MPI, gauge, absorption-only and widths checks; the frequency-grid check;
  the Tamm-Dancoff branch; the bubble check against `ndb.Chi`; the `BSEChiKoutW` call; descriptors; timing.
- `mod_Chi.F` (171 -> 60): everything not used here (G-route data, occupations, flags, descriptors).

Edits to upstream files:

- `K_kernel.F`
  - all `ndb.Chi` guards: optical q0, G basis, exchange G, `Lkind`, `BSEGinplane`, QP corrections, occupations,
    band energies, serial, number of transitions, first frequency, Hermiticity of fxc(0);
  - the restart errors, the MPI, spin and anti-resonant-block errors, the section titles and messages;
  - what stays: reading the header and fragment of `ndb.Chi`, the static/dynamic branch, the vertex
    collection, the projection, the frees.
- `K_solvers.F`: the change of the `K_diago_driver` argument list (it feeds the binding report, left out).
- `K_dipoles.F`: the guard on the G=0 head.
- The collision files: spin-polarization and magnon errors, the dipole-loaded check.

Not included at all (outside this path):

- the Green-function route: `Chi_G_bubble.F`, `Chi_fxc.F`, `mod_Chi_Green.F`, the `GRetarded` edits in
  `src/qp`, `io_QP_and_GF.F`, `mod_QP.F`;
- `ndb.W`: `W_potential.F`, `io_W.F`, `mod_Wpot.F`, `W_Sigma_check.F`;
- file I/O: `io_Chi.F` (its calls stay in `K_Chi_export.F` and `K_kernel.F`), `io_fragment*.F`;
- the text dump `K_Chi_Kdump.F` and the `[BSE]` binding report in `K_diago_driver.F`;
- input plumbing: `INIT*.F`, `options_yambo.c`, `driver/yambo.F`, `SET_defaults.F`, `mod_drivers.F`,
  `K_driver_init.F`, `X_dielectric_matrix.F`, `X_pre_setup.F`, `TDDFT_do_X_W_typs.F`, the output-name edits,
  the `.objects` lists, `cutoff_driver.F`, the two `K_correlation_collisions_*` files (a compiler-bug fix);
- scripts, documentation, tests.

## How "same numerics" was checked

1. **The linear algebra, on random data.** `TDDFT_Chi_dynamic.F` of this branch was compiled next to the
   fork's original and both were called on random complex data (14 transitions, 28 rows, 5 plane waves,
   both values of `exch_head`): bubble, G-resolved BSE response, kernel, polarization and dynamic response.
   Maximum difference: 0, bit for bit.
2. **Types and argument lists, against the real modules.** All 12 changed files were preprocessed and
   compiled with `gfortran -fsyntax-only` against the module files of the fork's build, with the reduced
   `Chi_m` and `BS` modules first in the search path. All compile. The check can fail: the fork's own
   `K_kernel.F` is rejected against the reduced `TDDFT_Chi_Fourier` module (it has another argument list),
   and a deliberate argument error is caught.
3. **End to end, on the WS2 test.** A Yambo binary was built with nine of these reading copies in place of the
   fork's files (`TDDFT_Chi_dynamic`, `TDDFT_Chi_Fourier`, `K_Chi_export`, `K_Chi_dynamic`, `K_kernel`,
   `K_dipoles`, `K_exchange_collisions` and the two collision files; the fork's `mod_BS`, `mod_Chi`,
   `io_Chi` and `K_solvers` stayed). On monolayer WS2 with spin-orbit coupling (bands 23-30, 21 G,
   `BSEGinplane`, QP scissor), `BSEChiOut = "exc"` at q = 0, then the static and the dynamic Casida run:

   | compared with the fork's results | max difference |
   | --- | --- |
   | frequencies, chi0, P and fxc of `ndb.Chi` (322 x 21 x 21) | 0, bit for bit |
   | BSE spectrum of the export run | 0 |
   | static Casida spectrum | 0 |
   | dynamic Casida spectrum | 0 |

Not tested: `"full"`, finite q, partly filled bands, and the removed failure paths (the guards stop a run; the
reading copies simply do not). The reduced `mod_BS.F`/`mod_Chi.F` were type-checked, not linked into the run.
The branch as a whole does not build: files outside this path are missing on purpose.

## Build branch (`numerics-reading-build`): what is added beyond the reading copy

This branch is `numerics-reading` plus only the infrastructure needed to **build and switch on** the same path.
The numerics files are the reading copies, with the guards still removed.

| added | purpose |
| --- | --- |
| `src/io/io_Chi.F` (fork, unchanged) and the full fork `src/modules/mod_Chi.F` | `ndb.Chi` read/write; the reduced module lacks what `io_Chi` needs |
| `"Chi"` in `src/Yio/io_fragment*.F` | fragment I/O for `ndb.Chi` |
| `.objects` in `src/bse`, `src/io`, `src/modules`, `src/tddft` | the new files are compiled |
| `l_chi_bs_fxc` (`mod_drivers.F`, `INIT*.F`), `BS_K_is_Chi` default (`SET_defaults.F`) | `BSKmod = "CHI"`, `-k chi` |
| `BSEChiOut`, `BSEChiDyn`, `BSEGinplane` (`INIT_load.F`, `INIT_activate.F`, `K_driver_init.F`) | the input variables |
| `TDDFT_do_X_W_typs.F`, `K_OUTPUT.F`, `X_OUTPUT_driver.F` | skip the local TDDFT setup; output name `ChiDyn` |

Still left out: `BSEChiKoutW/B` (text dump), the Green-function route, `ndb.W`, all guards and the `Lkind` checks.
Without the guards a mismatched Casida run is not stopped.

### Check: built and run on the WS2 test

The branch was configured and built in its own tree with the fork's configure options (gfortran, MPI + OpenMP, NetCDF/HDF5).
On the same WS2 inputs as above (`BSEChiOut = "exc"` at q = 0, then static and dynamic Casida), compared with the fork's binary:

| compared with the fork's results | result |
| --- | --- |
| `CHI_FREQ`, `CHI0`, `P`, `FXC` of `ndb.Chi` (322 x 21 x 21) | identical, bit for bit |
| BSE spectrum of the export run | identical |
| static Casida spectrum | identical |
| dynamic Casida spectrum | identical |

`ndb.Chi` written by this branch lacks two metadata fields the fork writes (`CHI_NTRANS`, `CHI_SERIAL`); they feed only the guards.
Not tested: `"full"`, finite q, partly filled bands.
