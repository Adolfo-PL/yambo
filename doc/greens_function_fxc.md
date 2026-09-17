# Green-function bubble route to fxc

This change adds an initial diagonal, frozen-screening GW route to the
`ChiFxc` export. The existing level-correction route remains the default:
`ChiGMode="LEVELS"`. The new modes are `G0` and `DYSON`.

Base: Adolfo-PL/yambo, branch `claude/nbd-file-generation-path-yidcr2`,
commit `2d64938d8328b184ffbf7ce82972d737a8e5c12a`.

## Implemented approximation

The orbitals, KS reference, chemical potential and screening are fixed. The
existing PPA GW machinery evaluates the diagonal correlation self-energy on a
real energy grid. With the new `GRetarded` flag it uses retarded denominators
for both occupied and empty intermediate states. The exported propagator is

```
DeltaSigma_nk(z) = Sigma_x,nk + Sigma_c,nk(z) - Vxc_KS,nk
G_nk(z) = 1 / [z - epsilon_KS,nk - DeltaSigma_nk(z)]
A_nk(E) = -Im G_nk(E + i eta) / pi
```

All energies here are absolute in Yambo's KS energy reference. The chemical
potential enters the occupations. Equivalently, energies measured relative
to the chemical potential give `G(z)=1/[z+mu-epsilon_KS-DeltaSigma]`.
Positive `GDmRnge` controls eta; `GDamping` must be zero so Sigma and G are
evaluated at the same complex energy.

For each finite q, the response loops over the complete requested band-pair
range and BZ k points, using Yambo's `qindx_X` mapping k -> k-q and density
vertices from `scatter_Bamp`. The frequency factor for a pair a,k and b,k-q is

```
integral dE dE' A_a,k(E) A_b,k-q(E')
    * [f(E') - f(E)] / [z - E + E']
```

It multiplies `rho*(G) rho(G')/(Nk_BZ V_cell)`, including the appropriate
spin factor. Spectral quadrature uses nonuniform trapezoid weights. It does
not divide spectra by their integrated weight. Missing spectral tails or
poor resolution trigger an error when `ChiGNormTol` is exceeded.

The independent KS bubble uses the analytic G0 pole limit of this expression.
Before export, it must agree with the native `X_irredux` KS response within
`ChiGCheckTol`. The export contains

```
ChiGMode="G0":    chi0=P_G0; fxc=0
ChiGMode="DYSON": chi0=P_G0; P=P_G; fxc=chi0^-1-P_G^-1
```

Inversions use double-precision, Coulomb-symmetrized matrices with LAPACK
factorization, condition estimation and inversion checks. No pseudoinverse
or regularization is applied. An ill-conditioned basis is rejected.

This is a **bare-vertex dressed-bubble kernel**. It does not implement
`delta Sigma/delta G`, `delta W/delta G`, a Bethe-Salpeter vertex, or the full
physical GW density response. It should not be interpreted as the final
excitonic kernel from the background derivation. No extra `-v` enters the
inverse difference because P is irreducible here. KS G0 is evaluated
analytically; no large real-space G0 database is constructed.

## Build and first run

Apply the patch to the base checkout, then rebuild Yambo using your usual
Linux/HPC configuration, including NetCDF and LAPACK. Both new source files
are registered in their directories' `.objects` lists. A clean rebuild is
recommended because the QP and Chi module types changed.

Start with one CPU MPI rank, a small insulating SAVE database containing
finite q points, a small response band range and a small density basis.
The spectral convolution costs approximately
`Nk * Nb^2 * NE^2 * N_response_frequencies`, before the density-matrix
accumulation. Gamma-only SAVE data cannot run this finite-q prototype.

### 1. Export retarded GW spectra

Generate a PPA GW input using the Green-function Dyson solver:

```sh
yambo -p p -g g -F gw_g.in -J gw_g
```

Retain your converged screening and GW cutoffs, and set/add the following
entries. These energy bounds and sample counts are illustrative: choose a
range that covers the poles and satellites of every requested state, then
refine until the spectral norm and the response converge.

```text
GRetarded
DysSolver="g"
GEnMode="absolute"
GEnSteps=1001
% GEnRnge
 -40.0 | 40.0 | eV
%
% GDmRnge
 0.2 | 0.2 | eV
%
GDamping=0.0 eV
GreenFTresh=0.0
GTermKind="none"
```

Remove the `GreenF2QP` flag. Do not apply `GfnQPdb` energy shifts. Set
`QPkrange` to cover **every k point, band and spin** that will enter the
response, with diagonal states only. Use the generated input's syntax for
your spin configuration. Set the internal `GbndRnge` independently to the
converged self-energy band range.

Run the edited input:

```sh
yambo -F gw_g.in -J gw_g
```

The database is normally `gw_g/ndb.G`. This step consumes the native PPA
screening databases (`ndb.pp` and associated `ndb.em1s`/fragments) underlying
your W. Keep these companion databases and use the same screening settings
to reuse them. **It does not import the standalone `ndb.W` export.** If the
native screening databases are missing or incompatible, Yambo's usual
screening machinery recalculates them. Do not enable `ChiFxc` in this GW
step; it forces response recalculation and the new bubble path excludes q=0.

The old `ndb.G` convention is not accepted by `ChiGDb`; rebuild it with
`GRetarded`. The new convention is also rejected by the historical
`XfnQPdb` Green-function remapper: use `ChiGDb` for this method.

### 2. Validate the independent KS bubble

Generate a dynamic dielectric input:

```sh
yambo -d d -F chi_g0.in -J chi_g0
```

Use the same SAVE and KS reference. Add:

```text
ChiFxc
ChiGMode="G0"
ChiGAxis="IMAG"
ChiGNuSteps=16
% ChiGNuRnge
 0.0 | 20.0 | eV
%
ChiGCheckTol=0.0001
ChiGNormTol=0.05
ChiRcondMin=1.0e-12
XTermKind="none"
% EhEngyXd
 0.0 | 0.0 | eV
%
```

Choose `QpntsRXd` entirely after q index 1; for an initial test, use one
existing finite q. Set `BndsRnXd` to the intended response bands. Disable
double-grid interpolation and external X energy, width, residue and
Green-function corrections. The implementation selects retarded ordering
and zero coarse-grid spacing for the KS baseline.

```sh
yambo -F chi_g0.in -J chi_g0
```

The run checks the native KS response and exports zero fxc. `IMAG` means a
continuous zero-temperature bosonic grid z=i nu, including nu=0 when
requested. This export grid does not change the native dielectric frequency
grid or replace the matrix used to calculate screening. `ChiGAxis="REAL"`
instead exports on the native response grid in the upper half plane; use
positive response damping to avoid real-axis poles.

### 3. Build the Dyson bubble kernel

Copy the validated response input to `chi_dyson.in` and change/add:

```text
ChiGMode="DYSON"
ChiGDb="gw_g/ndb.G"
```

```sh
yambo -F chi_dyson.in -J chi_dyson
```

The loader checks the k mesh and ordering, band/spin coverage, diagonal
state table, KS starting energies, chemical potential, positive damping,
causality, spectral norms, and consistency of exported G with its Dyson
self-energy. A partial QPkrange fails explicitly. `ChiGNormTol=0.05` is an
initial diagnostic threshold, not a convergence claim: tighten it while
expanding/refining the spectral grid and reducing eta.

## Database additions

The existing variables and per-q fragment layout of `ndb.Chi` are retained:

- `CHI_FREQ_Q_<iq>`: actual complex export frequencies in Hartree.
- `CHI0_Q_<iq>`, `P_Q_<iq>`, `FXC_Q_<iq>`: unsymmetrized PW matrices.
- `CHI_RCOND_Q_<iq>`: reciprocal condition estimates for chi0 and P,
  shape (2, Nfreq).
- `CHI_G0_ERROR_Q_<iq>`: relative G0/native-KS comparison error; -1 means
  the comparison was not performed in LEVELS mode.
- Header: `CHI_G_MODE`, `CHI_G_AXIS`, `CHI_G_DB`, and
  `CHI_G_PARAMETERS` = [norm tolerance, KS comparison tolerance, rcond
  threshold, chemical potential, kBT].

Older Chi files remain readable without the new optional metadata; absent
conditioning diagnostics read as zero and absent G0 diagnostics as -1.

Retarded `ndb.G` additionally stores `G_RETARDED`,
`G_CHEMICAL_POTENTIAL`, `G_REFERENCE_ENERGIES`, `G_SIGMA_X` and
`G_VXC_REFERENCE`. Existing `SE_Operator` stores DeltaSigma. Recover
Sigma_c as `SE_Operator - G_SIGMA_X + G_VXC_REFERENCE`, broadcasting the
static terms over frequencies. `Green_Functions_Energies` includes eta.

## Current scope and validation

The response path explicitly rejects MPI distribution, GPU execution,
q=0, finite temperature, metallic KS references, double grids,
transition-energy filtering and response terminators. The retarded GW
export requires PPA, diagonal Sigma, unshifted KS starting energies and no
GW terminator, Green-function zoom, GreenF2QP, self-consistent GW or mixed
electron-phonon/photon self-energies. Off-diagonal Sigma, optical-limit
vertices and the GW response vertex remain separate extensions.

Portable tests compile the production numerical module:

```sh
python tests/chi_green/run_tests.py --fc gfortran
```

To compile and exercise the production bubble, inversion and Dyson routines
with actual Reference LAPACK and synthetic material/IO fixtures:

```sh
python tests/chi_green/run_tests.py --fc gfortran --lapack-source /path/to/lapack
python tests/chi_green/run_tests.py --fc gfortran --single-precision --lapack-source /path/to/lapack
```

The tests check analytic poles and satellites, an independent fermionic
Matsubara sum, the G0 limit, causality, quadrature, XC subtraction, the
imaginary export grid, inverse response differences, screening restoration,
and rejection of missing states, legacy spectra, reference mismatches and
ill-conditioned inversions. They passed with GNU Fortran 16.2.0 and
Reference LAPACK 3.12.1 in both host precision settings.

These tests use mock material and database IO. They do not establish a
complete Yambo build, NetCDF round-trip, or a material benchmark. Those
checks require your Linux/HPC build and SAVE/screening data. Begin with the
G0 run, then converge spectral coverage/resolution, eta, bands, k mesh and
density basis before drawing physical conclusions from the Dyson kernel.
