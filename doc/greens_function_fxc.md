# Green-function bubble route to fxc

This change adds an initial diagonal, frozen-screening GW route to the
`ChiFxc` export. The existing level-correction route remains the default:
`ChiGMode="LEVELS"`. The G modes are `G0`, `COHSEX` and `DYSON`.
`COHSEX` works with native static screening; `DYSON` uses PPA dynamic screening.

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
evaluated at the same complex energy. The retarded PPA path preserves this
imaginary part when preparing its self-energy sampling grid.

For static W, Yambo's native COHSEX/Newton solver produces diagonal energies
`E_COHSEX = epsilon_KS + Sigma_SEX + Sigma_COH - Vxc_KS` and unit residues.
`ChiGMode="COHSEX"` reads its `ndb.QP` and constructs

```
G_nk(z) = 1 / [z - E_COHSEX,nk]
A_nk(E) = delta(E - E_COHSEX,nk)
```

The bubble is evaluated analytically from these poles. No spectral energy
grid, numerical delta broadening, fitted Z, or second subtraction of Vxc is
introduced. With fixed orbitals this static diagonal approximation gives
the same bubble as using the corrected energies with unit residues. It
does not reproduce dynamical GW satellites or lifetimes.

For each q, the response loops over the complete requested band-pair
range and BZ k points, using Yambo's `qindx_X` mapping k -> k-q and density
vertices from `scatter_Bamp`. At the optical q=0 point, the G=0 density
vertex instead uses Yambo's position dipoles and its numerical q0 norm;
the remaining G components use the native k-star rotation. The frequency
factor for a pair a,k and b,k-q is

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
ChiGMode="COHSEX":chi0=P_G0; P=P_G_COHSEX; fxc=chi0^-1-P_G_COHSEX^-1
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

Start with one CPU MPI rank, a small insulating SAVE database, a small
response band range and a small density basis. Include q index 1 for the
optical test.
The COHSEX route has no spectral convolution and needs only pole-pair sums.
The dynamical spectral convolution costs approximately
`Nk * Nb^2 * NE^2 * N_response_frequencies`, before the density-matrix
accumulation. Validate the G0 baseline against native `X_irredux` before
using the COHSEX kernel.

### Static screening: COHSEX route

Use this route when you have only a static W calculation. Retain SAVE and
the native `ndb.em1s` database and its fragments. The standalone `ndb.W`
export is not an input to the native COHSEX solver. If the native screening
databases are absent or incompatible, Yambo recalculates static screening.

Generate a COHSEX input with the Newton solver:

```sh
yambo -p c -g n -F cohsex.in -J cohsex
```

Keep your converged static screening settings and choose `QPkrange` to
cover **every response band, k point and spin**, with diagonal states only.
Keep the KS starting energies and orbitals unmodified; do not enable
external energy/width/residue corrections, self-consistent iterations or
mixed self-energies. Converge the native COHSEX band and screening cutoffs.
Run the edited input to produce `cohsex/ndb.QP`:

```sh
yambo -F cohsex.in -J cohsex
yambo -d s -F chi_cohsex.in -J chi_cohsex
```

In the generated static response input, set/add:

```text
ChiFxc
ChiGMode="G0"
ChiGAxis="IMAG"
ChiGNuSteps=16
% ChiGNuRnge
 0.0 | 20.0 | eV
%
ChiGCheckTol=0.0001
ChiRcondMin=1.0e-12
XTermKind="none"
% EhEngyXs
 0.0 | 0.0 | eV
%
```

Set `BndsRnXs` to the response band range, `NGsBlkXs` to the density basis,
and include q index 1 in `QpntsRXs` for an optical calculation. Set
`LongDrXs` to the desired optical field direction and later use the same
direction in BSE's `BLongDir`. Disable double grids and external
X energy, width, residue and Green-function corrections. Run once in G0
mode in a fresh `chi_g0` job to validate the KS baseline. Then change/add:

```sh
yambo -F chi_cohsex.in -J chi_g0
```

```text
ChiGMode="COHSEX"
ChiGDb="cohsex/ndb.QP"
```

```sh
yambo -F chi_cohsex.in -J chi_cohsex
```

The imaginary export grid can contain many frequencies even though native
screening is static. Its response matrix and native frequency grid are
preserved. The loader rejects a non-COHSEX QP database, nonunit residues,
linewidths, off-diagonal or duplicate states, incomplete state coverage,
different KS reference energies, and corrected poles that cross the fixed
KS chemical potential. Such crossings need a separate treatment of the
chemical potential and occupations. `ChiGNormTol` is used only for the
dynamic spectral route; static poles have exact unit weight.

For an optical Casida calculation, create a matrix BSE input from the same
`SAVE`, set `BSKmod="CHI"`, and set `BSENGfxc` to the exact G-basis size stored
in `ndb.Chi`. Point the BSE run at the response job database and use KS
transition energies (`KfnQPdb` off). Set `BSEQptR` to q index 1,
`QShiftOrder=1`, and `BLongDir` to the response `LongDrXs` direction. The
loader checks the direction,
normalization, G basis, and zero-frequency Hermiticity before projection.
The BSE path currently requires one CPU rank, scalar unpolarized states,
and a matrix solver. Its static Casida matrix uses only the first (zero)
frequency of the exported kernel. Start by comparing `G0` against the
ordinary KS calculation; then compare the `COHSEX` result. These optical
steps still require an end-to-end check with a real Yambo calculation. Use
a new BSE job when changing the response database or optical direction;
Chi kernel restarts are rejected because the old kernel cannot be
verified against the new inputs.

### MPI for the response

G0, COHSEX and DYSON support k-point distribution on CPUs. For N ranks,
set the response layout explicitly, for example with four ranks:

```text
X_and_IO_CPU="1 1 4 1 1"
X_and_IO_ROLEs="q g k c v"
X_and_IO_nCPU_LinAlg_INV=1
```

```sh
mpirun -np 4 yambo -F chi_cohsex.in -J chi_cohsex
```

Use N no larger than the number of BZ k points. Automatic layout selection
may distribute bands and is therefore not sufficient. Each rank computes
its owned k points; Yambo's native wavefunction partition includes their
k-q partners. The full bubble is summed in double precision across the
response communicator. The density matrices, G data and inversion work
remain replicated, and only the master writes `ndb.Chi`. q, G-vector and
conduction/valence-band distribution are rejected. GPU bubble kernels
remain unimplemented; use a CPU build for the response. Native COHSEX/GW
solver capabilities are unchanged by this response restriction.

### PPA dynamic screening route

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
  threshold, chemical potential, kBT]. `CHI_OPTICAL_Q0` stores the optical
  field direction with Yambo's numerical q0 norm for q=0 projection.

Older Chi files remain readable without the new optional metadata; absent
conditioning diagnostics read as zero and absent G0 diagnostics as -1.

Retarded `ndb.G` additionally stores `G_RETARDED`,
`G_CHEMICAL_POTENTIAL`, `G_REFERENCE_ENERGIES`, `G_SIGMA_X` and
`G_VXC_REFERENCE`. Existing `SE_Operator` stores DeltaSigma. Recover
Sigma_c as `SE_Operator - G_SIGMA_X + G_VXC_REFERENCE`, broadcasting the
static terms over frequencies. `Green_Functions_Energies` includes eta.

## Current scope and validation

The response G0 and COHSEX paths support optical q=0, serial CPUs and k-only
MPI. The DYSON optical limit remains unsupported. They reject GPU execution,
finite temperature, metallic KS references, double grids,
transition-energy filtering and response terminators. The retarded GW
export requires PPA, diagonal Sigma, unshifted KS starting energies and no
GW terminator, Green-function zoom, GreenF2QP, self-consistent GW or mixed
electron-phonon/photon self-energies. Off-diagonal Sigma and the GW response
vertex remain separate extensions. The optical G0/COHSEX path has synthetic
tests but has not yet been validated on a real material or a full Yambo build.

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
static COHSEX pole loading and inverse differences, preservation of retarded
sampling damping, and disjoint k partitions whose sum matches the serial
bubble (including an empty synthetic rank). Negative tests cover invalid
COHSEX databases, occupation crossings, missing states, legacy spectra,
reference mismatches, unsupported MPI layouts and ill-conditioned inversions.
They passed with GNU Fortran 16.2.0 and Reference LAPACK 3.12.1 in both host
precision settings.

A separate test uses real MPI collectives around the production bubble
with synthetic material/IO. Run it on a machine with a GNU-compatible MPI
Fortran compiler and launcher; no LAPACK source is needed for this test:

```sh
python tests/chi_green/run_tests.py --fc gfortran --mpifc mpifort --mpiexec mpiexec
python tests/chi_green/run_tests.py --fc gfortran --single-precision --mpifc mpifort --mpiexec mpiexec
```

The default runs use two and three ranks for a two-k-point fixture, testing
an empty rank as well. This test has not been executed in the Windows
validation environment, which has no MPI runtime. It does not exercise
Yambo's native communicator creation or wavefunction distribution.

These tests use mock material and database IO. They do not establish a
 complete Yambo build, NetCDF round-trip, real MPI material run, or a material
 benchmark. Those
checks require your Linux/HPC build and SAVE/screening data. Begin with the
G0 run, then converge spectral coverage/resolution, eta, bands, k mesh and
density basis before drawing physical conclusions from the Dyson kernel.
