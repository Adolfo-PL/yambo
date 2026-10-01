# Green-function bubble route to fxc

This change adds an initial diagonal, frozen-screening GW route to the
`ChiFxc` export. The existing level-correction route remains the default:
`ChiGMode="LEVELS"`. The G modes are `G0`, `COHSEX`, `QP` and `DYSON`.
`COHSEX` works with native static screening; `QP` takes the energies of any
diagonal GW run (G0W0 or GW0, PPA/MPA or real axis) as unit-weight poles;
`DYSON` uses PPA dynamic screening.

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

All energies here are in Yambo's internal KS frame, which OCCUPATIONS_Fermi
shifts so that the Fermi level (the VBM of an insulator) is zero. Yambo keeps
the unshifted Fermi energy in `E_Fermi`; `ndb.G` stores it as
`G_CHEMICAL_POTENTIAL` only to identify the reference. The loaders never use
it as a chemical potential: occupations are those of the T=0 KS reference, and
the chemical potential of the response is placed in the middle of the QP gap
(between the highest occupied and the lowest empty QP pole or peak).
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
ChiGMode="QP":    chi0=P_G0; P=P_G_QP;     fxc=chi0^-1-P_G_QP^-1
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
different KS reference energies, and corrected poles that close the gap (an
occupied pole above an empty one). Such crossings need a separate treatment of
the chemical potential and occupations. `ChiGNormTol` is used only for the
dynamic spectral route; static poles have exact unit weight.

### GW energies as unit-weight poles: QP route

`ChiGMode="QP"` builds the same pole bubble from the energies of any diagonal
GW `ndb.QP` written by the Newton solver: COHSEX, G0W0 or GW0, with PPA/MPA or
real-axis screening. Each state becomes

```
A_nk(E) = delta(E - Re E_QP,nk)
```

The linearized energy `E_QP = e_KS + Z <Sigma(e_KS) - Vxc>` already contains Z,
but the pole keeps unit weight. The loader therefore drops the residue Z and
Im E, and reports both in the report file (`Z dropped from QP poles`, `max |Im E|
dropped from QP poles`). A width larger than a quarter of the state's QP shift
is a lifetime rather than damping noise; it is dropped with a warning that
points to `DYSON`. Weighting the poles by Z alone would break the spectral sum
rule, so it is not offered. All other checks are those of the COHSEX route:
diagonal states, full band/k/spin coverage, the KS reference energies and an
open QP gap. A Green-function `ndb.G` is rejected;
use the Newton `ndb.QP`. `ChiGMode="COHSEX"` keeps its strict checks and still
rejects a non-COHSEX database.

For G0W0 energies, run PPA GW with the Newton solver and a `QPkrange` that
covers every response band and k point, then set:

```text
ChiGMode= "QP"
ChiGDb= "gw/ndb.QP"
```

The `ndb.Chi` header records `QP poles from G0W0 (PPA/MPA)` (or the actual
self-energy kind) as its QP action.

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
ordinary KS calculation; then compare the `COHSEX` result. With coupling,
the BSE bands and the G basis of the `ndb.Chi`, the static Casida response
reproduces the `ndb.Chi` response at zero frequency exactly (checked on
monolayer WS2). Use a new BSE job when changing the response database or
optical direction; Chi kernel restarts are rejected because the old kernel
cannot be verified against the new inputs.

### Frequency-dependent Casida: `BSEChiDyn`

A static kernel reproduces the zero-frequency response but cannot move the
absorption peaks: the Casida poles stay close to the KS transitions and only
their weights change. `BSEChiDyn` uses the whole real-axis kernel instead.
With `BSKmod="CHI"`, `K_kernel` keeps fxc out of the BSE matrix and stores the
transition vertices; `K_Chi_dynamic` then solves, at every BSE frequency
z = w + i eta,

```
Resp(w)  = -Co B^T [ z - H0 - Kf(z) ]^-1 A
Kf_xy(z) = spin_occ/(Nk V) diag(sqrt f) V_x^H fxc(z) V_y diag(sqrt f)
```

`H0` is the ordinary BSE matrix (KS transition energies plus exchange), `V_x`
the vertices of the static projection (excitation vertices, and de-excitation
vertices for the anti-resonant pairs), `A` and `B` the residual vectors of the
diagonalization solver. All blocks see the same retarded fxc(z). Since
fxc = chi0^-1 - P^-1, the Woodbury identity makes the coupled response equal
to P plus local fields whenever the BSE transitions build the same chi0 as the
`ndb.Chi`: same bands, k grid, G basis, frequencies and damping. For a
z-independent Hermitian fxc it equals the diagonalization result.

```text
BSKmod= "CHI"
BSEChiDyn
BSEmod= "coupling"
BSENGexx= <G-basis size of ndb.Chi>  RL
BSENGfxc= <G-basis size of ndb.Chi>  RL
% BEnRange
 <first Re w> | <last Re w> | eV       # the ndb.Chi grid
%
BEnSteps= <number of ndb.Chi frequencies>
% BDmRange
 <eta> | <eta> | eV
%
BSEprop= "abs"
```

The `ndb.Chi` must be a real-axis, retarded export (`ChiGAxis="REAL"` in the G
modes, which force retarded ordering; `GrFnTpXd="R"` in LEVELS mode). fxc is
not interpolated in frequency: the BSE grid must be the `ndb.Chi` grid, and the
error message prints the `BEnRange`, `BEnSteps` and `BDmRange` to use. The
solver replaces the requested `BSSmod` and writes `o-<job>.eps_q1_chidyn_bse`
and `o-<job>.eel_q1_chidyn_bse`. Each frequency costs one dense LU solve of
size N (TDA, twice) or 2N (coupling); frequencies are distributed over OpenMP
threads, each holding about three (2N)^2 complex(8) matrices with coupling.
The restrictions of the static path apply (one MPI rank, scalar unpolarized
states), plus the length gauge and absorption only.

Three conditions matter in practice:

- **Use coupling.** fxc(z) inverts the full resonant plus anti-resonant
  bubbles. In TDA their mixing is dropped and the KS poles are not cancelled.
  In silicon, TDA makes the main peak 1.9 times too high and eps1(0) 13%
  too large.
  The solver warns when TDA is used.
- **The BSE band window must not split a degenerate multiplet.** Yambo warns
  `User bands ... break level degeneracy`. A split multiplet makes the
  transition basis depend on the wavefunction gauge. The zero-frequency
  response hardly notices, but near resonances the dynamic identity fails. In
  silicon, bands 3-6 (split Gamma multiplets) give up to 25% error at the
  peaks, whereas bands 1-8 agree to 3e-4.
- **More bands in `ndb.Chi` than in the BSE break the identity.** fxc(z) is the
  inverse difference of the full bubbles, while the BSE window contains only
  part of chi0. The exact reference is an `ndb.Chi` built with the BSE bands;
  a wider `ndb.Chi` gives an approximation whose error has to be converged.

Validation on bulk silicon (LDA Troullier-Martins pseudopotential, 4x4x4
Gamma-centred grid with 64 BZ points, symmorphic operations, 15 G vectors,
LEVELS mode with a 1 eV scissor on P, retarded ordering, 81 real frequencies
0-8 eV, eta = 0.1 eV, bands 1-8):

| Run | Reference | Max. relative difference |
| --- | --- | --- |
| Hartree BSE, coupling | KS RPA+LF of `ndb.Chi` | 5.6e-4 |
| `BSEChiDyn`, coupling | QP RPA+LF of `ndb.Chi` | 3.0e-4 |
| IP column of the same run | KS IP of `ndb.Chi` | 1.1e-4 |

The dynamic kernel moves the main eps2 peak from the KS position (4.0 eV) to
the QP position (5.0 eV) with the QP height; the static kernel cannot.

### Kernel from the BSE and exciton binding energies: `BSEChiOut`

A P built from G0W0 Green functions is a QP bubble: it has no electron-hole
attraction, so Casida with its fxc opens the gap but binds no exciton (in
monolayer WS2 the static onset stays at the KS gap). The attraction is the W
term of the BSE. `BSEChiOut` takes P from a coupled BSE instead. With the BSE
matrix M (transition energies, exchange vbar, and W when `BSKmod="SEX"`) and
the Chi transition vertices collected on the fly, `K_Chi_export` computes at
every frequency

```
chib(z) = Co/(4 pi) B^T [ z - M ]^-1 A        (G-resolved BSE response)
P^-1    = chib^-1 + vbar
"exc" : fxc = Pref^-1  - P^-1,  Pref = bubble of the BSE transition energies
"full": fxc = chi0^-1  - P^-1,  chi0 = KS bubble (Yambo's bare energies)
```

and writes `ndb.Chi` (chi0 column = Pref or chi0, P, fxc) into the BSE job
directory. By the Woodbury identity, Casida with the same reference energies,
the same vbar and this fxc(z) reproduces the BSE response exactly, its exciton
poles included. The binding energy of the frequency-dependent Casida is
therefore the BSE binding energy. The first frequency of the file is z = 0
exactly (static Casida); the remaining ones are the BSE grid (`BSEChiDyn`).

```text
BSKmod= "SEX"
BSEmod= "coupling"
BSSmod= "d"
BSEChiOut= "exc"                 # or "full"
BSENGexx= <n> RL                 # vbar and the Chi basis must coincide
BSENGfxc= <n> RL
BSENGBlk= <W block>
KfnQPdb= "E < <gw job>/ndb.QP"   # QP transition energies
```

Casida then uses a copy of the exported `ndb.Chi` (without the BSE kernel
database, which would be taken as a restart) with `BSKmod="CHI"`, coupling,
the same bands, `BSENGexx` and `BSENGfxc`, and:

- "exc": the same `KfnQPdb` (the kernel is measured from the QP bubble; the
  run stops if the QP corrections are missing);
- "full": KS energies (the kernel carries the QP shift too; the run stops if QP
  corrections are applied).

Add `BSEChiDyn` and the BSE grid for the dynamic solver; without it, the static
diagonalization uses fxc(0).

Every diagonalization now reports its lowest excitations, their relative
optical strength |R_left R_right|, the lowest bright one (strength above 1% of
the maximum), the lowest transition of the window and the binding energy
`gap - lowest bright excitation` (in the `r-` file, `[BSE]` lines). With QP
transitions the gap is the QP gap on the k grid. For a "full" kernel the
transitions are KS ones; take the gap from the exporting BSE.

Restrictions: optical q, one MPI rank (OpenMP threads parallelize the
frequencies), coupling with the diagonalization solver, length gauge, scalar
unpolarized states, no transition widths or Z factors, `BSENGexx = BSENGfxc`.
The export stops if the reference bubble or the BSE response is
ill-conditioned (`[Chi/BSK] minimum rcond` below `ChiRcondMin`); reduce the
G basis then. Monolayer WS2 (bands 11-16, 9 G, no Coulomb cutoff) passes with
rcond 2e-9: the eight pure-G_z plane waves are nearly dependent, the eigenvalues of
D chi0 D span 1e-8 to 2.4, and the body of fxc reaches 3e6 in units of v. The head
and the optical response are exact (the exported P gives the BSE eps to 3e-7, the
8-digit precision of the text output), but for the kernel itself a smaller basis is
cleaner. `plot_summary.py --chi` prints this diagnostic. The same run with one plane
wave (`BSENGexx = BSENGfxc = 1`) gives the same BSE (2.2342 eV, binding 0.309 eV), the
same spectra (5e-7) and the same head alpha(w) (4e-6), with rcond = 1. Its 13->14 fxc
block is then exactly rank one, F_tt' = -c alpha(z) d_t conj(d_t') with d_t the optical
dipoles: for in-plane optics the exact kernel is a frequency-dependent LRC kernel. Each frequency costs one LU solve of size 2N. The Coulomb factor
is D = sqrt(4 pi)/bare_qpg as a complex number, the form the BSE exchange
(1/bare_qpg**2) uses, so the identity also holds with a Coulomb cutoff where
v_cut(q+G) < 0 for some G.

`BSEChiDyn` also compares the bubble of the BSE transitions with the `ndb.Chi`
chi0 at its first frequency (`[Chi/BSK] BSE bubble vs ndb.Chi chi0`, in the
Coulomb-symmetrized basis) and warns above 1e-3: the kernel is then an
approximation, not the exact one (bands, energies or k grid differ). Silicon:
3e-8 with an exported kernel, 4e-6 with a native `ndb.Chi` of the same bands,
0.15 (warning) for bands 3-6 against a bands 1-8 `ndb.Chi`.

### The Casida matrix as text: `BSEChiKoutW`, `BSEChiKoutB`

With `BSEChiDyn`, the matrix the dynamic solver inverts can be written at the
`ndb.Chi` frequencies whose real part lies in a window:

```text
% BSEChiKoutW
 0.0 | 0.0 | eV        # static point z=0 (and z=0+i eta) of an exported ndb.Chi
%
% BSEChiKoutB
 13 | 14 |             # one valence | conduction pair; 0 | 0 = all transitions
%
```

`o-<job>.Ktt_q1_transitions` lists each transition t (index T in the BSE, BZ and
IBZ k, k in reduced units, v, c, energy E_t in eV, occupation factor f_t and
residual a_t = d_t sqrt(f_t)). `o-<job>.Ktt_q1_w<iw>` holds, per pair t t',
the blocks in eV

```
M(z) = [[ diag(E) + X_rr + F_rr ,           X_rc + F_rc ],
        [           X_cr + F_cr , -diag(E) + X_cc + F_cc ]]
```

X is the static part of the BSE matrix without the energies (exchange), F the
fxc(z) part as it enters M (F_rr = Kf_rr, F_rc = i Kf_rc, F_cr = i Kf_cr,
F_cc = -Kf_cc). The response is Resp(z) = -Co B^T [z - M(z)]^-1 A with
A = (a, i b), B = (b, i a), b = conj(a). In TDA only the rr columns are
written. With a band pair, the file holds those rows and columns of the full
matrix, not a separate Casida problem. Silicon check: the eigenvalues of M(0)
rebuilt from the files match the static Casida run to 4e-5 eV (single-precision
matrix, 8 printed digits); the response at z = 0.1i eV matches the dynamic
solver to 1e-6.

`sbin/chi_tools/ktt_casida.py` reads these files, builds M(z), and prints its
lowest eigenvalues with their oscillator strengths, the lowest bright one and
its binding energy; `--tda`, `--no-fxc` and `--spectrum` are optional:

```sh
python sbin/chi_tools/ktt_casida.py o-<job>.Ktt_q1_transitions o-<job>.Ktt_q1_w<iw>
```

The kernel is exact only on the whole transition space. A single band pair
loses the binding: silicon at z = 3.4 + 0.1i eV gives 3.452 eV with all
transitions (BSE exciton 3.428 eV) but 3.589 eV with the 4->5 pair alone.

A single matrix is not a binding energy: with the exact kernel the exciton
solves lambda(w) = w, where lambda is the lowest bright eigenvalue of M(w). A
kernel exported with `BDmRange 0 | 0` sits on the undamped real axis (accepted by
`BSEChiDyn` for BSE-exported kernels only). Silicon, full matrix on 3.00-3.50 eV:
lambda falls from 3.486 eV and crosses w at 3.4276 eV (BSE 3.4278 eV); at w = 0
it is 3.558 eV. A damped z = w + i eta, or one band pair, misses the root.

`sbin/chi_tools/sc_scan.py` automates this. It diagonalizes each w-file,
records lambda(w), and interpolates the crossing. It flags sign changes
where lambda jumps across a pole of fxc(w), because those are not solutions.
`--delete` removes each file once it is read, and `--append` merges the
batches of a scan into one table:

```sh
python sbin/chi_tools/sc_scan.py o-<job>.Ktt_q1_transitions o-<job>.Ktt_q1_w* \
       --out sc_scan.dat --exciton <BSE energy> [--delete] [--append]
```

Silicon, 7 points on 3.30-3.48 eV in three batches: the crossing is at
3.4272 eV. The finer 0.02 eV grid gives 3.4276 eV. The BSE value is 3.4278 eV.

`sbin/chi_tools/plot_summary.py` draws the figures of the chain. Each figure
is made only when its inputs are given:

| Figure | Contents |
| --- | --- |
| `energies.png` | GW corrections of the states near the gap; the vb->cb gap at every k, KS and QP, with the exciton and static-Casida levels |
| `spectra.png` | eps2 for the QP transitions, the static-fxc Casida and the fxc(w) Casida |
| `fxc_omega.png` | ndb.Chi kernel in units of v: head alpha(w), one diagonal body element per shell, and the matrix at w = 0 |
| `fkq.png` | F part of K_tt' as a matrix; its diagonal next to the exchange X_tt; the diagonal over the zone (`--hex` for hexagonal lattices) |
| `selfconsistency.png` | lambda(w) and the line lambda = w |

```sh
python sbin/chi_tools/plot_summary.py --qp o-g0w0.qp --vb 13 \
   --eps "QP transitions=o-cas_st_exc.eps_q1_diago_bse:ip" \
   --eps "static fxc(0)=o-cas_st_exc.eps_q1_diago_bse" \
   --eps "fxc(w)=o-cas_dyn_exc.eps_q1_chidyn_bse" \
   --chi kx_exc --ktt o-cas_dyn_exc.Ktt_q1_transitions o-cas_dyn_exc.Ktt_q1_w225 \
   --scan sc_scan.dat --bse-report r-bse_sex_exc_* --static-report r-cas_st_exc_*
```

The `fxc_omega.png` figure needs netCDF4 (the script reads `ndb.Chi` itself).

Silicon, same setup as above (bands 1-8, 15 G, 1 eV scissor through
`KfnQP_E`), static screening with 12 bands, `BSKmod="SEX"`, `BSENGBlk= 15 RL`:

| Run | Lowest bright exciton | Binding energy | Spectrum vs BSE |
| --- | --- | --- | --- |
| SEX BSE, coupling (reference) | 3.428 eV | 0.173 eV | - |
| exported P, QP RPA+LF (`chi_spectrum.py`) | - | - | 1.4e-4 |
| Casida, "exc", fxc(z), `BSEChiDyn` | peak 3.40 eV | = BSE | 1.4e-4 |
| Casida, "full", fxc(z), `BSEChiDyn` | peak 3.40 eV | = BSE | 1.4e-4 |
| Casida, "exc", static fxc(0), QP energies | 3.558 eV | 0.043 eV | 0.71 |
| Casida, "full", static fxc(0), KS energies | 2.611 eV | none (onset at KS gap) | 0.92 |

The window gap is 3.601 eV (2.601 eV KS gap plus the scissor). Both static
kernels reproduce eps1(0) (19.84) but not the exciton: the frequency
dependence of the exact kernel carries the binding. The spectra were compared
on the 0.1 eV grid with eta = 0.1 eV. In this build, gfortran 13 at -O3
miscompiled `K_correlation_collisions_std` (the SEX kernel crashed without
the export and gave NaN with it); the file compiled at -O1 runs correctly.

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
state table, KS starting energies, energy reference, an open QP gap, positive damping,
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

The response G0, COHSEX and QP paths support optical q=0, serial CPUs and
k-only MPI. The DYSON optical limit remains unsupported. They force retarded
ordering on the response, so they stop in a job that also runs GW or a BSE:
build `ndb.Chi` in a screening-only job. They reject GPU execution,
finite temperature, metallic KS references, double grids,
transition-energy filtering and response terminators. The retarded GW
export requires PPA, diagonal Sigma, unshifted KS starting energies and no
GW terminator, Green-function zoom, GreenF2QP, self-consistent GW or mixed
electron-phonon/photon self-energies. Off-diagonal Sigma and the GW response
vertex remain separate extensions. The optical COHSEX and QP kernels and the
static Casida projection have been run on monolayer WS2; the frequency-dependent
Casida solver and the BSE kernel export have been validated end to end on bulk
silicon (see above).

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

The frequency-dependent Casida numerics have their own test, which compiles
the production module against LAPACK:

```sh
python tests/chi_dynamic/run_tests.py --fc gfortran
python tests/chi_dynamic/run_tests.py --fc gfortran --openmp
```

It checks the exact pair-space identity against a G-space Dyson inversion for
arbitrary non-Hermitian fxc(z), with coupling and in TDA; the static limit
against the matrix assembled as in `K_stored_in_a_big_matrix` and the
Lorentzian sums of `K_diago_response_functions`; error codes; and a threaded
frequency loop against the serial one. For `BSEChiOut` it checks that the pair
response of diag(E,-E) is the pair bubble and that Casida with the exported
"exc" and "full" kernels reproduces the response of a random coupled BSE
(QP energies, vbar and a random W part) at real and complex frequencies.

These tests use mock material and database IO. They do not establish a
 complete Yambo build, NetCDF round-trip, real MPI material run, or a material
 benchmark. Those
checks require your Linux/HPC build and SAVE/screening data. Begin with the
G0 run, then converge spectral coverage/resolution, eta, bands, k mesh and
density basis before drawing physical conclusions from the Dyson kernel.
