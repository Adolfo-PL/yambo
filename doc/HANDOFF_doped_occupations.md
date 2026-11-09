# Handoff: occupations of the doped BSE (WS2, jellium doping)

From the scripts/results chat to the technical chat. One question to settle before the doped
exciton energies go into the abstract; the exact-kernel checks do not depend on it.

## Setup

- WS2 monolayer, no SOC, PBE; doping by `tot_charge` (jellium background), N-X = X extra electrons per
  cell: 0.01, 0.05, 0.1 e/cell (1.1e13, 5.7e13, 1.1e14 cm^-2). 0.2 e/cell was dropped: a vacuum
  (nearly-free-electron) state falls below the conduction band at Gamma (jellium artifact, Brumme et al.
  PRB 89, 245406); 0.01-0.1 were checked free of it.
- G0W0 (PPA, RIM-W, slab cutoff, 400 bands + BG) on 12x12, QP bands 13-15 at all 19 k, ElecTemp 300 K.
- BSE on 24x24 (`exact_kernel.sh`, doped mode): `KfnQPdb= "E < gwqp/ndb.QP ubi"`, bands 13-15, TDA
  reference + coupled export at q1, ElecTemp 300 K. The static screening (em1s) and the GW use the DFT
  occupations; the BSE gets the QP energies through `KfnQPdb`.
- K = k19 of the 12x12 grid; band 13 = VB, band 14 = CB (partly filled).

## Data

G0W0 at K (KS energies from the DFT Fermi level, o-gw.qp):

| doping | VB at K (QP) | CB at K (QP) | QP gap at K | CB at K, KS, from E_F | KS occupation of CB at K (300 K) |
|---|---|---|---|---|---|
| 0.01 | -2.2823 | +0.3902 | 2.6725 | +0.016 eV | ~0.35 |
| 0.05 | -2.3089 | -0.0632 | 2.2457 | -0.072 eV | ~0.94 |
| 0.1  | -2.3026 | -0.1167 | 2.1859 | -0.092 eV | ~0.97 |

Reference BSE at q1 (fork report: interband gap = lowest transition out of a full band present in
the BSE; A = lowest bright interband excitation):

| doping | interband gap | A | strength of A | gap - A | spectrum peak |
|---|---|---|---|---|---|
| 0     | 2.9238 | 2.4021 | -    | 0.522 | - |
| 0.01  | 2.6720 | 2.5009 | 0.81 | 0.171 | 2.502 |
| 0.05  | 2.6445 | 2.4353 | 0.28 | 0.209 | 2.436 |
| 0.1   | 2.1847 | 2.1590 | 1.00 | 0.026 | 2.161 |

The contradiction: by the KS occupations, K is more blocked at 0.1 (CB ~97% full) than at 0.05 (~94%).
Yet the BSE keeps the K transition at 0.1 (interband gap 2.1847 = K gap 2.1859) and drops it at 0.05
(interband gap 0.40 eV above the K gap). A follows: 2.40 -> 2.50 -> 2.44 -> 2.16 eV, and at 0.05 A lies
0.19 eV above the K gap.

Report lines (r-bse_ref): the BSE refits the Fermi level with the QP energies.

| doping | [X] Fermi level | [K+QP] Fermi level | [X] avg metallic occ. | [K+QP] avg metallic occ. | [K+QP] N of met el |
|---|---|---|---|---|---|
| 0.01 | 1.476341 | 0.439227 | 0.007694 | 0.008829 | 0.03531 |
| 0.05 | 2.106190 | 0.234043 | 0.018522 | 0.019746 | 0.07898 |
| 0.1  | 2.780561 | 0.091744 | 0.029977 | 0.030604 | 0.12241 |

Code path, as far as the scripts chat can tell: `src/qp_control/QP_apply.F` (about lines 174-176) calls
`OCCUPATIONS_Fermi(en,k,"K+QP")` after the QP corrections, unless `en%fo` is allocated or the run is
self-consistent. So the BSE occupations are a Fermi-Dirac fill of the interpolated QP bands at 300 K,
not the DFT ones. HANDOFF_code.md already notes "QP corrections recompute the BSE occupations,
QP_apply.F:175-176" as untested with the "exc" kernel.

## Questions

1. Is the refill the cause of the 0.05 anomaly, i.e. does the QP shift move the doping electrons
   between valleys (K, Q) differently at 0.05 and 0.1? The per-transition list that settles it
   (13->14 transitions near K and their f_t, from `prod/o-xs_q1.Ktt_q1_transitions`) is being produced
   by the current runs; the user will forward it.
2. Which occupations should the doped BSE use?
   - (a) DFT occupations (consistent with the GW and the screening, which were computed with them), or
   - (b) the QP refill (electrons in the lowest QP valley, but inconsistent with W and Sigma).
   Is there a literature convention for G0W0+BSE of doped TMDs (e.g. Liang-Yang, Druppel et al.)?
3. Can Yambo keep the DFT occupations while applying QP energies (input flag, or a small fork change
   around `QP_apply.F` / `en%fo`)? If yes, what is the input, and does `K_Transition_check.F` then keep
   the K transitions with their DFT f_t?
4. Does the 24x24 grid resolve this at all? At 0.01-0.1 e/cell the Fermi sea covers about 1-4 k-points
   per valley, so one k-point moving across E_F changes the onset by tenths of an eV. Is a denser BSE
   grid (or a smaller ElecTemp) needed before comparing dopings?

## What stays valid either way

The exact-kernel checks (Casida with fxc(q1, w) reproduces the export BSE's A to < 1 meV) hold for
whatever occupations the BSE uses. The doped mode of `exact_kernel.sh` takes the interband gap and A
from the report lines added on Oct 20 (`K_diago_driver.F`, binding_report), so it follows the BSE's own
transitions. Only the physical reading of A and Eb vs doping waits on these answers.

---

## Answers (technical chat)

Status words used below. **Code**: read in the source, line given. **Tested**: run in this session on a
small test (silicon, 6x6x6, 30 bands, undoped and with 0.06 extra electrons; not WS2). **Inferred**: my
reasoning from the numbers in this file; needs the cluster check of section A4. **Not verified**: I could
not check it.

Short version:

- The refill exists and works as the scripts chat described (code). With it, K is almost completely
  Pauli-blocked at both 0.05 and 0.1; whether the K transition is *in* or *out* of the BSE is decided by
  a fixed cutoff of 1e-5 on the hole occupation. By the numbers above the hole at K is about 1e-5 at
  0.05 (on the cutoff) and about 3e-4 at 0.1 (30 times above it). That is the 0.05 / 0.1 asymmetry
  (inferred; one cheap run confirms it).
- There is a second, larger problem underneath: the GW correction of a conduction state at K depends on
  that state's own occupation (about +0.59 eV if empty, about -0.03 eV if full). The refill amplifies
  it; keeping the DFT occupations does not remove it (inferred).
- New input flag `KfnQP_KeepOcc`: QP energies with the reference (KS) occupations. Off by default; with it
  off the results are identical to the previous code (tested).
- The 24x24 BSE grid and, more seriously, the 12x12 GW grid do not resolve these Fermi seas. Lowering
  ElecTemp would make it worse, not better.

### A1. Mechanism, through the code

Line numbers are those of this branch after the change (`QP_apply.F` moved by a few lines).

1. **Reference occupations.** `src/setup/setup.F:378`, `OCCUPATIONS_Fermi(Xen,Xk,"X")`: Yambo fills the KS
   energies itself, with a Fermi-Dirac function at `ElecTemp`, on the k-grid of the run. These are the
   `[X]` lines of every report. They are not read from QE: "DFT occupations" in Yambo means "Yambo's
   Fermi-Dirac filling of the KS energies". `setup.F:379` passes them to the GW levels (`en`) and
   `setup.F:403` copies them to the BSE levels: `E_duplicate(Xen,Ken)`. So the screening, the GW and, at
   this point, the BSE all hold the same occupations.
2. **QP energies.** `src/setup/External_corrections.F:79-80`: for a BSE run, `QP_apply(BS_bands,Ken,Xk,"K")`.
   With `KfnQPdb`, `QP_apply_interpolate_QP.F:205` sets `en%E = en%Eo + (interpolated correction)` for the
   bands that are both in the database and in `BSEBands` (here 13-15). All other bands keep their KS
   energies.
3. **The refill.** `src/qp_control/QP_apply.F:183-185`: unless `en%fo` is allocated (real-time runs) or
   the run is self-consistent, `OCCUPATIONS_Fermi(en,k,"K+QP")`. That routine
   (`src/common/OCCUPATIONS_Fermi.F`) sorts the levels of all bands (lines 62-79), finds the Fermi level
   that holds `nel` electrons at the temperature (131-236), overwrites the occupations with
   `spin_occ * Fermi_fnc(E - E_F, T)` (289-295) and shifts the energies so that the new Fermi level is
   zero (319). It then prints the `[K+QP]` block. The KS occupations of the BSE levels are gone.
4. **The cutoff.** `filled_tresh = 1e-5` (`src/modules/SET_defaults.F:312`). `Fermi_fnc`
   (`src/modules/mod_functions.F:120-131`) sets any value below it to exactly zero.
5. **Transitions.** `src/bse/K_Transition_check.F`: `f_eh = (f_v - f_c)/spin_occ` (line 64); the electron
   and hole room `f_e = f_v/spin_occ`, `f_h = (spin_occ - f_c)/spin_occ` (85-86); the transition is
   dropped if `f_e <= filled_tresh` or `f_h <= filled_tresh` (137). A conduction state that is full to
   within 1e-5 loses all transitions into it: this is how Pauli blocking enters, as a yes/no decision.
6. **Kernel.** The kept transitions carry `f_eh` (`K_Transitions_setup.F:180`); the Hamiltonian is
   H = D + sqrt(f) K sqrt(f) with `sqrt_f` at `K_kernel.F:1000-1001`. A transition with a small `f_eh`
   stays in the matrix but barely couples to the others, and its optical weight is proportional to `f_eh`.

Two side facts (code): the refill also changes the temperature. `OCCUPATIONS_Fermi.F:201` and `:232` raise
`Tel` by 0.1 K when the search stalls, which is why a `[K+QP]` block can print 300.2 K after `[X]` printed
300.1 K (seen in the test), and why job chains differ by 0.1 K. And "N of met el" in the report is not an
electron count: `src/communicate/REPORT_Occupations.F:119-120` sums the occupations over the IBZ k-points
without weights. The values 0.035 / 0.079 / 0.122 in the table above are therefore not expected to equal
the dopings.

**Is the refill what removes K at 0.05 and keeps it at 0.1?** With the numbers of this file, yes within
their accuracy (inferred). `[K+QP] Fermi level` is the shift of the Fermi level from the KS one, and the QP
energies in `o-gw.qp` are also measured from the KS Fermi level, so the conduction band at K lies

| doping | E_F(K+QP) - E_CB(K), QP | hole left at K, 1/(1+exp(that/kT)) | against the cutoff 1e-5 |
|---|---|---|---|
| 0.01 | 0.439 - 0.390 = 0.049 eV | 0.13 | kept, f_eh = 0.13 |
| 0.05 | 0.234 + 0.063 = 0.297 eV | 1.0e-5 | on the cutoff: dropped |
| 0.1  | 0.092 + 0.117 = 0.208 eV | 3.2e-4 | kept, f_eh = 3e-4 |

(kT = 0.025852 eV.) So under the refill K is more than 99.9 % blocked at both dopings. At 0.1 the
transition survives only as an almost empty entry. That fits the report: at 0.1 "A" sits 26 meV below the
K gap with relative strength 1.00, which is what a nearly uncoupled K transition looks like, not a bound
exciton; at 0.05 the same transition is removed and the lowest remaining one is 0.40 eV higher. The
non-monotonic A (2.44 then 2.16 eV) is then a comparison of two different objects.

Caveats: 0.297 eV gives a hole of 1.02e-5, and 1 meV changes it by 4 %, so the table cannot decide 0.05 by
itself. The energies also mix two grids (the Fermi shift is from the 24x24 run, the QP energies from the
12x12 one). `STAGE=occ` (A4) prints the actual occupations and settles it.

**Why the refill moves the Fermi level so far above K (inferred).** Take the three GW corrections of the
conduction band at K from the first table, E_QP - E_KS = +0.374, +0.009, -0.025 eV, against the KS
occupation of that state, 0.35, 0.94, 0.97. They lie on one line, correction = 0.59 - 0.62 f eV (two
points fix the line; the third is off by 14 meV). The correction of the state follows the state's own
occupation: an empty state gets about +0.59 eV, a full one about -0.03 eV.

In plain terms: in GW an occupied state is lowered by its exchange with the occupied states, itself
included. On a k-grid the "itself" part is the q -> 0 cell of the screened interaction times the state's
occupation. On the 12x12 grid the whole Fermi sea of a valley is about one k-point (0.4 to 4 states,
A5), so this single-point term carries the whole effect. My estimate of its size, bare e^2 q0 = 1.4 eV for
the radius q0 = 0.10 1/A of one 12x12 cell, divided by a 2D screening of 2 to 3 at that q, is 0.5 to
0.7 eV, the same size as the fitted 0.62 eV. In a converged calculation the states just outside the Fermi
sea are lowered almost as much, and the band stays smooth through k_F. Here the neighbours on the grid get
nothing, so the QP band has a dip of about 0.6 eV at the occupied points. A lowering of the occupied
pocket by several hundred meV is real physics (band-gap renormalization, A2); the sharp step to the
neighbours is the grid.

The refill then fills this dip: the Fermi level ends up somewhere inside the step, 0.2 to 0.3 eV above K,
the states in the dip become full to 1e-5 or better, and the 300 K smearing of the KS filling is lost.
Where exactly the level lands depends on how many 24x24 points the interpolated dip holds and on the Q
valleys, so 0.297 against 0.208 eV is not a trend in the doping.

To confirm, send for each doping: `prod/r-occ_*` from `STAGE=occ` (A4), the 12x12 `prod_gw/o-gw.qp`, and
the existing `prod/r-bse_ref_*`. With `o-gw.qp` I can plot the correction of band 14 against the KS
occupation at all 19 k-points; if they fall on the line above, the step is established.

### A2. Which occupations should a doped G0W0+BSE use?

**Literature.** I could verify these three papers only at the level of their abstracts and publisher
records: from this session the full texts (arXiv, APS, ACS, eScholarship) are blocked. They exist and say
the following.

- Y. Liang and L. Yang, "Carrier Plasmon Induced Nonlinear Band Gap Renormalization in Two-Dimensional
  Semiconductors", Phys. Rev. Lett. 114, 063001 (2015), arXiv:1408.5815. A plasmon-pole theory for the
  coupling of the carrier plasmon to the quasiparticles; gap renormalization above 400 meV in doped
  MoS2 and MoSe2.
- S. Gao, Y. Liang, C. D. Spataru and L. Yang, "Dynamical Excitonic Effects in Doped Two-Dimensional
  Semiconductors", Nano Lett. 16, 5568 (2016). First-principles BSE for doped MoS2 with dynamical screening
  from a generalized plasmon-pole model. At moderate doping the binding energy changes by a few hundred
  meV while the exciton peak stays nearly constant, because the gap renormalization cancels it; at higher
  density the peak rises linearly and merges into a Fermi-edge singularity.
- A. Champagne, J. B. Haber, S. Pokawanvit, D. Y. Qiu, S. Biswas, H. A. Atwater, F. H. da Jornada and
  J. B. Neaton, "Quasiparticle and Optical Properties of Carrier-Doped Monolayer MoTe2 from First
  Principles", Nano Lett. (2023), arXiv:2303.12122. Gap renormalization of several hundred meV, a similar
  drop of the binding energy, a nearly unchanged lowest exciton energy, and a plasmon-pole model with
  separate poles for the intrinsic and the carrier plasmon.

**Not verified:** how any of them sets the occupations in the BSE, their k-sampling, and whether the
carriers are in the DFT or added to the screening only. I therefore cannot tell you "the convention is
(a)" or "(b)" from the literature. What to look for in the methods and supplements: the words Fermi
level, occupation, Pauli blocking; whether the band structure is the doped or the undoped one; the
k-grid or patch near K. If you send the PDFs I will read those parts.

What the abstracts do establish is a check on our numbers: the exciton peak should move little at moderate
doping. A that moves by 0.3 eV between 0.05 and 0.1 is not that behaviour. They also all treat the carrier
screening with a dedicated dynamical model, where our runs use the standard single plasmon pole and static
W; the Yambo report itself warns "The system is a metal but Drude term not included".

**My answer, from the formalism (reasoned, not a citation).**

- In the exact theory the BSE is built on the QP Green's function, so its occupations are those of the QP
  levels around the QP Fermi level. In that sense a refill is the right idea: option (b).
- But then W and Sigma must be computed with the same occupations. A one-shot G0W0 computes them with
  the KS filling. If the refill changes the filling, the BSE mixes QP energies that were computed for one
  set of occupied states with another set. Option (a), the KS filling everywhere, is the consistent
  one-shot choice.
- The two agree whenever the QP correction keeps the order of the states near the Fermi level. The number
  of electrons is fixed, so if the order is kept the same states are filled (the Fermi sea keeps its
  size; the thermal tails rescale a little). They differ only when the correction reorders states: between
  valleys (K against Q), which can be physical, or through the step of A1, which is not.
- In our runs the difference comes from the step: the states are low because they were occupied in the
  GW, and the refill then fills them completely. This is counting the same effect twice, and it replaces
  a 300 K filling by an almost sharp one.

So for these runs I recommend (a), `KfnQP_KeepOcc`. Be exact about what it does and does not do: it makes
the BSE consistent with W and Sigma and removes the cutoff accident; it does not repair the QP energies. The
K gap of the first table (2.67, 2.25, 2.19 eV) already contains the occupation-dependent shift, so the
doped energies stay unpublishable until the GW side is checked (A5).

### A3. The fix: `KfnQP_KeepOcc`

A flag in the BSE input (no value). With it, `QP_apply` applies the QP energies and skips the refill, so
`K_Transition_check` sees the KS occupations of step 1. It uses the switch Yambo already has for this
(the condition on `en%fo` at `QP_apply.F:183`), extended by one logical.

- `src/qp_control/QP_apply.F:175-183`: reads the flag (BSE levels only), skips `OCCUPATIONS_Fermi`.
- `src/qp_control/QP_apply_keep_occupations.F` (new): does the refill on a copy, restores the
  temperatures, and writes to the report what the refill would have changed. The BSE levels are not
  touched.
- `src/interface/INIT_QP_ctl_load.F`: declares the flag. It is not activated, so generated inputs do not
  contain it.

Report lines with the flag (from the silicon test):

```
[K+QP] KfnQP_KeepOcc: occupations of the reference levels kept, no refill of the QP levels
[K+QP] Fermi level of a QP refill, from the reference one [eV]:  0.590436
[K+QP] Electrons per cell a QP refill would move:  0.001753
[K+QP] Largest occupation change of a QP refill (f/spin_occ): -0.031793
[K+QP] ... at band, IBZ k-point:   5  13
[K+QP] States a QP refill changes by more than 0.001:  3
[K+QP] band, IBZ k, spin | E ref, E QP [eV] | f kept, f QP refill (f/spin_occ) | weight | full: kept, refill
[K+QP]      5    10     1    0.00778    0.59491  0.425390  0.456844  0.027778   n   n
...
[K+QP] Of these states, full with the kept occupations, with a QP refill: 0 0
```

"full: y" marks a state with no room left by the 1e-5 cutoff: the BSE drops every transition into it.
"Electrons a QP refill would move" is half the weighted sum of |f_refill - f_kept|. A warning is printed
when any state changes by more than 0.001.

With the flag the energies stay measured from the KS Fermi level (the refill would have shifted them);
the BSE uses only differences. The `[K+QP] Indirect Gap` block of the default path is not printed.

Tested here (gfortran 13.3, serial build, silicon):

- Compiles. Flag off, patched against unpatched binary, undoped and doped: the `o-*.eps` and `o-*.eel`
  data are byte-identical, and `ndb.BS_diago_Q1` (energies, residuals, parameters) is identical at
  17 digits except for the random run identifier that Yambo stores (`PARS1`, first entry).
- Flag on, undoped: not identical, and it should not be. At 300 K the KS filling has thermal tails up to
  1.6e-5 that the larger QP gap removes; the spectrum changes by at most 0.006 on a scale of 90. For
  undoped runs leave the flag off.
- Flag on, doped: the BSE runs with the KS occupations; a refill would have moved 0.0018 electrons,
  the lowest bright interband excitation moves by 4e-6 eV. The silicon GW shift is nearly rigid, so this
  test does not reproduce a WS2-sized effect: that is what A4 is for.

Not tested: the kernel export and Casida with the flag (the occupation guard of `K_kernel.F` compares
them, and `exact_kernel.sh` puts the flag in every BSE input, so a mismatch stops with an error);
`BSEChiOut= "full"` still stops with QP corrections and partly filled bands (`K_driver_init.F:78-83`),
although with the flag its reason no longer applies; MPI runs.

### A4. Check on the cluster

`exact_kernel.sh` has two new options and one new stage (other stages are unchanged when they are not
set):

- `STAGE=occ`: one independent-particle run, no kernel, about a minute. Prints the table above.
- `KEEPOCC=1`: adds `KfnQP_KeepOcc` to every BSE input.
- `PROD=<name>`: folder of the run (default `prod`). Runs with and without the flag must not share a
  folder; the script stops if they would. The screening is linked from `prod` (it does not depend on this
  choice) and the GW is read from where it is.

Fail-fast checks: the source has the flag and the `yambo` binary is newer than it (before anything
runs); the folder's existing reference BSE has the same kind of occupations; after the reference BSE with
`KEEPOCC=1` the report must contain the `KfnQP_KeepOcc` line, or the script stops.

```
cd ~/yambo-cpu && git pull && make yambo          # then, in each doped 24x24 folder, with your usual GWDIR=...

# 1. a minute per doping: what the refill does (no kernel is built)
GWDIR=<gw folder> STAGE=occ sbatch --time=00:15:00 --ntasks-per-node=1 ~/yambo-cpu/sbin/chi_tools/jobs/exact_kernel.sh

# 2. reference BSE with the KS occupations, in a new folder (reuses the GW and prod/w)
GWDIR=<gw folder> KEEPOCC=1 PROD=prod_ko STAGE=ref sbatch ~/yambo-cpu/sbin/chi_tools/jobs/exact_kernel.sh
```

Step 1, what to look for in `log/xkernel.o*` (or `prod/r-occ_*`), band 14, the k-point with the lowest
`E ref` (K):

- 0.05: "f QP refill" = 1.000000 and "full: kept n, refill y". That confirms A1.
- 0.1: "f QP refill" about 0.9997 and "refill n".
- The `E QP` column of band 14 shows the dip: the listed (occupied) points several tenths of an eV below
  their neighbours.
- "Electrons a QP refill would move" of a few 0.01: the two fillings differ by a sizeable part of the doping.

If instead the table is nearly empty (as in silicon), the refill is not the cause and A1 is wrong: stop
and send me the report.

Step 2, what the fix should give, in `prod_ko/r-bse_ref_*`:

- the `[K+QP] KfnQP_KeepOcc` lines;
- "interband gap, lowest transition from a full band" equal to the K gap of the 24x24 grid at every
  doping, so about 2.67, 2.25, 2.19 eV: at 0.05 it comes down from 2.6445, at 0.01 and 0.1 it stays;
- the K transition present with f_eh about 0.65, 0.06, 0.03 (visible later in
  `o-xs_q1.Ktt_q1_transitions`, last column).

What I cannot predict is A itself. With f_eh of a few percent the K transitions are still nearly
decoupled, so the lowest bright interband state may again be an almost bare K transition just under the
K gap, at 0.05 as well as at 0.1. If so, A would become monotonic (about 2.50, 2.2, 2.16 eV) for the
wrong reason: it would follow the occupation-dependent K gap. Check the absolute height of the peak in
`o-bse_ref.alpha_q1_*`, not only the relative strength: a real exciton has a height comparable to the
undoped A, a nearly blocked transition is smaller by its f_eh.

### A5. Does the grid resolve these dopings?

Count the states. The conduction band holds two electrons per k-point, so X extra electrons fill
X/2 x N_k k-points in total, shared by the two K valleys and whatever the Q valleys take.

| doping | k-points in the sea, 24x24 BSE grid (576) | per K valley, at most | same, 12x12 GW grid (144) |
|---|---|---|---|
| 0.01 | 2.9 | 1.4 | 0.36 |
| 0.05 | 14.4 | 7.2 (K and its 6 neighbours) | 1.8 |
| 0.1 | 28.8 | 14.4 | 3.6 |

Energy scales, for a conduction mass of about 0.3 m_e: the first shell of k-points around K on the 24x24
grid is 0.095 1/A away, 0.115 eV above K, more than four times kT = 0.026 eV. At the Fermi wave vector of
0.05 e/cell (0.13 1/A) neighbouring shells are about 0.3 eV apart.

Consequences:

- **The 24x24 BSE does not resolve the Fermi edge.** The filling changes in jumps of one shell of
  k-points, the Fermi level is a staircase in the doping, and each doping samples the sea with a
  different handful of points. Differences between 0.01, 0.05 and 0.1 of the order of 0.1 eV are within
  what one k-point crossing the edge does. I would not compare dopings quantitatively on this grid.
- **Do not lower ElecTemp.** The smearing is already much smaller than the level spacing; a lower
  temperature makes the occupations sharper, pushes more states over the 1e-5 cutoff and makes the
  staircase worse. A larger temperature would smooth it, but then it is a numerical smearing, not 300 K,
  and it enters the screening and the GW as well.
- **The GW grid is the more serious limit.** On 12x12 the sea of a valley is one k-point or less. The
  occupation-dependent step of A1 comes from there, and it sets the K gap that every later number uses.
- **What converged would need.** A level spacing at the Fermi edge near kT means a spacing of about
  0.008 1/A, a 300x300 grid: not reachable with a uniform BSE grid. I have not checked whether Yambo's
  BSE can use a denser patch around K only; a double grid changes the independent-particle part, not the
  kernel, so I do not count on it.

What I would do, in this order, each step cheap enough to stop after:

1. `STAGE=occ` at the three dopings (minutes).
2. From the three `o-gw.qp`: correction of band 14 against KS occupation, all k. If it is the line of A1,
   redo the GW of one doping (0.05) on 18x18 and 24x24 for bands 13-15 only and see whether the step
   shrinks. If it does not converge, the doped GW needs the carrier screening treated as in the papers of
   A2 before any BSE.
3. Reference BSE of 0.05 with `KEEPOCC=1` on 24x24, 36x36 and 48x48 (three bands: 1152, 2592 and 4608
   transitions). Compare the interband gap and A. Only if they settle to within the accuracy you want to
   quote is a comparison between dopings meaningful.
4. Lower dopings are harder, not easier: at 1e13 1/cm^2 the 24x24 sea is about one k-point per valley.

Open after this: the three papers' treatment of the occupations (A2); the export and Casida with the
flag; the doped screening at q -> 0 (already listed in `HANDOFF_code.md`).
