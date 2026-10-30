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
