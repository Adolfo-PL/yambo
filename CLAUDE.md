# CLAUDE.md

Yambo 5.3.0 fork (`Adolfo-PL/yambo`) that builds the exact TDDFT kernel f_xc = χ₀⁻¹ − P⁻¹ from GW/BSE,
so that TDDFT (Casida) reproduces the BSE. Physics goal: excitons, and later trions, in n-doped monolayer
WS2 (also WSe2). The advisor's TDDFT trion code will consume these kernels.

## Read first

- `doc/HANDOFF_code.md`: theory in one page, the code map, and the working rules of the Fortran work.
- `doc/greens_function_fxc.md`: the full reference (GW → G → χ → f_xc).
- `doc/casida_solvers.md`: the Casida solvers' equations mapped to code lines.
- `doc/NUMERICS_PROVENANCE.md`: where every numerical piece comes from, and every source file of the fork.
- `doc/casida_comparison.md`: our Casida solver vs Volodymyr's k-space program, with the math.
- `doc/HANDOFF_doped_occupations.md`: the open doped-occupations problem.

Fork code: `src/bse/K_Chi_*.F`, `src/tddft/TDDFT_Chi_dynamic.F`, `src/modules/mod_Chi.F`, `src/io/io_Chi.F`.
Analysis scripts: `sbin/chi_tools/*.py`. Cluster job scripts: `sbin/chi_tools/jobs/*.sh` (the main
chain is `exact_kernel.sh`, with `STAGE=` options; GW: `gw_converge.sh`, `gw_run.sh`).

## Git

- Work and push only on `claude/nbd-file-generation-path-yidcr2`: `git push -u origin claude/nbd-file-generation-path-yidcr2`,
  retrying network failures after 2, 4, 8 and 16 s. If the push is rejected as divergent (another session
  pushes to the same branch), `git pull --no-rebase --no-edit` and push again.
- No pull requests unless the user asks.
- No model names or identifiers in commits, code or documents.
- End commit messages with the attribution lines the session's system reminder gives.
- Run `rm -rf sbin/chi_tools/__pycache__ sbin/chi_tools/jobs/__pycache__` before every commit; never commit them.

## Code rules

- New source files carry `Copyright (C) 2026 Adolfo Partida` and
  `Authors: Adolfo Partida, with AI assistance (Claude Code)`. Upstream files keep their own authors.
- Fortran: match the surrounding Yambo style (the `!` spacer lines, `ONLY:` lists, short comments).
- Never modify Volodymyr's program (`eigenvalues_SPN.f90`, outside this repo). Changed versions are
  separate copies; the comparison is made by converting our data to his input (`ktt_to_kq.py`).
- Every SBATCH script header includes
  ```
  #SBATCH --mail-type=ALL,TIME_LIMIT_80
  #SBATCH --mail-user=ad531114@ucf.edu
  ```
- Cluster time is expensive: job scripts must check their inputs and fail fast before the expensive step.

## Cluster (UCF Stokes)

- The fork lives at `~/yambo-cpu` there; modules `cuda_hpc_sdk/24.5`, `hdf5/hdf5-1.14.3-cuda_hpc_sdk_24.5.0`,
  `yambo/5.3-cpu`, then `hpcx_load` (see `jobs/gw_run.sh`).
- Typical node: 24 MPI tasks × 2 OpenMP threads, 150G, partition `normal`. In Yambo inputs every
  `*_CPU` product must equal the total number of MPI tasks.
- Production WS2: PseudoDojo NC fully relativistic PBE (identical in the standard and stringent sets;
  W 5s5p5d6s, S 3s3p: 13 occupied bands), 100/400 Ry, slab cutoff, RIM-W, PPA, SOC off.
  Bands are numbered from 1 with the VBM at 13 (20 with SG15, whose W includes 4f).

## Working with the user

- PhD student, new to this formalism. Explain plainly, with the math when asked.
- Back claims with sources (a Yambo source line, a paper); do not guess. Be exact: the user corrects
  overstatements.
- When the user asks for a file, send the file.
- Two sessions share this branch: one handles cluster runs, results and explanations; the other the theory,
  Fortran, tests and new features. Pull before editing.
