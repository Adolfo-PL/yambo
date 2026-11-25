#!/bin/bash
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --job-name=dyson_q0
#SBATCH --output=log/dyson_q0.o%A
#SBATCH --error=log/dyson_q0.e%A
#SBATCH --mail-type=ALL,TIME_LIMIT_80
#SBATCH --mail-user=ad531114@ucf.edu
#
# ChiGMode="DYSON" at q=0 (optical limit): chi0 (KS), P = chi^[G1] from the retarded PPA G0W0
# Green functions, and fxc = chi0^-1 - P^-1, in the BSE band window.
#   1. gw_g      : PPA G0W0 with GRetarded -> gw_g/ndb.G (G on a real grid, every IBZ k, window bands)
#   2. chi_dyson : screening-only em1d run at q=0 with ChiFxc -> chi_dyson/ndb.Chi
# Both run here, on this folder's SAVE (no GWDIR interpolation: ndb.G must be on the response k grid),
# at zero electronic temperature (GRetarded and the Chi/G routes stop at finite T).
# Run from the folder with SAVE, undoped, e.g. the 12x12 GW folder:
#   mkdir -p log && STAGE=test sbatch dyson_q0.sh     # one QP state: times the self-energy, builds the T=0 screening
#   mkdir -p log && sbatch dyson_q0.sh                # full: gw_g (reusing gw_g_test's screening) then chi_dyson
#   STAGE=chi sbatch dyson_q0.sh                      # only step 2, on an existing gw_g/ndb.G
#   STAGE=validate sbatch dyson_q0.sh                 # GRetarded vs Newton, with and without the BG terminator
# Terminator: GTERM=BG (default, 400 bands: converged for WS2 with BG) or GTERM= (none, 800 bands).
# STAGE=validate runs four GW on the states VSTATES (default: the window bands at k 1), on one T=0
# screening: Newton BG, Newton none (ExtendOut, GDamping = GDAMP) and GRetarded BG, GRetarded none.
# gf_check.py then requires, for each GRetarded ndb.G: no negative spectral weight, |norm - 1| <= NORMTOL,
# Re Sigma_c(Eo) on the grid = Newton's Sc|Eo (same PPA sum at Eo + i*GDAMP; the real part is the same
# for retarded and time-ordered G), and the BG shift of Re Sigma_c(Eo) on the grid = the Newton BG shift.
# Options: BANDS=NV:NC (window, default 1:2 = VB..VB+2), GW_BANDS, NP_CHI (tasks of step 2: with several
# nodes the window dipoles may not fit all tasks, e.g. NP_CHI=24 for bands 13-15 on 18x18).
# The self-energy is evaluated at GSTEPS frequencies per state (Newton uses 2): from the test run,
# full self-energy time ~ (test self-energy time) x (number of QP states printed in the plan).

STAGE=${STAGE:-all}                  # test | all | gw | chi | validate
GTERM=${GTERM-BG}                    # BG terminator (GRetarded accepts BG or none) or empty for none
GW_EXX=${GW_EXX:-40}; GW_BLOCK=${GW_BLOCK:-12}; RIMW_G=${RIMW_G:-15}
XP_BANDS=${XP_BANDS:-400}            # screening bands (production)
PPA=${PPA:-'PPAPntXp= 27.21138 eV'}
GRANGE=${GRANGE:-40}                 # eV: G grid is E_KS +- GRANGE (holds the QP peak and the satellites)
GDAMP=${GDAMP:-0.136}                # eV: broadening of G; the grid step is set equal to it
NGX=${NGX:-"2 Ry"}                   # response / kernel basis
WMAX=${WMAX:-10.88}; WSTEPS=${WSTEPS:-800}; ETA=${ETA:-0.136}   # eV: response grid (exciting 0-0.4 Ha, 800, 0.005 Ha)
NORMTOL=${NORMTOL:-0.05}; CHECKTOL=${CHECKTOL:-0.0001}; RCOND=${RCOND:-1.0e-12}
CUT=${CUT:-'CUTGeo= "slab z"'}
RIM=${RIM:-$'RandQpts= 1000000\nRandGvec= 100 RL'}
YSRC=${YSRC:-$HOME/yambo-cpu}

if type module >/dev/null 2>&1 && [ -n "$SLURM_JOB_ID" ]; then
  module purge
  module use "$HOME/modulefiles"
  module load cuda_hpc_sdk/24.5
  module load hdf5/hdf5-1.14.3-cuda_hpc_sdk_24.5.0
  module load yambo/5.3-cpu
  source "$NVHPC_ROOT/comm_libs/12.4/hpcx/latest/hpcx-init.sh" && hpcx_load
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}
NP=${NP:-${SLURM_NTASKS:-1}}
NP_CHI=${NP_CHI:-$NP}                # tasks of step 2 (k-only; the window dipoles need k*c*v = NP_CHI)
RUN="srun -n $NP --cpus-per-task=$OMP_NUM_THREADS"
[ "$NP" -gt 1 ] || RUN=""
RUN_CHI="srun -n $NP_CHI --cpus-per-task=$OMP_NUM_THREADS"
[ "$NP_CHI" -gt 1 ] || RUN_CHI=""

case $STAGE in test|all|gw|chi|validate) ;; *) echo "!!! STAGE=$STAGE: use test, all, gw, chi or validate"; exit 1;; esac
case "$GTERM" in BG|"") ;; *) echo "!!! GTERM=$GTERM: GRetarded accepts BG or none (empty)"; exit 1;; esac

# --- the fork: q=0 DYSON and GRetarded, and a binary built after them ---
grep -q "optical q=0 DYSON" $YSRC/src/pol_function/X_pre_setup.F 2>/dev/null &&
grep -q "GRetarded" $YSRC/src/qp/QP_driver.F 2>/dev/null &&
grep -q 'GRetarded supports GTermKind="BG"' $YSRC/src/qp/QP_driver.F 2>/dev/null || {
  echo "!!! $YSRC lacks q=0 DYSON or GRetarded with BG: git pull there (branch claude/nbd-file-generation-path-yidcr2)"; exit 1; }
[ "$(which yambo)" -nt $YSRC/src/qp/QP_driver.F ] || {
  echo "!!! $(which yambo) is older than $YSRC/src/qp/QP_driver.F: make yambo"; exit 1; }
GFCHECK=$YSRC/sbin/chi_tools/gf_check.py
python3 -c "import netCDF4, numpy" 2>/dev/null || { echo "!!! python3 lacks netCDF4/numpy (needed by gf_check.py)"; exit 1; }
[ "$(which yambo)" -nt $YSRC/src/pol_function/Chi_G_bubble.F ] || {
  echo "!!! $(which yambo) is older than $YSRC/src/pol_function/Chi_G_bubble.F: make yambo"; exit 1; }

# --- the system ---
[ -f SAVE/ns.db1 ] || { echo "!!! SAVE/ns.db1 not found: run from the folder with SAVE"; exit 1; }
[ -f r_setup ] || $RUN yambo > log_setup.out 2>&1
NB=$(awk -F: '/^ *Bands  /{print $2+0; exit}' r_setup 2>/dev/null)
NK=$(awk -F: '/IBZ K-points/{print $2+0; exit}' r_setup 2>/dev/null)
NKBZ=$(awk -F: '/BZ +K-points/ && !/IBZ/{print $2+0; exit}' r_setup 2>/dev/null)
NQ=$(awk -F: '/IBZ Q-points/{print $2+0; exit}' r_setup 2>/dev/null)
VB=$(awk -F: '/Filled Bands/{print $2+0; exit}' r_setup 2>/dev/null)
NSP=$(awk -F: '/Spinor components/{print $2+0; exit}' r_setup 2>/dev/null)
[ -n "$NB" ] && [ -n "$NK" ] && [ -n "$NKBZ" ] && [ -n "$NQ" ] && [ -n "$VB" ] && [ -n "$NSP" ] || {
  echo "!!! cannot read bands/k/q/spinors from r_setup (yambo setup failed?)"; tail -5 log_setup.out 2>/dev/null; exit 1; }
grep -q "Metallic Bands" r_setup && { echo "!!! metallic (doped) SAVE: the Chi/G routes need an undoped insulator"; exit 1; }
SF=$NSP
GBASE=800; [ -n "$GTERM" ] && GBASE=400            # converged for WS2: 400 with BG, or 800
GW_BANDS=${GW_BANDS:-$((GBASE * SF))}
BANDS=${BANDS:-$([ "$SF" = 2 ] && echo 2:4 || echo 1:2)}
NV=${BANDS%%:*}; NC=${BANDS##*:}
BLO=$((VB - NV + 1)); BHI=$((VB + NC))
NQP=$((NK * (BHI - BLO + 1)))
GSTEPS=$(awk -v r=$GRANGE -v d=$GDAMP 'BEGIN{printf "%d", 2*r/d + 1.5}')

msg=""
[ "$GW_BANDS" -le "$NB" ] || msg="$msg\n!!! GW_BANDS=$GW_BANDS > the $NB bands of the SAVE (WS2 converged: 400 with GTERM=BG, 800 without)"
[ "$XP_BANDS" -le "$GW_BANDS" ] || msg="$msg\n!!! XP_BANDS=$XP_BANDS > GW_BANDS=$GW_BANDS"
[ "$NV" -le "$VB" ] || msg="$msg\n!!! $NV valence bands requested, there are $VB"
[ "$NP_CHI" -le "$NKBZ" ] || msg="$msg\n!!! NP_CHI=$NP_CHI tasks > $NKBZ BZ k: the response uses k-only MPI"
[ "$NP_CHI" -le "$NP" ] || msg="$msg\n!!! NP_CHI=$NP_CHI > the $NP tasks of the job"
[ -z "$msg" ] || { echo -e "$msg"; exit 1; }

# --- MPI layouts (every *_CPU product = NP) ---
split () {   # ntasks cap1 cap2 ... -> tasks per role, greedy largest divisor; nothing if impossible
  local r=$1 out="" c d; shift
  for c in "$@"; do
    for ((d = c < r ? c : r; d > 1; d--)); do [ $((r % d)) -eq 0 ] && break; done
    [ "$d" -lt 1 ] && d=1
    out="$out $d"; r=$((r / d))
  done
  [ "$r" -eq 1 ] && echo $out
}
par_gw () {   # nqp -> DIP, X, SE lines of the GW
  [ "$NP" -gt 1 ] || return 0
  local nc=$((GW_BANDS - VB)) d x s
  d=($(split $NP $nc $NK)); [ ${#d[@]} -eq 2 ] || { echo "!!! $NP tasks do not fit the GW dipoles" >&2; return 1; }
  x=($(split $NP $nc $NKBZ $NQ)); [ ${#x[@]} -eq 3 ] || { echo "!!! $NP tasks do not fit the GW screening" >&2; return 1; }
  s=($(split $NP $BHI $NQ $1)); [ ${#s[@]} -eq 3 ] || { echo "!!! $NP tasks do not fit the self-energy" >&2; return 1; }
  printf 'DIP_CPU= "%d %d 1"\nDIP_ROLEs= "k c v"\n' ${d[1]} ${d[0]}
  printf 'X_and_IO_CPU= "%d 1 %d %d 1"\nX_and_IO_ROLEs= "q g k c v"\n' ${x[2]} ${x[1]} ${x[0]}
  printf 'SE_CPU= "%d %d %d"\nSE_ROLEs= "q qp b"\n' ${s[1]} ${s[2]} ${s[0]}
}
par_chi () {  # DIP (window bands) and k-only X lines of the response
  [ "$NP_CHI" -gt 1 ] || return 0
  local d; d=($(split $NP_CHI $NK $NC $NV)); [ ${#d[@]} -eq 3 ] || {
    echo "!!! NP_CHI=$NP_CHI tasks do not fit the dipoles of bands $BLO-$BHI: k <= $NK, c <= $NC, v <= $NV (set NP_CHI, e.g. 24)" >&2; return 1; }
  printf 'DIP_CPU= "%d %d %d"\nDIP_ROLEs= "k c v"\n' ${d[0]} ${d[1]} ${d[2]}
  printf 'X_and_IO_CPU= "1 1 %d 1 1"\nX_and_IO_ROLEs= "q g k c v"\n' $NP_CHI
}

solver_lines () {   # g | n
  if [ "$1" = g ]; then
    printf 'GRetarded\nDysSolver= "g"\nGsampType= "ra"\nGEnMode= "centered"\nGEnSteps= %s\n%% GEnRnge\n -%s | %s | eV\n%%\n%% GDmRnge\n %s | %s | eV\n%%\nGDamping= 0 eV\nGreenFTresh= 0\n' \
      $GSTEPS $GRANGE $GRANGE $GDAMP $GDAMP
  else
    printf 'DysSolver= "n"\nGDamping= %s eV\nExtendOut\n' $GDAMP
  fi
}
gw_in () {   # file qpkrange-line nqp [solver g|n] [terminator BG|none]
  local par; par=$(par_gw $3) || return 1
  cat > $1 <<EOF
$par
rim_w
RandGvecW= $RIMW_G RL
HF_and_locXC
gw0
ppa
el_el_corr
dyson
rim_cut
em1d
ElecTemp= 0 eV
BoseTemp= 0 eV
X_Threads= $OMP_NUM_THREADS
DIP_Threads= $OMP_NUM_THREADS
SE_Threads= $OMP_NUM_THREADS
$CUT
$RIM
EXXRLvcs= $GW_EXX Ry
VXCRLvcs= $GW_EXX Ry
Chimod= "HARTREE"
% BndsRnXp
 1 | $XP_BANDS |
%
NGsBlkXp= $GW_BLOCK Ry
% LongDrXp
 1.0 | 1.0 | 0.0 |
%
$PPA
% GbndRnge
 1 | $GW_BANDS |
%
GTermKind= "${5:-${GTERM:-none}}"
$(solver_lines ${4:-g})
%QPkrange
 $2
%
EOF
}

chi_in () {
  local par; par=$(par_chi) || return 1
  cat > $1 <<EOF
$par
em1d
rim_cut
ChiFxc
ChiGMode= "DYSON"
ChiGDb= "gw_g/ndb.G"
ChiGAxis= "REAL"
ChiGNormTol= $NORMTOL
ChiGCheckTol= $CHECKTOL
ChiRcondMin= $RCOND
ElecTemp= 0 eV
BoseTemp= 0 eV
X_Threads= $OMP_NUM_THREADS
DIP_Threads= $OMP_NUM_THREADS
$CUT
$RIM
Chimod= "HARTREE"
XTermKind= "none"
% QpntsRXd
 1 | 1 |
%
% BndsRnXd
 $BLO | $BHI |
%
NGsBlkXd= $NGX
% EnRngeXd
 0.0 | $WMAX | eV
%
% DmRngeXd
 $ETA | $ETA | eV
%
ETStpsXd= $WSTEPS
% LongDrXd
 1.0 | 0.0 | 0.0 |
%
% EhEngyXd
 -1.0 | -1.0 | eV
%
EOF
}

# both layouts must fit before anything runs
VSTATES=${VSTATES:-"1|1|$BLO|$BHI|"}
NV_QP=$(echo "$VSTATES" | awk -F'|' '{print ($2-$1+1)*($4-$3+1)}')
[ "$STAGE" = chi ] || par_gw $NQP > /dev/null || exit 1
[ "$STAGE" = validate ] && { par_gw $NV_QP > /dev/null || exit 1; }
[ "$STAGE" = test ] || [ "$STAGE" = validate ] || par_chi > /dev/null || exit 1

fresh () { rm -rf $1 o-$1.* r-$1_* LOG/l-$1_*; }
run () {   # job -J-string done-test [runner]
  local t0=$SECONDS r=${4-$RUN}
  $r yambo -F $1.in -J "$2" > $1.out 2>&1
  if eval "$3"; then echo "  $1 done in $(( (SECONDS - t0) / 60 )) min ($(date +%H:%M))"; return 0; fi
  echo "!!! $1 failed: see $1.out and LOG/l-$1_* (or l-$1_*)"
  grep -h "ERROR" LOG/l-$1_* l-$1_* 2>/dev/null | sed 's/.*\[ERROR\]/  [ERROR]/' | sort -u | head -5
  return 1
}
timing () { awk '/Timing/{p=1} p' r-$1_* 2>/dev/null | head -40; }

echo "=== $(basename "$PWD"): $NB bands, $NK IBZ k, $NKBZ BZ k, VB $VB, $([ $SF = 2 ] && echo spinors || echo no SOC); $NP MPI x $OMP_NUM_THREADS threads ==="
echo "    GW : PPA, T=0, GRetarded, $GW_BANDS bands, terminator ${GTERM:-none}, screening $XP_BANDS bands $GW_BLOCK Ry"
echo "         G grid E_KS +- $GRANGE eV, $GSTEPS steps, broadening $GDAMP eV; QP states: bands $BLO-$BHI at $NK k = $NQP"
echo "    chi: $NP_CHI MPI, q=0, bands $BLO-$BHI, $NGX, 0-$WMAX eV, $WSTEPS steps, eta $ETA eV, ChiGNormTol $NORMTOL"

if [ "$STAGE" = test ]; then
  fresh gw_g_test; gw_in gw_g_test.in "1|1|$VB|$VB|" 1 || exit 1
  run gw_g_test "$(scr gw_g_test)" "ls gw_g_test/ndb.G* >/dev/null 2>&1" || exit 1
  timing gw_g_test
  echo "  self-energy of the full run ~ (self-energy time above) x $NQP; the screening (gw_g_test/ndb.pp) is reused"
  exit 0
fi

scr () {   # job -> -J string reusing an existing T=0 screening (only ndb.pp*, linked in scr_T0)
  local d; for d in gw_g_test v_nbg; do
    [ -e $d/ndb.pp ] || continue
    mkdir -p scr_T0; local f; for f in $d/ndb.pp*; do ln -sf ../$f scr_T0/; done; echo "$1,scr_T0"; return
  done; echo $1
}

if [ "$STAGE" = validate ]; then
  echo "  validate on QPkrange $VSTATES ($NV_QP states), $GW_BANDS bands"
  for v in "nbg n BG" "nno n none" "gbg g BG" "gno g none"; do set -- $v
    j=v_$1
    if [ $2 = n ]; then done_t="ls o-$j.qp >/dev/null 2>&1"; else done_t="ls $j/ndb.G* >/dev/null 2>&1"; fi
    if eval "$done_t"; then echo "  $j exists: kept"; continue; fi
    fresh $j; gw_in $j.in "$VSTATES" $NV_QP $2 $3 || exit 1
    run $j "$(scr $j)" "$done_t" || exit 1
  done
  ok=0
  echo "--- GRetarded + BG (vs Newton + BG; BG shift vs the Newton BG shift)"
  python3 $GFCHECK v_gbg/ndb.G --newton o-v_nbg.qp --ref v_gno/ndb.G --ref-newton o-v_nno.qp --norm-tol $NORMTOL || ok=1
  echo "--- GRetarded, no terminator (vs Newton, no terminator)"
  python3 $GFCHECK v_gno/ndb.G --newton o-v_nno.qp --norm-tol $NORMTOL || ok=1
  [ $ok = 0 ] && echo "  VALIDATED: GRetarded with ${GTERM:-no} terminator can be used (sbatch dyson_q0.sh)" ||
    echo "!!! validation failed: see above"
  exit $ok
fi

if [ "$STAGE" = all ] || [ "$STAGE" = gw ]; then
  if ls gw_g/ndb.G* >/dev/null 2>&1; then
    echo "  gw_g/ndb.G exists: kept (rm -rf gw_g to redo)"
  else
    fresh gw_g; gw_in gw_g.in "1|$NK|$BLO|$BHI|" $NQP || exit 1
    run gw_g "$(scr gw_g)" "ls gw_g/ndb.G* >/dev/null 2>&1" || exit 1
    python3 $GFCHECK gw_g/ndb.G --norm-tol $NORMTOL | tail -3
    timing gw_g
  fi
fi

if [ "$STAGE" = all ] || [ "$STAGE" = chi ]; then
  ls gw_g/ndb.G* >/dev/null 2>&1 || { echo "!!! gw_g/ndb.G not found: run STAGE=gw (or all) first"; exit 1; }
  fresh chi_dyson; chi_in chi_dyson.in || exit 1
  run chi_dyson chi_dyson "ls chi_dyson/ndb.Chi* >/dev/null 2>&1" "$RUN_CHI" || exit 1
  grep -h "\[Chi/G\]\|\[Chi\]" r-chi_dyson_* LOG/l-chi_dyson_* 2>/dev/null | sort -u | head -20
  echo "  chi_dyson/ndb.Chi: chi0, P and fxc at q=0"
fi
