#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=normal
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --job-name=fkq_simple
#SBATCH --output=fkq_simple.o%A

# F_kq for ONE q: Yambo reads the kernel f_xc(G,G',q,w) from KERNEL/ndb.Chi, projects it on
# the transitions,
#     F_tt'(w) = c sum_GG' rho_t(q+G) f_xc(G,G',q,w) rho_t'(q+G')*
# and writes the Casida matrix pieces (X and F) at every kernel frequency between W1 and W2:
#     o-fkq_q<Q>.Ktt_q<Q>_transitions    one line per transition: k, v, c, E_t, dipole
#     o-fkq_q<Q>.Ktt_q<Q>_antiresonant   (q > 1 only) the anti-resonant rows
#     o-fkq_q<Q>.Ktt_q<Q>_w<iw>          one file per frequency, ~1.1 GB for 2025 transitions
# Then:  ./casida_simple o-fkq_q<Q>.Ktt_q<Q>_transitions o-fkq_q<Q>.Ktt_q<Q>_w*
# Run from the directory with SAVE/ and g0w0/.

Q=1                  # q index (1 = optical limit)
W1=0.0; W2=0.0       # frequency window [eV]; 0 0 = static kernel (z = 0)

# the frequency grid MUST be the one the kernel was exported on:
KERNEL=kx_fq; WLO=0.0;  WHI=4.0;  NSTEP=201; ETA=0.05    # fq_export.sh
#KERNEL=kx_xq; WLO=2.15; WHI=2.54; NSTEP=40;  ETA=0.0    # xq_dispersion.sh (undamped, for the scan)

module purge
module use "$HOME/modulefiles"
module load cuda_hpc_sdk/24.5
module load hdf5/hdf5-1.14.3-cuda_hpc_sdk_24.5.0
module load yambo/5.3-cpu
source "$NVHPC_ROOT/comm_libs/12.4/hpcx/latest/hpcx-init.sh" && hpcx_load
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "${SLURM_SUBMIT_DIR:-.}"

cat > fkq_q$Q.in <<EOF
optics
bss
bse
dipoles
BSKmod= "CHI"
BSEmod= "coupling"
BSSmod= "d"
BSEChiDyn
BSEGinplane
BSENGexx= 21 RL
BSENGfxc= 21 RL
K_Threads= $OMP_NUM_THREADS
DIP_Threads= $OMP_NUM_THREADS
% BSEQptR
 $Q | $Q |
%
% BSEBands
 11 | 16 |
%
% BEnRange
 $WLO | $WHI | eV
%
% BDmRange
 $ETA | $ETA | eV
%
BEnSteps= $NSTEP
% BLongDir
 1.0 | 0.0 | 0.0 |
%
BSEprop= "abs"
% BSEChiKoutW
 $W1 | $W2 | eV
%
% BSEChiKoutB
 0 | 0 |
%
KfnQPdb= "E < g0w0/ndb.QP"
EOF

rm -rf fkq_q$Q
yambo -nompi -F fkq_q$Q.in -J "fkq_q$Q,$KERNEL" > fkq_q$Q.log 2>&1
ls -lh o-fkq_q$Q.Ktt_*
