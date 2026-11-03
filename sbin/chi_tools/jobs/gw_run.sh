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
#SBATCH --time=24:00:00
#SBATCH --job-name=gw
#SBATCH --output=gw.o%A
#SBATCH --error=gw.e%A
#SBATCH --mail-type=ALL,TIME_LIMIT_80
#SBATCH --mail-user=ad531114@ucf.edu

# A plain G0W0 plasmon-pole run (yambo -gw0 p), in the folder with SAVE:
#   1st time:  bash gw_run.sh     -> writes gw.in (Yambo's defaults, no editor) and stops: edit it
#   then:      sbatch gw_run.sh   -> runs yambo -F gw.in -J gw
# IN= sets another input name (e.g. IN=gw_test.in), JOB= the job name.

IN=${IN:-gw.in}; JOB=${JOB:-gw}
cd "${SLURM_SUBMIT_DIR:-.}"
if type module >/dev/null 2>&1; then
  module purge
  module use "$HOME/modulefiles"
  module load cuda_hpc_sdk/24.5
  module load hdf5/hdf5-1.14.3-cuda_hpc_sdk_24.5.0
  module load yambo/5.3-cpu
  source "$NVHPC_ROOT/comm_libs/12.4/hpcx/latest/hpcx-init.sh" && hpcx_load
fi
[ -f SAVE/ns.db1 ] || { echo "!!! no SAVE here"; exit 1; }
if [ ! -f "$IN" ]; then
  yambo -nompi -F "$IN" -gw0 p -hf -dyson n -rim_cut -Q > /dev/null 2>&1 || yambo -F "$IN" -gw0 p -hf -dyson n -rim_cut -Q
  echo "wrote $IN: edit it (bands, cutoffs, QPkrange, CUTGeo), then sbatch this script"
  exit 0
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
srun yambo -F "$IN" -J "$JOB"
ls o-$JOB.qp* 2>/dev/null && tail -25 o-$JOB.qp* | grep -v "^#"
