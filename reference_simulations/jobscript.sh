#!/bin/bash

#SBATCH -J pyuvsim
#SBATCH --mem=100G
#SBATCH -t 29:00:00
#SBATCH -n 50
# SBATCH -n 25
# SBATCH --cpus-per-task=4
#SBATCH --cpus-per-task=2

echo JOBID ${SLURM_ARRAY_JOB_ID}
echo TASKID ${SLURM_ARRAY_TASK_ID}

task=${SLURM_ARRAY_TASK_ID}
jobid=${SLURM_ARRAY_JOB_ID}

obsparam_configs=("$@")
obsparam=${obsparam_configs[$task]}

echo -e `date`'\t'$task'\t'$obsparam >> RUN_LOG_$jobid

dir=$(pwd -P)

scriptdir=$(python -c "import pyuvsim; print('/'.join(pyuvsim.__file__.split('/')[:-2])+ '/scripts')")

srun --kill-on-bad-exit --mpi=pmi2 python $scriptdir/run_param_pyuvsim.py --profile='profiling/tprof_'$jobid'_'$task'.out' $dir"/"$obsparam
