#!/bin/bash
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

#SBATCH -J pyuvsim
#SBATCH --array=0-0

echo JOBID ${SLURM_ARRAY_JOB_ID}
echo TASKID ${SLURM_ARRAY_TASK_ID}

ntasks=${SLURM_NTASKS}
nnodes=${SLURM_JOB_NUM_NODES}
task=${SLURM_ARRAY_TASK_ID}
jobid=${SLURM_ARRAY_JOB_ID}

if [ "$#" -ne 5 ]; then
    echo 'Usage: sbatch batch_pyuvsim.sh <Nsrcs> <Ntimes> <Nfreqs> <Nbls> <beam>'
    exit
fi

branch=`git branch | grep \* | cut -d ' ' -f2`

nsrcs=$1
ntimes=$2
nfreqs=$3
nbls=$4
beam=$5

dir1=$branch'_profiling/sim_'$nsrcs'src_'$nfreqs'freq_'$ntimes'time_'$nbls'bls_'$beam'beam_'$nnodes'nodes_'$ntasks'_cpus'

if [ ! -d "$dir1" ]; then
    mkdir -p $dir1
fi

if [ "$task" -eq '0' ]; then
    srun --mpi=pmi2 kernprof -l -v run_profile_pyuvsim.py --Nsrcs $nsrcs --Ntimes $ntimes --Nfreqs $nfreqs --Nbls $nbls --beam $beam > $dir1/time_profile.out
fi

## Try to clean up the scripts directory
ofilename='slurm-'$jobid'_'$task'.out'
mv $ofilename $dir1/$ofilename
