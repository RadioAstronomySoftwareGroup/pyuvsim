#!/bin/bash
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

#SBATCH -J pyuvsim
#SBATCH --array=0-1
#SBATCH --mem=15G
#SBATCH --nodes=2
#SBATCH --tasks-per-node=24
#SBATCH -t 10:00:00

export PYUVSIM_BATCH_JOB=1

echo JOBID ${SLURM_ARRAY_JOB_ID}
echo TASKID ${SLURM_ARRAY_TASK_ID}


ncorepernode=${SLURM_NTASKS_PER_NODE}
nnodes=${SLURM_JOB_NUM_NODES}
task=${SLURM_ARRAY_TASK_ID}
jobid=${SLURM_ARRAY_JOB_ID}

if [ "$#" -ne 4 ]; then
    echo 'Usage: sbatch batch_pyuvsim.sh <Nsrcs> <Nchan> <Ntimes> <Instrument>'
    exit
fi

nsrc=$1
nchan=$2
ntime=$3
array=$4

dir1=$array'_sim_'$nsrc'src_'$nchan'chan_'$ntime'time_'$ncorepernode'cores_per_'$nnodes'nodes'

if [ ! -d "$dir1" ]; then
    mkdir $dir1
fi

if [ "$array" == "triangle" ]; then
    infile="../pyuvsim/data/28m_triangle_10time_10chan.uvfits"
fi
if [ "$array" == "mwa" ]; then
    infile="/users/alanman/data/alanman/MWA128/GoldenSet/1061315448.uvfits"
fi

echo "File: "$infile
echo "$nsrc sources, $nchan chan, $ntime time"

if [ "$task" -eq '0' ]; then
    srun --mpi=pmi2 python -m memory_profiler run_pyuvsim.py --outdir $dir1 --mock_arrangement zenith --Nsrcs $nsrc --Nchans $nchan --Ntimes $ntime $infile > $dir1/mem_profile.out
fi

#pids+=($!)

if [ "$task" -eq '1' ]; then
    srun --mpi=pmi2 kernprof -l -v run_pyuvsim.py --outdir $dir1 --mock_arrangement zenith --Nsrcs $nsrc --Nchans $nchan --Ntimes $ntime $infile > $dir1/time_profile.out
fi

#pids+=($!)

#echo ${pids[@]}

#wait ${pids[@]}

export PYUVSIM_BATCH_JOB=0

## Try to clean up the scripts directory
ofilename='slurm-'$jobid'_'$task'.out'
mv $ofilename $dir1/$ofilename
