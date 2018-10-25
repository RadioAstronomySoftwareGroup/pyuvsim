#!/bin/bash
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

#SBATCH -J pyuvsim_profile
#SBATCH --array=0-0

echo JOBID ${SLURM_ARRAY_JOB_ID}
echo TASKID ${SLURM_ARRAY_TASK_ID}

ntasks=${SLURM_NTASKS}
nnodes=${SLURM_JOB_NUM_NODES}
task=${SLURM_ARRAY_TASK_ID}
jobid=${SLURM_ARRAY_JOB_ID}

branch=`git branch | grep \* | cut -d ' ' -f2`

_IFS=$IFS
IFS=','
read -r -a Nsrcs <<< "$1"
read -r -a Ntimes <<< "$2"
read -r -a Nfreqs <<< "$3"
read -r -a Nbls <<< "$4"
read -r -a beams <<< "$5"
IFS=$_IFS

echo ${Ntimes[@]}

dir1=$branch'_profiling/sim_'$nnodes'nodes_'$ntasks'cpus'
#
if [ ! -d "$dir1" ]; then
    mkdir -p $dir1
fi

slids_out="slurm_ids_"$branch".out"
if [ ! -f $slids_out ]; then
   echo 'Npus, Nnodes, Nsrcs, Ntimes, Nfreqs, Nbls, beam, MaxMem, Elapsed' > $slids_out
fi

for nsrcs in "${Nsrcs[@]}"; do
    for ntimes in "${Ntimes[@]}"; do
        for nfreqs in "${Nfreqs[@]}"; do
            for nbls in "${Nbls[@]}"; do
                for beam in "${beams[@]}"; do
                    START=$(date +%s)   # Timing
                    if [ "$task" -eq '0' ]; then
                        srun --mpi=pmi2 python run_profile_pyuvsim.py --Nsrcs $nsrcs --Ntimes $ntimes --Nfreqs $nfreqs --Nbls $nbls --beam $beam \
                                --prof_out $dir1/time_profile.out --mem_out $dir1/memory_usage.out
                    fi
                    END=$(date +%s)
                    DIFF=$(( $END - $START ))
                    mem_used=$(<$dir1'/memory_usage.out')
                    echo $ntasks", "$nnodes', '$nsrcs', '$ntimes', '$nfreqs', '$nbls', '$beam', '$mem_used', '$DIFF
                done
            done
        done
    done
done

rm $dir1'/memory_usage.out'

## Try to clean up the scripts directory
ofilename='slurm-'$jobid'_'$task'.out'
mv $ofilename $dir1/$ofilename
