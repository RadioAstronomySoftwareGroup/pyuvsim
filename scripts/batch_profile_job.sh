#!/bin/bash
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

#SBATCH -J pyuvsim_profile
#SBATCH --array=0-0
# SBATCH -q jpober-condo
# SBATCH -A jpober-condo

echo JOBID ${SLURM_ARRAY_JOB_ID}
echo TASKID ${SLURM_ARRAY_TASK_ID}

npus=${SLURM_NTASKS}
nnodes=${SLURM_JOB_NUM_NODES}
task=${SLURM_ARRAY_TASK_ID}
jobid=${SLURM_ARRAY_JOB_ID}
cpuspertask=${SLURM_CPUS_PER_TASK}

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

dir1=$branch'_profiling/sim_'$nnodes'nodes_'$npus'pus_'$cpuspertask'cpuspertask'
#
if [ ! -d "$dir1" ]; then
    mkdir -p $dir1
fi

slids_out="prof_data_"$branch".out"
if [ ! -f $slids_out ]; then
   echo 'JobID,Start,MaxRSS (GB),NNodes,NProcs,Ncpus_per_task,Nbls,Ntimes,Nchan,Nsrc,Beam,Ntasks,Runtime_Seconds' > $slids_out
#   echo 'Npus, Nnodes, Nsrcs, Ntimes, Nfreqs, Nbls, beam, MaxMemGB, ElapsedSec' > $slids_out
fi

mem_out=$dir1"/memory_usage_"$task".out"
time_out=$dir1"/timing_"$task".out"
function do_run {
    nsrcs=$1
    ntimes=$2
    nfreqs=$3
    nbls=$4
    beam=$5

    START=$(date +%s)   # Timing
    start_str=$(date)
    srun --kill-on-bad-exit --mpi=pmi2 python run_profile_pyuvsim.py --Nsrcs $nsrcs --Ntimes $ntimes --Nfreqs $nfreqs --Nbls $nbls --beam $beam \
                --prof_out $dir1"/time_profile_"$nsrcs"src_"$ntimes"t_"$nfreqs"f_"$nbls"bl_"$beam".out" --mem_out $dir1/memory_usage.out
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    mem_used=$(<$dir1'/memory_usage.out')
                --prof_out $dir1"/time_profile_"$nsrcs"src_"$ntimes"t_"$nfreqs"f_"$nbls"bl_"$beam".out" --mem_out $mem_out --time_out $time_out
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    mem_used=$(<$mem_out)
    time_used=$(<$time_out)
    ntasks=$(( $nsrcs * $ntimes * $nfreqs * nbls ))
    echo $jobid'_'$task','$start_str','$mem_used", "$nnodes','$npus', '$cpuspertask', '$nbls', '$ntimes', '$nfreqs', '$nsrcs', '$beam','$ntasks','$time_used >> $slids_out

}

counter=-1
Nother=15
for beam in "${beams[@]}"; do
    for nsrcs in "${Nsrcs[@]}"; do
        counter=$(( $counter + 1 ))
        if [ $counter -lt $min_config ] || [ $counter -gt $max_config ]
        then
                continue
        fi
        do_run $nsrcs $Nother $Nother $Nother $beam 
        echo $counter
    done
    for ntimes in "${Ntimes[@]}"; do
        counter=$(( $counter + 1 ))
        if [ $counter -lt $min_config ] || [ $counter -gt $max_config ]
        then
                continue
        fi
        do_run $Nother $ntimes $Nother $Nother $beam 
        echo $counter
    done
    for nfreqs in "${Nfreqs[@]}"; do
        counter=$(( $counter + 1 ))
        if [ $counter -lt $min_config ] || [ $counter -gt $max_config ]
        then
                continue
        fi
        do_run $Nother $Nother $nfreqs $Nother $beam
        echo $counter
    done
    for nbls in "${Nbls[@]}"; do
        counter=$(( $counter + 1 ))
        if [ $counter -lt $min_config ] || [ $counter -gt $max_config ]
        then
                continue
        fi
        do_run $Nother $Nother $Nother $nbls $beam
        echo $counter
    done
done

rm $mem_out
rm $time_out

## Try to clean up the scripts directory
ofilename='slurm-'$jobid'_'$task'.out'
mv $ofilename $dir1/$ofilename
