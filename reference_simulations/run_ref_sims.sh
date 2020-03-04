#!/bin/bash


# Run reference simulations and save slurmIDs

obsparam_files=( "$@" )
nfiles="$#"

echo ${obsparam_files[@]}
echo $nfiles

dirname=$(date +%Y-%m-%d)   # Name directory with today's date.
mkdir -p $dirname"/profiling"
mkdir -p $dirname"/slurm_out"
mkdir -p $dirname"/slurm_out"

prevdir=$(pwd -P)

ln -sf $prevdir/telescope_config $dirname/telescope_config
ln -sf $prevdir/catalog_files $dirname/catalog_files

for fn in "${obsparam_files[@]}"
do
        fn=$(basename $fn)
        echo $fn
        cp $fn $dirname/$fn
done


cd $dirname

sbatch -o slurm_out/pyuvsim-%A_%a.out --array=0-$(( $nfiles - 1  )) $prevdir/ref_sim_job.sh ${obsparam_files[@]}
