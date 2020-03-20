#!/bin/bash

# Run reference simulations and save slurmIDs

obsparam_files=( "$@" )
nfiles="$#"

v1files=()
v2files=()

for fn in "${obsparam_files[@]}"
do
        fdir=$(dirname $fn)
        fnme=$(basename $fn)
        if [[ $fdir = "first_reference_simulations" ]]
        then
                v1files+=($fnme)
        fi
        if [[ $fdir = "second_reference_simulations" ]]
        then
                v2files+=($fnme)
        fi
        if [[ "$fnme" == *"1.4"*  && ! -f 'first_reference_simulations/catalog_files/gleam.vot' ]]
        then
            echo -e "Need to download the GLEAM catalog to run refsim 1.4."\
                    "\nRun \"python get_gleam.py\" to do this automatically, and then rerun run_ref_sims.sh"
            exit
        fi
done

prevdir=$(pwd -P)
dirname=$(date +%Y-%m-%d)       # Name directory with today's date.

if [[ ${#v1files[@]} -gt 0 ]]
then
    fulldir=$dirname"_v1"
    echo $fulldir
    sourcedir=$prevdir"/first_reference_simulations"
    mkdir -p $fulldir"/profiling"
    mkdir -p $fulldir"/slurm_out"
    ln -sf $sourcedir/telescope_config $fulldir/telescope_config
    ln -sf $sourcedir/catalog_files $fulldir/catalog_files
    for fn in "${v1files[@]}"
    do
            cp 'first_reference_simulations/'$fn $fulldir"/"$(basename $fn)
    done
    cd $fulldir
    sbatch -o slurm_out/pyuvsim-%A_%a.out --array=0-$(( ${#v1files[@]} - 1  )) $prevdir/jobscript.sh ${v1files[@]}
    cd $prevdir
fi

if [[ ${#v2files[@]} -gt 0 ]]
then
    fulldir=$dirname"_v2"
    echo $fulldir
    sourcedir=$prevdir"/second_reference_simulations"
    mkdir -p $fulldir"/profiling"
    mkdir -p $fulldir"/slurm_out"
    ln -sf $sourcedir/telescope_config $fulldir/telescope_config
    ln -sf $sourcedir/catalog_files $fulldir/catalog_files
    for fn in "${v1files[@]}"
    do
            cp 'second_reference_simulations/'$fn $fulldir"/"$(basename $fn)
    done
    cd $fulldir
    sbatch -o slurm_out/pyuvsim-%A_%a.out --array=0-$(( ${#v2files[@]} - 1  )) $prevdir/jobscript.sh ${v2files[@]}
    cd $prevdir
fi
