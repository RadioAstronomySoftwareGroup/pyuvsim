#!/bin/bash

#SBATCH -J pyuvsim_benchmark
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1-1
#SBATCH --ntasks=10
#SBATCH -m cyclic

srun --mpi=pmi2 python ../scripts/run_param_pyuvsim.py configdir/obsparam_benchmark.yaml --profile='profdata/time_profile' --raw_profile