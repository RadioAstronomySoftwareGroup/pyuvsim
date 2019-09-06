#!/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

# Run profiling jobs for multiple configurations
# and with different numbers of cores.

import subprocess

import numpy as np

# Nsrcs = [5, 10, 20, 30]
# Ntimes = [1, 5, 10, 15]
# Nfreqs = [1, 5, 10, 15]
# Nbls = [10, 20, 50]
# beam = ['hera', 'uniform']

Nsrcs = np.logspace(4, 5, 10).astype(int).tolist()
#Ntimes = np.logspace(np.log10(500), np.log10(4000), 10).astype(int).tolist()
Ntimes = np.linspace(200, 500, 10).astype(int).tolist()
Nfreqs = np.linspace(200, 500, 10).astype(int).tolist()
Nbls = np.logspace(3, np.log10(3900), 10).astype(int).tolist()
beam = ['uniform']#, 'hera']
#beam = ['hera']

#Nsrcs = []
#Nsrcs = [] #np.logspace(4,6,20).astype(int).tolist()[10:]
#Ntimes = []
#Nfreqs = [1024]*5
#Nbls = []
#beam = ['uniform']#, 'hera']

meshgrid = False

#Npus = [8, 16, 32, 64, 128]
#Npus = [40, 60, 80, 100]
#Npus = [20,40]

#Ncpus = [40, 60, 80]
Ncpus = [20, 40, 80]

Ncpus.reverse()

# Submit largest jobs earlier
Nsrcs.reverse()
Ntimes.reverse()
Nfreqs.reverse()
Nbls.reverse()

if meshgrid:
    Nconfigs = len(Nsrcs) * len(Ntimes) * len(Nfreqs) * len(Nbls) * len(beam)
else:
    Nconfigs = (len(Nsrcs) + len(Ntimes) + len(Nfreqs) + len(Nbls)) * len(beam)


Nsrcs = ','.join(map(str, Nsrcs))
Ntimes = ','.join(map(str, Ntimes))
Nfreqs = ','.join(map(str, Nfreqs))
Nbls = ','.join(map(str, Nbls))
beam = ','.join(map(str, beam))

mem = '20G'
time = '48:00:00'
print("Meshgrid?: ", meshgrid)
array_size = np.min([Nconfigs, 100])
cpus_per_task=4
Npus = [int(n/cpus_per_task) for n in Ncpus]
for n in Npus:
    script = 'batch_profile_job.sh'
    if meshgrid:
        script = 'batch_profile_meshgrid_job.sh'
    cmd = ['sbatch', '-n ' + str(n), '--array=0-' + str(array_size), '--cpus-per-task='+str(cpus_per_task), '--mem=' + mem, '--time=' + time,
           script, Nsrcs, Ntimes, Nfreqs, Nbls, beam, str(Nconfigs)]
    print(" ".join(cmd))
    results = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    print(results)
