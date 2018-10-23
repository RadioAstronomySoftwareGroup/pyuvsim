#!/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

# Run profiling jobs for multiple configurations
# and with different numbers of cores.

from __future__ import absolute_import, division, print_function

import numpy as np
import subprocess
<<<<<<< HEAD
from pyuvsim.utils import check_file_exists_and_increment
=======
>>>>>>> 8f8d834... Loop within bash jobscript, need memory_profiler

Nsrcs = [5, 10, 20]
Ntimes = [1, 5, 10, 15]
Nfreqs = [1, 5, 10, 15]
Nbls = [10, 20, 50]
beam = ['uniform', 'hera']

Ncores = [8, 16, 32, 64]

Nsrcs = ','.join(map(str, Nsrcs))
Ntimes = ','.join(map(str, Ntimes))
Nfreqs = ','.join(map(str, Nfreqs))
Nbls = ','.join(map(str, Nbls))
beam = ','.join(map(str, beam))

mem = '40G'
time = '48:00:00'

output = subprocess.check_output('git branch | grep \\* | cut -d \' \' -f2', shell=True)
git_branch = output.strip()

fname = check_file_exists_and_increment(git_branch + '_slurm_ids.out')

sids_out = open(fname, 'w')
sids_out.write('Nsrcs, Ntimes, Nfreqs, Nbls, beam, slurm_id\n')
for n in Ncores:
    #cmd = ['sbatch', '-n ' + str(n), '--cpus-per-task=1', '--mem=' + mem, '--time=' + time,
    #       'batch_profile_job.sh', Nsrcs, Ntimes, Nfreqs, Nbls, beam]
    cmd = ['./batch_profile_job.sh', Nsrcs, Ntimes, Nfreqs, Nbls, beam]
    print(" ".join(cmd))
    results = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    print(results)
#    slurm_id = results.strip().split(' ')[-1]
#    parms = [
#        str(Nsrcs[i]),
#        str(Ntimes[i]),
#        str(Nfreqs[i]),
#        str(Nbls[i]),
#        str(beam[i])]
#    sids_out.write(','.join(parms) + ',' + slurm_id + '\n')
