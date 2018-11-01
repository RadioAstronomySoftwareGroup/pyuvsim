#!/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

# Run profiling jobs for multiple configurations
# and with different numbers of cores.

from __future__ import absolute_import, division, print_function

import numpy as np
import subprocess

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

for n in Ncores:
    cmd = ['sbatch', '-n ' + str(n), '--cpus-per-task=1', '--mem=' + mem, '--time=' + time,
           'batch_profile_job.sh', Nsrcs, Ntimes, Nfreqs, Nbls, beam]
    print(" ".join(cmd))
    results = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    print(results)
