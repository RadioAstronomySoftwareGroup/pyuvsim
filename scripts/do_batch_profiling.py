#!/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

# Run profiling jobs for multiple configurations
# and with different numbers of cores.

from __future__ import absolute_import, division, print_function

import numpy as np
import subprocess

# Nsrcs = [5, 10, 20, 30]
# Ntimes = [1, 5, 10, 15]
# Nfreqs = [1, 5, 10, 15]
# Nbls = [10, 20, 50]
# beam = ['hera', 'uniform']

Nsrcs = np.logspace(2, 4, 20).astype(int).tolist()
Ntimes = np.linspace(1, 99, 20).astype(int).tolist()
Nfreqs = [1]
Nbls = np.logspace(2, np.log10(3000), 20).astype(int).tolist()
beam = ['uniform', 'hera']

meshgrid = False

Ncores = [8, 16, 32, 64, 128]

Nsrcs = ','.join(map(str, Nsrcs))
Ntimes = ','.join(map(str, Ntimes))
Nfreqs = ','.join(map(str, Nfreqs))
Nbls = ','.join(map(str, Nbls))
beam = ','.join(map(str, beam))

mem = '80G'
time = '48:00:00'
print("Meshgrid?: ", meshgrid)
for n in Ncores:
    script = 'batch_profile_job.sh'
    if meshgrid:
        script = 'batch_profile_meshgrid_job.sh'
    cmd = ['sbatch', '-n ' + str(n), '--cpus-per-task=1', '--mem=' + mem, '--time=' + time,
           script, Nsrcs, Ntimes, Nfreqs, Nbls, beam]
    print(" ".join(cmd))
    results = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    print(results)
