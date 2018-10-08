#!/bin/env python

# Run profiling jobs for multiple configurations
# and with different numbers of cores.

import numpy as np
import subprocess


Nsrcs = [5, 10, 20]
Ntimes = [1, 5, 10]
Nfreqs = [1, 5, 10]
Nbls = [3, 10]
beam = ['uniform', 'hera']

Nsrcs, Ntimes, Nfreqs, Nbls, beam = map(np.ndarray.flatten, np.meshgrid(Nsrcs, Ntimes, Nfreqs, Nbls, beam))

Nconfigs = Nsrcs.size

Ncores = [8, 16, 32, 64]

mem='40G'
#time='48:00:00'
time='00:02:00'

sids_out = open('slurm_ids.out', 'a+')

Nconfigs = 1

for n in Ncores:
    for i in range(Nconfigs):
        cmd = ['sbatch', '-n '+str(n), '--cpus-per-task=1', '--mem='+mem, '--time='+time,
               '-J pyuvsim_profile',
               'batch_profile_job.sh',
                str(Nsrcs[i]),
                str(Ntimes[i]),
                str(Nfreqs[i]),
                str(Nbls[i]),
                str(beam[i])]
        print " ".join(cmd)
        results = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        print results
        slurm_id = results.strip().split(' ')[-1]
        sids_out.write(slurm_id + '\n')

        
         

