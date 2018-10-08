# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Summarize profiling results in a table, given a set of directories"""

from __future__ import absolute_import, division, print_function

import os
import sys
import glob
import re
import numpy as np
from numpy.lib.recfunctions import append_fields
import subprocess
from matplotlib.mlab import rec2csv

dirs = sys.argv[1:]

slurmids = []
Nall = []  # Each entry is a list of Nbl, Ntimes, Nchan, Nsrcs

diruse = []
for dr in dirs:
    # Get slurm jobID
    if not os.path.isdir(dr):
        continue
    try:
        slurmfile = glob.glob(dr + "/slurm-*_0.out")[0]
    except IndexError:
        continue

    diruse.append(dr)
    slurmid = re.search('-.*_', slurmfile).group(0)[1:-1]
    slurmids.append(slurmid + "_0.0")

# Future --- Set the comment to contain info on Nsrcs, etc.

p = subprocess.Popen('sacct --jobs=\"' + ",".join(slurmids) + '\" --format=\"JobID, Start, Elapsed, MaxRSS, NNodes, NTasks, NCPUS\"', shell=True, stdout=subprocess.PIPE)

table = np.genfromtxt(p.stdout, dtype=None, names=True, comments='--', encoding=None)
table['MaxRSS'] = map(lambda x: float(x[:-1]) * 1e3 / 1e9, table['MaxRSS'])

dt = table.dtype.descr
ind = dt.index(('MaxRSS', table.dtype['MaxRSS']))
dt[ind] = ('MaxRSS (GB)', 'f')
ind = dt.index(('NTasks', table.dtype['NTasks']))
dt[ind] = ('NCores', 'i')    # So there's no confusion with slurm tasks vs. UVTasks.
dt = np.dtype(dt)
table = table.astype(dt)

njobs = len(slurmids)
Nsrcs, Nchans, Ntimes, Nbls = np.zeros((4, njobs)).astype(int)
# They got resorted in the sacct task.
for di, dr in enumerate(diruse):
    sid = slurmids[di]
    ind = np.where(table['JobID'] == sid)
    for k in dr.split("_"):
        if k.endswith("src"):
            Nsrcs[ind] = int(k[:-3])
        if k.endswith("time"):
            Ntimes[ind] = int(k[:-4])
        if k.endswith('chan'):
            Nchans[ind] = int(k[:-4])
        if k == 'mwa':
            Nbls[ind] = 8128
        if k == 'triangle':
            Nbls[ind] = 3


def hms2sec(hms):
    h, m, s = map(float, hms.split(":"))
    return h * 60.**2 + m * 60. + s


runtime_sec = np.array(map(hms2sec, table['Elapsed']))
cores_per_node = table['NCores'] / table['NNodes']
ntasks = Nsrcs * Ntimes * Nbls * Nchans
timepertask = runtime_sec / ntasks


table = append_fields(table, ['CoresPerNode', 'Nbls', 'Ntimes', 'Nchan', 'Nsrc', 'Ntasks', 'Runtime_Seconds', 'RuntimePerTask_seconds'], [cores_per_node, Nbls, Ntimes, Nchans, Nsrcs, ntasks, runtime_sec, timepertask])
print table
print table.dtype

rec2csv(table, 'profiling_results_table.csv')
