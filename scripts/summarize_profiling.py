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
import subprocess

slurmids = []
Nall = []  # Each entry is a list of Nbl, Ntimes, Nchan, Nsrcs

fname = sys.argv[1]

with open(fname, 'a+') as fhandle:
    header = fhandle.readline()

header = [h.strip().upper() for h in header.split(',')]
dt = np.format_parser(['i4', 'i4', 'i4', 'i4', 'U8', 'i4'],
                      ['Nsrcs', 'Ntimes', 'Nfreqs', 'Nbls', 'beam', 'slurm_id'], header)

filedat = np.genfromtxt(fname, autostrip=True, skip_header=1,
                        delimiter=',', dtype=dt.dtype)

slurmids = filedat['slurm_id'].astype(str)
slurmids = [sid + '_0.0' for sid in slurmids]

p = subprocess.Popen('sacct --jobs=\"' + ",".join(slurmids) + '\" --format=\"JobID, Start, Elapsed, MaxRSS, NNodes, NTasks, NCPUS\"', shell=True, stdout=subprocess.PIPE)

table = np.genfromtxt(p.stdout, dtype=None, names=True, comments='--', encoding=None)
table['MaxRSS'] = map(lambda x: float(x[:-1]) * 1e3 / 1e9, table['MaxRSS'])

dt = table.dtype.descr
ind = dt.index(('MaxRSS', table.dtype['MaxRSS']))
dt[ind] = ('MaxRSS (GB)', 'f')
ind = dt.index(('NTasks', table.dtype['NTasks']))
dt[ind] = ('NProcs', 'i')    # So there's no confusion with slurm tasks vs. UVTasks.
dt = np.dtype(dt)
table = table.astype(dt)

njobs = len(slurmids)
Nsrcs, Nchans, Ntimes, Nbls = np.zeros((4, njobs)).astype(int)
beam_type = []
# They got resorted in the sacct task.
for ent in filedat:
    sid = str(ent['slurm_id']) + "_0.0"
    ind = np.where(table['JobID'] == sid)
    Nsrcs[ind] = ent['Nsrcs']
    Nchans[ind] = ent['Nfreqs']
    Ntimes[ind] = ent['Ntimes']
    Nbls[ind] = ent['Nbls']
    beam_type.append(ent['beam'])


def hms2sec(hms):
    h, m, s = map(float, hms.split(":"))
    return h * 60.**2 + m * 60. + s


runtime_sec = np.array(map(hms2sec, table['Elapsed']))
cores_per_node = table['NProcs'] / table['NNodes']
ntasks = Nsrcs * Ntimes * Nbls * Nchans

table = append_fields(table, ['CoresPerNode', 'Nbls', 'Ntimes', 'Nchan', 'Nsrc', 'Beam', 'Ntasks', 'Runtime_Seconds'], [cores_per_node, Nbls, Ntimes, Nchans, Nsrcs, beam_type, ntasks, runtime_sec], usemask=False)

rec2csv(table, 'profiling_results_table.csv')
