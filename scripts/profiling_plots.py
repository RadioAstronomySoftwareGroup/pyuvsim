#!/bin/env python

# Make plots of profiling results under different constraints:

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

dat = np.genfromtxt(sys.argv[1], names=True, max_rows=2, delimiter=',')


fields = ['JobID', 'Start', 'MaxRSS_GB', 'NNodes', 'NProcs', 'Nbls', 'Ntimes', 'Nchan', 'Nsrc', 'Beam', 'Ntasks', 'Runtime_Second']
titles = [f.lower() for f in fields]
fmts = ['U10', 'U10', 'f8', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'U10', 'i4', 'f8']
dt = np.format_parser(fmts, dat.dtype.names, titles)

dat = np.genfromtxt(sys.argv[1], autostrip=True, dtype=dt.dtype, delimiter=',', skip_header=1)

Ncpus = np.unique(dat['NProcs'])
beams = np.unique(dat['Beam'])
Ntasks, inv, counts = np.unique(dat['Ntasks'], return_inverse=True, return_counts=True)
NNodes = np.unique(dat['NNodes'])
markers = ['.', 'o', '+', '<', '>', '*', 'v', '^', 'h', 'p']
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(len(Ncpus))]


params = ['Nbls', 'Ntimes', 'Nchan', 'Nsrc']
fig1, axes = plt.subplots(nrows=1, ncols=len(beams))
handles, labels = [], []
ymin, ymax = np.min(dat['Runtime_Seconds']), np.max(dat['Runtime_Seconds'])
ymax *= 1.20
for nni in range(len(NNodes)):
    condN = (dat['NNodes'] == NNodes[nni])
    for bi in range(len(beams)):
        condB = (dat['Beam'] == beams[bi])
        for nc in range(len(Ncpus)):
            inds = np.where(condN & condB & (dat['NProcs'] == Ncpus[nc]))
            if len(inds[0]) == 0:
                continue
            axes[bi].scatter(dat['Ntasks'][inds], dat['Runtime_Seconds'][inds], label="{:d} cores, {:d} Nodes".format(Ncpus[nc], NNodes[nni]), color=colors[nc], marker=markers[nni])

        axes[bi].set_title("{} beam".format(beams[bi]))
        axes[bi].set_ylabel("Runtime (seconds)")
        axes[bi].set_xlabel("Ntasks")
        axes[bi].set_ylim([ymin, ymax])


def f(m, c):
    return plt.plot([], [], marker=m, color=c, ls="none")[0]


colhandles = [f("s", colors[i]) for i in range(len(Ncpus))]
markhandles = [f(markers[i], "k") for i in range(len(NNodes))]

collabels = map('{:d} cpus'.format, Ncpus.tolist())
marklabels = map('{:d} nodes'.format, NNodes.tolist())
plt.figlegend(labels=collabels, handles=colhandles, loc='center right')
leg2 = plt.figlegend(labels=marklabels, handles=markhandles, loc='lower right')
plt.gca().add_artist(leg2)

fig2, axes = plt.subplots(nrows=1, ncols=len(beams))
ymin, ymax = np.min(dat['MaxRSS_GB']), np.max(dat['MaxRSS_GB'])
ymax *= 1.20
for nni in range(len(NNodes)):
    condN = (dat['NNodes'] == NNodes[nni])
    for bi in range(len(beams)):
        condB = (dat['Beam'] == beams[bi])
        for nc in range(len(Ncpus)):
            inds = np.where(condN & condB & (dat['NProcs'] == Ncpus[nc]))
            if len(inds[0]) == 0:
                continue
            axes[bi].scatter(dat['Ntasks'][inds], dat['MaxRSS_GB'][inds], label="{:d} cores, {:d} Nodes".format(Ncpus[nc], NNodes[nni]), color=colors[nc], marker=markers[nni])

        axes[bi].set_title("{} beam".format(beams[bi]))
        axes[bi].set_ylabel("MaxRSS (GB)")
        axes[bi].set_xlabel("Ntasks")
        axes[bi].set_ylim([ymin, ymax])
plt.figlegend(labels=collabels, handles=colhandles, loc='center right')
leg2 = plt.figlegend(labels=marklabels, handles=markhandles, loc='lower right')
plt.gca().add_artist(leg2)

plt.show()
