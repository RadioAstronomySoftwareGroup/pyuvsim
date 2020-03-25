
import numpy as np
import line_profiler as lp
import sys
import os

# Estimate scalings on different axes from profiling data.

# NB -- This script cannot tell the difference between lines that are
#          hit once and axes with length 1.
#       It is only useful for jobs where the number of MPI processes is smaller
#       than each simulation axis.


def _func_times(timings, Nlist, dt=1e-6):
    outarr = np.zeros(len(Nlist))
    for key, values in timings.items():
        if key[2] not in funcs_to_check:
            continue
        for lnum, nhits, time in values:
            if nhits in Nlist:
                ind = Nlist.index(nhits)
                outarr[ind] += dt * time / nhits    # Time per hit.
    return outarr


# Only include functions that are called in loops.
funcs_to_check = ['interp', 'get_beam_jones', 'apply_beam', 'make_visibility',
                  'uvdata_to_task_iter', 'update_positions', 'coherency_calc']
profname = 'profdata/time_profile.lprof'
axesname = 'profdata/time_profile_axes.npz'

# Prepend path
basepath = sys.argv[1]
profname = os.path.join(basepath, profname)
axesname = os.path.join(basepath, axesname)

lstat = lp.load_stats(profname)
axes_npz = np.load(axesname)
axes = {k: axes_npz[k][0] for k in axes_npz.keys()}

# Set up combinations of axes.
Naxes = len(axes)

combos = [
    ('Ntimes',),
    ('Ntimes', 'Nfreqs'),
    ('Ntimes', 'Nfreqs', 'Nbls')
]

if not axes['Nsrcs_loc'] == 1:
    combos.append(('Ntimes', 'Nfreqs', 'Nbls', 'Nsrcs_loc'))

Nlist = []
for comb in combos:
    length = 1
    for key in comb:
        length *= axes[key]
    Nlist.append(length)

results = _func_times(lstat.timings, Nlist, dt=lstat.unit)
print(Nlist)
print(dict(zip(combos, results)))
print(np.sum(np.array(Nlist) * np.array(results)))
