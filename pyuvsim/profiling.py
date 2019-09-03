# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
Use the line profiler when requested.
"""

from __future__ import absolute_import, division, print_function

import atexit
from inspect import isclass, isfunction

import six

import pyuvsim as _pyuvsim
try:
    from . import mpi
except ImportError:
    mpi = None

try:
    from line_profiler import LineProfiler
except ImportError:  # pragma: no cover
    def LineProfiler():
        return None

default_profile_funcs = ['interp', 'get_beam_jones', 'initialize_uvdata_from_params',
                         'coherency_calc', 'update_positions', 'apply_beam', 'make_visibility',
                         'uvdata_to_task_iter', 'run_uvsim']


def set_profiler(func_list=default_profile_funcs, rank=0, outfile_name='time_profile.out',
                 dump_raw=False):
    """
    Applies a line profiler to the listed functions, wherever they appear in pyuvsim.

    Places a LineProfiler object in the builtins list, and registers its dumping/printing
    functions to run at the end.

    Args:
        func_list: list of function names (strings). Defaults to the list above.
        rank: int, optional
            Which rank process should write out to file? (only one rank at a time will).
        outfile_name: Filename for printing profiling results.
        dump_raw: Write out a pickled LineStats object to <outfile_name>.lprof (Default False)
    """
    if mpi is None:
        raise ImportError("You need mpi4py to use the profiling module. "
                          "Install it by running pip install pyuvsim[sim] "
                          "or pip install pyuvsim[all] if you also want h5py "
                          "and line_profiler installed.")

    mpi.start_mpi()
    global prof
    prof = LineProfiler()
    if prof is None:  # pragma: no cover
        raise ImportError("You need line_profiler to use the profiling module. "
                          "Install it by running pip install pyuvsim[all]. ")

    # Add module functions to profiler.
    for mod_it in _pyuvsim.__dict__.values():
        if isfunction(mod_it):
            if mod_it.__name__ in func_list:
                prof.add_function(mod_it)
        if isclass(mod_it):
            for item in mod_it.__dict__.values():
                if isfunction(item):
                    if item.__name__ in func_list:
                        prof.add_function(item)

    # Write out profiling report to file.
    if mpi.get_rank() == rank:
        ofile = open(outfile_name, 'w')
        atexit.register(ofile.close)
        atexit.register(prof.print_stats, stream=ofile)
        if dump_raw:
            outfile_raw_name = outfile_name + ".lprof"
            if isinstance(dump_raw, six.string_types):
                outfile_raw_name = dump_raw
            atexit.register(prof.dump_stats, outfile_raw_name)

        prof.enable_by_count()


def get_profiler():
    return prof
