# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
Use the line profiler when requested.
"""

from __future__ import absolute_import, division, print_function

from .mpi import start_mpi, get_rank
from inspect import isclass, isfunction
import atexit
import six

import pyuvsim as _pyuvsim

try:
    from line_profiler import LineProfiler
except ImportError:   # pragma: no cover
    def LineProfiler():
        return None


default_profile_funcs = ['interp', 'get_beam_jones', 'initialize_uvdata_from_params', 'coherency_calc', 'alt_az_calc', 'apply_beam', 'make_visibility', 'uvdata_to_task_iter', 'init_uvdata_out', 'run_uvsim']


def set_profiler(func_list=default_profile_funcs, rank=0, outfile_name='time_profile.out', dump_raw=False):
    """
    Applies a line profiler to the listed functions, wherever they appear in pyuvsim.

    Places a LineProfiler object in the builtins list, and registers its dumping/printing functions to run
    at the end.

    Args:
        func_list: list of function names (strings). Defaults to the list above.
        rank: (int) Which rank process should write out to file? (only one rank at a time will). Default 0
        outfile_name: Filename for printing profiling results.
        dump_raw: Write out a pickled LineStats object to <outfile_name>.lprof (Default False)
    """
    start_mpi()
    global prof
    prof = LineProfiler()
    if prof is None:   # pragma: no cover
        raise ImportError("line_profiler module required to use profiling tools.")

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
    if get_rank() == rank:
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
