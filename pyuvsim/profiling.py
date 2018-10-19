# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
Use the line profiler when requested.
"""

from __future__ import absolute_import, division, print_function

from .mpi import start_mpi, get_rank
from inspect import isclass, isfunction
import sys
import atexit

import pyuvsim as _pyuvsim

try:
    from line_profiler import LineProfiler
except ImportError:   # pragma: no cover
    def LineProfiler():
        return None

PY3 = sys.version_info[0] == 3

if PY3:
    import builtins
else:
    import __builtin__ as builtins

prof = LineProfiler()

default_profile_funcs = ['interp', 'get_beam_jones', 'initialize_uvdata_from_params', 'uvdata_to_telescope_config', 'uvdata_to_config_file', 'coherency_calc', 'alt_az_calc', 'apply_beam', 'make_visibility', 'uvdata_to_task_list', 'initialize_uvdata', 'run_uvsim']


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
    if prof is None:   # pragma: no cover
        raise ImportError("line_profiler module required to use profiling tools.")

    # Add module functions to profiler.
    for mod_it in _pyuvsim.__dict__.values():
        if isfunction(mod_it):
            if mod_it.func_name in func_list:
                prof.add_function(mod_it)
        if isclass(mod_it):
            for item in mod_it.__dict__.values():
                if isfunction(item):
                    if item.func_name in func_list:
                        prof.add_function(item)

    builtins.__dict__['time_profiler'] = prof       # So it can accessed by tests
    # Write out profiling report to file.
        if get_rank() == rank:
            ofile = file(outfile_name, 'w')
            atexit.register(ofile.close)
            atexit.register(prof.print_stats, stream=ofile)
            if dump_raw:
                outfile_raw_name = outfile_name+".lprof"
                atexit.register(prof.dump_stats, outfile_raw_name)

            prof.enable_by_count()
