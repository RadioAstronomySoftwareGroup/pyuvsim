# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
Use the line profiler when requested.

If the set_profiler method is called, the profile decorators
throughout the code will do line profiling on all instance methods.
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


def set_profiler(func_list=default_profile_funcs, rank=0, outfile_name='time_profile.out'):
    start_mpi()
    if prof is None:
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

    # Write out profiling report to file.
        if get_rank() == rank:
            ofile = file(outfile_name, 'w')
            builtins.__dict__['profiling_vars'] = [prof, ofile]
            atexit.register(ofile.close)
            atexit.register(prof.print_stats, stream=ofile)
            prof.enable_by_count()
