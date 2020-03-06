# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
Use the line profiler when requested.
"""

import atexit
from inspect import isclass, isfunction
from itertools import chain

import pyuvsim as _pyuvsim
import pyradiosky as _pyradiosky
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
                         'apply_beam', 'make_visibility', 'update_positions', 'coherency_calc',
                         'uvdata_to_task_iter', 'run_uvdata_uvsim', 'run_uvsim']

prof = None

def set_profiler(func_list=default_profile_funcs, rank=0, outfile_prefix='time_profile.out',
                 dump_raw=False):
    """
    Applies a line profiler to the listed functions, wherever they appear in pyuvsim.

    Places a LineProfiler object in the module namespace, and registers its dumping/printing
    functions to run at the end.

    Parameters
    ----------
    func_list: list
        List of function names (strings) to profile.
        Defaults to ``profiling.default_profile_funcs``.
    rank: int, optional
        Which rank process should write out to file? (only one rank at a time will).
    outfile_prefix: str
        Filename prefix for printing profiling results.
            Human-readable line by line profiling goes to <outfile_prefix>.out
            LineStats data goes to <outfile_prefix>.lprof (if dump_raw)
            Axis sizes go to <outfile_prefix>_axes.npz
    dump_raw: bool
        Write out a pickled LineStats object to <outfile_name>.lprof (Default False)

    Sets
    ----
    prof : LineProfiler
        An instance of LineProfiler in the module namespace.
    exit functions:
        When the Python environment closes, the profiler functions
        print_stats (and dump_stats, if dump_raw is True) will execute,
        saving profiler data to file.

    """
    global prof
    prof = LineProfiler()
    if mpi is None or prof is None:  # pragma: no cover
        raise ImportError("You need mpi4py and line_profiler to use the "
                          "profiling module. Install them both by running pip "
                          "install pyuvsim[all].")

    mpi.start_mpi()

    # Add module functions to profiler.
    mod_iter = chain(_pyuvsim.__dict__.values(), _pyradiosky.__dict__.values())
    for mod_it in mod_iter:
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
        outfile_name = outfile_prefix
        if not outfile_name.endswith(".out"):
            outfile_name += '.out'
        else:
            outfile_prefix = outfile_prefix[:-3]    # Strip extension

        ofile = open(outfile_name, 'w')
        atexit.register(ofile.close)
        atexit.register(prof.print_stats, stream=ofile)
        if dump_raw:
            outfile_raw_name = outfile_prefix + ".lprof"
            atexit.register(prof.dump_stats, outfile_raw_name)
        setattr(prof, 'rank', rank)     # Add "rank" as an attribute to the profiler.
        setattr(prof, 'meta_file', outfile_prefix + '_meta.out')
        prof.enable_by_count()


def get_profiler():
    return prof
