# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Use the line profiler when requested."""

import atexit
import warnings
from inspect import isclass, isfunction
from itertools import chain

import pyradiosky as _pyradiosky

import pyuvsim as _pyuvsim

try:
    from . import mpi
except ImportError:
    mpi = None

try:
    from line_profiler import LineProfiler
except ImportError:

    def LineProfiler():  # noqa
        """Mock to fix imports."""
        return None


default_profile_funcs = [
    "interp",
    "get_beam_jones",
    "initialize_uvdata_from_params",
    "apply_beam",
    "make_visibility",
    "update_positions",
    "coherency_calc",
    "uvdata_to_task_iter",
    "run_uvdata_uvsim",
    "run_uvsim",
]

prof = None


# Note that enabling the profiler interferes with pytest-cov, so several lines here are
# marked with "nocover" though they are hit by tests.
# These "nocover" comments will need to remain until issue 179 on line_profiler is resolved.
# https://github.com/rkern/line_profiler/issues/179


def set_profiler(
    func_list=default_profile_funcs,
    rank=0,
    outfile_prefix="time_profile.out",
    dump_raw=False,
):
    """
    Apply a line profiler to the listed functions, wherever they appear in pyuvsim.

    Places a LineProfiler object in the module namespace, and registers its
    dumping/printing functions to run at the end. When the Python environment closes,
    the profiler functions print_stats (and dump_stats, if dump_raw is True) will
    execute, saving profiler data to file.

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

    """
    global prof

    if outfile_prefix.endswith(".out"):
        outfile_prefix = outfile_prefix[:-4]  # Strip extension

    outfile_name = outfile_prefix + ".out"

    # Can only set up profiling once per Python session.
    if prof is not None:  # pragma: nocover
        warnings.warn("Profiler already set. Returning now.")
        return

    prof = LineProfiler()
    if mpi is None or prof is None:
        raise ImportError(
            "You need mpi4py and line_profiler to use the "
            "profiling module. Install them both by running pip "
            "install pyuvsim[all]."
        )

    mpi.start_mpi()

    # Add module functions to profiler.
    mod_iter = chain(_pyuvsim.__dict__.values(), _pyradiosky.__dict__.values())
    for mod_it in mod_iter:
        if isfunction(mod_it) and mod_it.__name__ in func_list:
            prof.add_function(mod_it)
        if isclass(mod_it):
            for item in mod_it.__dict__.values():
                if isfunction(item) and item.__name__ in func_list:
                    prof.add_function(item)

    # Write out profiling report to file.
    if mpi.get_rank() == rank:
        # don't use a context manager here because it the file needs to stay
        # open for the various jobs to write to it. File closing is handled
        # explicitly by the `atexit.register` call.
        ofile = open(outfile_name, "w")  # noqa
        atexit.register(ofile.close)
        atexit.register(prof.print_stats, stream=ofile)
        if dump_raw:
            outfile_raw_name = outfile_prefix + ".lprof"
            atexit.register(prof.dump_stats, outfile_raw_name)
        prof.rank = rank  # Add "rank" as an attribute to the profiler.
        prof.meta_file = outfile_prefix + "_meta.out"

        prof.enable_by_count()


def get_profiler():  # pragma: nocover
    """Get the profiler."""
    return prof
