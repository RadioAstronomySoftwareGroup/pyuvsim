# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import shutil
import atexit

import numpy as np

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = False


def profdata_dir_setup():
    # The pytest temporary test directory doesn't work with the profiling outputs
    # because they are only created after the exit functions are run, and pytest
    # deletes its temporary directory before then.
    # As a workaround, make a temporary directory for profiling data and register
    # its deletion in the exit functions.
    outpath = os.path.join(SIM_DATA_PATH, 'temporary_profiling_data')
    os.mkdir(outpath)
    atexit.register(shutil.rmtree, outpath)
    return outpath


def test_profiler():
    if LineProfiler:
        outpath = profdata_dir_setup()
        testprof_fname = os.path.join(outpath, 'time_profile')
        pyuvsim.profiling.set_profiler(outfile_prefix=testprof_fname, dump_raw=True)
        param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
        pyuvsim.uvsim.run_uvsim(param_filename, return_uv=True)

        time_profiler = pyuvsim.profiling.get_profiler()
        assert isinstance(time_profiler, LineProfiler)
        lstats = time_profiler.get_stats()

        assert len(lstats.timings) != 0
        func_names = [k[2] for k in lstats.timings.keys()]
        assert np.unique(func_names).tolist() == sorted(pyuvsim.profiling.default_profile_funcs)
        time_profiler.disable()
