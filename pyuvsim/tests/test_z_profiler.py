# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os

import numpy as np

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = False

testprof_fname = '/dev/null'


def test_profiler():
    if LineProfiler:
        pyuvsim.profiling.set_profiler(outfile_name=testprof_fname, dump_raw='/dev/null')
        param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
        pyuvsim.uvsim.run_uvsim(param_filename, return_uv=True)

        time_profiler = pyuvsim.profiling.get_profiler()
        assert isinstance(time_profiler, LineProfiler)
        lstats = time_profiler.get_stats()

        assert len(lstats.timings) != 0
        func_names = [k[2] for k in lstats.timings.keys()]
        assert np.unique(func_names).tolist() == sorted(pyuvsim.profiling.default_profile_funcs)
        time_profiler.disable()
