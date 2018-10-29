# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import sys
import yaml
import numpy as np
import nose.tools as nt

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = False

PY3 = sys.version_info[0] == 3

if PY3:
    import builtins
    import _pickle as pkl
else:
    import __builtin__ as builtins
    import cPickle as pkl

testprof_fname = '/dev/null'


def test_profiling():
    if LineProfiler:
        pyuvsim.profiling.set_profiler(outfile_name=testprof_fname)
        param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
        with open(param_filename, 'r') as pfile:
            params_dict = yaml.safe_load(pfile)
        uv_in, beam_list, beam_dict, beam_ids = pyuvsim.simsetup.initialize_uvdata_from_params(param_filename)
        beam_list[0] = pyuvsim.analyticbeam.AnalyticBeam('uniform')  # Replace the one that's a HERA beam
        catalog = os.path.join(SIM_DATA_PATH, params_dict['sources']['catalog'])
        uv_out = pyuvsim.run_uvsim(uv_in, beam_list, catalog_file=catalog, beam_dict=beam_dict)

        # The profiler is in the builtins as time_profiler. Check results.
        lstats = time_profiler.get_stats()
        del builtins.__dict__['time_profiler']  # Remove from builtins

        nt.assert_false(len(lstats.timings) == 0)
        func_names = [k[2] for k in lstats.timings.keys()]
        nt.assert_equal(np.unique(func_names).tolist(), sorted(pyuvsim.profiling.default_profile_funcs))
