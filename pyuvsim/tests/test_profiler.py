# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import pyuvsim
try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler=False

def test_profiling():
    if LineProfiler:
        pyuvsim.profiling.set_profiler()
