# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import pyuvsim


def test_profiling():
    pyuvsim.profiling.set_profiler()
