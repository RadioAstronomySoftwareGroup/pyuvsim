# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import pyuvsim
import nose.tools as nt


def test_profiling():
    pyuvsim.profiling.set_profiler()
