# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2024 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import copy

import numpy as np

import pyuvsim


def test_baseline_comparison():
    antenna1 = pyuvsim.Antenna("ant1", 1, np.array([0, 10, 0]), 1)
    antenna2 = pyuvsim.Antenna("ant2", 2, np.array([0, 20, 0]), 1)

    baseline12 = pyuvsim.Baseline(antenna1, antenna2)
    baseline21 = pyuvsim.Baseline(antenna2, antenna1)

    bl12_copy = copy.deepcopy(baseline12)

    assert baseline12 < baseline21
    assert baseline12 <= baseline21
    assert baseline12 <= baseline12
    assert baseline12 == bl12_copy

    assert baseline21 > baseline12
    assert baseline21 >= baseline21
    assert baseline21 >= baseline12
