#!/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

from pyuvdata import UVData
import sys
import matplotlib.pyplot as plt
from pyuvsim.simsetup import _parse_layout_csv as parse_csv

antpos = parse_csv(sys.argv[1])

enu = antpos[['e', 'n', 'u']]
antnames = antpos['number']

plt.scatter(enu['e'], enu['n'])
for i, t in enumerate(enu):
    plt.annotate(antnames[i], (enu['e'][i], enu['n'][i]))

plt.xlabel('East [m]')
plt.ylabel('North [m]')
plt.show()
