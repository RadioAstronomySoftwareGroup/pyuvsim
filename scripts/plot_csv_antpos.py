#!/bin/env python
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Plot antenna positions from a layout csv."""

import sys

import matplotlib.pyplot as plt

from pyuvsim.simsetup import _parse_layout_csv as parse_csv

antpos = parse_csv(sys.argv[1])

enu = antpos[["e", "n", "u"]]
antnames = antpos["number"]

plt.scatter(enu["e"], enu["n"])
for i, _ in enumerate(enu):
    plt.annotate(antnames[i], (enu["e"][i], enu["n"][i]))

plt.xlabel("East [m]")
plt.ylabel("North [m]")
plt.show()
