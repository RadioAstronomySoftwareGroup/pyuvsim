# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np


class Baseline:
    """
    Defines a single pair of Antenna objects, the means for calculating their
    uvw coordinates, and comparing them
    """

    def __init__(self, antenna1, antenna2):
        self.antenna1 = antenna1
        self.antenna2 = antenna2
        self.enu = antenna2.pos_enu - antenna1.pos_enu
        # we're using the local alt/az frame so uvw is just enu
        self.uvw = self.enu

    def __eq__(self, other):
        return ((self.antenna1 == other.antenna1)
                and (self.antenna2 == other.antenna2)
                and np.allclose(self.enu.to('m').value, other.enu.to('m').value, atol=1e-3)
                and np.allclose(self.uvw.to('m').value, other.uvw.to('m').value, atol=1e-3))

    def __gt__(self, other):
        if self.antenna1 == other.antenna1:
            return self.antenna2 > other.antenna2
        return self.antenna1 > other.antenna1

    def __ge__(self, other):
        if self.antenna1 == other.antenna1:
            return self.antenna2 >= other.antenna2
        return self.antenna1 >= other.antenna1

    def __lt__(self, other):
        return not self.__ge__(other)

    def __le__(self, other):
        return not self.__gt__(other)
