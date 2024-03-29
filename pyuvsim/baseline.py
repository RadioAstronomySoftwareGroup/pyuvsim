# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Definition of Baseline objects, to describe a pair of :class:`~Antenna` objects."""

import numpy as np


class Baseline:
    """
    Defines a single pair of :class:`~Antenna` objects with useful methods.

    Defines the means for calculating their uvw coordinates, and comparing them.

    Parameters
    ----------
    antenna1, antenna2 : :class:`pyuvsim.Antenna`
        The two antennas that make up this baseline.

    Attributes
    ----------
    antenna1, antenna2 : :class:`pyuvsim.Antenna`
        The two antennas that make up this baseline.
    enu : array_like of float
        The vector pointing from antenna1 to antenna2 in ENU coordinates
        (relative to the telescope location.)
    uvw : array_like of float
        The uvw coordinate of this baseline. Identical to `enu` since we work in the
        local alt/az frame.

    """

    def __init__(self, antenna1, antenna2):
        self.antenna1 = antenna1
        self.antenna2 = antenna2
        self.enu = antenna2.pos_enu - antenna1.pos_enu
        # we're using the local alt/az frame so uvw is just enu
        self.uvw = self.enu

    def __eq__(self, other):
        """Define baseline equality."""
        return (
            (self.antenna1 == other.antenna1)
            and (self.antenna2 == other.antenna2)
            and np.allclose(self.enu.to("m").value, other.enu.to("m").value, atol=1e-3)
            and np.allclose(self.uvw.to("m").value, other.uvw.to("m").value, atol=1e-3)
        )

    def __gt__(self, other):
        """Check if baseline is greater than another based on antenna numbers."""
        if self.antenna1 == other.antenna1:
            return self.antenna2 > other.antenna2
        return self.antenna1 > other.antenna1

    def __ge__(self, other):
        """Check if baseline is greater than or equal to another based on antenna numbers."""
        if self.antenna1 == other.antenna1:
            return self.antenna2 >= other.antenna2
        return self.antenna1 >= other.antenna1

    def __lt__(self, other):
        """Check if baseline is less than another based on antenna numbers."""
        return not self.__ge__(other)

    def __le__(self, other):
        """Check if baseline is less than or equal to another based on antenna numbers."""
        return not self.__gt__(other)
