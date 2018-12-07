# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np


class Telescope(object):
    """
    Defines the location and name of the observing site, and holds the list
    of beam objects used by the array.
    """

    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.location = telescope_location
        self.name = telescope_name

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list

    def __eq__(self, other):
        this_vector_loc = np.array([self.location.x.to('m').value,
                                    self.location.y.to('m').value,
                                    self.location.z.to('m').value])
        other_vector_loc = np.array([other.location.x.to('m').value,
                                     other.location.y.to('m').value,
                                     other.location.z.to('m').value])
        return ((np.allclose(this_vector_loc, other_vector_loc, atol=1e-3))
                and (self.beam_list == other.beam_list)
                and (self.name == other.name))
