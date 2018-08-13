# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import numpy as np

from . import profiling


@profile
class Telescope(object):
    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.telescope_location = telescope_location
        self.telescope_name = telescope_name

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list

    def __eq__(self, other):
        return ((np.allclose(self.telescope_location.to('m').value, other.telescope_location.to("m").value, atol=1e-3))
                and (self.beam_list == other.beam_list)
                and (self.telescope_name == other.telescope_name))
