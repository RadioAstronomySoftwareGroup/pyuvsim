# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle

import pyuvsim


def test_tee_ra_loop():
    time = Time(2457458.1739, scale='utc', format='jd')
    tee_ra = Angle(np.pi / 4., unit='rad')  # rad
    cirs_ra = pyuvsim.utils.tee_to_cirs_ra(tee_ra, time)
    new_tee_ra = pyuvsim.utils.cirs_to_tee_ra(cirs_ra, time)
    nt.assert_equal(new_tee_ra, tee_ra)
