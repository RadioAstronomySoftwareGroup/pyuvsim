# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle

import pyuvsim.utils as simutils


def test_tee_ra_loop():
    time = Time(2457458.1739, scale='utc', format='jd')
    tee_ra = Angle(np.pi / 4., unit='rad')  # rad
    cirs_ra = simutils.tee_to_cirs_ra(tee_ra, time)
    new_tee_ra = simutils.cirs_to_tee_ra(cirs_ra, time)
    nt.assert_equal(new_tee_ra, tee_ra)


def test_altaz_to_za_az():
    # 15 degrees off zenith in the North direction
    alt_az_1 = (75, 0)
    za_az_1 = (15, 90)

    calc_za_az1 = simutils.altaz_to_zenithangle_azimuth(*np.deg2rad(np.array(alt_az_1)))
    nt.assert_true(np.allclose(calc_za_az1, np.deg2rad(np.array(za_az_1))))


def test_altaz_to_za_az_below_horizon():
    # below the horizon in the North direction
    alt_az_3 = (-10, 0)
    za_az_3 = (100, 90)

    calc_za_az3 = simutils.altaz_to_zenithangle_azimuth(*np.deg2rad(np.array(alt_az_3)))
    nt.assert_true(np.allclose(calc_za_az3, np.deg2rad(np.array(za_az_3))))


def test_za_az_to_altaz():
    # 5 degrees off zenith in the East direction
    za_az_2 = (5, 0)
    alt_az_2 = (85, 90)

    calc_alt_az_2 = simutils.zenithangle_azimuth_to_altaz(*np.deg2rad(np.array(za_az_2)))
    nt.assert_true(np.allclose(calc_alt_az_2, np.deg2rad(np.array(alt_az_2))))
