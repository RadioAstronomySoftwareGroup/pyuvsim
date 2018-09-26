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
    alts = [90, 75, 30, 0, -10, -45, -90]
    # 0=North, 90=East
    azs = [0, 45, 90, 135, 180, 270, 350]

    zas = [0, 15, 60, 90, 100, 135, 180]
    # 0=East, 90=North
    beam_azs = [90, 45, 0, 315, 270, 180, 100]

    calc_za, calc_az = simutils.altaz_to_zenithangle_azimuth(np.deg2rad(alts),
                                                             np.deg2rad(azs))
    nt.assert_true(np.allclose(calc_za, np.deg2rad(zas)))
    nt.assert_true(np.allclose(calc_az, np.deg2rad(beam_azs)))


def test_single_altaz_to_za_az():
    alts = 90
    # 0=North, 90=East
    azs = 135

    zas = 0
    # 0=East, 90=North
    beam_azs = 315

    calc_za, calc_az = simutils.altaz_to_zenithangle_azimuth(np.deg2rad(alts),
                                                             np.deg2rad(azs))
    nt.assert_true(np.isclose(calc_za, np.deg2rad(zas)))
    nt.assert_true(np.isclose(calc_az, np.deg2rad(beam_azs)))


def test_za_az_to_altaz():
    # 5 degrees off zenith in the East direction
    zas = [0, 5, 45, 90, 120, 150, 180]
    # 0=East, 90=North
    azs = [0, 45, 90, 135, 180, 270, 350]

    alts = [90, 85, 45, 0, -30, -60, -90]
    # 0=North, 90=East
    astropy_azs = [90, 45, 0, 315, 270, 180, 100]

    calc_alt, calc_az = simutils.zenithangle_azimuth_to_altaz(np.deg2rad(zas),
                                                              np.deg2rad(azs))
    nt.assert_true(np.allclose(calc_alt, np.deg2rad(alts)))
    nt.assert_true(np.allclose(calc_az, np.deg2rad(astropy_azs)))


def test_za_az_to_altaz():
    # 5 degrees off zenith in the East direction
    zas = 5
    # 0=East, 90=North
    azs = 180

    alts = 85
    # 0=North, 90=East
    astropy_azs = 270

    calc_alt, calc_az = simutils.zenithangle_azimuth_to_altaz(np.deg2rad(zas),
                                                              np.deg2rad(azs))
    nt.assert_true(np.isclose(calc_alt, np.deg2rad(alts)))
    nt.assert_true(np.isclose(calc_az, np.deg2rad(astropy_azs)))


def test_altaz_za_az_errors():
    nt.assert_raises(ValueError, simutils.altaz_to_zenithangle_azimuth, 0, [0, np.pi / 2])
    nt.assert_raises(ValueError, simutils.zenithangle_azimuth_to_altaz, 0, [0, np.pi / 2])
