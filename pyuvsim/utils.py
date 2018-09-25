# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import time as pytime
import sys
import numpy as np
from astropy import _erfa as erfa
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.coordinates.builtin_frames.utils import get_jd12

from . import version as simversion


def get_version_string():
    version_string = ('Simulated with pyuvsim version: ' + simversion.version + '.')
    if simversion.git_hash is not '':
        version_string += ('  Git origin: ' + simversion.git_origin
                           + '.  Git hash: ' + simversion.git_hash
                           + '.  Git branch: ' + simversion.git_branch
                           + '.  Git description: ' + simversion.git_description + '.')
    return version_string


class progsteps:
    """
        Similar to progress bar. Prints a percentage of completion.
        For when running in batch and progress bar doesn't work well.
    """

    def __init__(self, maxval=None):
        self.t0 = pytime.time()
        if maxval is None:
            raise ValueError("Maximum value is needed.")
        self.maxval = float(maxval)
        step = self.maxval * 0.01
        if step < 1.0:
            step = 1
        self.step = step

    def update(self, count):
        if count % self.step == 0:
            print("{:0.2f}% completed. {:0.3f} minutes elapsed \n".format(
                  (count / self.maxval) * 100., (pytime.time() - self.t0) / 60.))
            sys.stdout.flush()

    def finish(self):
        self.update(self.maxval)


# The frame radio astronomers call the apparent or current epoch is the
# "true equator & equinox" frame, notated E_upsilon in the USNO circular
# astropy doesn't have this frame but it's pretty easy to adapt the CIRS frame
# by modifying the ra to reflect the difference between
# GAST (Grenwich Apparent Sidereal Time) and the earth rotation angle (theta)
def tee_to_cirs_ra(tee_ra, time):
    era = erfa.era00(*get_jd12(time, 'ut1'))
    theta_earth = Angle(era, unit='rad')

    assert(isinstance(time, Time))
    assert(isinstance(tee_ra, Angle))
    gast = time.sidereal_time('apparent', longitude=0)
    cirs_ra = tee_ra - (gast - theta_earth)
    return cirs_ra


def cirs_to_tee_ra(cirs_ra, time):
    era = erfa.era00(*get_jd12(time, 'ut1'))
    theta_earth = Angle(era, unit='rad')

    assert(isinstance(time, Time))
    assert(isinstance(cirs_ra, Angle))
    gast = time.sidereal_time('apparent', longitude=0)
    tee_ra = cirs_ra + (gast - theta_earth)
    return tee_ra


def altaz_to_zenithangle_azimuth(altitude, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Args:
        altitude: in radians
        azimuth: in radians in astropy convention: East of North (N=0, E=90 degrees)

    Returns:
        zenith_angle in radians
        azimuth in radians in uvbeam convention: North of East(East=0, North=90 degrees)
    """
    zenith_angle = np.pi / 2 - np.array(altitude)
    new_azimuth = np.pi / 2 - np.array(azimuth)

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < 0)
        if wh_neg[0].size > 0:
            new_azimuth[wh_neg] = new_azimuth + np.pi * 2
    else:
        if new_azimuth < 0:
            new_azimuth = new_azimuth + np.pi * 2

    return zenith_angle, new_azimuth


def zenithangle_azimuth_to_altaz(zenith_angle, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Args:
        zenith_angle: in radians
        azimuth: in radians in uvbeam convention: North of East(East=0, North=90)

    Returns:
        altitude in radians
        azimuth in radians in astropy convention: East of North (N=0, E=90 degrees)
    """
    altitude = np.pi / 2 - np.array(zenith_angle)
    new_azimuth = np.pi / 2 - np.array(azimuth)

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < 0)
        if wh_neg[0].size > 0:
            new_azimuth[wh_neg] = new_azimuth + np.pi * 2
    else:
        if new_azimuth < 0:
            new_azimuth = new_azimuth + np.pi * 2

    return altitude, new_azimuth
