# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import time as pytime
import sys
import os
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
        Similar to a progress bar, this prints a percentage of task completion.
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
        self.curval = -1

    def update(self, count):
        if count % self.step == 0:
            doprint = False
            if not self.curval == count:
                doprint = True
                self.curval = count
            if doprint:
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
        azimuth: in radians in astropy convention: East of North (N=0, E=pi/2)

    Returns:
        zenith_angle in radians
        azimuth in radians in uvbeam convention: North of East(East=0, North=pi/2)
    """
    input_alt = np.array(altitude)
    input_az = np.array(azimuth)
    if input_alt.size != input_az.size:
        raise ValueError('number of altitude and azimuth values must match.')

    zenith_angle = np.pi / 2 - input_alt
    new_azimuth = np.pi / 2 - input_az

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < -1e-9)
        if wh_neg[0].size > 0:
            new_azimuth[wh_neg] = new_azimuth[wh_neg] + np.pi * 2
    else:
        if new_azimuth < -1e-9:
            new_azimuth = new_azimuth + np.pi * 2

    return zenith_angle, new_azimuth


def zenithangle_azimuth_to_altaz(zenith_angle, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Args:
        zenith_angle: in radians
        azimuth: in radians in uvbeam convention: North of East(East=0, North=pi/2)

    Returns:
        altitude in radians
        azimuth in radians in astropy convention: East of North (N=0, E=pi/2)
    """
    input_za = np.array(zenith_angle)
    input_az = np.array(azimuth)
    if input_za.size != input_az.size:
        raise ValueError('number of zenith_angle and azimuth values must match.')

    altitude = np.pi / 2 - input_za
    new_azimuth = np.pi / 2 - input_az

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < -1e-9)
        if wh_neg[0].size > -1e-9:
            new_azimuth[wh_neg] = new_azimuth[wh_neg] + np.pi * 2
    else:
        if new_azimuth < -1e-9:
            new_azimuth = new_azimuth + np.pi * 2

    return altitude, new_azimuth


def strip_extension(filepath):
    """ Remove extension from file. """
    if '.' not in filepath:
        return filepath, ''
    file_list = filepath.split('.')
    ext = file_list[-1]
    # miriad files might not have an extension
    # support miriad, uvfits, and uvh5
    # yaml is included for the tests
    if ext not in ['uvfits', 'uvh5', 'yaml']:
        return filepath, ''
    return ".".join(file_list[:-1]), '.' + file_list[-1]


def check_file_exists_and_increment(filepath):
    """
        Given filepath (path + filename), check if it exists. If so, add a _1
        at the end, if that exists add a _2, and so on.
    """
    base_filepath, ext = strip_extension(filepath)
    bf_list = base_filepath.split('_')
    if bf_list[-1].isdigit():
        base_filepath = '_'.join(bf_list[:-1])
    n = 0
    while os.path.exists(filepath):
        filepath = "{}_{}".format(base_filepath, n) + ext
        n += 1
    return filepath
