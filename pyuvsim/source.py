# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import copy
import sys

from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from . import utils as simutils


def _skymodel_basesize():
    """
    Estimate the memory footprint of a SkyModel with a single source.

    Sum the sizes of the data types that go into SkyModel
    """
    attrs = ['',
             Quantity(1.0, 'Hz'),
             [0.0] * 4,
             Angle(np.pi, 'rad'),
             Angle(np.pi, 'rad'),
             [1.5] * 4,
             [0.3] * 2,
             [0.0] * 3]
    return np.sum([sys.getsizeof(a) for a in attrs])


class SkyModel(object):
    """
    Defines a set of point source components at given ICRS ra/dec coordinates, with a
    flux densities defined by stokes parameters.

    Used to calculate local coherency matrix for source brightness in AltAz
    frame at a specified time.
    """

    Ncomponents = None      # Number of point source components represented here.

    _basesize = _skymodel_basesize()

    _Ncomp_attrs = ['ra', 'dec', 'coherency_radec', 'coherency_local',
                    'stokes', 'alt_az', 'rise_lst', 'set_lst', 'freq',
                    'pos_lmn', 'name', 'horizon_mask']
    _scalar_attrs = ['Ncomponents', 'time', 'pos_tol']

    _member_funcs = ['coherency_calc', 'update_positions']

    def __init__(self, name, ra, dec, freq, stokes,
                 rise_lst=None, set_lst=None, pos_tol=np.finfo(float).eps):
        """
        Initialize from source catalog

        Args:
            name: Unique identifier for the source.
            ra: astropy Angle object, shape (Ncomponents,)
                source RA in J2000 (or ICRS) coordinates
            dec: astropy Angle object, shape (Ncomponents,)
                source Dec in J2000 (or ICRS) coordinates
            stokes:
                shape (Ncomponents,4)
                4 element vector giving the source [I, Q, U, V]
            freq: astropy quantity, shape (Ncomponents,)
                reference frequencies of flux values
            rise_lst: (float), shape (Ncomponents,)
                Approximate lst (radians) when the source rises. Set by coarse horizon cut in simsetup.
                Default is nan, meaning the source never rises.
            set_lst: (float), shape (Ncomponents,)
                Approximate lst (radians) when the source sets.
                Default is None, meaning the source never sets.
            pos_tol: float, defaults to minimum float in numpy
                position tolerance in degrees
        """

        if not isinstance(ra, Angle):
            raise ValueError('ra must be an astropy Angle object. '
                             'value was: {ra}'.format(ra=ra))

        if not isinstance(dec, Angle):
            raise ValueError('dec must be an astropy Angle object. '
                             'value was: {dec}'.format(dec=dec))

        if not isinstance(freq, Quantity):
            raise ValueError('freq must be an astropy Quantity object. '
                             'value was: {f}'.format(f=freq))

        self.Ncomponents = ra.size

        self.name = np.atleast_1d(np.asarray(name))
        self.freq = np.atleast_1d(freq)
        self.stokes = np.asarray(stokes)
        self.ra = np.atleast_1d(ra)
        self.dec = np.atleast_1d(dec)
        self.pos_tol = pos_tol

        self.has_rise_set_lsts = False
        if (rise_lst is not None) and (set_lst is not None):
            self.rise_lst = np.asarray(rise_lst)
            self.set_lst = np.asarray(set_lst)
            self.has_rise_set_lsts = True

        self.alt_az = np.zeros((2, self.Ncomponents), dtype=float)
        self.pos_lmn = np.zeros((3, self.Ncomponents), dtype=float)

        self.horizon_mask = np.zeros(self.Ncomponents).astype(bool)     # If true, source component is below horizon.

        if self.Ncomponents == 1:
            self.stokes = self.stokes.reshape(4, 1)

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        # Shape = (2,2,Ncomponents)
        self.coherency_radec = .5 * np.array([[self.stokes[0, :] + self.stokes[1, :],
                                               self.stokes[2, :] - 1j * self.stokes[3, :]],
                                              [self.stokes[2, :] + 1j * self.stokes[3, :],
                                               self.stokes[0, :] - self.stokes[1, :]]])

        self.time = None

        assert np.all([self.Ncomponents == l for l in
                       [self.ra.size, self.dec.size, self.freq.size, self.stokes.shape[1]]]), 'Inconsistent quantity dimensions.'

    def coherency_calc(self, telescope_location):
        """
        Calculate the local coherency in alt/az basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local alt/az.

        Args:
            telescope_location: astropy EarthLocation object

        Returns:
            local coherency in alt/az basis
        """
        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        Ionly_mask = np.sum(self.stokes[1:, :], axis=0) == 0.0
        NstokesI = np.sum(Ionly_mask)   # Number of unpolarized sources

        # For unpolarized sources, there's no need to rotate the coherency matrix.
        coherency_local = self.coherency_radec.copy()

        if NstokesI < self.Ncomponents:
            # If there are any polarized sources, do rotation.

            # First need to calculate the sin & cos of the parallactic angle
            # See Meeus's astronomical algorithms eq 14.1
            # also see Astroplan.observer.parallactic_angle method
            polarized_sources = np.where(~Ionly_mask)
            sinX = np.sin(self.hour_angle)
            cosX = np.tan(telescope_location.lat) * np.cos(self.dec) - np.sin(self.dec) * np.cos(self.hour_angle)
    
            rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]]).astype(float)       # (2, 2, Ncomponents)
            rotation_matrix = rotation_matrix[..., polarized_sources]

            rotation_matrix_T = np.swapaxes(rotation_matrix, 0, 1)
            coherency_local[:,:,polarized_sources] = np.einsum('abx,bcx,cdx->adx', rotation_matrix_T,
                                                                self.coherency_radec[:,:,polarized_sources],
                                                                rotation_matrix)
        # Zero coherency on sources below horizon.
        coherency_local[:, :, self.horizon_mask] *= 0.0

        return coherency_local

    def update_positions(self, time, telescope_location):
        """
        Calculate the altitude/azimuth positions for these sources
        From alt/az, calculate direction cosines (lmn)

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Sets:
            self.pos_lmn: (3, Ncomponents)
            self.alt_az: (2, Ncomponents)
            self.time: (1,) Time object
        """

        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        if self.time == time:  # Don't repeat calculations
            return

        skycoord_use = SkyCoord(self.ra, self.dec, frame='icrs')
        source_altaz = skycoord_use.transform_to(AltAz(obstime=time, location=telescope_location))

        time.location = telescope_location
        lst = time.sidereal_time('apparent')

        cirs_source_coord = skycoord_use.transform_to('cirs')
        tee_ra = simutils.cirs_to_tee_ra(cirs_source_coord.ra, time)

        self.hour_angle = None
        self.hour_angle = (lst - tee_ra).rad

        self.time = time
        alt_az = np.array([source_altaz.alt.rad, source_altaz.az.rad])

        self.alt_az = alt_az

        pos_l = np.sin(alt_az[1, :]) * np.cos(alt_az[0, :])
        pos_m = np.cos(alt_az[1, :]) * np.cos(alt_az[0, :])
        pos_n = np.sin(alt_az[0, :])

        self.pos_lmn[0, :] = pos_l
        self.pos_lmn[1, :] = pos_m
        self.pos_lmn[2, :] = pos_n

        # Horizon mask:
        self.horizon_mask = self.alt_az[0, :] < 0.0

    def __eq__(self, other):
        return (np.allclose(self.ra.deg, other.ra.deg, atol=self.pos_tol)
                and np.allclose(self.stokes, other.stokes)
                and np.all(self.name == other.name)
                and np.all(self.freq == other.freq))

    def get_size(self):
        """
        Estimate the SkyModel memory footprint in bytes
        """
        return self.Ncomponents * self._basesize
