# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import copy
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from . import utils as simutils


class Source(object):
    """
    Does nothing. Placeholder while I fix other things

    """

    def __init__(self):
        self.nothing = None


class _view_source(object):
    """
    Interface for accessing a subset of sources in a source list.

    This is returned by Sources.__getitem__ and should not be used directly.
    """

    def __init__(self, sourcelist, index):
        self.slice = index
        self.sourcelist = sourcelist

    def __getattr__(self, key):

        if key in self.sourcelist._Ncomp_attrs:
            val = getattr(self.sourcelist, key)
            if val is None:
                return val
            val = val[..., self.slice]
            if isinstance(val, np.ndarray):
                # Turn length 0 arrays into scalars
                if val.size == 1:
                    val = val.flatten()[0]
            return val
        elif key in self.sourcelist._scalar_attrs:
            return getattr(self.sourcelist, key)
        elif key in self.sourcelist._member_funcs:
            func = getattr(self.sourcelist, key)
            return lambda x: func(x, src_inds=self.slice)

    def copy(self):
        """
        Copy this selected segment as a new Sources object.
        """
        src_obj = Sources(None, None, None, None, None)
        attrs = self.sourcelist._Ncomp_attrs + self.sourcelist._scalar_attrs
        for key in attrs:
            if hasattr(self.sourcelist, key):
                val = copy.copy(self.__getattr__(key))
                setattr(src_obj, key, val)

        src_obj.Ncomponents = src_obj.ra.size
        return src_obj


class Sources(object):
    """
    Defines a set of point source components at given ICRS ra/dec coordinates, with a
    flux densities defined by stokes parameters.

    Used to calculate local coherency matrix for source brightness in AltAz
    frame at a specified time.
    """

    Ncomponents = None      # Number of point source components represented here.

    name = None
    freq = None
    stokes = None
    ra = None
    dec = None
    coherency_radec = None
    alt_az = None
    pos_lmn = None


    _Ncomp_attrs = ['ra', 'dec', 'coherency_radec', 'coherency_local',
                       'stokes', 'alt_az', 'rise_lst', 'set_lst', 'freq',
                       'pos_lmn', 'name']
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

        if stokes is None:
            if np.all([None] * 4 == [name, ra, dec, freq]):
                # Not initializing here.
                return

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

        self.name = name
        self.freq = freq
        self.stokes = stokes
        self.ra = ra
        self.dec = dec
        self.pos_tol = pos_tol
        self.time = None
        self.rise_lst = rise_lst
        self.set_lst = set_lst

        # If rise/set lsts are not passed in they should be Ncomponent arrays of NaN
        if self.rise_lst is None:
            self.rise_lst = np.array([np.nan] * self.Ncomponents)
        if self.set_lst is None:
            self.set_lst = np.array([np.nan] * self.Ncomponents)


        # Get a list of all member functions
#        method_list = [func for func in dir(self) if callable(getattr(self, func)) and not func.startswith("__")]

        assert np.all([self.Ncomponents == l for l in
                       [self.ra.size, self.dec.size, self.freq.size, self.stokes.shape[1]]]), 'Inconsistent quantity dimensions.'

        if self.Ncomponents == 1:
            self.stokes = self.stokes.reshape(1, 1)

        self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        # Shape = (2,2,Ncomponents)
        self.coherency_radec = .5 * np.array([[self.stokes[0, :] + self.stokes[1, :],
                                               self.stokes[2, :] - 1j * self.stokes[3, :]],
                                              [self.stokes[2, :] + 1j * self.stokes[3, :],
                                               self.stokes[0, :] - self.stokes[1, :]]])

    def select(self, s):
        """
        Subselect from the Ncomponents
        """

        for k in self._Ncomp_attrs:
            if not hasattr(self, k):
                continue
            val = getattr(self, k)
            if val is not None:
                val = val[..., s]
            if isinstance(val, np.ndarray):
                if val.size == 1:
                    val = val.flatten()[0]
            setattr(self, k, val)

        self.Ncomponents = self.ra.size


    def __getitem__(self, i):
        # i = valid index object
        return _view_source(self, i)

    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        if self._counter < self.Ncomponents:
            i = self._counter
            self._counter += 1
            return self[i]
        else:
            raise StopIteration

    def __len__(self):
        return self.Ncomponents

    def next(self):
         # For python 2 compatibility
        return self.__next__()

    def coherency_calc(self, time, telescope_location, src_inds=slice(None)):
        """
        Calculate the local coherency in alt/az basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local alt/az.

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object
            src_inds: Optionally, only calculate for a subset of sources.

        Returns:
            local coherency in alt/az basis
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

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

        return coherency_local

    def update_positions(self, time, telescope_location):
        """
        Calculate the altitude/azimuth positions for these sources
        From alt/az, calculate direction cosines (lmn)

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object
            src_inds: Optionally, only calculate for a subset of sources.

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

        source_altaz = self.skycoord[src_inds].transform_to(AltAz(obstime=time, location=telescope_location))

        self.time = time
        alt_az = np.array([source_altaz.alt.rad, source_altaz.az.rad])
        self.alt_az = alt_az

        pos_l = np.sin(alt_az[1,:]) * np.cos(alt_az[0,:])
        pos_m = np.cos(alt_az[1,:]) * np.cos(alt_az[0,:])
        pos_n = np.sin(alt_az[0,:])

        if self.pos_lmn is None:
            self.pos_lmn = np.array([np.zeros(self.Ncomponents)] * 3)

        self.pos_lmn[0,src_inds] = pos_l
        self.pos_lmn[1,src_inds] = pos_m
        self.pos_lmn[2,src_inds] = pos_n

        # Horizon mask:
        self.horizon_mask = self.alt_az[:,0] < 0.0

    def __eq__(self, other):
        return (np.isclose(self.ra.deg, other.ra.deg, atol=self.pos_tol)
                and np.isclose(self.dec.deg, other.dec.deg, atol=self.pos_tol)
                and np.all(self.stokes == other.stokes)
                and (self.name == other.name)
                and (self.freq == other.freq))
