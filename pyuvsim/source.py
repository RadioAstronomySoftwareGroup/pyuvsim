# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz


class Source(object):
    """
    Defines a single point source at a given ICRS ra/dec coordinate, with a
    flux density defined by stokes parameters.

    Used to calculate local coherency matrix for source brightness in AltAz
    frame at a specified time.
    """
    name = None
    freq = None
    stokes = None
    ra = None
    dec = None
    coherency_radec = None
    alt_az = None

    def __init__(self, name, ra, dec, freq, stokes, rise_lst=None, set_lst=None, pos_tol=np.finfo(float).eps):
        """
        Initialize from source catalog

        Args:
            name: Unique identifier for the source.
            ra: astropy Angle object
                source RA in J2000 (or ICRS) coordinates
            dec: astropy Angle object
                source Dec in J2000 (or ICRS) coordinates
            stokes:
                4 element vector giving the source [I, Q, U, V]
            freq: astropy quantity
                frequency of source catalog value
            rise_lst: (float)
                Approximate lst (radians) when the source rises. Set by coarse horizon cut in simsetup.
                Default is None, meaning the source never rises.
            set_lst: (float)
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

        self.name = name
        self.freq = freq
        self.stokes = stokes
        self.ra = ra
        self.dec = dec
        self.pos_tol = pos_tol
        self.time = None
        self.rise_lst = rise_lst
        self.set_lst = set_lst

        self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        self.coherency_radec = .5 * np.array([[self.stokes[0] + self.stokes[1],
                                               self.stokes[2] - 1j * self.stokes[3]],
                                              [self.stokes[2] + 1j * self.stokes[3],
                                               self.stokes[0] - self.stokes[1]]])

    def coherency_calc(self, time, telescope_location):
        """
        Calculate the local coherency in alt/az basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local alt/az.

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            local coherency in alt/az basis
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        if np.sum(np.abs(self.stokes[1:])) == 0:
            rotation_matrix = np.array([[1, 0], [0, 1]])
        else:
            # First need to calculate the sin & cos of the parallactic angle
            # See Meeus's astronomical algorithms eq 14.1
            # also see Astroplan.observer.parallactic_angle method
            time.location = telescope_location
            lst = time.sidereal_time('apparent')

            cirs_source_coord = self.skycoord.transform_to('cirs')
            tee_ra = simutils.cirs_to_tee_ra(cirs_source_coord.ra, time)

            hour_angle = (lst - tee_ra).rad
            sinX = np.sin(hour_angle)
            cosX = np.tan(telescope_location.lat) * np.cos(self.dec) - np.sin(self.dec) * np.cos(hour_angle)

            rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]])

        coherency_local = np.einsum('ab,bc,cd->ad', rotation_matrix.T,
                                    self.coherency_radec, rotation_matrix)

        return coherency_local

    def alt_az_calc(self, time, telescope_location):
        """
        calculate the altitude & azimuth for this source at a time & location

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            (altitude, azimuth) in radians
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        source_altaz = self.skycoord.transform_to(AltAz(obstime=time, location=telescope_location))

        alt_az = (source_altaz.alt.rad, source_altaz.az.rad)
        self.time = time
        self.alt_az = alt_az
        return alt_az

    def get_alt_az(self, time, telescope_location):
        """ Reuse alt_az if already calculated """
        if (self.alt_az is None) or (not time == self.time):
            return self.alt_az_calc(time, telescope_location)
        else:
            return self.alt_az

    def pos_lmn(self, time, telescope_location):
        """
        calculate the direction cosines of this source at a time & location

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            (l, m, n) direction cosine values
        """
        # calculate direction cosines of source at current time and array location
        # Will only do the calculation if time has changed
        alt_az = self.get_alt_az(time, telescope_location)

        # Horizon Mask
        if alt_az[0] < 0:
            return None

        pos_l = np.sin(alt_az[1]) * np.cos(alt_az[0])
        pos_m = np.cos(alt_az[1]) * np.cos(alt_az[0])
        pos_n = np.sin(alt_az[0])

        return (pos_l, pos_m, pos_n)

    def __eq__(self, other):
        return (np.isclose(self.ra.deg, other.ra.deg, atol=self.pos_tol)
                and np.isclose(self.dec.deg, other.dec.deg, atol=self.pos_tol)
                and np.all(self.stokes == other.stokes)
                and (self.name == other.name)
                and (self.freq == other.freq))
