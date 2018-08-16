# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import numpy as np
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz


class Source(object):
    name = None
    freq = None
    stokes = None
    ra = None
    dec = None
    coherency_radec = None
    az_za = None

    def __init__(self, name, ra, dec, freq, stokes):
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

        self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        self.coherency_radec = .5 * np.array([[self.stokes[0] + self.stokes[1],
                                               self.stokes[2] - 1j * self.stokes[3]],
                                              [self.stokes[2] + 1j * self.stokes[3],
                                               self.stokes[0] - self.stokes[1]]])

    def __eq__(self, other):
        return ((self.ra == other.ra)
                and (self.dec == other.dec)
                and np.all(self.stokes == other.stokes)
                and (self.name == other.name)
                and (self.freq == other.freq))

    @profile
    def coherency_calc(self, time, telescope_location):
        """
        Calculate the local coherency in az/za basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local az/za.

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            local coherency in az/za basis
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

    @profile
    def az_za_calc(self, time, telescope_location):
        """
        calculate the azimuth & zenith angle for this source at a time & location

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            (azimuth, zenith_angle) in radians
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        source_altaz = self.skycoord.transform_to(AltAz(obstime=time, location=telescope_location))

        az_za = (source_altaz.az.rad, source_altaz.zen.rad)
        self.az_za = az_za
        return az_za

    @profile
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
        if self.az_za is None:
            self.az_za_calc(time, telescope_location)

        # Need a horizon mask, for now using pi/2
        if self.az_za[1] > (np.pi / 2.):
            return None

        pos_l = np.sin(self.az_za[0]) * np.sin(self.az_za[1])
        pos_m = np.cos(self.az_za[0]) * np.sin(self.az_za[1])
        pos_n = np.cos(self.az_za[1])
        return (pos_l, pos_m, pos_n)
