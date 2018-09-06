# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from .spherical_coordinates_basis_transformation import spherical_basis_transformation_components


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
            basis_transformation_matrix = np.array([[1, 0], [0, 1]])
        else:
            # Calculate the coherency transformation matrix

            # unit vectors to be transformed by astropy
            x_c = np.array([1., 0, 0])
            y_c = np.array([0, 1., 0])
            z_c = np.array([0, 0, 1.])

            ''' We are using GCRS rather than ICRS to explicitly neglect the effect of aberration, which is a position-dependent
                effect, and obtain a proper, orthgonal rotation matrix between the RA/Dec and Alt/Az coordinate systems.  It is
                not clear (to JEA) what the correct thing for the transformation of Stokes parameters is in the presence of
                aberration '''
            axes_icrs = SkyCoord(x=x_c, y=y_c, z=z_c,
                                 obstime=time,
                                 location=telescope_location,
                                 frame='gcrs',
                                 equinox=self.epoch,
                                 representation='cartesian')
            axes_altaz = axes_icrs.transform_to('altaz')
            axes_altaz.representation = 'cartesian'

            # Test that this actually is a rotation matrix (i.e., orthogonal)
            R = np.array(axes_altaz.cartesian.xyz)  # the 3D rotation matrix that defines the mapping (RA,Dec) <--> (Alt,Az)

            # This calculation is for a single point on the sphere.  The rotation_matrix below is different for every point
            # on the sphere.
            # Check np.power(cosX,2)+np.power(sinX,2) = 1
            cosX, sinX = spherical_basis_transformation_components(self.dec.rad, self.ra.rad, R)
            rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]])

            alt_to_za = np.array([[-1., 0], [0, 1]]) # ??? Supposedly fixes theta_hat points south for za and north for alt
            basis_transformation_matrix = np.einsum('ab,bc->ac', alt_to_za, rotation_matrix)

        coherency_local = np.einsum('ab,bc,cd->ad', basis_transformation_matrix.T, self.coherency_radec, basis_transformation_matrix)

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
