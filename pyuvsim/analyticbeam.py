# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Definition of Analytic Beam objects."""

import warnings

import numpy as np
import pyuvdata.utils as uvutils
from scipy.special import j1
from astropy.constants import c as speed_of_light

c_ms = speed_of_light.to('m/s').value


def diameter_to_sigma(diam, freqs):
    """
    Find the sigma that gives a beam width similar to an Airy disk.

    Find the stddev of a gaussian with fwhm equal to that of
    an Airy disk's main lobe for a given diameter.

    Parameters
    ----------
    diam : float
        Antenna diameter in meters
    freqs : array
        Frequencies in Hz

    Returns
    -------
    sigma : float
        The standard deviation in zenith angle radians for a Gaussian beam
        with FWHM equal to that of an Airy disk's main lobe for an aperture
        with the given diameter.
    """
    wavelengths = c_ms / freqs

    scalar = 2.2150894  # Found by fitting a Gaussian to an Airy disk function

    sigma = np.arcsin(scalar * wavelengths / (np.pi * diam)) * 2 / 2.355

    return sigma


class AnalyticBeam:
    """
    Calculate Jones matrices from analytic functions.

    This provides similar functionality to :class:`pyuvdata.UVBeam`, but the :meth:`~.interp` method
    evaluates the function at given azimuths and zenith angles, instead of interpolating
    from data.

    Supported types include:

        * Uniform beam: Unit response from all directions.
        * Airy: An Airy disk pattern (the 2D Fourier transform of a circular aperture of
          width given by `diameter`)
        * Gaussian: A peak-normalized gaussian function.
            * If given a `diameter`, then this makes a chromatic beam with FWHMs
              matching an equivalent Airy disk beam at each frequency.
            * If given a `sigma`, this makes an achromatic beam with standard deviation
              set to `sigma`
            * If given a `sigma`, `ref_freq`, and `spectral_index`, then this will make
              a chromatic beam with standard deviation defined by a power law:
              `stddev(f) = sigma * (f/ref_freq)**(spectral_index)`

    Parameters
    ----------
    type : str, {'uniform', 'airy', 'gaussian'}
        Beam type to use.
    sigma : float
        standard deviation [radians] for gaussian beam.
        When spectral index is set, this represents the FWHM at the ref_freq.
    spectral_index : float
        Scale gaussian beam width as a power law with frequency.
    ref_freq : float
        If set, this sets the reference frequency for the beam width power law.

    """

    supported_types = ['uniform', 'gaussian', 'airy']

    def __init__(self, type_, sigma=None, diameter=None, spectral_index=0.0, ref_freq=None):
        if type_ in self.supported_types:
            self.type = type_
        else:
            raise ValueError('type not recognized')

        self.sigma = sigma
        if self.type == 'gaussian' and self.sigma is not None:
            warnings.warn("Achromatic gaussian beams will not be supported in the future. "
                          + "Define your gaussian beam by a dish diameter from now on.",
                          PendingDeprecationWarning)

        if (spectral_index != 0.0) and (ref_freq is None):
            raise ValueError("ref_freq must be set for nonzero gaussian beam spectral index")
        elif ref_freq is None:
            ref_freq = 1.0
        self.ref_freq = ref_freq
        self.spectral_index = spectral_index
        self.diameter = diameter
        self.data_normalization = 'peak'
        self.freq_interp_kind = 'linear'
        self.beam_type = 'efield'

    def peak_normalize(self):
        """Do nothing, mocks the :meth:`pyuvdata.UVBeam.peak_normalize` method."""
        pass

    def efield_to_power(self):
        """Tell :meth:`~.interp` to return values corresponding with a power beam."""
        self.beam_type = 'power'
        pol_strings = ['XX', 'XY', 'YX', 'YY']
        self.polarization_array = np.array([uvutils.polstr2num(ps.upper()) for ps in pol_strings])

    def interp(self, az_array, za_array, freq_array, reuse_spline=None, spline_opts=None):
        """
        Evaluate the primary beam at given sky coordinates and frequencies.

        (mocks :meth:`pyuvdata.UVBeam.interp`, but these are analytic, so no interpolation is done.)

        Parameters
        ----------
        az_array : array-like of float
            Azimuth values to evaluate at in radians. Should be a 1D array with the same
            length as `za_array`. The azimuth here has the :class:`pyuvdata.UVBeam` convention:
            North of East (East=0, North=pi/2)
        za_array : array-like of float
            Zenith angle values to evaluate at in radians. Should be a 1D array with the
            same length as `az_array`.
        freq_array : array-like of float
            Frequency values to evaluate at in Hz. Should be a 1D array.
        reuse_spline : bool
            Unused. Here for compatibility with :meth:`pyuvdata.UVBeam.interp`.
        spline_opts : dict
            Unused. Here for compatibility with :meth:`pyuvdata.UVBeam.interp`.

        Returns
        -------
        beam_vals : array-like of float
            Array of beam values, shape (Naxes_vec, Nspws, Nfeeds or Npols,
                Nfreqs or freq_array.size if freq_array is passed,
                Npixels/(Naxis1, Naxis2) or az_array.size if az/za_arrays are passed)
        interp_basis_vectors : None
            Currently returns None. In :meth:`pyuvdata.UVBeam.interp`, this is the set
            of basis vectors for the electric field component values.

        Notes
        -----
        See :meth:`pyuvdata.UVBeam.interp` documentation for more details on the returned data.

        """
        if az_array.ndim > 1 or za_array.ndim > 1 or freq_array.ndim > 1:
            raise ValueError("az_array, za_array and freq_array must all be one dimensional.")

        if az_array.shape != za_array.shape:
            raise ValueError("az_array and za_array must have the same shape.")

        if self.type == 'uniform':
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=float)
            interp_data[1, 0, 0, :, :] = 1
            interp_data[0, 0, 1, :, :] = 1
            interp_data[1, 0, 1, :, :] = 0
            interp_data[0, 0, 0, :, :] = 0

            # If beam_type == "power", we need to know the shape of "values"
            values = interp_data[0, 0, 0]

            interp_basis_vector = None
        elif self.type == 'gaussian':
            if (self.diameter is None) and (self.sigma is None):
                raise ValueError("Dish diameter needed for gaussian beam -- units: meters")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=float)
            # gaussian beam only depends on Zenith Angle (symmetric is azimuth)
            # standard deviation of sigma is referring to the standard deviation of e-field beam!
            # copy along freq. axis
            if self.diameter is not None:
                sigmas = diameter_to_sigma(self.diameter, freq_array)
            elif self.sigma is not None:
                sigmas = self.sigma * (freq_array / self.ref_freq) ** self.spectral_index
            values = np.exp(-(za_array[np.newaxis, ...] ** 2) / (2 * sigmas[:, np.newaxis] ** 2))
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_basis_vector = None
        elif self.type == 'airy':
            if self.diameter is None:
                raise ValueError("Dish diameter needed for airy beam -- units: meters")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=float)
            za_grid, f_grid = np.meshgrid(za_array, freq_array)
            xvals = self.diameter / 2. * np.sin(za_grid) * 2. * np.pi * f_grid / c_ms
            values = np.zeros_like(xvals)
            nz = xvals != 0.
            ze = xvals == 0.
            values[nz] = 2. * j1(xvals[nz]) / xvals[nz]
            values[ze] = 1.
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_basis_vector = None
        else:
            raise ValueError('no interp for this type: {}'.format(self.type))

        if self.beam_type == 'power':
            # Cross-multiplying feeds, adding vector components
            pairs = [(i, j) for i in range(2) for j in range(2)]
            power_data = np.zeros((1, 1, 4) + values.shape, dtype=float)
            for pol_i, pair in enumerate(pairs):
                power_data[:, :, pol_i] = ((interp_data[0, :, pair[0]]
                                           * np.conj(interp_data[0, :, pair[1]]))
                                           + (interp_data[1, :, pair[0]]
                                           * np.conj(interp_data[1, :, pair[1]])))
            interp_data = power_data

        return interp_data, interp_basis_vector

    def __eq__(self, other):
        """Define equality for Analytic Beams."""
        if not isinstance(other, self.__class__):
            return False
        if self.type == 'gaussian':
            return ((self.type == other.type)
                    and (self.sigma == other.sigma))
        elif self.type == 'uniform':
            return other.type == 'uniform'
        elif self.type == 'airy':
            return ((self.type == other.type)
                    and (self.diameter == other.diameter))
        else:
            return False
