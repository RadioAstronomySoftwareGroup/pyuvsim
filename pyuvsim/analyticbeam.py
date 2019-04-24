# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import warnings
from scipy.special import j1

import pyuvdata.utils as uvutils


def diameter_to_sigma(diam, freqs):
    """
    Find the stddev of a gaussian with fwhm equal to that of
    an Airy disk's main lobe for a given diameter.

    Args:
        diam: Antenna diameter in meters
        freqs: Array of frequencies, in Hz
    Returns:
        sigma: The standard deviation in zenith angle radians for a Gaussian beam
               with FWHM equal to that of an Airy disk's main lobe for an aperture
               with the given diameter.
    """

    c_ms = 299792458.
    wavelengths = c_ms / freqs

    scalar = 2.2150894        # Found by fitting a Gaussian to an Airy disk function

    sigma = np.arcsin(scalar * wavelengths / (np.pi * diam)) * 2 / 2.355

    return sigma


class AnalyticBeam(object):
    """
    Defines an object with similar functionality to pyuvdata.UVBeam

    Directly calculates jones matrices at given azimuths and zenith angles
    from analytic functions.

    Supports uniform (unit response in all directions), gaussian, and Airy
    function beam types.

    Supported types:
        Uniform beam : Unit response from all directions.
        Airy: An Airy disk pattern (the 2D Fourier transform of a circular aperture of width given by `diameter`)
        Gaussian : A peak-normalized gaussian function.
                   If given a `diameter`, then this makes a chromatic beam with FWHMs matching an equivalent Airy disk beam at each frequency.
                   If given a `sigma`, this makes an achromatic beam with standard deviation set to `sigma`
                   If given a `sigma`, `ref_freq`, and `spectral_index`, then this will make a chromatic beam
                        with standard deviation defined by a power law:
                            stddev(f) = sigma * (f/ref_freq)**(spectral_index)

    Args:
        type: (string)
            Beam type to use
        sigma: (float)
            standard deviation [radians] for gaussian beam
            When spectral index is set, this represents the FWHM at the ref_freq.
        spectral_index : (float, optional)
            Scale gaussian beam width as a power law with frequency.
        ref_freq : (float, optional)
            If set, this sets the reference frequency for the beam width power law.
    """

    supported_types = ['uniform', 'gaussian', 'airy']

    def __init__(self, type, sigma=None, diameter=None, spectral_index=0.0, ref_freq=None):
        if type in self.supported_types:
            self.type = type
        else:
            raise ValueError('type not recognized')

        self.sigma = sigma
        if self.type == 'gaussian' and self.sigma is not None:
            warnings.warn("Achromatic gaussian beams will not be supported in the future."
                          + "Define your gaussian beam by a dish diameter from now on.", PendingDeprecationWarning)

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
        pass

    def efield_to_power(self):
        """
        Tell interp to return values corresponding with a power beam.
        """

        self.beam_type = 'power'
        pol_strings = ['XX', 'XY', 'YX', 'YY']
        self.polarization_array = np.array([uvutils.polstr2num(ps.upper()) for ps in pol_strings])

    def interp(self, az_array, za_array, freq_array, reuse_spline=None):
        """
        Evaluate the primary beam at given az, za locations (in radians).

        (similar to UVBeam.interp)

        Args:
            az_array: az values to evaluate at in radians (same length as za_array)
                The azimuth here has the UVBeam convention: North of East(East=0, North=pi/2)
            za_array: za values to evaluate at in radians (same length as az_array)
            freq_array: frequency values to evaluate at
            reuse_spline: Does nothing for analytic beams. Here for compatibility with UVBeam.

        Returns:
            an array of beam values, shape (Naxes_vec, Nspws, Nfeeds or Npols,
                Nfreqs or freq_array.size if freq_array is passed,
                Npixels/(Naxis1, Naxis2) or az_array.size if az/za_arrays are passed)
            an array of interpolated basis vectors (or self.basis_vector_array
                if az/za_arrays are not passed), shape: (Naxes_vec, Ncomponents_vec,
                Npixels/(Naxis1, Naxis2) or az_array.size if az/za_arrays are passed)
        """

        if self.type == 'uniform':
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            interp_data[1, 0, 0, :, :] = 1
            interp_data[0, 0, 1, :, :] = 1
            interp_data[1, 0, 1, :, :] = 1
            interp_data[0, 0, 0, :, :] = 1
            interp_basis_vector = None
        elif self.type == 'gaussian':
            if (self.diameter is None) and (self.sigma is None):
                raise ValueError("Dish diameter needed for gaussian beam -- units: meters")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            # gaussian beam only depends on Zenith Angle (symmetric is azimuth)
            # standard deviation of sigma is referring to the standard deviation of e-field beam!
            # copy along freq. axis
            if self.diameter is not None:
                sigmas = diameter_to_sigma(self.diameter, freq_array)
            elif self.sigma is not None:
                sigmas = self.sigma * (freq_array / self.ref_freq)**(self.spectral_index)
            values = np.exp(-(za_array[np.newaxis, ...]**2) / (2 * sigmas[:, np.newaxis]**2))
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_data[1, 0, 1, :, :] = values
            interp_data[0, 0, 0, :, :] = values
            interp_basis_vector = None
        elif self.type == 'airy':
            if self.diameter is None:
                raise ValueError("Dish diameter needed for airy beam -- units: meters")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            za_grid, f_grid = np.meshgrid(za_array, freq_array)
            xvals = self.diameter / 2. * np.sin(za_grid) * 2. * np.pi * f_grid / 3e8
            values = np.zeros_like(xvals)
            nz = xvals != 0.
            ze = xvals == 0.
            values[nz] = 2. * j1(xvals[nz]) / xvals[nz]
            values[ze] = 1.
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_data[1, 0, 1, :, :] = values
            interp_data[0, 0, 0, :, :] = values
            interp_basis_vector = None
        else:
            raise ValueError('no interp for this type: ', self.type)

        if self.beam_type == 'power':
            # Cross-multiplying feeds, adding vector components
            pairs = [(i, j) for i in range(2) for j in range(2)]
            power_data = np.zeros((1, 1, 4) + values.shape, dtype=np.float)
            for pol_i, pair in enumerate(pairs):
                power_data[:, :, pol_i] = ((interp_data[0, :, pair[0]]
                                            * np.conj(interp_data[0, :, pair[1]]))
                                           + (interp_data[1, :, pair[0]]
                                              * np.conj(interp_data[1, :, pair[1]])))
            interp_data = power_data

        return interp_data, interp_basis_vector

    def __eq__(self, other):
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
