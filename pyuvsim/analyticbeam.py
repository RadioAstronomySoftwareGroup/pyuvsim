# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import numpy as np
from scipy.special import spherical_jn as jn


class AnalyticBeam(object):

    supported_types = ['uniform', 'gaussian', 'airy']

    @profile
    def __init__(self, type, sigma=None, diameter=None):
        if type in self.supported_types:
            self.type = type
        else:
            raise ValueError('type not recognized')

        self.sigma = sigma
        self.diameter = diameter
        self.data_normalization = 'peak'

    def peak_normalize(self):
        pass

    @profile
    def interp(self, az_array, za_array, freq_array):
        # (Naxes_vec, Nspws, Nfeeds or Npols, freq_array.size, az_array.size)

        if self.type == 'uniform':
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            interp_data[1, 0, 0, :, :] = 1
            interp_data[0, 0, 1, :, :] = 1
            interp_data[1, 0, 1, :, :] = 1
            interp_data[0, 0, 0, :, :] = 1
            interp_basis_vector = None
        elif self.type == 'gaussian':
            if self.sigma is None:
                raise ValueError("Sigma needed for gaussian beam -- units: radians")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            # gaussian beam only depends on Zenith Angle (symmetric is azimuth)
            # standard deviation of sigma is refereing to the standard deviation of e-field beam!
            values = np.exp(-(za_array**2) / (2 * self.sigma**2))
            # copy along freq. axis
            values = np.broadcast_to(values, (freq_array.size, az_array.size))
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_data[1, 0, 1, :, :] = values
            interp_data[0, 0, 0, :, :] = values
            interp_basis_vector = None
        elif self.type == 'airy':
            if self.diameter is None:
                raise ValueError("Diameter needed for airy beam -- units: meters")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            za_grid, f_grid = np.meshgrid(za_array, freq_array)
            xvals = self.diameter / 2. * np.sin(za_grid) * 2. * np.pi * f_grid / 3e8
            values = np.zeros_like(xvals)
            values[xvals > 0.] = 2. * jn(1, xvals[xvals > 0.]) / xvals[xvals > 0.]
            values[xvals == 0.] = 1.
            interp_data[1, 0, 0, :, :] = values
            interp_data[0, 0, 1, :, :] = values
            interp_data[1, 0, 1, :, :] = values
            interp_data[0, 0, 0, :, :] = values
            interp_basis_vector = None
        else:
            raise ValueError('no interp for this type: ', self.type)

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
