# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import astropy.units as units

from . import utils as simutils


class Antenna(object):
    """
    Defines an object that can return a Jones matrix and has a specified
    location, name, and number.

    One of these defined for each antenna in the array.
    """

    def __init__(self, name, number, enu_position, beam_id):
        self.name = name
        self.number = number
        # ENU position in meters relative to the telescope_location
        self.pos_enu = enu_position * units.m
        # index of beam for this antenna from array.beam_list
        self.beam_id = beam_id

    def get_beam_jones(self, array, source_alt_az, frequency, reuse_spline=True):
        """
        2x2 array of Efield vectors in Az/Alt

        Args:
            array: az values to evaluate (same length as za_array)
            source_az_za: Tuple or list (azimuth, zenith angle) in radians
            frequency: (float) frequency in Hz
            reuse_spline: (bool) Reuse interpolation fits in beam objects.

        Returns:
            jones_matrix, A 2x2 float array. The first axis is feed, the
                          second axis is vector component on the sky in az/za.
        """
        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in alt/az

        # convert to UVBeam az/za convention
        source_za, source_az = simutils.altaz_to_zenithangle_azimuth(source_alt_az[0], source_alt_az[1])
        source_za = np.array([source_za])
        source_az = np.array([source_az])

        freq = np.array([frequency.to('Hz').value])

        if array.beam_list[self.beam_id].data_normalization != 'peak':
            array.beam_list[self.beam_id].peak_normalize()
        array.beam_list[self.beam_id].interpolation_function = 'az_za_simple'

        interp_data, interp_basis_vector = \
            array.beam_list[self.beam_id].interp(az_array=source_az,
                                                 za_array=source_za,
                                                 freq_array=freq,
                                                 reuse_spline=reuse_spline)

        # interp_data has shape: (Naxes_vec, Nspws, Nfeeds, 1 (freq), 1 (source position))
        jones_matrix = np.zeros((2, 2), dtype=np.complex)
        # first axis is feed, second axis is theta, phi (opposite order of beam!)
        jones_matrix[0, 0] = interp_data[1, 0, 0, 0, 0]
        jones_matrix[1, 1] = interp_data[0, 0, 1, 0, 0]
        jones_matrix[0, 1] = interp_data[0, 0, 0, 0, 0]
        jones_matrix[1, 0] = interp_data[1, 0, 1, 0, 0]

        return jones_matrix

    def __eq__(self, other):
        return ((self.name == other.name)
                and np.allclose(self.pos_enu.to('m').value, other.pos_enu.to('m').value, atol=1e-3)
                and (self.beam_id == other.beam_id))

    def __gt__(self, other):
        return (self.number > other.number)

    def __ge__(self, other):
        return (self.number >= other.number)

    def __lt__(self, other):
        return not self.__ge__(other)

    def __le__(self, other):
        return not self.__gt__(other)
