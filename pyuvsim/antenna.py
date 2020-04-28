# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import astropy.units as units
import numpy as np

from . import utils as simutils
from .telescope import BeamList


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

    def get_beam_jones(self, array, source_alt_az, frequency, reuse_spline=True,
                       interpolation_function='az_za_simple', freq_interp_kind=None):
        """
        2x2 array of Efield vectors in Az/Alt

        Parameters
        ----------
        array : Telescope object
            Provides the beam list.
        source_alt_az : array_like
            Positions to evaluate in alt/az, where
            source_alt_az[0] gives list of alts
            soruce_alt_az[1] gives list of corresponding az
        frequency : float or Quantity
            Frequency. Assumed to be Hz if float.
        reuse_spline : bool
            Option to keep and reuse interpolation splines in UVBeam.
        interpolation_function: str
            Set the angular interpolation function on the UVBeam.
            See UVBeam.interp for options.
            Defaults to az_za_simple.
        freq_interp_kind : str
            Interpolation method for frequencies.
            Note -- This overrides whatever method may be set on the
            UVBeam objects.
        Returns
        -------

        jones_matrix : (2,2) ndarray, dtype float
            The first axis is feed, the second axis is vector component
            on the sky in az/za.
        """

        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in alt/az

        # convert to UVBeam az/za convention
        source_za, source_az = simutils.altaz_to_zenithangle_azimuth(
            source_alt_az[0], source_alt_az[1]
        )

        if isinstance(frequency, units.Quantity):
            freq = np.array([frequency.to('Hz').value])
        else:
            freq = np.array([frequency])

        if array.beam_list[self.beam_id].data_normalization != 'peak':
            array.beam_list[self.beam_id].peak_normalize()

        if freq_interp_kind is not None:
            array.beam_list[self.beam_id].freq_interp_kind = freq_interp_kind

        if interpolation_function is not None:
            array.beam_list[self.beam_id].interpolation_function = interpolation_function

        spline_opts = None
        if isinstance(array.beam_list, BeamList):
            spline_opts = array.beam_list.spline_interp_opts

        interp_kwargs = {'az_array' : source_az, 'za_array' : source_za,
                         'freq_array' : freq, 'reuse_spline' : reuse_spline}

        if spline_opts is not None:
            interp_kwargs['spline_opts'] = spline_opts

        try:
            interp_data, interp_basis_vector = \
                array.beam_list[self.beam_id].interp(**interp_kwargs)
        except TypeError as err:   # pragma: nocover
            raise TypeError(
                "pyuvdata version >=2.0.1 required to use spline_interp_opts"
            ) from err

        Ncomponents = source_za.shape[-1]

        # interp_data has shape:
        #   (Naxes_vec, Nspws, Nfeeds, 1 (freq),  Ncomponents (source positions))
        jones_matrix = np.zeros((2, 2, Ncomponents), dtype=np.complex)
        # first axis is feed, second axis is theta, phi (opposite order of beam!)
        jones_matrix[0, 0] = interp_data[1, 0, 0, 0, :]
        jones_matrix[1, 1] = interp_data[0, 0, 1, 0, :]
        jones_matrix[0, 1] = interp_data[0, 0, 0, 0, :]
        jones_matrix[1, 0] = interp_data[1, 0, 1, 0, :]

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
