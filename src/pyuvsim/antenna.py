# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Definition of Antenna objects, to describe a single interferometric element."""

import astropy.units as units
import numpy as np
import numpy.typing as npt

from . import utils as simutils


class Antenna:
    """
    Describe a single interferometric element.

    One of these defined for each antenna in the array.

    Parameters
    ----------
    name :  str
        Name of this antenna.
    number : int
        Number of this antenna.
    enu_position : array_like of float
        Position of this antenna in meters in the East, North, Up frame centered on the
        telescope location.
    beam_id : int
        Index of the beam for this antenna from :class:`~BeamList`.

    """

    def __init__(
        self, name: str, number: int, enu_position: npt.NDArray[float], beam_id: int
    ):
        assert isinstance(name, str), "name must be a string"
        assert isinstance(number, int | np.int32 | np.int64), "number must be an int"
        assert isinstance(beam_id, int | np.int32 | np.int64), "beam_id must be an int"
        self.name = name
        self.number = number
        self.pos_enu = enu_position * units.m
        self.beam_id = beam_id

    def get_beam_jones(
        self,
        array,
        source_alt_az,
        frequency,
        reuse_spline=True,
        interpolation_function=None,
        freq_interp_kind=None,
        beam_interp_check=True,
    ):
        """
        Calculate the jones matrix for this antenna in the direction of sources.

        A Nfeeds x 2 x Nsources array of Efield vectors in Az/Alt.

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
            Option to keep and reuse interpolation splines in :class:`pyuvdata.UVBeam`.
        interpolation_function: str
            Set the angular interpolation function on the :class:`pyuvdata.UVBeam`.
            See :meth:`pyuvdata.UVBeam.interp` for options.
        freq_interp_kind : str
            Interpolation method for frequencies. See :meth:`pyuvdata.UVBeam.interp`
            for options.
        beam_interp_check : bool
            Option to check whether the beam covers all the skymodel az/za values
            to ensure that they are covered by the intrinsic data array. Checking
            them can be quite computationally expensive. Defaults to True to run
            the check, can be set to False if the beam coverage is known.

        Returns
        -------
        jones_matrix : array_like of float
            Jones matricies for each source location, shape (Nfeeds,2, Ncomponents). The
            first axis is feed, the second axis is vector component on the sky in az/za.

        """
        # convert to UVBeam az/za convention
        source_za, source_az = simutils.altaz_to_zenithangle_azimuth(
            source_alt_az[0], source_alt_az[1]
        )

        if isinstance(frequency, units.Quantity):
            freq = np.array([frequency.to("Hz").value])
        else:
            freq = np.array([frequency])

        if self.beam_id not in range(len(array.beam_list)):
            raise ValueError(
                f"This antenna beam_id is {self.beam_id}, which is too large "
                f"for the beam_list, which has length {len(array.beam_list)}."
            )
        if array.beam_list.beam_type != "efield":
            raise ValueError("Beam type must be efield!")

        if (
            array.beam_list.data_normalization is not None
            and array.beam_list.data_normalization != "peak"
        ):
            raise ValueError("UVBeams must be peak normalized.")

        bi = array.beam_list[self.beam_id]

        interp_kwargs = {
            "az_array": source_az,
            "za_array": source_za,
            "freq_array": freq,
            "reuse_spline": reuse_spline,
            "check_azza_domain": beam_interp_check,
        }

        if interpolation_function is not None:
            interp_kwargs["interpolation_function"] = interpolation_function

        spline_opts = array.beam_list.spline_interp_opts
        if spline_opts is not None:
            interp_kwargs["spline_opts"] = spline_opts
        bl_freq_interp_kind = array.beam_list.freq_interp_kind

        if freq_interp_kind is None:
            freq_interp_kind = bl_freq_interp_kind

        if freq_interp_kind is not None:
            interp_kwargs["freq_interp_kind"] = freq_interp_kind

        interp_data = bi.compute_response(**interp_kwargs)

        Ncomponents = source_za.shape[-1]
        Nfeeds = interp_data.shape[1]
        # interp_data has shape:
        #   (Naxes_vec, Nfeeds, 1 (freq),  Ncomponents (source positions))
        jones_matrix = np.zeros((Nfeeds, 2, Ncomponents), dtype=complex)

        # first axis is feed, second axis is theta, phi (opposite order of beam!)
        jones_matrix[:, 0] = interp_data[1, :, 0, :]
        jones_matrix[:, 1] = interp_data[0, :, 0, :]
        return jones_matrix

    def __eq__(self, other):
        """Test for equality of objects."""
        return (
            (self.name == other.name)
            and np.allclose(
                self.pos_enu.to("m").value, other.pos_enu.to("m").value, atol=1e-3
            )
            and (self.beam_id == other.beam_id)
        )

    def __gt__(self, other):
        """Test for larger antenna number."""
        return self.number > other.number

    def __ge__(self, other):
        """Test for larger or equal antenna number."""
        return self.number >= other.number

    def __lt__(self, other):
        """Test for smaller antenna number."""
        return not self.__ge__(other)

    def __le__(self, other):
        """Test for smaller or equal antenna number."""
        return not self.__gt__(other)
