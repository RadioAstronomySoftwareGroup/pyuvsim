# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Definition of Telescope objects, for metadata common to all antennas in an array."""

from __future__ import annotations

import warnings
from dataclasses import KW_ONLY, InitVar, dataclass
from typing import Literal

import numpy as np
from astropy.coordinates import EarthLocation
from pyuvdata import BeamInterface, UVBeam, parameter, utils as uvutils
from pyuvdata.analytic_beam import AnalyticBeam

try:
    from . import mpi
except ImportError:
    mpi = None

try:
    from lunarsky import MoonLocation

    hasmoon = True
except ImportError:
    hasmoon = False


class BeamConsistencyError(Exception):
    """An Exception class to mark inconsistencies among beams in a BeamList."""

    pass


@dataclass(frozen=True)
class BeamList:
    """
    A container for the set of beam models and related parameters.

    When used with MPI, this can be initialized simultaneously on all
    processes such that sky_in is provided only on the root process.
    The data can then be shared across all processes by running the `share`
    method.

    This behaves just as a normal Python list in terms of indexing but it does
    not allow changes to the list.

    Attributes
    ----------
    beam_list : list of pyuvdata.BeamInterface
        The list of BeamInterface objects.
    spline_interp_opts : dict
        Degrees of the bivariate spline for the angular interpolation (passed
        directly to numpy.RectBivariateSpline, e.g. use {'kx' : 2, 'ky' : 2}
        for 2nd degree splines in each direction.) Default is 3rd order in
        each direction.
    freq_interp_kind : str
        The type of frequency interpolation, anything accepted by
        scipy.interpolate.interp1d. Default is "cubic".

    Parameters
    ----------
    beam_list : list pyuvdata.UVBeam or pyuvdata.AnalyticBeam
        A list of pyuvdata UVBeam or AnalyticBeam objects. These will be converted
        into BeamInterface objects during initialization.
    beam_type : str
        Desired beam type, either "efield" or "power", defaults to "efield".
    spline_interp_opts : dict
        Provide options to numpy.RectBivariateSpline. This includes spline
        order parameters `kx` and `ky`, and smoothing parameter `s`.
        Only applies UVBeams and for `az_za_simple` interpolation.
    freq_interp_kind : str
        Interpolation method to use for frequencies. See scipy.interpolate.interp1d
        for details. Defaults to "cubic".
    peak_normalize : bool
        Option to peak normalize UVBeams in the beam list. This is required to
        be True for the pyuvsim simulator. Defaults to True.

    """

    beam_list: [UVBeam | AnalyticBeam]
    _: KW_ONLY
    beam_type: Literal["efield", "power"] | None = "efield"
    spline_interp_opts: dict[str, int] | None = None
    freq_interp_kind: str = "cubic"
    peak_normalize: InitVar(bool) = True

    def __post_init__(self, peak_normalize):
        """
        Post-initialization validation and conversions, define the beam_list attribute.

        Parameters
        ----------
        peak_normalize : bool
            Option to peak normalize UVBeams in the beam list. This is required to
            be True for the pyuvsim simulator. Defaults to True.

        """
        beam_list = []
        for beam in self.beam_list:
            if peak_normalize and isinstance(beam, UVBeam):
                beam = beam.copy()
                beam.peak_normalize()
            # turn them into BeamInterface objects
            beam_list.append(BeamInterface(beam, beam_type=self.beam_type))
        object.__setattr__(self, "beam_list", beam_list)

        self._check_consistency()

    def _get_beam_basis_type(self):
        beam_basis = {}
        for index, bi in enumerate(self):
            if bi._isuvbeam:
                if bi.beam.pixel_coordinate_system not in ["az_za", "healpix"]:
                    raise ValueError(
                        "pyuvsim currently only supports UVBeams with 'az_za' or "
                        "'healpix' pixel coordinate systems."
                    )
                if (
                    np.all(bi.beam.basis_vector_array[0, 0] == 1)
                    and np.all(bi.beam.basis_vector_array[0, 1] == 0)
                    and np.all(bi.beam.basis_vector_array[1, 0] == 0)
                    and np.all(bi.beam.basis_vector_array[1, 1] == 1)
                ):
                    beam_basis[index] = "az_za"
                else:
                    raise ValueError(
                        "pyuvsim currently only supports beams with basis vectors that"
                        "are aligned with the azimuth and zenith angle in each pixel."
                        "Work is in progress to add other basis vector systems."
                    )
            else:
                # AnalyticBeams all have the az_za basis by construction
                beam_basis[index] = "az_za"
        return beam_basis

    def _check_consistency(self):
        """
        Check the consistency of all beams in the list.

        This checks basic parameters of the objects for consistency, eg. the
        ``beam_type``. It is called in the __post_init__.

        Raises
        ------
        BeamConsistencyError
            If any beam is inconsistent with the rest of the beams.

        """
        if len(self) == 0:
            return

        def check_thing(item):
            if item == "basis":
                items = self._get_beam_basis_type()
            elif item == "beam_type":
                items = {i: bi.beam_type for i, bi in enumerate(self)}
            elif item == "feed_array":
                items = {}
                for i, bi in enumerate(self):
                    if "e" in bi.beam.feed_array or "n" in bi.beam.feed_array:
                        feed_map = uvutils.x_orientation_pol_map(bi.beam.x_orientation)
                        inv_feed_map = {value: key for key, value in feed_map.items()}
                        items[i] = [inv_feed_map[feed] for feed in bi.beam.feed_array]
                    else:
                        items[i] = bi.beam.feed_array
            else:
                items = {
                    i: getattr(bi.beam, item)
                    for i, bi in enumerate(self)
                    if hasattr(bi.beam, item)
                }

            if not items:
                # If none of the beams has this item, move on.
                return

            min_index = min(items.keys())

            for indx, thing in items.items():
                if np.any(thing != items[min_index]):
                    raise BeamConsistencyError(
                        f"{item} of beam {indx + 1} is not consistent with beam "
                        f"{min_index + 1}: {thing} vs. {items[min_index]}"
                    )

        check_thing("beam_type")
        check_thing("x_orientation")
        check_thing("data_normalization")

        if self[0].beam_type == "efield":
            check_thing("Nfeeds")
            check_thing("feed_array")
            check_thing("basis")

        elif self[0].beam_type == "power":
            check_thing("Npols")
            check_thing("polarization_array")

    def check_all_azza_beams_full_sky(self):
        """
        Check if all az_za UVBeams cover the full sky.

        Used to decide whether to turn off checking that the beam covers the
        interpolation location.
        """
        all_azza_full_sky = True
        for bi in self:
            if bi._isuvbeam and bi.beam.pixel_coordinate_system == "az_za":
                # regular gridding is enforced in UVBeam checking, so can use first diff
                axis1_diff = np.diff(bi.beam.axis1_array)[0]
                axis2_diff = np.diff(bi.beam.axis2_array)[0]
                max_axis_diff = np.max([axis1_diff, axis2_diff])
                if not (
                    np.max(bi.beam.axis2_array) >= np.pi / 2.0 - max_axis_diff * 2.0
                    and np.min(bi.beam.axis1_array) <= max_axis_diff * 2.0
                    and np.max(bi.beam.axis1_array) >= 2.0 * np.pi - max_axis_diff * 2.0
                ):
                    all_azza_full_sky = False
                    break
        return all_azza_full_sky

    @property
    def x_orientation(self):
        """Return the x_orientation of all beams in list."""
        if len(self) == 0:
            return None

        for bi in self:
            if hasattr(bi.beam, "x_orientation"):
                xorient = bi.beam.x_orientation
                break
        else:
            # None of the constituent beams has defined x_orientation.
            # Here we return None, instead of raising an attribute error, to maintain
            # backwards compatibility (pyuvsim access the beamlist.x_orientation without
            # checking for its existence).
            return None

        if xorient is None:
            warnings.warn(
                "All polarized beams have x_orientation set to None. This will make it "
                "hard to interpret the polarizations of the simulated visibilities."
            )

        return xorient

    @property
    def data_normalization(self):
        """Return the data_normalization of all beams in list."""
        if len(self) == 0:
            return None

        for bi in self:
            if hasattr(bi.beam, "data_normalization"):
                data_norm = bi.beam.data_normalization
                break
        else:
            # None of the constituent beams has defined data_normalization.
            return None

        return data_norm

    def __len__(self):
        """Define the length of a BeamList."""
        return len(self.beam_list)

    def __iter__(self):
        """Get the list as an iterable."""
        return iter(self.beam_list)

    def __getitem__(self, ind):
        """Get a particular beam from the list."""
        return self.beam_list[ind]

    def _share_uvbeam(self, bi, root):
        mpi.start_mpi()

        # Get list of attributes that are set.
        rank = mpi.world_comm.Get_rank()

        set_values = []
        if rank == root:
            set_values = [
                key
                for key, value in bi.beam.__dict__.items()
                if value is not None and isinstance(value, parameter.UVParameter)
            ]

        mpi.world_comm.Barrier()

        set_values = mpi.world_comm.bcast(set_values, root=root)

        # share these attributes to the other ranks
        for key in set_values:
            if rank == root:
                attr = getattr(bi.beam, key)
            else:
                attr = parameter.UVParameter(key)

            # for the data array on a UVBeam
            # share the actual response pattern in shared memory
            # this could in theory be quite large
            if key == "_data_array":
                attr.value = mpi.shared_mem_bcast(attr.value, root=root)

                for metadata in attr.__dict__:
                    if metadata == "value":
                        continue
                    meta_to_assign = getattr(attr, metadata)
                    meta_to_assign = mpi.world_comm.bcast(meta_to_assign, root=root)
                    setattr(attr, metadata, meta_to_assign)
            else:
                attr = mpi.world_comm.bcast(attr, root=root)

            #  assign the initialized UVParameter to the new beam
            if rank != root:
                setattr(bi.beam, key, attr)

        return bi

    def share(self, root=0):
        """
        Share across MPI processes (requires mpi4py to use).

        All attributes are put in shared memory.

        Parameters
        ----------
        root : int
            Root rank on COMM_WORLD, from which data will be broadcast.

        """
        if mpi is None:
            raise ImportError(
                "You need mpi4py to use this method. "
                "Install it by running pip install pyuvsim[sim] "
                "or pip install pyuvsim[all] if you also want the "
                "line_profiler installed."
            )

        mpi.start_mpi()
        rank = mpi.get_rank()

        beam_opts = {
            "beam_type": None,
            "spline_interp_opts": None,
            "freq_interp_kind": None,
            "peak_normalize": None,
        }
        if rank == root:
            beam_opts = {
                "beam_type": self.beam_type,
                "spline_interp_opts": self.spline_interp_opts,
                "freq_interp_kind": self.freq_interp_kind,
                "peak_normalize": self.peak_normalize,
            }

        beam_opts = mpi.world_comm.bcast(beam_opts, root=root)

        if rank != root:
            for key, val in beam_opts.items():
                object.__setattr__(self, key, val)

        # Get list of beam indices that are set.
        sharelist = None
        if rank == root:
            sharelist = [
                (bi._isuvbeam, bi.beam_type, bi.include_cross_pols) for bi in self
            ]
        sharelist = mpi.world_comm.bcast(sharelist, root=root)

        if rank != root:
            beam_list = []

        for cnt, (uvbeam, beam_type, cross_pols) in enumerate(sharelist):
            mpi.world_comm.Barrier()

            if rank == root:
                bi = self[cnt]
            else:
                if not uvbeam:
                    bi = None
                else:
                    tmp_beam = UVBeam()
                    tmp_beam.beam_type = beam_type

                    bi = BeamInterface(
                        beam=tmp_beam,
                        beam_type=beam_type,
                        include_cross_pols=cross_pols,
                    )

            if not uvbeam:
                bi = mpi.world_comm.bcast(bi, root=root)
            else:
                bi = self._share_uvbeam(bi, root)

            if rank != root:
                beam_list.append(bi)

        if rank != root:
            object.__setattr__(self, "beam_list", beam_list)

        mpi.world_comm.Barrier()


class Telescope:
    """
    Container for data common to all antennas in the array.

    Defines the location and name of the observing site, and holds the list
    of beam objects used by the array.

    Parameters
    ----------
    telescope_name : str
        Name of the telescope.
    telescope_location : :class:`astropy.coordinates.EarthLocation`, :class:`lunarsky.MoonLocation`
        Location of the telescope.
    beam_list : :class:`pyuvsim.BeamList`
        BeamList carrying beam models.

    """

    def __init__(
        self,
        telescope_name: str,
        telescope_location: EarthLocation,
        beam_list: BeamList,
    ):
        allowed_location_types = [EarthLocation]
        if hasmoon:
            allowed_location_types.append(MoonLocation)

        if not isinstance(telescope_location, tuple(allowed_location_types)):
            raise ValueError(
                "location must be an EarthLocation object or a "
                "lunarsky.MoonLocation object"
            )
        self.location = telescope_location
        self.name = telescope_name

        # BeamList object, length of number of unique beams
        if not isinstance(beam_list, BeamList):
            raise ValueError("beam_list must be a BeamList object")

        self.beam_list = beam_list

    def __eq__(self, other):
        """Define Telescope equality."""
        this_vector_loc = np.array(
            [
                self.location.x.to("m").value,
                self.location.y.to("m").value,
                self.location.z.to("m").value,
            ]
        )
        other_vector_loc = np.array(
            [
                other.location.x.to("m").value,
                other.location.y.to("m").value,
                other.location.z.to("m").value,
            ]
        )
        Nbeams_self = len(self.beam_list)
        Nbeams_other = len(other.beam_list)
        if Nbeams_self == Nbeams_other:
            return (
                np.allclose(this_vector_loc, other_vector_loc, atol=1e-3)
                and np.all(
                    [
                        self.beam_list[bi] == other.beam_list[bi]
                        for bi in range(Nbeams_self)
                    ]
                )
                and self.name == other.name
            )
        return False
