# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Definition of Telescope objects, for metadata common to all antennas in an array."""

from __future__ import annotations

import warnings

import numpy as np
from pyuvdata import UVBeam, parameter

try:
    from . import mpi
except ImportError:
    mpi = None

from .analyticbeam import AnalyticBeam

weird_beamfits_extension_warning = (
    "This beamfits file does not have a '.fits' or '.beamfits' extension, so "
    "UVBeam does not recognize it as a beamfits file. Either change the file "
    "extension or specify the beam_type. This is currently handled with a "
    "try/except clause in pyuvsim, but this will become an error in version 1.4"
)


class BeamConsistencyError(Exception):
    """An Exception class to mark inconsistencies among beams in a BeamList."""

    pass


class BeamList:
    """
    A container for the set of beam models and related parameters.

    Each rank in the simulation gets a copy of the set of beam objects
    used for calculating Jones matrices. Rather than broadcasting the objects
    themselves, the list is composed of strings which provide either a complete
    description of an analytic beam or the path to a UVBeam readable file. After
    the broadcast, the beams are initialized from these strings.

    This class provides methods to transform between the string and object
    representations, and also carries parameters used globally by all UVBeams,
    including their frequency / angular interpolation options.

    This behaves just as a normal Python list in terms of indexing, appending, etc.

    Attributes
    ----------
    string_mode : bool
        Is True if the beams are represented as strings, False if they're objects
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
    beam_list : list (optional)
        A list of UVBeam/AnalyticBeam objects or strings describing beams. If
        beam_list consists of objects, then the BeamList will be initialized in
        object mode with string_mode == False. If beam_list consists of strings,
        the BeamList will be initialized in string mode.
        Passing in a mixture of strings and objects will error.
    uvb_read_kwargs : dict (optional)
        A nested dictionary that can contain parameters to pass to the UVBeam
        read method for each UVBeam in the beam dict. The top-level keys are
        beam_id integers, the value for each beam id is a dict with the kwargs
        specified for each beam. This can include ``file_type`` or any other
        parameter accepted by the UVBeam.read method. Note that this can in
        principal overlap with entries in select_params. Entries in uvb_read_kwargs
        will supercede anything that is also in select_params.
    select_params : dict (optional)
        A dictionary that can contain parameters for selecting parts of the beam
        to read. Example keys include ``freq_range`` and ``za_range``. Note that
        these will only be used for UVBeam readable files and they apply to all
        such beams unless they are superceded by uvb_read_kwargs.
    check : bool
        Whether to perform a consistency check on the beams (i.e. asserting that several
        of their defining parameters are the same for all beams in the list).
    force_check : bool
        The consistency check is only possible for object-beams (not string mode). If
        `force_check` is True, it will convert beams to object-mode to force the check
        to run (then convert back to string mode).

    Raises
    ------
    ValueError
        For an invalid beam_list (mix of strings and objects).

    Notes
    -----
    If an object beam is added while in string mode, it will be converted to a string
    before being inserted, and vice versa. ***In order to be converted to a string, a UVBeam
    needs to have the entry 'beam_path' in its extra_keywords, giving the path to the beam
    file it came from.***

    """

    _float_params = {
        "sig": "sigma",
        "diam": "diameter",
        "reff": "ref_freq",
        "ind": "spectral_index",
    }

    string_mode = True

    def __init__(
        self,
        beam_list=None,
        uvb_read_kwargs: dict[str : tuple[float, float]] | None = None,
        select_params: dict[str : tuple[float, float]] | None = None,
        spline_interp_opts: dict[str:int] | None = None,
        freq_interp_kind: str = "cubic",
        check: bool = True,
        force_check: bool = False,
    ):
        self.spline_interp_opts = spline_interp_opts
        self.freq_interp_kind = freq_interp_kind
        self._str_beam_list = []
        self._obj_beam_list = []
        if beam_list is not None:
            if all(isinstance(b, str) for b in beam_list):
                self._str_beam_list[:] = beam_list[:]
                self.string_mode = True
            elif all(not isinstance(b, str) for b in beam_list):
                self._obj_beam_list[:] = beam_list[:]
                for beam in self._obj_beam_list:
                    # always use future shapes
                    if (
                        isinstance(beam, UVBeam)
                        and hasattr(beam, "use_current_array_shapes")
                        and not beam.future_array_shapes
                    ):
                        beam.use_future_array_shapes()
                self.string_mode = False
            else:
                raise ValueError("Invalid beam list: " + str(beam_list))

        self.uvb_read_kwargs = uvb_read_kwargs or {}
        self.select_params = select_params or {}
        self.is_consistent = False
        if check and (not self.string_mode or force_check):
            # only call the check if it will do something. This avoids an annoying warning.
            self.check_consistency(force=force_check)

    def _get_beam_basis_type(self):
        if self.string_mode:
            raise ValueError("Cannot get beam basis type from a string-mode BeamList.")

        beam_basis = {}
        for index, bm in enumerate(self):
            if isinstance(bm, UVBeam):
                if bm.pixel_coordinate_system not in ["az_za", "healpix"]:
                    raise ValueError(
                        "pyuvsim currently only supports UVBeams with 'az_za' or "
                        "'healpix' pixel coordinate systems."
                    )
                if (
                    np.all(bm.basis_vector_array[0, 0] == 1)
                    and np.all(bm.basis_vector_array[0, 1] == 0)
                    and np.all(bm.basis_vector_array[1, 0] == 0)
                    and np.all(bm.basis_vector_array[1, 1] == 1)
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

    def check_consistency(self, force: bool = False):
        """
        Check the consistency of all beams in the list.

        This checks basic parameters of the objects for consistency, eg. the ``beam_type``.
        It is meant to be manually called by the user.

        Parameters
        ----------
        force : bool
            Whether to force the consistency check even if the object is in string
            mode. This wil convert to to object mode, perform the check, then
            convert back.

        Raises
        ------
        BeamConsistencyError
            If any beam is inconsistent with the rest of the beams.

        """
        if self.string_mode and len(self._str_beam_list) > 0:
            if not force:
                warnings.warn(
                    "Cannot check consistency of a string-mode BeamList! Set force=True"
                    " to force consistency checking."
                )
                return
            else:
                self.set_obj_mode(check=True)
                self.set_str_mode()
                return

        if len(self) == 0:
            return

        def check_thing(item):
            if item == "basis":
                items = self._get_beam_basis_type()
            else:
                items = {
                    i: getattr(b, item) for i, b in enumerate(self) if hasattr(b, item)
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

        if self[0].beam_type == "efield":
            check_thing("Nfeeds")
            check_thing("feed_array")
            check_thing("basis")

        elif self[0].beam_type == "power":
            check_thing("Npols")
            check_thing("polarization_array")

        self.is_consistent = True

    def check_all_azza_beams_full_sky(self):
        """
        Check if all az_za UVBeams cover the full sky.

        Used to decide whether to turn off checking that the beam covers the
        interpolation location.
        """
        all_azza_full_sky = True
        for beam in self:
            if isinstance(beam, UVBeam) and beam.pixel_coordinate_system == "az_za":
                # regular gridding is enforced in UVBeam checking, so can use first diff
                axis1_diff = np.diff(beam.axis1_array)[0]
                axis2_diff = np.diff(beam.axis2_array)[0]
                max_axis_diff = np.max([axis1_diff, axis2_diff])
                if not (
                    np.max(beam.axis2_array) >= np.pi / 2.0 - max_axis_diff * 2.0
                    and np.min(beam.axis1_array) <= max_axis_diff * 2.0
                    and np.max(beam.axis1_array) >= 2.0 * np.pi - max_axis_diff * 2.0
                ):
                    all_azza_full_sky = False
                    break
        return all_azza_full_sky

    @property
    def x_orientation(self):
        """
        Return the x_orientation of all beams in list (if consistent).

        Raises
        ------
        BeamConsistencyError
            If the beams in the list are not consistent with each other.

        """
        if not self.is_consistent:
            self.check_consistency(force=True)

        if len(self) == 0:
            return None

        if self.string_mode:
            # Can't get xorientation from beams in string mode.
            # Convert to obj mode, get xorient, and then convert back.
            self.set_obj_mode()
            xorient = self.x_orientation
            self.set_str_mode()
            return xorient

        for b in self:
            if hasattr(b, "x_orientation"):
                xorient = b.x_orientation
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
    def beam_type(self):
        """
        Return the beam_type of all beams in list (if consistent).

        Raises
        ------
        BeamConsistencyError
            If the beams in the list are not consistent with each other.

        """
        if not self.is_consistent:
            self.check_consistency(force=True)

        if len(self) == 0:
            return None

        return getattr(self._obj_beam_list[0], "beam_type", None)

    def __len__(self):
        """Define the length of a BeamList."""
        # Note that only one of these lists has nonzero length at a given time.
        return len(self._obj_beam_list) + len(self._str_beam_list)

    def __iter__(self):
        """Get the list as an iterable."""
        lst = self._str_beam_list if self.string_mode else self._obj_beam_list
        return iter(lst)

    def __getitem__(self, ind):
        """Get a particular beam from the list."""
        if self.string_mode:
            return self._str_beam_list[ind]
        return self._obj_beam_list[ind]

    def __setitem__(self, ind, value, uvb_read_kwargs=None):
        """
        Insert into the beam list.

        Converts objects to strings if in object mode,
        or vice versa if in string mode.

        """
        if self.string_mode:
            self._str_beam_list[ind] = self._obj_to_str(value)
        else:
            if uvb_read_kwargs is not None:
                self.uvb_read_kwargs[ind] = uvb_read_kwargs
            try:
                value = self._str_to_obj(ind, value)
            except FileNotFoundError as err:
                raise ValueError(f"Invalid file path: {value}") from err
            self._obj_beam_list[ind] = value

    def __eq__(self, other):
        """Define BeamList equality."""
        if self.string_mode:
            return self._str_beam_list == other._str_beam_list
        return self._obj_beam_list == other._obj_beam_list

    def append(self, value, uvb_read_kwargs=None):
        """
        Append to the beam list.

        Converts objects to strings if in object mode,
        or vice versa if in string mode.

        """
        if self.string_mode:
            self._str_beam_list.append("")
        else:
            self._obj_beam_list.append("")
        try:
            self.__setitem__(-1, value, uvb_read_kwargs=uvb_read_kwargs)
        except ValueError as err:
            if self.string_mode:
                # This line is hard to cover with a test. I've verified that I
                # can get here, but I when I wrap it in  pytest.raises it stops
                # on the earlier error.
                self._str_beam_list.pop()  # pragma: no cover
            else:
                self._obj_beam_list.pop()
            raise err

    def _str_to_obj(self, beam_id, beam_model, use_shared_mem=False):
        """Convert beam strings to objects."""
        if isinstance(beam_model, AnalyticBeam | UVBeam):
            return beam_model
        if beam_model.startswith("analytic"):
            bspl = beam_model.split("_")
            model = bspl[1]

            to_set = {}
            for extra in bspl[2:]:
                par, val = extra.split("=")
                full = self._float_params[par]
                to_set[full] = float(val)

            return AnalyticBeam(model, **to_set)

        if use_shared_mem and mpi is None:
            raise ImportError(
                "You need mpi4py to use shared memory. Either call this method with "
                "use_shared_mem=False or install mpi4py. You can install it by running "
                "pip install pyuvsim[sim] or pip install pyuvsim[all] if you also want "
                "the line_profiler installed."
            )

        path = beam_model  # beam_model = path to UVBeam readable file
        uvb = UVBeam()
        read_kwargs = {**self.select_params, **self.uvb_read_kwargs.get(beam_id, {})}

        # always use future shapes if it's an option
        if hasattr(uvb, "use_current_array_shapes"):
            read_kwargs["use_future_array_shapes"] = True

        if (
            (use_shared_mem and mpi.world_comm is not None and mpi.rank == 0)
            or not use_shared_mem
            or mpi is None
            or mpi.world_comm is None
        ):
            try:
                uvb.read(path, **read_kwargs)
            except ValueError:
                # If file type is not recognized, assume beamfits,
                # which was originally the only option.
                uvb.read_beamfits(path, **read_kwargs)
                warnings.warn(weird_beamfits_extension_warning, DeprecationWarning)
            uvb.peak_normalize()

        if use_shared_mem and (mpi.world_comm is not None):
            for key, attr in uvb.__dict__.items():
                if not isinstance(attr, parameter.UVParameter):
                    continue
                if key == "_data_array":
                    uvb.__dict__[key].value = mpi.shared_mem_bcast(attr.value, root=0)
                else:
                    uvb.__dict__[key].value = mpi.world_comm.bcast(attr.value, root=0)
            mpi.world_comm.Barrier()
        uvb.extra_keywords["beam_path"] = path
        return uvb

    def _obj_to_str(self, beam_model):
        """Convert beam objects to strings that may generate them."""
        if isinstance(beam_model, str):
            return beam_model
        if isinstance(beam_model, AnalyticBeam):
            btype = beam_model.type

            bm_str = "analytic_" + btype
            for abbrv, full in self._float_params.items():
                val = getattr(beam_model, full)
                if val is not None:
                    bm_str += "_" + abbrv + "=" + str(val)
            return bm_str

        # If not AnalyticBeam, it's UVBeam.
        try:
            path = beam_model.extra_keywords["beam_path"]
        except KeyError as ke:
            raise ValueError(
                "Need to set 'beam_path' key in extra_keywords for UVBeam objects."
            ) from ke

        return path

    def set_str_mode(self):
        """
        Convert object beams to their strings representations.

        Convert AnalyticBeam objects to a string representation
        that can be used to generate them. Replaces UVBeam objects with
        the path to the UVBeam readable file they originate from.

        Any UVBeams in the list must have the path to their beam files
        stored as 'beam_path' in the `extra_keywords` dictionary.

        Sets :attr:`~.string_mode` to True.

        Raises
        ------
        ValueError:
            If UVBeam objects in the list have inconsistent UVBeam parameters.
        ValueError:
            If any UVBeam objects lack a "beam_path" keyword in the "extra_keywords"
            attribute, specifying the path to the UVBeam readable file that generates it.

        """
        if self._obj_beam_list != []:
            # Convert object beams to string definitions
            self._str_beam_list = [
                self._obj_to_str(bobj) for bobj in self._obj_beam_list
            ]
        self._obj_beam_list = []
        self.string_mode = True

    def set_obj_mode(self, use_shared_mem: bool = False, check: bool = False):
        """
        Initialize AnalyticBeam and UVBeam objects from string representations.

        Sets :attr:`~.string_mode` to False.

        Parameters
        ----------
        use_shared_mem : bool
            Whether to use shared memory for the beam data.
        check : bool
            Whether to perform consistency checks on all the beams in the list after
            conversion to object mode.

        """
        if self._str_beam_list != []:
            beam_list = []
            for bid, bstr in enumerate(self._str_beam_list):
                beam_list.append(
                    self._str_to_obj(bid, bstr, use_shared_mem=use_shared_mem)
                )
            self._obj_beam_list = beam_list
        self._str_beam_list = []
        self.string_mode = False

        if check:
            self.check_consistency()


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

    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.location = telescope_location
        self.name = telescope_name

        # list of UVBeam objects, length of number of unique beams
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
