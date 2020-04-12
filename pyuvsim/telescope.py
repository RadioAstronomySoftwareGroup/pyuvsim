# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
import warnings
from pyuvdata import UVBeam, parameter

from .analyticbeam import AnalyticBeam
try:
    from . import mpi
except ImportError:
    mpi = None


class Telescope(object):
    """
    Container for data common to all antennas in the array.

    Defines the location and name of the observing site, and holds the list
    of beam objects used by the array.
    """

    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.location = telescope_location
        self.name = telescope_name

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list

    def __eq__(self, other):
        this_vector_loc = np.array([self.location.x.to('m').value,
                                    self.location.y.to('m').value,
                                    self.location.z.to('m').value])
        other_vector_loc = np.array([other.location.x.to('m').value,
                                     other.location.y.to('m').value,
                                     other.location.z.to('m').value])
        Nbeams_self = len(self.beam_list)
        Nbeams_other = len(other.beam_list)
        if Nbeams_self == Nbeams_other:
            return (
                np.allclose(this_vector_loc, other_vector_loc, atol=1e-3)
                and np.all(
                    [self.beam_list[bi] == other.beam_list[bi] for bi in range(Nbeams_self)]
                )
                and self.name == other.name
            )
        return False


class BeamList(object):
    """
    A container for the set of beam models and related parameters.

    Each rank in the simulation gets a copy of the set of beam objects
    used for calculating Jones matrices. Rather than broadcasting the objects
    themselves, the list is composed of strings which provide either a complete
    description of an analytic beam or the path to a beamfits file. After the
    broadcast, the beams are initialized from these strings.

    This class provides methods to transform between the string and object
    representations, and also carries parameters used globally by all UVBeams,
    including their frequency / angular interpolation options.

    This behaves just as a normal Python list in terms of indexing, appending, etc.

    Attributes
    ----------
    string_mode : bool
        Is True if the beams are represented as strings, False if they're objects

    uvb_params : dict
        Set of additional attributes to set on UVBeam objects.

    Parameters
    ----------
    beam_list : list (optional)
        A list of UVBeam/AnalyticBeam objects or strings describing beams. If
        beam_list consists of objects, then the BeamList will be initialized in
        object mode with string_mode == False. If beam_list consists of strings,
        the BeamList will be initialized in string mode.

        Passing in a mixture of strings and objects will error.

    uvb_params : dict (optional)
        Options to set uvb_params, overriding settings from passed-in UVBeam objects.

    Raises
    ------
    ValueError :
        For an invalid beam_list (mix of strings and objects).

    Notes
    -----
    If an object beam is added while in string mode, it will be converted to a string
    before being inserted, and vice versa. ***In order to be converted to a string, a UVBeam
    needs to have the entry 'beam_path' in its extra_keywords, giving the path to the beamfits
    file it came from.***

    """

    _float_params = {'sig': 'sigma', 'diam': 'diameter',
                     'reff': 'ref_freq', 'ind': 'spectral_index'}

    string_mode = True

    def __init__(self, beam_list=None, uvb_params=None):

        self.uvb_params = {'freq_interp_kind': 'cubic',
                           'interpolation_function': 'az_za_simple'}
        self.spline_interp_opts = None
        self._str_beam_list = []
        self._obj_beam_list = []
        if beam_list is not None:
            if all(isinstance(b, str) for b in beam_list):
                self._str_beam_list[:] = beam_list[:]
                self.string_mode = True
            elif all(not isinstance(b, str) for b in beam_list):
                self._obj_beam_list[:] = beam_list[:]
                self.string_mode = False
            else:
                raise ValueError("Invalid beam list: " + str(beam_list))

        # If any UVBeam objects are passed in, update uvb_params:
        list_uvb_params = {}
        if uvb_params is None:
            uvb_params = {}
            if not self.string_mode:
                list_uvb_params = self._scrape_uvb_params(self._obj_beam_list, strict=False)
        self.uvb_params.update(list_uvb_params)
        self.uvb_params.update(uvb_params)    # Optional parameters for the UVBeam objects

    def _scrape_uvb_params(self, beam_objs, strict=True):
        """
        Collect uvb_params from the set of beam objects.

        Parameters
        ----------
        beam_objs : list
            If any UVBeams are found in this list, they will be scraped
            for relevant parameters.

        strict : bool
            If True, will raise an error if the UVBeams in beam_objs
            have conflicting parameters.

        Notes
        -----
        If beam_objs has strings in it, this will simply return an empty dict.
        """

        new_pars = {}
        for bm in beam_objs:
            if isinstance(bm, UVBeam):
                for key in self.uvb_params.keys():
                    val = getattr(bm, key)
                    if key not in new_pars:
                        new_pars[key] = []
                    if val is not None:
                        new_pars[key].append(val)
        for key, vals in new_pars.items():
            svals = set(vals)
            if len(svals) == 1:
                new_pars[key] = svals.pop()
            elif len(svals) > 1:
                try:
                    raise ValueError('Conflicting settings for {}: {}'.format(key, str(svals)))
                except ValueError as err:
                    if strict:
                        raise err
                    else:
                        warnings.warn(err.args[0])  # Use exception message as warning message.

        new_pars = {k: v for k, v in new_pars.items() if not v == []}

        return new_pars

    def _set_params_on_uvbeams(self, beam_objs):
        """Set parameters on UVBeam objects using uvb_params."""
        for bm in np.atleast_1d(beam_objs):
            if isinstance(bm, UVBeam):
                for key, val in self.uvb_params.items():
                    setattr(bm, key, val)

    def __len__(self):
        # Note that only one of these lists has nonzero length at a given time.
        return len(self._obj_beam_list) + len(self._str_beam_list)

    def __getitem__(self, ind):
        if self.string_mode:
            return self._str_beam_list[ind]
        return self._obj_beam_list[ind]

    def __setitem__(self, ind, value):
        """
        Insert into the beam list.

        Converts objects to strings if in object mode,
        or vice versa if in string mode.
        """
        if self.string_mode:
            # If value is a string, then _scrape_uvb_params returns an empty dictionary.
            newobj_params = self._scrape_uvb_params([value], strict=False)
            mismatches = []
            for key, val in newobj_params.items():
                if key in self.uvb_params:
                    if not val == self.uvb_params[key]:
                        mismatches.append(key)
            if len(mismatches) > 0:
                compare = "\n\t".join([
                    "current = {}, new = {}".format(self.uvb_params[key], newobj_params[key])
                    for key in mismatches]
                )
                raise ValueError('UVBeam parameters do not match '
                                 'those currently set: \n\t' + compare)

            self._str_beam_list[ind] = self._obj_to_str(value)
        else:
            value = self._str_to_obj(value)
            self._obj_beam_list[ind] = value
            self._scrape_uvb_params(self._obj_beam_list, strict=False)

    def append(self, value):
        """
        Append to the beam list.

        Converts objects to strings if in object mode,
        or vice versa if in string mode.
        """
        if self.string_mode:
            self._str_beam_list.append('')
        else:
            self._obj_beam_list.append('')
        try:
            self.__setitem__(-1, value)
        except ValueError as err:
            if self.string_mode:
                self._str_beam_list.pop()
            else:
                self._obj_beam_list.pop()
            raise err

    def _str_to_obj(self, beam_model, use_shared_mem=False):
        # Convert beam strings to objects.
        if isinstance(beam_model, (AnalyticBeam, UVBeam)):
            return beam_model
        if beam_model.startswith('analytic'):
            bspl = beam_model.split('_')
            model = bspl[1]

            to_set = {}
            for extra in bspl[2:]:
                par, val = extra.split('=')
                full = self._float_params[par]
                to_set[full] = float(val)

            return AnalyticBeam(model, **to_set)

        path = beam_model  # beam_model = path to beamfits
        uvb = UVBeam()
        if use_shared_mem and (mpi is not None):
            if mpi.world_comm is None:
                mpi.start_mpi()
            if mpi.rank == 0:
                uvb.read_beamfits(path)
            for key, attr in uvb.__dict__.items():
                if not isinstance(attr, parameter.UVParameter):
                    continue
                if key == '_data_array':
                    uvb.__dict__[key].value = mpi.shared_mem_bcast(attr.value, root=0)
                else:
                    uvb.__dict__[key].value = mpi.world_comm.bcast(attr.value, root=0)
            mpi.world_comm.Barrier()
        else:
            uvb.read_beamfits(path)
        for key, val in self.uvb_params.items():
            setattr(uvb, key, val)
        uvb.extra_keywords['beam_path'] = path
        return uvb

    def _obj_to_str(self, beam_model):
        # Convert beam objects to strings that may generate them.

        if isinstance(beam_model, str):
            return beam_model
        if isinstance(beam_model, AnalyticBeam):
            btype = beam_model.type

            bm_str = 'analytic_' + btype
            for abbrv, full in self._float_params.items():
                val = getattr(beam_model, full)
                if val is not None:
                    bm_str += '_' + abbrv + '=' + str(val)
            return bm_str

        # If not AnalyticBeam, it's UVBeam.
        try:
            path = beam_model.extra_keywords['beam_path']
        except KeyError:
            raise ValueError("Need to set 'beam_path' key in extra_keywords for UVBeam objects.")

        return path

    def set_str_mode(self):
        """
        Convert object beams to their strings representations.

        Convert AnalyticBeam objects to a string representation
        that can be used to generate them. Replaces UVBeam objects with
        the path to the beamfits file they originate from. Additional
        UVBeam attributes are stored in the BeamList.uvb_params dictionary.

        Any UVBeams in the list must have the path to their beamfits files
        stored as 'beam_path' in the `extra_keywords` dictionary.

        Sets
        ----
        string_mode : bool
            Sets to True

        Raises
        ------
        ValueError:
            If UVBeam objects in the list have inconsistent UVBeam parameters.
        ValueError:
            If any UVBeam objects lack a "beam_path" keyword in the "extra_keywords"
            attribute, specifying the path to the beamfits file that generates it.
        """
        if not self._obj_beam_list == []:
            # Convert object beams to string definitions
            new_uvb_par = self._scrape_uvb_params(self._obj_beam_list, strict=True)
            self.uvb_params.update(new_uvb_par)
            self._str_beam_list = [self._obj_to_str(bobj) for bobj in self._obj_beam_list]
        self._obj_beam_list = []
        self.string_mode = True

    def set_obj_mode(self, use_shared_mem=False):
        """
        Initialize AnalyticBeam and UVBeam objects from string representations.

        For any UVBeams, additional attributes will be set by the BeamList.uvb_params
        dictionary. This overwrites any settings that may have been read from file.

        Sets
        ----
        string_mode : bool
            Sets to False
        """
        if not self._str_beam_list == []:
            self._obj_beam_list = [self._str_to_obj(bstr, use_shared_mem=use_shared_mem)
                                   for bstr in self._str_beam_list]
        self._set_params_on_uvbeams(self._obj_beam_list)
        self._str_beam_list = []
        self.string_mode = False
