# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
from pyuvdata import UVBeam

from .analyticbeam import AnalyticBeam


class Telescope(object):
    """
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
        For invalid beam_list, or if UVBeams in the input beam_list have
        conflicting parameters (e.g., one has freq_interp_kind == 'cubic' and
        another has freq_interp_kind == 'linear').

    Notes
    -----
    When a UVBeam is appended to the beam list, its parameters override any in
    uvb_params.

    If an object beam is added while in string mode, it will be converted to a string
    before being inserted, and vice versa. ***In order to be converted to a string, a UVBeam
    needs to have the entry 'beam_path' in its extra_keywords, giving the path to the beamfits
    file it came from.***

    """

    _float_params = {'sig': 'sigma', 'diam': 'diameter',
                     'reff': 'ref_freq', 'ind': 'spectral_index'}

    string_mode = True

    def __init__(self, beam_list=None, uvb_params={}):

        self.uvb_params = {'freq_interp_kind': 'cubic',
                           'interpolation_function': 'az_za_simple'}

        self._str_beam_list = []
        self._obj_beam_list = []
        if beam_list is not None:
            if all([isinstance(b, str) for b in beam_list]):
                self._str_beam_list[:] = beam_list[:]
                self.string_mode = True
            elif all([not isinstance(b, str) for b in beam_list]):
                self._obj_beam_list[:] = beam_list[:]
                self.string_mode = False
            else:
                raise ValueError("Invalid beam list: " + str(beam_list))

        # If any UVBeam objects are passed in, update uvb_params:
        self._update_uvb_params()
        self.uvb_params.update(uvb_params)    # Optional parameters for the UVBeam objects

    def _update_uvb_params(self, ind=slice(None)):
        """
        Updates the uvb_params according to the parameters
        of any UVBeams in the list.

        This will check to ensure that the parameters are consistent
        among UVBeams.

        Parameters
        ----------
        ind : int or slice (optional)
            Choose a beam object, or set of them, to read for parameters.

        """

        if self.string_mode:
            return

        new_pars = {k: [] for k in self.uvb_params.keys()}
        for bm in np.atleast_1d(self._obj_beam_list[ind]):
            if isinstance(bm, UVBeam):
                for key in self.uvb_params.keys():
                    val = getattr(bm, key)
                    new_pars[key].append(val)
        for key, vals in new_pars.items():
            svals = set(vals)
            if len(svals) == 1:
                self.uvb_params[key] = svals.pop()
            elif len(svals) > 1:
                raise ValueError('Conflicting settings for {}: {}'.format(key, str(svals)))

    def append(self, value):
        """
        Append to the beam list, converting objects to strings if in object mode,
        or vice versa if in string mode.
        """
        if self.string_mode:
            self._str_beam_list.append(self._obj_to_str(value))
            self._update_uvb_params(ind=-1)
        self._obj_beam_list.append(self._str_to_obj(value))

    def __len__(self):
        # Note that only one of these lists has nonzero length at a given time.
        return len(self._obj_beam_list) + len(self._str_beam_list)

    def __getitem__(self, ind):
        if self.string_mode:
            return self._str_beam_list[ind]
        return self._obj_beam_list[ind]

    def __setitem__(self, ind, value):
        """
        Insert into the beam list, converting objects to strings if in object mode,
        or vice versa if in string mode.
        """
        if self.string_mode:
            self._str_beam_list[ind] = self._obj_to_str(value)
        else:
            self._obj_beam_list[ind] = self._str_to_obj(value)
            self._update_uvb_params(ind=ind)

    def _str_to_obj(self, beam_model):
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
        path = beam_model.extra_keywords['beam_path']

        for key, val in self.uvb_params.items():
            self.uvb_params[key] = getattr(beam_model, key)
        return path

    def set_str_mode(self):
        """
        Set string_mode = True, converting object beams to their
        string representations and putting UVBeam parameters into the
        uvb_params dictionary.
        """
        if not self._obj_beam_list == []:
            # Convert object beams to string definitions
            self._str_beam_list = [self._obj_to_str(bobj) for bobj in self._obj_beam_list]
        self._obj_beam_list = []
        self.string_mode = True

    def set_obj_mode(self):
        """
        Set string_mode = False, initializing beam objects from the strings.
        Additional parameters on UVBeam objects will be set from the uvb_params dictionary.
        """
        if not self._str_beam_list == []:
            self._obj_beam_list = [self._str_to_obj(bstr) for bstr in self._str_beam_list]
        self._update_uvb_params()
        self._str_beam_list = []
        self.string_mode = False
