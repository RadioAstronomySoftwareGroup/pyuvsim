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
    """
    float_params = {'sig': 'sigma', 'diam': 'diameter',
                    'reff': 'ref_freq', 'ind': 'spectral_index'}

    uvb_params = {'freq_interp_kind': 'cubic',
                  'interpolation_function': 'az_za_simple'}

    string_mode = True

    def __init__(self, beam_list=None, uvb_params={}):
        self._str_beam_list = []
        self._obj_beam_list = []
        if beam_list is not None:
            if all([isinstance(b, str) for b in beam_list]):
                self._str_beam_list = beam_list
                self.string_mode = True
            elif all([not isinstance(b, str) for b in beam_list]):
                self._obj_beam_list = beam_list
                self.string_mode = False
            else:
                raise ValueError("Invalid beam list: " + str(beam_list))

        self.uvb_params.update(uvb_params)    # Optional parameters for the UVBeam objects

    def append(self, value):
        if self.string_mode:
            self._str_beam_list.append(self._obj_to_str(value))
        self._obj_beam_list.append(self._str_to_obj(value))

    def __len__(self):
        return len(self._obj_beam_list) + len(self._str_beam_list)

    def __getitem__(self, ind):
        if self.string_mode:
            return self._str_beam_list[ind]
        return self._obj_beam_list[ind]

    def __setitem__(self, ind, value):
        if self.string_mode:
            self._str_beam_list[ind] = self._obj_to_str(value)
        self._obj_beam_list[ind] = self._str_to_obj(value)

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
                full = self.float_params[par]
                if (par == 'sig') and ('diameter' in to_set.keys()):
                    continue    # Do not overwrite if diameter is already set.
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
            for abbrv, full in self.float_params.items():
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
        if not self._obj_beam_list == []:
            # Convert object beams to string definitions
            self._str_beam_list = [self._obj_to_str(bobj) for bobj in self._obj_beam_list]
        self._obj_beam_list = []
        self.string_mode = True

    def set_obj_mode(self):
        if not self._str_beam_list == []:
            self._obj_beam_list = [self._str_to_obj(bstr) for bstr in self._str_beam_list]
        self._str_beam_list = []
        self.string_mode = False
