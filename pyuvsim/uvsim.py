# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import os
import sys
import copy
import six
import yaml
from six.moves import map, range, zip
import warnings
import astropy.constants as const
import astropy.units as units
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import EarthLocation

from pyuvdata import UVData, UVBeam
import pyuvdata.utils as uvutils

from . import profiling
from .antenna import Antenna
from .baseline import Baseline
from .telescope import Telescope
from . import utils as simutils
from . import simsetup
from . import mpi

__all__ = ['UVTask', 'UVEngine', 'uvdata_to_task_iter', 'run_uvsim', 'run_uvdata_uvsim', 'init_uvdata_out', 'serial_gather']


class UVTask(object):
    # holds all the information necessary to calculate a single src, t, f, bl, array
    # need the array because we need an array location for mapping to local alt/az

    def __init__(self, source, time, freq, baseline, telescope):
        self.time = time
        self.freq = freq
        self.source = source
        self.baseline = baseline
        self.telescope = telescope
        self.visibility_vector = None
        self.uvdata_index = None        # Where to add the visibility in the uvdata object.

        if isinstance(self.time, float):
            self.time = Time(self.time, format='jd')
        if isinstance(self.freq, float):
            self.freq = self.freq * units.Hz

    def __eq__(self, other):
        return (np.isclose(self.time.jd, other.time.jd, atol=1e-4)
                and np.isclose(self.freq.value, other.freq.value, atol=1e-4)
                and (self.source == other.source)
                and (self.baseline == other.baseline)
                and (self.visibility_vector == other.visibility_vector)
                and (self.uvdata_index == other.uvdata_index)
                and (self.telescope == other.telescope))

    def __gt__(self, other):
        blti0, _, fi0 = self.uvdata_index
        blti1, _, fi1 = other.uvdata_index
        if self.baseline == other.baseline:
            if fi0 == fi1:
                return blti0 > blti1
            return fi0 > fi1
        return self.baseline > other.baseline

    def __ge__(self, other):
        blti0, _, fi0 = self.uvdata_index
        blti1, _, fi1 = other.uvdata_index
        if self.baseline == other.baseline:
            if fi0 == fi1:
                return blti0 >= blti1
            return fi0 >= fi1
        return self.baseline >= other.baseline

    def __lt__(self, other):
        return not self.__ge__(other)

    def __le__(self, other):
        return not self.__gt__(other)


class UVEngine(object):

    def __init__(self, task=None, reuse_spline=True):
        self.reuse_spline = reuse_spline  # Reuse spline fits in beam interpolation
        if task is not None:
            self.set_task(task)

    def set_task(self, task):
        self.task = task

    def apply_beam(self):
        """ Get apparent coherency from jones matrices and source coherency. """
        baseline = self.task.baseline
        source = self.task.source
        # coherency is a 2x2 matrix
        # [ |Ex|^2, Ex* Ey, Ey* Ex |Ey|^2 ]
        # where x and y vectors along the local alt/az axes.

        # Apparent coherency gives the direction and polarization dependent baseline response to a source.
        beam1_jones = baseline.antenna1.get_beam_jones(self.task.telescope,
                                                       source.get_alt_az(self.task.time,
                                                                         self.task.telescope.location),
                                                       self.task.freq, reuse_spline=self.reuse_spline)
        beam2_jones = baseline.antenna2.get_beam_jones(self.task.telescope,
                                                       source.get_alt_az(self.task.time,
                                                                         self.task.telescope.location),
                                                       self.task.freq, reuse_spline=self.reuse_spline)
        this_apparent_coherency = np.dot(beam1_jones,
                                         source.coherency_calc(self.task.time,
                                                               self.task.telescope.location))
        this_apparent_coherency = np.dot(this_apparent_coherency,
                                         (beam2_jones.conj().T))

        self.apparent_coherency = this_apparent_coherency

    def make_visibility(self):
        """ Visibility contribution from a single source """
        assert(isinstance(self.task.freq, Quantity))

        pos_lmn = self.task.source.pos_lmn(self.task.time, self.task.telescope.location)
        if pos_lmn is None:
            return np.array([0., 0., 0., 0.], dtype=np.complex128)

        self.apply_beam()

        # need to convert uvws from meters to wavelengths
        uvw_wavelength = self.task.baseline.uvw / const.c * self.task.freq.to('1/s')
        fringe = np.exp(2j * np.pi * np.dot(uvw_wavelength, pos_lmn))
        vij = self.apparent_coherency * fringe

        # Reshape to be [xx, yy, xy, yx]
        vis_vector = [vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]]
        return np.array(vis_vector)


def uvdata_to_task_iter(task_ids, input_uv, catalog, beam_list, beam_dict):
    """
    Generate local tasks, reusing quantities where possible.

    Args:
        task_ids (numpy.ndarray of ints): Task index in the full flattened meshgrid of parameters.
        input_uv (UVData): UVData object to use
        catalog: array of Source objects
        beam_list: (list of UVBeam or AnalyticBeam objects
        beam_dict (dict, optional): dict mapping antenna number to beam index in beam_list

    Yields:
        Iterable of task objects to be done on current rank.
    """

    # Loops, outer to inner: times, sources, frequencies, baselines
    # The task_ids refer to tasks on the flattened meshgrid.
    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    if not isinstance(catalog, np.ndarray):
        raise TypeError("sources must be a numpy array")

    # There will always be relatively few antennas, so just build the full list.
    antenna_names = input_uv.antenna_names
    antennas = []
    antpos_enu, antnums = input_uv.get_ENU_antpos()
    for num, antname in enumerate(antenna_names):
        if beam_dict is None:
            beam_id = 0
        else:
            beam_id = beam_dict[antname]
        antennas.append(Antenna(antname, num, antpos_enu[num], beam_id))

    baselines = {}
    Ntimes = input_uv.Ntimes
    Nfreqs = input_uv.Nfreqs
    Nsrcs = len(catalog)
    Nbls = input_uv.Nbls

    prev_src_ind = None
    prev_time_ind = None
    prev_freq_ind = None

    telescope = Telescope(input_uv.telescope_name,
                          EarthLocation.from_geocentric(*input_uv.telescope_location, unit='m'),
                          beam_list)

    freq_array = input_uv.freq_array * units.Hz
    time_array = Time(input_uv.time_array, scale='utc', format='jd', location=telescope.location)

    # Shape indicates slowest to fastest index. (time is slowest, baselines is fastest).
    for task_index in task_ids:
        time_i, src_i, freq_i, bl_i = np.unravel_index(task_index, (Ntimes, Nsrcs, Nfreqs, Nbls))
        blti = bl_i + time_i * Nbls   # baseline is the fast axis

        # We reuse a lot of baseline info, so make the baseline list on the first go and reuse.
        if bl_i not in baselines.keys():
            antnum1 = input_uv.ant_1_array[blti]
            antnum2 = input_uv.ant_2_array[blti]
            index1 = np.where(input_uv.antenna_numbers == antnum1)[0][0]
            index2 = np.where(input_uv.antenna_numbers == antnum2)[0][0]
            baselines[bl_i] = Baseline(antennas[index1], antennas[index2])

        source = catalog[src_i]

        if source.rise_lst and source.set_lst:
            r_lst = source.rise_lst
            s_lst = source.set_lst
            now_lst = input_uv.lst_array[blti]

            # Compare time since rise to time between rise and set,
            # properly accounting for phase wrap.
            dt0 = now_lst - r_lst       # radians since rise time.
            if dt0 < 0:
                dt0 += 2 * np.pi

            dt1 = s_lst - r_lst         # radians between rise and set
            if dt1 < 0:
                dt1 += 2 * np.pi

            if dt1 < dt0:
                continue

        time = time_array[blti]
        bl = baselines[bl_i]
        freq = freq_array[0, freq_i]  # 0 = spw axis

        task = UVTask(source, time, freq, bl, telescope)
        task.uvdata_index = (blti, 0, freq_i)    # 0 = spectral window index

        prev_src_ind = src_i
        prev_time_ind = time_i
        prev_freq_ind = freq_i

        yield task


def init_uvdata_out(uv_in, source_list_name,
                    obs_param_file=None, telescope_config_file=None,
                    antenna_location_file=None):
    """
    Initialize an empty uvdata object to fill with simulated data.
    Args:
        uv_in: The input uvdata object.
               This is usually an incomplete object, containing only metadata.
        source_list_name: Name of source list file or mock catalog.
        obs_param_file: Name of observation parameter config file
        telescope_config_file: Name of telescope config file
        antenna_location_file: Name of antenna location file
    """
    if not isinstance(source_list_name, str):
        raise ValueError('source_list_name must be a string')

    if not isinstance(obs_param_file, str):
        raise ValueError('obs_param_file must be a string')
    if not isinstance(telescope_config_file, str):
        raise ValueError('telescope_config_file must be a string')
    if not isinstance(antenna_location_file, str):
        raise ValueError('antenna_location_file must be a string')

    # Version string to add to history
    history = simutils.get_version_string()

    history += ' Sources from source list: ' + source_list_name + '.'

    history += (' Based on config files: ' + obs_param_file + ', '
                + telescope_config_file + ', ' + antenna_location_file)

    history += ' Npus = ' + str(mpi.Npus) + '.'

    uv_obj = copy.deepcopy(uv_in)

    uv_obj.set_drift()
    uv_obj.vis_units = 'Jy'
    uv_obj.polarization_array = np.array([-5, -6, -7, -8])
    uv_obj.instrument = uv_obj.telescope_name
    uv_obj.set_lsts_from_time_array()
    uv_obj.spw_array = np.array([0])
    if uv_obj.Nfreqs == 1:
        uv_obj.channel_width = 1.  # Hz
    else:
        uv_obj.channel_width = np.diff(uv_obj.freq_array[0])[0]
    uv_obj.set_uvws_from_antenna_positions()
    if uv_obj.Ntimes == 1:
        uv_obj.integration_time = np.ones_like(uv_obj.time_array, dtype=np.float64)  # Second
    else:
        # Note: currently only support a constant spacing of times
        uv_obj.integration_time = (np.ones_like(uv_obj.time_array, dtype=np.float64)
                                   * np.diff(np.unique(uv_obj.time_array))[0] * (24. * 60**2))  # Seconds
    # add pyuvdata version info
    history += uv_obj.pyuvdata_version_str

    # Clear existing data, if any
    uv_obj.data_array = np.zeros((uv_obj.Nblts, uv_obj.Nspws, uv_obj.Nfreqs, uv_obj.Npols), dtype=np.complex)
    uv_obj.flag_array = np.zeros((uv_obj.Nblts, uv_obj.Nspws, uv_obj.Nfreqs, uv_obj.Npols), dtype=bool)
    uv_obj.nsample_array = np.ones_like(uv_obj.data_array, dtype=float)
    uv_obj.history = history

    uv_obj.extra_keywords = {}
    uv_obj.check()

    return uv_obj


def serial_gather(uvtask_list, uv_out):
    """
        Loop over uvtask list, acquire visibilities and add to uvdata object.
    """
    for task in uvtask_list:
        blt_ind, spw_ind, freq_ind = task.uvdata_index
        uv_out.data_array[blt_ind, spw_ind, freq_ind, :] += task.visibility_vector

    return uv_out


def run_uvdata_uvsim(input_uv, beam_list, beam_dict=None, catalog=None, source_list_name=None,
                     obs_param_file=None,
                     telescope_config_file=None, antenna_location_file=None):
    """
    Run uvsim from UVData object.

    Arguments:
        input_uv: An input UVData object, containing baseline/time/frequency information.
        beam_list: A list of UVBeam and/or AnalyticBeam identifier strings.

    Keywords:
        beam_dict: Dictionary of {antenna_name : beam_ID}, where beam_id is an index in
                   the beam_list. This assigns beams to antennas.
                   Default: All antennas get the 0th beam in the beam_list.
        source_list_name: Catalog identifier string for file metadata. Only required on rank 0
        catalog: array of source.Source objects
        obs_param_file: Parameter filename if running from config files.
        telescope_config_file: Telescope configuration file if running from config files.
        antenna_location_file: antenna_location file if running from config files.
    """
    mpi.start_mpi()
    rank = mpi.get_rank()
    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")
    # The Head node will initialize our simulation
    # Read input file and make uvtask list
    if rank == 0:
        print('Nbls:', input_uv.Nbls)
        print('Ntimes:', input_uv.Ntimes)
        print('Nfreqs:', input_uv.Nfreqs)
        print('Nsrcs:', len(catalog))

        if 'obs_param_file' in input_uv.extra_keywords:
            obs_param_file = input_uv.extra_keywords['obs_param_file']
            telescope_config_file = input_uv.extra_keywords['telescope_config_file']
            antenna_location_file = input_uv.extra_keywords['antenna_location_file']
        else:
            obs_param_file = ''
            telescope_config_file = ''
            antenna_location_file = ''

        uv_container = init_uvdata_out(input_uv, source_list_name,
                                       obs_param_file=obs_param_file,
                                       telescope_config_file=telescope_config_file,
                                       antenna_location_file=antenna_location_file)

    comm = mpi.get_comm()
    Npus = mpi.get_Npus()

    # Ensure all ranks have finished setup.
    # comm.Barrier()

    Nblts = input_uv.Nblts
    Nbls = input_uv.Nbls
    Ntimes = input_uv.Ntimes
    Nfreqs = input_uv.Nfreqs
    Nsrcs = len(catalog)

    Ntasks = Nblts * Nfreqs * Nsrcs

    Neach_section, extras = divmod(Ntasks, Npus)
    if rank < extras:
        length = Neach_section + 1
        start = rank * (length)
        end = start + length
    else:
        length = Neach_section
        start = extras * (Neach_section + 1) + (rank - extras) * length
        end = start + length
    task_ids = six.moves.range(start, end)

    # Construct beam objects from strings
    beam_models = [simsetup.beam_string_to_object(bm) for bm in beam_list]

    Ntasks_local = (end - start)

    local_task_iter = uvdata_to_task_iter(task_ids, input_uv, catalog, beam_models, beam_dict)

    summed_task_dict = {}

    if rank == 0:
        tot = Ntasks
        print("Tasks: ", tot)
        sys.stdout.flush()
        pbar = simutils.progsteps(maxval=tot)

    engine = UVEngine()
    for count, task in enumerate(local_task_iter):
        engine.set_task(task)
        if task.uvdata_index not in summed_task_dict.keys():
            summed_task_dict[task.uvdata_index] = task
        if summed_task_dict[task.uvdata_index].visibility_vector is None:
            summed_task_dict[task.uvdata_index].visibility_vector = engine.make_visibility()
        else:
            summed_task_dict[task.uvdata_index].visibility_vector += engine.make_visibility()
        if rank == 0:
            pbar.update(count * mpi.get_Npus())
    if rank == 0:
        pbar.finish()

    if rank == 0:
        print("Calculations Complete.")

    # All the sources in this summed list are foobar-ed
    # Source are summed over but only have 1 name
    # Some source may be correct
    summed_local_task_list = list(summed_task_dict.values())
    # Tasks contain attributes that are not pickle-able.
    # Remove everything except uvdata_index and visibility_vector
    for task in summed_local_task_list:
        del task.time
        del task.freq
        del task.source
        del task.baseline
        del task.telescope

    # gather all the finished local tasks into a list of list of len NPUs
    # gather is a blocking communication, have to wait for all PUs
    full_tasklist = comm.gather(summed_local_task_list, root=0)
    localtasks_count = comm.gather(Ntasks_local, root=0)
    # Concatenate the list of lists into a flat list of tasks
    if rank == 0:
        localtasks_count = np.sum(localtasks_count)
        uvtask_list = sum(full_tasklist, [])
        uvdata_out = serial_gather(uvtask_list, uv_container)

        return uvdata_out


def run_uvsim(params, return_uv=False):
    """
    Run a simulation off of an obsparam yaml file.

    Args:
        params (string): Path to a parameter yaml file.
        return_uv (bool): If true, do not write results to file (default False) and return uv_out

    Returns:
        uv_out (UVData): Finished simulation results. (if return_uv is True)
    """

    mpi.start_mpi()
    rank = mpi.get_rank()
    comm = mpi.get_comm()

    input_uv = UVData()
    beam_list = None
    beam_dict = None
    catalog = None
    source_list_name = None

    if rank == 0:
        input_uv, beam_list, beam_dict = simsetup.initialize_uvdata_from_params(params)
        catalog, source_list_name = simsetup.initialize_catalog_from_params(params, input_uv)

    input_uv = comm.bcast(input_uv, root=0)
    beam_list = comm.bcast(beam_list, root=0)
    beam_dict = comm.bcast(beam_dict, root=0)
    catalog = comm.bcast(catalog, root=0)

    uv_out = run_uvdata_uvsim(input_uv, beam_list, beam_dict=beam_dict, catalog=catalog, source_list_name=source_list_name)

    if rank == 0:
        with open(params, 'r') as pfile:
            param_dict = yaml.safe_load(pfile)
        simutils.write_uvdata(uv_out, param_dict, dryrun=return_uv)
    if return_uv:
        return uv_out
    comm.Barrier()
