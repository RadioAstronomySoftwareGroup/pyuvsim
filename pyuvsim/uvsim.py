# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""
The primary machinery of pyuvsim.

The :class:`~UVTask` and :class:`~UVEngine` classes and the functions that actually run
the simulation.
"""

import numpy as np
import yaml
import warnings
import astropy.units as units
from astropy.coordinates import EarthLocation
from astropy.units import Quantity
from astropy.constants import c as speed_of_light
from pyuvdata import UVData

from . import mpi
from . import simsetup
from . import utils as simutils
from .antenna import Antenna
from .baseline import Baseline
from .telescope import Telescope
from .simsetup import SkyModelData

from .astropy_interface import MoonLocation, hasmoon, Time


__all__ = ['UVTask', 'UVEngine', 'uvdata_to_task_iter', 'run_uvsim', 'run_uvdata_uvsim']


class UVTask:
    """
    An object to hold all the information to calculate a visibility for a set of sources.

    Parameters
    ----------
    sources : :class:`pyuvsim.simsetup.SkyModelData`
        The sources to include in the visibility.
    time : :class:`astropy.time.Time` object or float
        Time at which to calculate the visibility, either an astropy Time object or a
        float in Julian Day (will be converted to a Time object using `format='jd'`).
    freq : :class:`astropy.units.Quantity` or float
        Frequency at which to calculate the visibility, either an astropy Quantity with
        units compatible with Hz or a float in Hz (will be converted to an astropy
        Quantity with units of Hz)
    baseline : :class:`pyuvsim.Baseline`
        The baseline to calculate the visibility for.
    telescope : :class:`pyuvsim.Telescope`
        The telescope object this visibility belongs to, which carries the telescope
        location and beam information.
    freq_i : int
        Frequency index for this visibility's source coherency.

    Attributes
    ----------
    sources : :class:`pyuvsim.simsetup.SkyModelData`
        The sources to include in the visibility.
    time : :class:`astropy.time.Time`
        Time at which to calculate the visibility.
    freq : :class:`astropy.units.Quantity`
        Frequency at which to calculate the visibility.
    baseline : :class:`pyuvsim.Baseline`
        The baseline to calculate the visibility for.
    telescope : :class:`pyuvsim.Telescope`
        The telescope object this visibility belongs to.
    freq_i : int
        Frequency index for this visibility's source coherency. Always `0` for flat
        spectrum sources, otherwise should be the index of `freq` in the full frequency
        array for this simulation.
    uvdata_index : tuple of int
        Indexes for this visibility location in the uvdata data_array, length 3 giving
        (blt_index, 0, frequency_index), where the zero index is for the old spw axis.
    visibility_vector : ndarray of complex
        The calculated visibility, shape (4,) ordered as [xx, yy, xy, yx].

    """

    def __init__(self, sources, time, freq, baseline, telescope, freq_i=0):
        self.time = time
        self.freq = freq
        self.sources = sources  # SkyModel object
        self.baseline = baseline
        self.telescope = telescope
        self.freq_i = freq_i
        self.visibility_vector = None
        self.uvdata_index = None  # Where to add the visibility in the uvdata object.

        if isinstance(self.time, float):
            self.time = Time(self.time, format='jd')
        if isinstance(self.freq, float):
            self.freq = self.freq * units.Hz
        if sources.spectral_type == 'flat':
            self.freq_i = 0

    def __eq__(self, other):
        """
        Test for equality between two UVTask objects.

        Parameters
        ----------
        other : :class:`UVTask`
            The other

        """
        return (np.isclose(self.time.jd, other.time.jd, atol=1e-4)
                and np.isclose(self.freq.value, other.freq.value, atol=1e-4)
                and (self.sources == other.sources)
                and (self.baseline == other.baseline)
                and (self.visibility_vector == other.visibility_vector)
                and (self.uvdata_index == other.uvdata_index)
                and (self.telescope == other.telescope))

    def __gt__(self, other):
        """
        Check if this UVTask has a larger `uvdata_index` than another one.

        Parameters
        ----------
        other : :class:`UVTask`
            The other

        """
        blti0, fi0 = self.uvdata_index
        blti1, fi1 = other.uvdata_index
        if self.baseline == other.baseline:
            if fi0 == fi1:
                return blti0 > blti1
            return fi0 > fi1
        return self.baseline > other.baseline

    def __ge__(self, other):
        """
        Check if this UVTask has a larger or equal `uvdata_index` than another one.

        Parameters
        ----------
        other : :class:`UVTask`
            The other

        """
        blti0, fi0 = self.uvdata_index
        blti1, fi1 = other.uvdata_index
        if self.baseline == other.baseline:
            if fi0 == fi1:
                return blti0 >= blti1
            return fi0 >= fi1
        return self.baseline >= other.baseline

    def __lt__(self, other):
        """
        Check if this UVTask has a smaller `uvdata_index` than another one.

        Parameters
        ----------
        other : :class:`UVTask`
            The other

        """
        return not self.__ge__(other)

    def __le__(self, other):
        """
        Check if this UVTask has a smaller or equal `uvdata_index` than another one.

        Parameters
        ----------
        other : :class:`UVTask`
            The other

        """
        return not self.__gt__(other)


class UVEngine:
    """
    The object where the visibility calculations actually take place.

    Parameters
    ----------
    task : :class:`~UVTask`
        The task to do the calculations for.
    reuse_spline : bool
        Option to reuse the spline in the beam interpolation to save time.
    update_positions : bool
        Flag indicating that positions need to be updated (when time or source changes).
    update beams : bool
        Flag indicating that beams need to be updated (when time, sources, frequency,
        or beam pair changes).

    Attributes
    ----------
    reuse_spline : bool
        Option to reuse the spline in the beam interpolation to save time.
    update_positions : bool
        Flag indicating that positions need to be updated (when time or source changes).
    update_local_coherency : bool
        Flag indicating that local coherencies need to be updated (when time or source changes).
    update beams : bool
        Flag indicating that beams need to be updated (when time, sources, frequency,
        or beam pair changes).
    sources : :class:`pyuvsim.simsetup.SkyModelData`
        The sources to include in the visibility.
    current_time : :class:`astropy.time.Time`
        Time for the current calculation.
    current_freq : :class:`astropy.units.Quantity`
        Frequency for the current calculation, units compatible with Hz.
    current_beam_pair : 2-tuple of ints
        Tuple containing the beam ids for the current calculation.
    beam1_jones : array_like of float
        Jones matricies for the first beam for each source location, shape
        (2,2, Ncomponents). The first axis is feed, the second axis is vector component
        on the sky in az/za.
    beam2_jones : array_like of float
        Jones matricies for the second beam for each source location, shape
        (2,2, Ncomponents). The first axis is feed, the second axis is vector component
        on the sky in az/za.
    local_coherency : array_like of float
        Source local coherencies in alt/az basis, shape (2, 2, Nfreqs, Ncomponents).
    apparent_coherency : array_like of float
        Source apparent coherencies (including the beam response), shape
        (2, 2, Nfreqs, Ncomponents).
    task :  :class:`~UVTask`
        The task currently being calculated.
    """

    def __init__(self, task=None, update_positions=True, update_beams=True, reuse_spline=True):
        self.reuse_spline = reuse_spline  # Reuse spline fits in beam interpolation
        self.update_positions = update_positions
        self.update_beams = update_beams

        self.sources = None
        self.current_time = None
        self.current_freq = None
        self.current_beam_pair = None
        self.beam1_jones = None
        self.beam2_jones = None
        self.local_coherency = None
        self.apparent_coherency = None

        if task is not None:
            self.set_task(task)

    def set_task(self, task):
        """
        Set the task for the Engine.

        Parameters
        ----------
        task : :class:`UVTask`
            The task to do the calculations for.

        """
        self.task = task

        # These flags control whether quantities can be re-used from the last task:
        #   Update local coherency for all frequencies when time or sources changes.
        #   Update source positions when time or sources changes.
        #   Update saved beam jones matrices when time, sources, frequency, or beam pair changes.

        baseline = self.task.baseline
        beam_pair = (baseline.antenna1.beam_id, baseline.antenna2.beam_id)

        if (not self.current_time == task.time.jd) or (self.sources is not task.sources):
            self.update_positions = True
            self.update_local_coherency = True
            self.update_beams = True
        else:
            self.update_beams = False
            self.update_positions = False
            self.update_local_coherency = False

        if not self.current_freq == task.freq.to('Hz').value:
            self.update_beams = True

        if not self.current_beam_pair == beam_pair:
            self.current_beam_pair = beam_pair
            self.update_beams = True

        self.current_time = task.time.jd
        self.current_freq = task.freq.to("Hz").value
        self.sources = task.sources

    def apply_beam(self, beam_interp_check=True):
        """
        Set apparent coherency from jones matrices and source coherency.

        beam_interp_check :  bool
            Option to enable checking that the source positions are within the area
            covered by the beam. If all the az/za beams cover the full sky horizon to
            horizon this checking is turned off by default in `run_uvdata_uvsim`,
            otherwise it is turned on.
            Setting this to False can speed up simulations but if sources are simulated
            outside the beam area the response will be incorrect.
            This keyword only applies to beams that are regularly gridded in azimuth and
            zenith angle.
        """
        beam1_id, beam2_id = self.current_beam_pair

        sources = self.task.sources
        baseline = self.task.baseline

        if sources.alt_az is None:
            sources.update_positions(self.task.time, self.task.telescope.location)

        if self.update_local_coherency:
            self.local_coherency = sources.coherency_calc()

        self.beam1_jones = baseline.antenna1.get_beam_jones(
            self.task.telescope,
            sources.alt_az[..., sources.above_horizon],
            self.task.freq,
            reuse_spline=self.reuse_spline,
            beam_interp_check=beam_interp_check,
        )

        if beam1_id == beam2_id:
            self.beam2_jones = np.copy(self.beam1_jones)
        else:
            self.beam2_jones = baseline.antenna2.get_beam_jones(
                self.task.telescope,
                sources.alt_az[..., sources.above_horizon],
                self.task.freq,
                reuse_spline=self.reuse_spline,
                beam_interp_check=beam_interp_check,
            )

        # coherency is a 2x2 matrix
        # [ |Ex|^2, Ex* Ey, Ey* Ex |Ey|^2 ]
        # where x and y vectors along the local alt/az axes.

        # Apparent coherency gives the direction and polarization dependent baseline response to
        # a source.

        coherency = self.local_coherency[:, :, self.task.freq_i, :]

        self.beam2_jones = np.swapaxes(self.beam2_jones, 0, 1).conj()  # Transpose at each component

        self.apparent_coherency = np.einsum(
            "abz,bcz,cdz->adz", self.beam1_jones, coherency, self.beam2_jones
        )

    def make_visibility(self):
        """Calculate visibility contribution from a set of source components."""
        assert (isinstance(self.task.freq, Quantity))
        srcs = self.task.sources
        time = self.task.time
        location = self.task.telescope.location

        if self.update_positions:
            srcs.update_positions(time, location)

        if self.update_beams:
            self.apply_beam()

        pos_lmn = srcs.pos_lmn[..., srcs.above_horizon]

        # need to convert uvws from meters to wavelengths
        uvw_wavelength = self.task.baseline.uvw / speed_of_light * self.task.freq.to('1/s')
        fringe = np.exp(2j * np.pi * np.dot(uvw_wavelength, pos_lmn))
        vij = self.apparent_coherency * fringe

        # Sum over source component axis:
        vij = np.sum(vij, axis=2)

        # Reshape to be [xx, yy, xy, yx]
        vis_vector = np.asarray([vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]])
        return vis_vector


def _make_task_inds(Nblts, Nfreqs, Nsrcs, rank, Npus):
    """
    Make iterators defining task and sources computed on rank.

    The splitting depends on the relative sizes of the axes and the number of processing
    units. There are three possible ways the splitting is done:

        - if Npus < Nbltf -- Split by Nbltf, split sources in the task loop for memory's sake.
        - if Nbltf < Npus and Nsrcs > Npus -- Split by Nsrcs only
        - if (Nsrcs, Nbltf) < Npus -- Split by Nbltf

    Parameters
    ----------
    Nblts : int
        Number of blts included in the sim.
    Nfreqs :  int
        Number of frequencies included in the sim.
    Nsrcs : int
        Number of sources included in the sim.
    rank : int
        mpi rank.
    Npus : int
        Number of processing units.

    """
    Nbltf = Nblts * Nfreqs

    split_srcs = False

    if (Nbltf < Npus) and (Npus < Nsrcs):
        split_srcs = True

    if split_srcs:
        src_inds, Nsrcs_local = simutils.iter_array_split(rank, Nsrcs, Npus)
        task_inds = range(Nbltf)
        Ntasks_local = Nbltf
    else:
        task_inds, Ntasks_local = simutils.iter_array_split(rank, Nbltf, Npus)
        src_inds = range(Nsrcs)
        Nsrcs_local = Nsrcs

    return task_inds, src_inds, Ntasks_local, Nsrcs_local


def uvdata_to_task_iter(task_ids, input_uv, catalog, beam_list, beam_dict, Nsky_parts=1):
    """
    Generate UVTask objects.

    Parameters
    ----------
    task_ids : range
        Task indices in the full flattened meshgrid of parameters.
    input_uv : :class:`pyuvdata.UVData`
        UVData object to be filled with data.
    catalog : :class:`pyuvsim.simsetup.SkyModelData`
        Source components.
    beam_list : :class:`pyuvsim.BeamList`
        BeamList carrying beam model (in object mode).
    beam_dict : dict
        Map of antenna numbers to index in beam_list.

    Yields
    ------
        Iterable of UVTask objects.

    """
    # The task_ids refer to tasks on the flattened meshgrid.
    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object.")

    #   Skymodel will now be passed in as a catalog array.
    if not isinstance(catalog, SkyModelData):
        raise TypeError("catalog must be a SkyModelData object.")

    # use future array shapes
    if not input_uv.future_array_shapes:
        input_uv.use_future_array_shapes()

    # Splitting the catalog for memory's sake.
    Nsrcs_total = catalog.Ncomponents
    if Nsky_parts > 1:
        Nsky_parts = int(Nsky_parts)
        src_iter = [simutils.iter_array_split(s, Nsrcs_total, Nsky_parts)[0]
                    for s in range(Nsky_parts)]
    else:
        src_iter = [range(Nsrcs_total)]
    # Build the antenna list.
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
    Nfreqs = input_uv.Nfreqs
    Nblts = input_uv.Nblts

    tasks_shape = (Nblts, Nfreqs)
    blt_ax, freq_ax = range(2)

    # we need to get slice or indices to index into the flat iterators
    # it does not accept range arguments as indices.
    if isinstance(task_ids, range):
        task_slice = slice(task_ids.start, task_ids.stop, task_ids.step)
    else:
        # There is a test where the input ids is an np.arange object
        # direct indexing should be fine.
        task_slice = task_ids

    # broadcast_to returns a view into the array
    # flat returns an iterator object.
    # no overhead from these calls.
    _times = np.broadcast_to(
        input_uv.time_array.reshape(-1, 1),
        (input_uv.Nblts, input_uv.Nfreqs),
    ).flat

    _bls = np.broadcast_to(
        input_uv.baseline_array.reshape(-1, 1),
        (input_uv.Nblts, input_uv.Nfreqs),
    ).flat

    _freqs = np.broadcast_to(
        input_uv.freq_array.flatten().reshape(1, -1),
        (input_uv.Nblts, input_uv.Nfreqs),
    ).flat

    # indexing with the slice should return a view,
    # but the iterator has to be collected into an array.
    # This should incure 64 * Ntasks_local bits per array
    # some additional overhead for numpy arrays as well
    # usually seeing [0.00001, 0.00003] MiB / task memory required
    # based on the reference simulations.
    order = np.lexsort((_bls[task_slice], _freqs[task_slice], _times[task_slice])).flat

    tloc = [np.float64(x) for x in input_uv.telescope_location]

    world = input_uv.extra_keywords.get('world', 'earth')

    if world.lower() == 'earth':
        location = EarthLocation.from_geocentric(*tloc, unit='m')
    elif world.lower() == 'moon':
        if not hasmoon:
            raise ValueError("Need lunarsky module to simulate an array on the Moon.")
        location = MoonLocation.from_selenocentric(*tloc, unit='m')
    else:
        raise ValueError("If world keyword is set, it must be either 'moon' or 'earth'.")
    telescope = Telescope(input_uv.telescope_name, location, beam_list)
    freq_array = input_uv.freq_array * units.Hz
    time_array = Time(input_uv.time_array, scale='utc', format='jd', location=telescope.location)
    for src_i in src_iter:
        sky = catalog.get_skymodel(src_i)
        if (
            sky.spectral_type == 'flat'
            and sky.freq_array is None
            and sky.reference_frequency is None
        ):
            sky.freq_array = freq_array
        if sky.component_type == 'healpix' and hasattr(sky, 'healpix_to_point'):
            sky.healpix_to_point()
        if sky.spectral_type != 'flat':
            sky.at_frequencies(freq_array)

        for index in order:
            task_index = task_ids[index]
            # Shape indicates slowest to fastest index.
            if not isinstance(task_index, tuple):
                task_index = np.unravel_index(task_index, tasks_shape)
            blt_i = task_index[blt_ax]
            freq_i = task_index[freq_ax]

            bl_num = input_uv.baseline_array[blt_i]

            # We reuse a lot of baseline info, so make the baseline list on the first go and reuse.
            if bl_num not in baselines.keys():
                antnum1 = input_uv.ant_1_array[blt_i]
                antnum2 = input_uv.ant_2_array[blt_i]
                index1 = np.where(input_uv.antenna_numbers == antnum1)[0][0]
                index2 = np.where(input_uv.antenna_numbers == antnum2)[0][0]
                baselines[bl_num] = Baseline(antennas[index1], antennas[index2])

            time = time_array[blt_i]
            bl = baselines[bl_num]
            freq = freq_array[freq_i]  # 0 = spw axis

            task = UVTask(sky, time, freq, bl, telescope, freq_i)
            task.uvdata_index = (blt_i, freq_i)    # 0 = spectral window index

            yield task
        del sky


def _check_ntasks_valid(Ntasks_tot):
    """
    Check that the size of the task array won't overflow the gather.

    Parameters
    ----------
    Ntasks_tot : int
        Number of total tasks.

    """
    # Found experimentally, this is the limit on the number of tasks
    # that can be gathered with MPI.
    MAX_NTASKS_GATHER = 10226018

    if Ntasks_tot >= MAX_NTASKS_GATHER:
        raise ValueError(
            f"Too many tasks for MPI to gather successfully: {Ntasks_tot}\n"
            "\t Consider splitting your simulation into multiple smaller jobs "
            "and joining them together with UVData.\n"
            "\t This error will go away when issue #289 is closed."
        )


def run_uvdata_uvsim(
    input_uv,
    beam_list,
    beam_dict=None,
    catalog=None,
    beam_interp_check=None,
    quiet=False,
    block_nonroot_stdout=True
):
    """
    Run uvsim from UVData object.

    Parameters
    ----------
    input_uv : :class:`pyuvdata.UVData` instance
        Provides baseline/time/frequency information.
    beam_list : :class:`pyuvsim.BeamList`
        BeamList carrying beam models.
    beam_dict : dictionary, optional
        {`antenna_name` : `beam_id`}, where `beam_id` is an index in the beam_list.
        This is used to assign beams to antennas. Default: All antennas get the 0th
        beam in the `beam_list`.
    catalog : :class:`pyuvsim.simsetup.SkyModelData`
        Immutable source parameters.
    beam_interp_check :  bool
        Option to enable checking that the source positions are within the area covered
        by the beam. If the beam covers the full sky horizon to horizon this checking
        is turned off by default, otherwise it is turned on.
        Setting this to False can speed up simulations but if sources are simulated
        outside the beam area the response will be incorrect.
        This keyword only applies to beams that are regularly gridded in azimuth and
        zenith angle.
    quiet : bool
        Do not print anything.
    block_nonroot_stdout : bool
        Redirect stdout on nonzero ranks to /dev/null, for cleaner output.

    Returns
    -------
    :class:`pyuvdata.UVData`
        instance containing simulated visibilities.

    """
    if mpi is None:
        raise ImportError("You need mpi4py to use the uvsim module. "
                          "Install it by running pip install pyuvsim[sim] "
                          "or pip install pyuvsim[all] if you also want the "
                          "line_profiler installed.")

    mpi.start_mpi(block_nonroot_stdout=block_nonroot_stdout)
    rank = mpi.get_rank()
    comm = mpi.get_comm()
    Npus = mpi.get_Npus()

    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    if not ((input_uv.Npols == 4) and (input_uv.polarization_array.tolist() == [-5, -6, -7, -8])):
        raise ValueError("input_uv must have XX,YY,XY,YX polarization")

    if not input_uv.future_array_shapes:
        input_uv.use_future_array_shapes()

    input_order = input_uv.blt_order
    if input_order != ("time", "baseline"):
        input_uv.reorder_blts(order="time", minor_order="baseline")

    # The root node will initialize our simulation
    # Read input file and make uvtask list
    if rank == 0 and not quiet:
        print('Nbls:', input_uv.Nbls, flush=True)
        print('Ntimes:', input_uv.Ntimes, flush=True)
        print('Nfreqs:', input_uv.Nfreqs, flush=True)
        print('Nsrcs:', catalog.Ncomponents, flush=True)
    if rank == 0:
        uv_container = simsetup._complete_uvdata(input_uv, inplace=False)
        if 'world' in input_uv.extra_keywords:
            uv_container.extra_keywords['world'] = input_uv.extra_keywords['world']
        vis_data = mpi.MPI.Win.Create(uv_container._data_array.value, comm=mpi.world_comm)
    else:
        vis_data = mpi.MPI.Win.Create(None, comm=mpi.world_comm)

    Nbls = input_uv.Nbls
    Nblts = input_uv.Nblts
    Nfreqs = input_uv.Nfreqs
    Nsrcs = catalog.Ncomponents

    task_inds, src_inds, Ntasks_local, Nsrcs_local = _make_task_inds(
        Nblts, Nfreqs, Nsrcs, rank, Npus
    )

    # Construct beam objects from strings
    beam_list.set_obj_mode(use_shared_mem=True)

    # In case the user created the beam list without checking consistency:
    beam_list.check_consistency()

    if beam_list.beam_type != 'efield':
        raise ValueError("Beam type must be efield!")

    if beam_interp_check is None:
        # check that all the beams cover the full sky
        if beam_list.check_all_azza_beams_full_sky():
            beam_interp_check = False
        else:
            beam_interp_check = True

    # Estimating required memory to decide how to split source array.
    mem_avail = (simutils.get_avail_memory()
                 - mpi.get_max_node_rss(return_per_node=True) * 2**30)

    Npus_node = mpi.node_comm.Get_size()
    skymodel_mem_footprint = (
        simutils.estimate_skymodel_memory_usage(Nsrcs, catalog.Nfreqs) * Npus_node
    )

    # Allow up to 50% of available memory for SkyModel data.
    skymodel_mem_max = 0.5 * mem_avail

    Nsky_parts = np.ceil(skymodel_mem_footprint / float(skymodel_mem_max))
    Nsky_parts = max(Nsky_parts, 1)
    if Nsky_parts > Nsrcs:
        raise ValueError("Insufficient memory for simulation.")

    local_task_iter = uvdata_to_task_iter(
        task_inds, input_uv, catalog.subselect(src_inds),
        beam_list, beam_dict, Nsky_parts=Nsky_parts
    )

    Ntasks_tot = Ntasks_local * Nsky_parts
    # Sum all the tasks across each node
    Nsky_parts = comm.reduce(Nsky_parts, op=mpi.MPI.MAX, root=0)
    Ntasks_tot = comm.reduce(Ntasks_tot, op=mpi.MPI.SUM, root=0)
    if rank == 0 and not quiet:
        if Nsky_parts > 1:
            print(
                "The source list has been split into Nsky_parts"
                f" <= {Nsky_parts} chunks on some or all nodes"
                " due to memory limitations.",
                flush=True,
            )
        print("Tasks: ", Ntasks_tot, flush=True)
        pbar = simutils.progsteps(maxval=Ntasks_tot)

    engine = UVEngine()
    count = mpi.Counter()
    size_complex = np.ones(1, dtype=complex).nbytes
    data_array_shape = (Nblts, Nfreqs, 4)
    uvdata_indices = []

    for task in local_task_iter:
        engine.set_task(task)
        vis = engine.make_visibility()

        blti, freq_ind = task.uvdata_index

        uvdata_indices.append(task.uvdata_index)

        flat_ind = np.ravel_multi_index(
            (blti, freq_ind, 0), data_array_shape
        )
        offset = flat_ind * size_complex

        vis_data.Lock(0)
        vis_data.Accumulate(vis, 0, target=offset, op=mpi.MPI.SUM)
        vis_data.Unlock(0)

        cval = count.next()
        if rank == 0 and not quiet:
            pbar.update(cval)

    request = comm.Ibarrier()

    while not request.Test():
        if rank == 0 and not quiet:
            cval = count.current_value()
            pbar.update(cval)

    count.free()
    if rank == 0 and not quiet:
        pbar.finish()

    if rank == 0 and not quiet:
        print("Calculations Complete.", flush=True)

    # If profiling is active, save meta data:
    from .profiling import prof     # noqa
    if hasattr(prof, 'meta_file'):  # pragma: nocover
        # Saving axis sizes on current rank (local) and for the whole job (global).
        # These lines are affected by issue 179 of line_profiler, so the nocover
        # above will need to stay until this issue is resolved (see profiling.py).
        task_inds = np.array(uvdata_indices)
        bl_inds = task_inds[:, 0] % Nbls
        time_inds = (task_inds[:, 0] - bl_inds) // Nbls
        Ntimes_loc = np.unique(time_inds).size
        Nbls_loc = np.unique(bl_inds).size
        Nfreqs_loc = np.unique(task_inds[:, 1]).size
        axes_dict = {
            'Ntimes_loc': Ntimes_loc,
            'Nbls_loc': Nbls_loc,
            'Nfreqs_loc': Nfreqs_loc,
            'Nsrcs_loc': Nsky_parts,
            'prof_rank': prof.rank
        }

        with open(prof.meta_file, 'w') as afile:
            for k, v in axes_dict.items():
                afile.write("{} \t {:d}\n".format(k, int(v)))

    vis_data.Free()

    if rank == 0:
        if input_order is not None:
            if len(input_order) < 2:
                input_order = (input_order[0], None)
            # Don't call the function if we're already expecting (time, baseline)
            if input_order != ("time", "baseline"):
                uv_container.reorder_blts(
                    order=input_order[0],
                    minor_order=input_order[1],
                )
        else:
            warnings.warn(
                "The parameter `blt_order` could not be identified for input_uv "
                " so the ordering cannot be restored."
                "The output UVData object will have (time, baseline) ordering."
            )

        # Updating file history.
        history = simutils.get_version_string()
        if catalog.filename is not None:
            history += ' Sources from source list(s): [' + ', '.join(catalog.filename) + '].'
        if 'obsparam' in input_uv.extra_keywords:
            obs_param_file = input_uv.extra_keywords['obsparam']
            telescope_config_file = input_uv.extra_keywords['telecfg']
            antenna_location_file = input_uv.extra_keywords['layout']
            history += (' Based on config files: ' + obs_param_file + ', '
                        + telescope_config_file + ', ' + antenna_location_file)
        elif hasattr(uv_container, "filename") and uv_container.filename is not None:
            history += ' Based on uvdata file(s): [' + ', '.join(uv_container.filename) + '].'
        history += ' Npus = ' + str(mpi.Npus) + '.'

        # add pyradiosky version
        history += catalog.pyradiosky_version_str

        # add pyuvdata version info
        history += uv_container.pyuvdata_version_str

        uv_container.history = history

        return uv_container


def run_uvsim(
    params, return_uv=False, beam_interp_check=None, quiet=False, block_nonroot_stdout=True,
):
    """
    Run a simulation off of an obsparam yaml file.

    Parameters
    ----------
    params : str
        Path to a parameter yaml file.
    return_uv : bool
        If true, do not write results to file and return uv_out. (Default False)
    beam_interp_check :  bool
        Option to enable checking that the source positions are within the area covered
        by the beam. If the beam covers the full sky horizon to horizon this checking
        is turned off by default, otherwise it is turned on.
        Setting this to False can speed up simulations but if sources are simulated
        outside the beam area the response will be incorrect.
        This keyword only applies to beams that are regularly gridded in azimuth and
        zenith angle.
    quiet : bool
        If True, do not print anything to stdout. (Default False)
    block_nonroot_stdout : bool
        Redirect stdout on nonzero ranks to /dev/null, for cleaner output.

    Returns
    -------
    uv_out : `pyuvdata.UVData`
        Finished simulation results.
        Returned only if return_uv is True.

    """
    if mpi is None:
        raise ImportError("You need mpi4py to use the uvsim module. "
                          "Install it by running pip install pyuvsim[sim] "
                          "or pip install pyuvsim[all] if you also want the "
                          "line_profiler installed.")

    mpi.start_mpi(block_nonroot_stdout=block_nonroot_stdout)
    rank = mpi.get_rank()
    comm = mpi.get_comm()

    input_uv = UVData()
    beam_list = None
    beam_dict = None
    skydata = SkyModelData()

    if rank == 0:
        start = Time.now()
        input_uv, beam_list, beam_dict = simsetup.initialize_uvdata_from_params(
            params, return_beams=True
        )
        skydata = simsetup.initialize_catalog_from_params(
            params, input_uv, return_recarray=False, return_catname=False
        )
        if not quiet:
            print(f"UVData initialization took {(Time.now() - start).to('minute'):.3f}")
        start = Time.now()
        skydata = simsetup.SkyModelData(skydata)
        if not quiet:
            print(f"Skymodel setup took {(Time.now() - start).to('minute'):.3f}")

    input_uv = comm.bcast(input_uv, root=0)
    beam_list = comm.bcast(beam_list, root=0)
    beam_dict = comm.bcast(beam_dict, root=0)
    skydata.share(root=0)

    start = Time.now()
    uv_out = run_uvdata_uvsim(
        input_uv,
        beam_list,
        beam_dict=beam_dict,
        catalog=skydata,
        quiet=quiet,
        beam_interp_check=beam_interp_check,
    )
    if rank == 0 and not quiet:
        print(f"Run uvdata uvsim took {(Time.now() - start).to('minute'):.3f}")

    if rank == 0:
        if isinstance(params, str):
            with open(params, 'r') as pfile:
                param_dict = yaml.safe_load(pfile)
        else:
            param_dict = params

        simutils.write_uvdata(uv_out, param_dict, dryrun=return_uv)

    if return_uv:
        return uv_out
