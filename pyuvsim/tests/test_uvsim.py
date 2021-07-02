# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
import copy
import itertools
import os
import warnings

import astropy.constants as const
import pyuvdata.utils as uvutils
from astropy import units
from astropy.coordinates import Angle, SkyCoord, EarthLocation, Longitude, Latitude
import pytest
from pyuvdata import UVData, UVBeam
import pyradiosky

import pyuvsim
import pyuvsim.utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.astropy_interface import Time

EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
EW_uvfits_10time10chan = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')
herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')


def multi_beams():
    beam0 = UVBeam()
    beam0.read_beamfits(herabeam_default)
    beam0.freq_interp_kind = 'cubic'
    beam0.interpolation_function = 'az_za_simple'
    beam1 = pyuvsim.AnalyticBeam('uniform')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.6)
    beams = [beam0, beam1, beam2, beam3]

    try:
        beam4 = beam0.copy()
        beam4.to_healpix(nside=8)
        beam4.interpolation_function = 'healpix_simple'
        beams.append(beam4)
    except (ImportError, ModuleNotFoundError):
        pass

    return beams


multi_beams = multi_beams()


@pytest.fixture(scope='module')
def triangle_pos():
    hera_uv = UVData()
    hera_uv.read_uvfits(longbl_uvfits_file,
                        ant_str='cross')  # consists of a right triangle of baselines with w term
    hera_uv.unphase_to_drift(use_ant_pos=True)

    enu = hera_uv.get_ENU_antpos()[0]
    uvw = hera_uv.uvw_array[:hera_uv.Nbls]

    return enu, uvw


@pytest.fixture
def uvobj_beams_srcs():
    # A uvdata object, beam list, beam dict, and source array.
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_hex37_14.6m.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    param_dict['select'] = {'redundant_threshold': 0.1}
    uv_obj, beam_list, beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)

    # Add more beams to the list.
    # Don't use the uniform beam (need to see coherency change with positions).
    ref_freq, alpha = 100e6, -0.5
    beam_list[0] = pyuvsim.AnalyticBeam('airy', diameter=13.0)
    beam_list.append(pyuvsim.AnalyticBeam('airy', diameter=14.6))
    beam_list.append(pyuvsim.AnalyticBeam('gaussian', sigma=0.03))
    beam_list.append(pyuvsim.AnalyticBeam('gaussian', sigma=0.03,
                     ref_freq=ref_freq, spectral_index=alpha))

    # Assign the last few antennas to use these other beams.
    beam_dict['ANT36'] = 1
    beam_dict['ANT35'] = 2
    beam_dict['ANT34'] = 3

    time = Time(uv_obj.time_array[0], format='jd', scale='utc')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, arrangement='long-line', Nsrcs=30, return_data=True
    )

    sources.polarized = np.array([15, 16, 17])  # 15 -- 17th sources have polarization.
    sources.stokes_Q = np.array([0.1] * 3)[None, :]
    sources.stokes_U = np.array([2.0] * 3)[None, :]
    sources.stokes_V = np.array([0.0] * 3)[None, :]

    return uv_obj, beam_list, beam_dict, sources


def test_visibility_single_zenith_source(cst_beam, hera_loc):
    """Test single zenith source."""

    beam0 = cst_beam.copy()
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=np.radians(10.0))
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.0)

    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = hera_loc

    time.location = array_location

    freq = (150e6 * units.Hz)
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    source.update_positions(time, array_location)

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    for beam in [beam0, beam1, beam2, beam3]:
        beam_list = pyuvsim.BeamList([beam])
        array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

        task = pyuvsim.UVTask(source, time, freq, baseline, array)

        engine = pyuvsim.UVEngine(task)

        visibility = engine.make_visibility()
        assert np.allclose(visibility, np.array([.5, .5, 0, 0]), atol=5e-3)


def test_visibility_source_below_horizon(cst_beam, hera_loc):
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = hera_loc
    time.location = array_location

    freq = (150e6 * units.Hz)

    src_alt = Angle('-40d')

    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = cst_beam.copy()

    beam_list = pyuvsim.BeamList([beam])
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source_arr, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    assert np.allclose(visibility, np.array([0, 0, 0, 0]))


def test_visibility_source_below_horizon_radec(cst_beam, hera_loc):
    # redo with RA/Dec defined source
    time = Time(2458098.27471265, format='jd')

    array_location = hera_loc

    time.location = array_location
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(ra=Angle('13h20m'), dec=Angle('-30d43m17.5s'),
                            obstime=time, frame='icrs', location=array_location)

    source = pyradiosky.SkyModel('src_down', source_coord.ra, source_coord.dec,
                                 np.array([1.0, 0, 0, 0]).reshape(4, 1) * units.Jy, 'flat')

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = cst_beam

    beam_list = pyuvsim.BeamList([beam])
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    assert np.allclose(visibility, np.array([0, 0, 0, 0]))


def test_redundant_baselines(cst_beam, hera_loc):
    """Check that two perfectly redundant baselines are truly redundant. """

    time = Time(2458098.27471265, format='jd')
    array_location = hera_loc
    time.location = array_location
    freq = (150e6 * units.Hz)
    src_alt = Angle('85.0d')

    # Set up antenna positions in ENU:
    antpos = np.array([[0, 0, 0], [28, 0, 0]], dtype=float)

    en_shift = [5., 5., 0]
    antenna1 = pyuvsim.Antenna('ant1', 1, antpos[0, :], 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, antpos[1, :], 0)
    antenna3 = pyuvsim.Antenna('ant3', 3, antpos[0, :] + en_shift, 0)
    antenna4 = pyuvsim.Antenna('ant4', 4, antpos[1, :] + en_shift, 0)

    # make a source off zenith
    time.location = array_location
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    beam_list = pyuvsim.BeamList([cst_beam])

    baseline1 = pyuvsim.Baseline(antenna1, antenna2)
    baseline2 = pyuvsim.Baseline(antenna3, antenna4)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task1 = pyuvsim.UVTask(source, time, freq, baseline1, array)
    engine = pyuvsim.UVEngine(task1)

    visibility1 = engine.make_visibility()

    task2 = pyuvsim.UVTask(source, time, freq, baseline2, array)
    engine = pyuvsim.UVEngine(task2)

    visibility2 = engine.make_visibility()

    assert np.allclose(visibility1, visibility2)


@pytest.mark.parametrize('beam', multi_beams)
def test_single_offzenith_source(beam, hera_loc):
    """Test single off-zenith source."""

    time = Time(2458098.27471265, format='jd')
    array_location = hera_loc
    time.location = array_location
    freq = (123e6 * units.Hz)

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    antpos = np.array([[0, 0, 0], [28, 0, 0]], dtype=float)

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    # create_mock_catalog uses azimuth of 90
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    source.update_positions(time, array_location)
    src_alt_az = source.alt_az
    assert np.isclose(src_alt_az[0], src_alt.rad)
    assert np.isclose(src_alt_az[1], src_az.rad)

    src_lmn = source.pos_lmn
    assert np.isclose(src_lmn[0], src_l)
    assert np.isclose(src_lmn[1], src_m)
    assert np.isclose(src_lmn[2], src_n)

    beam_list = pyuvsim.BeamList([beam])

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    # analytically calculate visibility
    beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_alt.rad, src_az.rad)
    beam_za2, beam_az2 = simutils.altaz_to_zenithangle_azimuth(src_alt_az[0], src_alt_az[1])

    assert np.isclose(beam_za, beam_za2)
    assert np.isclose(beam_az, beam_az2)

    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array([beam_az]), za_array=np.array([beam_za]),
        freq_array=np.array([freq.to_value('Hz')]))

    jones = np.zeros((2, 2, 1), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    beam_jones = antenna1.get_beam_jones(array, src_alt_az, freq)
    assert np.allclose(beam_jones, jones)

    uvw_array = units.Quantity([[28., 0., 0.]], unit='m')
    uvw_wavelength_array = uvw_array / const.c * freq.to('1/s')
    # Remove source axis from jones matrix
    jones = jones.squeeze()
    vis_analytic = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(
        2j * np.pi * (
            uvw_wavelength_array[0, 0] * src_l
            + uvw_wavelength_array[0, 1] * src_m
            + uvw_wavelength_array[0, 2] * src_n
        )
    )
    vis_analytic = np.array(
        [vis_analytic[0, 0], vis_analytic[1, 1], vis_analytic[0, 1], vis_analytic[1, 0]]
    )

    assert np.allclose(baseline.uvw.to_value('m'), uvw_array.value)
    assert np.allclose(visibility, vis_analytic)


@pytest.mark.parametrize('beam', multi_beams)
def test_offzenith_source_multibl(beam, hera_loc, triangle_pos):
    """Calculate visibilities for a baseline triangle of an off-zenith source."""

    enu_antpos, uvw_array = triangle_pos
    time = Time(2458098.27471265, format='jd')
    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    array_location = hera_loc
    freq = (123e6 * units.Hz)

    antenna1 = pyuvsim.Antenna('ant1', 0, np.array(enu_antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 1, np.array(enu_antpos[1, :]), 0)
    antenna3 = pyuvsim.Antenna('ant3', 2, np.array(enu_antpos[2, :]), 0)

    # make a source off zenith
    time.location = array_location
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    beam_list = pyuvsim.BeamList([beam])

    baselines = [pyuvsim.Baseline(antenna2, antenna1),
                 pyuvsim.Baseline(antenna3, antenna1),
                 pyuvsim.Baseline(antenna3, antenna2)]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    tasks = [pyuvsim.UVTask(source, time, freq, bl, array) for bl in baselines]
    visibilities = []
    uvws = []
    engine = pyuvsim.UVEngine()
    for t in tasks:
        engine.set_task(t)
        visibilities.append(engine.make_visibility())
        uvws.append(t.baseline.uvw)
    uvws = np.array(uvws)

    # analytically calculate visibilities
    beam.peak_normalize()
    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array([0.0]), za_array=np.array([src_za.rad]),
        freq_array=np.array([freq.to_value('Hz')])
    )
    jones = np.zeros((2, 2, 1), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, :]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, :]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, :]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, :]

    uvw_wavelength_array = uvw_array * units.m / const.c * freq.to('1/s')

    visibilities_analytic = []
    for u, v, w in uvw_wavelength_array:
        vis = 0.5 * np.dot(jones[..., 0], np.conj(jones[..., 0]).T) * np.exp(
            2j * np.pi * (u * src_l + v * src_m + w * src_n))
        visibilities_analytic.append(np.array([vis[0, 0], vis[1, 1], vis[0, 1], vis[1, 0]]))

    assert np.allclose(uvws, uvw_array)
    assert np.allclose(visibilities, visibilities_analytic)


@pytest.mark.filterwarnings("ignore:The uvw_array does not match the expected values")
def test_file_to_tasks(cst_beam):
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', Nsrcs=5, return_data=True)

    beam_list = pyuvsim.BeamList([cst_beam])

    Nblts = hera_uv.Nblts
    Nfreqs = hera_uv.Nfreqs

    Ntasks = Nblts * Nfreqs
    beam_dict = None

    taskiter = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict
    )
    uvtask_list = list(taskiter)

    tlist = copy.deepcopy(uvtask_list)
    # Test task comparisons
    tlist.sort()
    assert np.all([tlist[i + 1] > tlist[i] for i in range(Ntasks - 1)])
    assert np.all([tlist[i + 1] >= tlist[i] for i in range(Ntasks - 1)])
    task0 = copy.deepcopy(tlist[0])
    task1 = copy.deepcopy(tlist[1])
    task0.baseline = task1.baseline
    task0.uvdata_index = (0, 0, 0)
    task1.uvdata_index = (1, 0, 0)
    assert task1 > task0
    assert task1 >= task0

    task0.uvdata_index = (0, 0, 0)
    task1.uvdata_index = (0, 0, 1)
    assert task1 > task0
    assert task1 >= task0
    assert task0 <= task1

    tel_loc = EarthLocation.from_geocentric(*hera_uv.telescope_location, unit='m')

    telescope = pyuvsim.Telescope(hera_uv.telescope_name, tel_loc, beam_list)

    ant_pos = hera_uv.antenna_positions + hera_uv.telescope_location
    ant_pos_enu = uvutils.ENU_from_ECEF(
        ant_pos, *hera_uv.telescope_location_lat_lon_alt
    )

    expected_task_list = []
    antenna_names = hera_uv.antenna_names
    antennas = []
    for num, antname in enumerate(antenna_names):
        beam_id = 0
        antennas.append(pyuvsim.Antenna(antname, num, ant_pos_enu[num], beam_id))

    antennas1 = []
    for antnum in hera_uv.ant_1_array:
        index = np.where(hera_uv.antenna_numbers == antnum)[0][0]
        antennas1.append(antennas[index])

    antennas2 = []
    for antnum in hera_uv.ant_2_array:
        index = np.where(hera_uv.antenna_numbers == antnum)[0][0]
        antennas2.append(antennas[index])

    sources = sources.get_skymodel()
    # set the frequency array to match what happens in uvdata_to_task_iter
    sources.freq_array = hera_uv.freq_array[0] * units.Hz

    for idx, antenna1 in enumerate(antennas1):
        antenna2 = antennas2[idx]
        baseline = pyuvsim.Baseline(antenna1, antenna2)
        task = pyuvsim.UVTask(sources, time.jd, hera_uv.freq_array[0, 0], baseline, telescope)
        task.uvdata_index = (idx, 0, 0)
        expected_task_list.append(task)

    expected_task_list.sort()
    for idx, task in enumerate(uvtask_list):
        exp_task = expected_task_list[idx]
        assert task == exp_task


@pytest.mark.filterwarnings("ignore:The uvw_array does not match the expected values")
def test_gather():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', return_data=True)

    beam_list = pyuvsim.BeamList([multi_beams[1]])

    Nblts = hera_uv.Nblts
    Nfreqs = hera_uv.Nfreqs

    Ntasks = Nblts * Nfreqs
    beam_dict = dict(zip(hera_uv.antenna_names, [0] * hera_uv.Nants_data))
    taskiter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, sources, beam_list,
                                           beam_dict)
    uvtask_list = list(taskiter)

    uv_out = pyuvsim.simsetup._complete_uvdata(hera_uv, inplace=False)

    size_complex = np.ones(1, dtype=complex).nbytes

    visbuf = bytearray(uv_out.data_array)
    for task in uvtask_list:
        engine = pyuvsim.UVEngine(task)
        vis = engine.make_visibility()

        blti, spw, freq_ind = task.uvdata_index

        flat_ind = np.ravel_multi_index(
            (blti, spw, freq_ind, 0), uv_out.data_array.shape
        )
        offset = flat_ind * size_complex
        val = np.frombuffer(visbuf[offset:offset + vis.nbytes], dtype=np.complex128)
        val += vis
        visbuf[offset:offset + vis.nbytes] = bytearray(val)[:]
    uv_out.data_array = np.frombuffer(visbuf, dtype=complex).reshape(uv_out.data_array.shape)

    assert np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3)


def test_local_task_gen():
    # Confirm I get the same results looping over the task list as I do with the generator function.
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_10time10chan)
    hera_uv.select(times=np.unique(hera_uv.time_array)[0:3], freq_chans=range(3))
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, arrangement='random', Nsrcs=5, return_data=True
    )

    beam_list = pyuvsim.BeamList([multi_beams[1]])

    Nblts = hera_uv.Nblts
    Nfreqs = hera_uv.Nfreqs
    Ntasks = Nblts * Nfreqs
    beam_dict = None

    # Check error conditions
    uv_iter0 = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), 'not_uvdata', sources, beam_list, beam_dict
    )
    with pytest.raises(TypeError, match='input_uv must be UVData object'):
        next(uv_iter0)
    uv_iter1 = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv,
                                           'not_skydata', beam_list, beam_dict)
    with pytest.raises(TypeError, match='catalog must be a SkyModelData'):
        next(uv_iter1)

    # Copy sources and beams so we don't accidentally reuse quantities.
    taskiter = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict
    )
    uvtask_list = list(taskiter)

    uvtask_iter = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), hera_uv, copy.deepcopy(sources),
        copy.deepcopy(beam_list), beam_dict
    )

    engine0 = pyuvsim.UVEngine(reuse_spline=False)

    for tki, task0 in enumerate(uvtask_iter):
        task1 = uvtask_list[tki]
        engine1 = pyuvsim.UVEngine(task1, reuse_spline=True)
        engine0.set_task(task0)
        assert np.allclose(engine1.make_visibility(), engine0.make_visibility())


def test_task_coverage():
    """
    Check that the task ids generated in different scenarios
    cover all combinations of baseline/time/frequency and source.
    """

    Npus_list = [1, 5, 13, 19]
    for Npus in Npus_list:
        # Case 1 -- (Npus < Nbltf)

        print(Npus)
        Nbls = 4
        Ntimes = 3
        Nfreqs = 1
        Nsrcs = 10

        Nbltf = Nbls * Ntimes * Nfreqs

        # List of pairs -- (bl/t/f index, source index)
        srci, bltfi = map(np.ndarray.flatten, np.meshgrid(np.arange(Nsrcs), np.arange(Nbltf)))

        tasks_expected = np.column_stack((bltfi, srci))
        tasks_all = []
        for rank in range(Npus):
            task_inds, src_inds, Ntasks_local, Nsrcs_local = pyuvsim.uvsim._make_task_inds(
                Nbls, Ntimes, Nfreqs, Nsrcs, rank, Npus
            )
            src_inds = np.arange(Nsrcs)[src_inds]   # Turn slice into iterator
            tasks = itertools.product(task_inds, src_inds)
            tasks_all.append(tasks)
        tasks_all = itertools.chain(*tasks_all)
        tasks = np.array(list(tasks_all))
        assert np.all(tasks == tasks_expected)

        # Case 2 -- (Nbltf < Npus and Nsrcs > Npus)

        if Npus == 1:
            continue  # case 2 won't work for 1 pu

        Nbls = 2
        Ntimes = 1
        Nfreqs = 2
        Nsrcs = 100

        Nbltf = Nbls * Ntimes * Nfreqs

        bltfi, srci = map(np.ndarray.flatten, np.meshgrid(np.arange(Nbltf), np.arange(Nsrcs)))

        tasks_expected = np.column_stack((bltfi, srci))
        tasks_all = []
        for rank in range(Npus):
            task_inds, src_inds, Ntasks_local, Nsrcs_local = pyuvsim.uvsim._make_task_inds(
                Nbls, Ntimes, Nfreqs, Nsrcs, rank, Npus
            )

            tasks = itertools.product(task_inds, src_inds)
            tasks_all.append(tasks)
        tasks_all = itertools.chain(*tasks_all)
        tasks = np.array(list(tasks_all))

        # Returned task indices are out of order, compared with the meshgrid.
        inds = np.lexsort((tasks[:, 0], tasks[:, 1]), axis=0)
        assert np.all(tasks[inds] == tasks_expected)


def test_source_splitting():
    # Check that if the available memory is less than the expected size of the source catalog,
    # then the task iterator will loop over chunks of the source array.
    pytest.importorskip('mpi4py')
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_10time10chan)
    hera_uv.select(times=np.unique(hera_uv.time_array)[0:3], freq_chans=range(3))
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, arrangement='random', Nsrcs=30, return_data=True
    )

    # Spoof environmental parameters.
    # Choose an absurdly large number of tasks per node and very small available memory
    # to trigger the source splitting condition in uvdata_to_task_iter
    # The alternative would be to make a very large source catalog, but that's not ideal in a test.
    Npus_node = 2000
    pyuvsim.mpi.Npus_node = Npus_node  # Absurdly large
    os.environ['SLURM_MEM_PER_NODE'] = str(400.0)  # Only 4MB of memory

    beam = pyuvsim.analyticbeam.AnalyticBeam('uniform')
    beam_list = pyuvsim.BeamList([beam])

    Nblts = hera_uv.Nblts
    Nfreqs = hera_uv.Nfreqs
    Nsrcs = sources.Ncomponents
    Ntasks = Nblts * Nfreqs
    beam_dict = None

    skymodel = sources.get_skymodel()
    skymodel_mem_footprint = (
        simutils.estimate_skymodel_memory_usage(
            skymodel.Ncomponents, skymodel.Nfreqs)
        * Npus_node
    )
    mem_avail = pyuvsim.utils.get_avail_memory()

    Nsky_parts = np.ceil(skymodel_mem_footprint / float(mem_avail))
    partsize = int(np.floor(Nsrcs / Nsky_parts))

    taskiter = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict, Nsky_parts=Nsky_parts
    )

    uvtask_list = list(taskiter)

    assert pyuvsim.estimate_skymodel_memory_usage(partsize, 1) * Npus_node < mem_avail

    # Normally, the number of tasks is Nbls * Ntimes * Nfreqs (partsize = Nsrcs)
    # If the source list is split within the task iterator, then it will be larger.
    assert len(uvtask_list) == Ntasks * Nsky_parts

    # Reset spoofed parameters.
    del os.environ['SLURM_MEM_PER_NODE']
    pyuvsim.mpi.Npus_node = 1


def test_quantity_reuse(uvobj_beams_srcs):
    # Check that the right quantities on the UVEngine are changed when
    # the time/frequency/antenna pair change.

    uv_obj, beam_list, beam_dict, sources = uvobj_beams_srcs

    beam_list.set_obj_mode()

    Ntasks = uv_obj.Nblts * uv_obj.Nfreqs

    taskiter = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), uv_obj, sources, beam_list, beam_dict
    )
    uvtask_list = list(taskiter)
    engine = pyuvsim.UVEngine()

    def allclose_or_none(first, second):
        if first is None or second is None:
            return False
        elif not first.shape == second.shape:
            return False
        return np.allclose(first, second)

    for ti, task in enumerate(uvtask_list):
        sky = task.sources
        prev_freq = engine.current_freq
        prev_time = engine.current_time
        prev_beam_pair = engine.current_beam_pair

        prev_local_coherency = copy.deepcopy(engine.local_coherency)
        prev_apparent_coherency = copy.deepcopy(engine.apparent_coherency)
        prev_jones1 = copy.deepcopy(engine.beam1_jones)
        prev_jones2 = copy.deepcopy(engine.beam2_jones)
        prev_source_pos = copy.deepcopy(sky.alt_az)

        engine.set_task(task)
        engine.make_visibility()

        apcoh_changed = not allclose_or_none(engine.apparent_coherency, prev_apparent_coherency)
        jones_changed = (not allclose_or_none(engine.beam1_jones, prev_jones1)
                         or not allclose_or_none(engine.beam2_jones, prev_jones2))
        locoh_changed = not allclose_or_none(engine.local_coherency, prev_local_coherency)
        srcpos_changed = not allclose_or_none(sky.alt_az, prev_source_pos)

        freq = task.freq.to_value("Hz")
        time = task.time.jd
        beampair = (task.baseline.antenna1.beam_id, task.baseline.antenna2.beam_id)

        # check that when time, freq, or beam pair changes, the relevant quantities also change.
        if (freq != prev_freq) or (time != prev_time) or (beampair != prev_beam_pair):
            assert apcoh_changed and jones_changed
        if time != prev_time:
            # Note -- local_coherency will only change if the sources are polarized.
            assert srcpos_changed and locoh_changed


def test_update_flags(uvobj_beams_srcs):
    # Ensure that the right update flags are set when certain
    # task attributes change.

    uv_obj, beam_list, beam_dict, sources = uvobj_beams_srcs

    Nsky_parts = 6
    Ntasks = uv_obj.Nblts * uv_obj.Nfreqs

    # Simulate the chunk of tasks covered on one rank.
    rank = 5
    Npus = 20
    task_inds, _ = pyuvsim.utils.iter_array_split(rank, Ntasks, Npus)

    taskiter = pyuvsim.uvdata_to_task_iter(
        task_inds, uv_obj, sources, beam_list, beam_dict, Nsky_parts=Nsky_parts
    )

    time, freq, beam_pair, src_chunk = [None] * 4
    engine = pyuvsim.UVEngine()

    axes_covered = [False] * 4
    for task in taskiter:
        engine.set_task(task)
        baseline = task.baseline
        task_beam_pair = (baseline.antenna1.beam_id, baseline.antenna2.beam_id)
        time_changed = task.time.jd != time
        freq_changed = task.freq != freq
        antpair_changed = task_beam_pair != beam_pair
        src_changed = src_chunk is not task.sources

        if time_changed or src_changed:
            axes_covered[0] = axes_covered[0] or time_changed
            axes_covered[1] = axes_covered[1] or src_changed

            assert (engine.update_positions
                    and engine.update_local_coherency
                    and engine.update_beams)

        if antpair_changed or freq_changed:
            axes_covered[2] = axes_covered[2] or antpair_changed
            axes_covered[3] = axes_covered[3] or freq_changed
            assert engine.update_beams

        time = task.time.jd
        freq = task.freq
        beam_pair = task_beam_pair
        src_chunk = task.sources

    # This is a test of the test itself, making sure that
    # each value changed at some point in the task iterator.
    assert all(axes_covered)


def test_overflow_check():
    # Ensure error before running sim for too many tasks.

    # This is a large number of tasks that has been run successfully.
    should_pass = 2004992
    pyuvsim.uvsim._check_ntasks_valid(should_pass)

    # This number of tasks is known to produce an overflow error
    # on bcast/gather operations in MPI.
    # This test just makes sure that the right error is raised
    # by the checker function for this value.
    should_fail = 22118400
    with pytest.raises(ValueError, match="Too many tasks"):
        pyuvsim.uvsim._check_ntasks_valid(should_fail)


def test_fullfreq_check(uvobj_beams_srcs):
    # Check that the task iter will error if 'spectral_type' is 'full'
    # and the frequencies on the catalog do not match the simulation's.

    uv_obj, beam_list, beam_dict, sources = uvobj_beams_srcs
    beam_list.set_obj_mode()

    Nsrcs = 30

    freqs0 = np.linspace(100, 130, uv_obj.Nfreqs) * 1e6 * units.Hz
    freqs1 = uv_obj.freq_array[0, :] * units.Hz

    stokes = np.zeros((4, uv_obj.Nfreqs, Nsrcs)) * units.Jy
    stokes[0, :, :] = 1.0 * units.Jy

    ra = Longitude(np.linspace(0, 2 * np.pi, Nsrcs), 'rad')
    dec = Latitude(np.linspace(-np.pi / 2, np.pi / 3, Nsrcs), 'rad')

    sky0 = pyradiosky.SkyModel(
        name=np.arange(Nsrcs).astype(str),
        ra=ra,
        dec=dec,
        stokes=stokes,
        spectral_type='full',
        freq_array=freqs0
    )

    sky1 = pyradiosky.SkyModel(
        name=np.arange(Nsrcs).astype(str),
        ra=ra,
        dec=dec,
        stokes=stokes,
        spectral_type='full',
        freq_array=freqs1
    )

    Ntasks = uv_obj.Nblts * uv_obj.Nfreqs

    sky0 = pyuvsim.simsetup.SkyModelData(sky0)
    sky1 = pyuvsim.simsetup.SkyModelData(sky1)

    taskiter0 = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), uv_obj, sky0, beam_list, beam_dict
    )
    taskiter1 = pyuvsim.uvdata_to_task_iter(
        np.arange(Ntasks), uv_obj, sky1, beam_list, beam_dict
    )

    with pytest.raises(ValueError, match="Some requested frequencies are not present"):
        next(taskiter0)

    next(taskiter1)


@pytest.mark.skipif('pyuvsim.astropy_interface.hasmoon')
def test_moonloc_error(uvobj_beams_srcs):
    # Break if the uvobj indicates that the sim is on the Moon, but lunarsky is not available.

    uv_obj, beam_list, beam_dict, sources = uvobj_beams_srcs

    uv_obj.extra_keywords['world'] = 'moon'

    with pytest.raises(ValueError, match='Need lunarsky module to simulate'):
        next(pyuvsim.uvsim.uvdata_to_task_iter(range(5), uv_obj, sources, beam_list, beam_dict))

    uv_obj.extra_keywords['world'] = 'gallifrey'

    with pytest.raises(ValueError, match="If world keyword is set, it must "):
        next(pyuvsim.uvsim.uvdata_to_task_iter(range(5), uv_obj, sources, beam_list, beam_dict))


def test_run_mpierr():
    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    )
    if pyuvsim.mpi is None:
        with pytest.raises(ImportError, match='You need mpi4py to use the uvsim module'):
            pyuvsim.run_uvsim(params, return_uv=True)

        with pytest.raises(ImportError, match='You need mpi4py to use the uvsim module'):
            pyuvsim.run_uvdata_uvsim(UVData(), ['beamlist'])
