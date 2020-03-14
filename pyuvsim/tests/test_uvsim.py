# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import copy
import itertools
import os

import astropy.constants as const
import numpy as np
import pyuvdata.utils as uvutils
from astropy import units
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy.time import Time
import pytest
from pyuvdata import UVData
from pyuvdata.data import DATA_PATH
import pyradiosky

import pyuvsim
import pyuvsim.tests as simtest
import pyuvsim.utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
EW_uvfits_10time10chan = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')


def test_visibility_single_zenith_source():
    """Test single zenith source."""

    beam0 = simtest.make_cst_beams()
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=np.radians(10.0))
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.0)

    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
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


def test_visibility_source_below_horizon():
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location

    freq = (150e6 * units.Hz)

    src_alt = Angle('-40d')

    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = simtest.make_cst_beams()

    beam_list = pyuvsim.BeamList([beam])
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source_arr, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    assert np.allclose(visibility, np.array([0, 0, 0, 0]))


def test_visibility_source_below_horizon_radec():
    # redo with RA/Dec defined source
    time = Time(2458098.27471265, format='jd')

    array_location = EarthLocation(
        lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.
    )
    time.location = array_location
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(ra=Angle('13h20m'), dec=Angle('-30d43m17.5s'),
                            obstime=time, frame='icrs', location=array_location)

    source = pyradiosky.SkyModel('src_down', source_coord.ra, source_coord.dec,
                                 np.array([1.0, 0, 0, 0]).reshape(4, 1), [1e8], 'flat')

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = simtest.make_cst_beams()

    beam_list = pyuvsim.BeamList([beam])
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    assert np.allclose(visibility, np.array([0, 0, 0, 0]))


def test_visibility_single_zenith_source_uvdata():
    """Test single zenith source using test uvdata file."""
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(
        *hera_uv.telescope_location, unit='m'
    )
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos = hera_uv.antenna_positions[0:2, :] + hera_uv.telescope_location
    lat, lon, alt = hera_uv.telescope_location_lat_lon_alt
    antpos = uvutils.ENU_from_ECEF(antpos, lat, lon, alt)

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source at zenith
    time.location = array_location
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = pyuvsim.BeamList([beam])

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    assert np.allclose(visibility, np.array([.5, .5, 0, 0]), atol=5e-3)


def test_redundant_baselines():
    """Check that two perfectly redundant baselines are truly redundant. """

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file, ant_str='cross')
    hera_uv.unphase_to_drift()

    src_alt = Angle('85.0d')

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(
        *hera_uv.telescope_location, unit='m'
    )
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos, _ = hera_uv.get_ENU_antpos()

    en_shift = [5., 5., 0]
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)
    antenna3 = pyuvsim.Antenna('ant3', 3, np.array(antpos[0, :]) + en_shift, 0)
    antenna4 = pyuvsim.Antenna('ant4', 4, np.array(antpos[1, :]) + en_shift, 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = pyuvsim.BeamList([beam])

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


def test_single_offzenith_source_uvfits():
    """Test single off-zenith source using test uvdata file."""
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file, ant_str='cross')
    hera_uv.unphase_to_drift()

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(
        *hera_uv.telescope_location, unit='m'
    )
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos, _ = hera_uv.get_ENU_antpos()

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

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = pyuvsim.BeamList([beam])

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    # analytically calculate visibility
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'
    beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_alt.rad, src_az.rad)
    beam_za2, beam_az2 = simutils.altaz_to_zenithangle_azimuth(src_alt_az[0], src_alt_az[1])

    assert np.isclose(beam_za, beam_za2)
    assert np.isclose(beam_az, beam_az2)

    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array([beam_az]), za_array=np.array([beam_za]),
        freq_array=np.array([freq.to('Hz').value]))

    jones = np.zeros((2, 2, 1), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    beam_jones = antenna1.get_beam_jones(array, src_alt_az, freq)
    assert np.allclose(beam_jones, jones)

    uvw_wavelength_array = hera_uv.uvw_array * units.m / const.c * freq.to('1/s')
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

    assert np.allclose(baseline.uvw.to('m').value, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4)
    assert np.allclose(visibility, vis_analytic, atol=1e-4)


def test_offzenith_source_multibl_uvfits():
    """Test single off-zenith source using test uvdata file.
        Calculate visibilities for a baseline triangle.
    """
    hera_uv = UVData()
    hera_uv.read_uvfits(longbl_uvfits_file,
                        ant_str='cross')  # consists of a right triangle of baselines with w term
    hera_uv.unphase_to_drift(use_ant_pos=True)

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(
        *hera_uv.telescope_location, unit='m'
    )
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos, ants = hera_uv.get_ENU_antpos()
    antenna1 = pyuvsim.Antenna('ant1', ants[0], np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', ants[1], np.array(antpos[1, :]), 0)
    antenna3 = pyuvsim.Antenna('ant3', ants[2], np.array(antpos[2, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    # create_mock_catalog uses azimuth of 90
    source, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)

    beam = pyuvsim.AnalyticBeam('uniform')

    beam_list = [beam]

    baselines = [pyuvsim.Baseline(antenna2, antenna1),
                 pyuvsim.Baseline(antenna3, antenna1),
                 pyuvsim.Baseline(antenna3, antenna2)]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    tasks = [pyuvsim.UVTask(source, time, freq, bl, array) for bl in baselines]
    visibilities = []
    uvws = []
    for t in tasks:
        engine = pyuvsim.UVEngine(t)
        visibilities.append(engine.make_visibility())
        uvws.append(t.baseline.uvw)

    uvws = np.array(uvws)
    # analytically calculate visibilities

    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'
    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array([src_az.rad]), za_array=np.array([src_za.rad]),
        freq_array=np.array([freq.to('Hz').value])
    )
    jones = np.zeros((2, 2), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    uvw_wavelength_array = hera_uv.uvw_array[0:hera_uv.Nbls] * units.m / const.c * freq.to('1/s')

    visibilities_analytic = []
    for u, v, w in uvw_wavelength_array:
        vis = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(
            2j * np.pi * (u * src_l + v * src_m + w * src_n))
        visibilities_analytic.append(np.array([vis[0, 0], vis[1, 1], vis[1, 0], vis[0, 1]]))

    # the file used different phasing code than the test uses -- increase the tolerance
    assert np.allclose(uvws, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4)

    # the file used different phasing code than the test uses -- increase the tolerance
    assert np.allclose(visibilities, visibilities_analytic, atol=1e-4)


def test_file_to_tasks():
    pytest.importorskip('mpi4py')
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', Nsrcs=5, return_table=True)

    beam = simtest.make_cst_beams()
    beam_list = pyuvsim.BeamList([beam])

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

    sources = pyradiosky.array_to_skymodel(sources)
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


def test_gather():
    pytest.importorskip('mpi4py')
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', return_table=True)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = pyuvsim.BeamList([beam])

    Nblts = hera_uv.Nblts
    Nfreqs = hera_uv.Nfreqs

    Ntasks = Nblts * Nfreqs
    beam_dict = dict(zip(hera_uv.antenna_names, [0] * hera_uv.Nants_data))
    taskiter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, sources, beam_list,
                                           beam_dict)
    uvtask_list = list(taskiter)

    uv_out = pyuvsim.simsetup._complete_uvdata(hera_uv, inplace=False)
    for task in uvtask_list:
        engine = pyuvsim.UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uv_out = pyuvsim.serial_gather(uvtask_list, uv_out)

    assert np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3)


def test_local_task_gen():
    # Confirm I get the same results looping over the task list as I do with the generator function.
    pytest.importorskip('mpi4py')
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_10time10chan)
    hera_uv.select(times=np.unique(hera_uv.time_array)[0:3], freq_chans=range(3))
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, arrangement='random', Nsrcs=5, return_table=True
    )

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = pyuvsim.BeamList([beam])

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
                                           'not_ndarray', beam_list, beam_dict)
    with pytest.raises(TypeError, match='catalog must be a record array'):
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


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    pytest.importorskip('mpi4py')

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    hera_uv.select(polarizations=['xx'])

    with pytest.raises(ValueError, match='input_uv must have XX,YY,XY,YX polarization'):
        pyuvsim.run_uvdata_uvsim(hera_uv, ['beamlist'])


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
            src_inds = range(Nsrcs)[src_inds]   # Turn slice into iterator
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
        time, arrangement='random', Nsrcs=30, return_table=True
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
    Nsrcs = len(sources)
    Ntasks = Nblts * Nfreqs
    beam_dict = None

    skymodel = pyradiosky.array_to_skymodel(sources)
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


def test_get_beam_jones():
    # check setting the interpolation method

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam.freq_interp_kind = None
    beam.interpolation_function = 'az_za_simple'
    beam_list = pyuvsim.BeamList([beam])
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 10, 0]), 0)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    source_altaz = np.array([[0.0], [np.pi / 4.]])
    freq = 100e6 * units.Hz

    with pytest.raises(ValueError, match='freq_interp_kind must be set'):
        antenna1.get_beam_jones(array, source_altaz, freq)

    jones = antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind='cubic')
    assert beam.freq_interp_kind == 'cubic'
    jones0 = antenna1.get_beam_jones(array, source_altaz, freq)
    jones1 = antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind='linear')
    assert beam.freq_interp_kind == 'linear'
    jones2 = antenna1.get_beam_jones(array, source_altaz, freq)

    assert (np.all(jones2 == jones0)
            and np.all(jones1 == jones)
            and np.all(jones1 == jones0))
