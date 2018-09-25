# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import numpy as np
import nose.tools as nt
import copy
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy import units
import astropy.constants as const

from pyuvdata import UVBeam, UVData
import pyuvdata.utils as uvutils
from pyuvdata.data import DATA_PATH
import pyuvdata.tests as uvtest

import pyuvsim
import pyuvsim.utils as simutils
import pyuvsim.tests as simtest
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
EW_uvfits_10time10chan = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')


def test_visibility_single_zenith_source():
    """Test single zenith source."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location

    freq = (150e6 * units.Hz)
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    source = source_arr[0]

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = simtest.make_cst_beams()

    beam_list = [beam]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([.5, .5, 0, 0]), atol=5e-3))


def test_visibility_source_below_horizon():
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location

    freq = (150e6 * units.Hz)

    src_alt = Angle('-40d')

    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)
    source = source_arr[0]

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = simtest.make_cst_beams()

    beam_list = [beam]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([0, 0, 0, 0])))


def test_visibility_source_below_horizon_radec():
    # redo with RA/Dec defined source
    time = Time(2458098.27471265, format='jd')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(ra=Angle('13h20m'), dec=Angle('-30d43m17.5s'),
                            obstime=time, frame='icrs', location=array_location)

    source = pyuvsim.Source('src_down', source_coord.ra, source_coord.dec, freq,
                            [1.0, 0, 0, 0])

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = simtest.make_cst_beams()

    beam_list = [beam]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([0, 0, 0, 0])))


def test_visibility_single_zenith_source_uvdata():
    """Test single zenith source using test uvdata file."""
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
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
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    source = source_arr[0]

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([.5, .5, 0, 0]), atol=5e-3))


def test_redundant_baselines():
    """Check that two perfectly redundant baselines are truly redundant. """

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file, ant_str='cross')
    hera_uv.unphase_to_drift()

    src_alt = Angle('85.0d')

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
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
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)
    source = source_arr[0]

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

    baseline1 = pyuvsim.Baseline(antenna1, antenna2)
    baseline2 = pyuvsim.Baseline(antenna3, antenna4)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task1 = pyuvsim.UVTask(source, time, freq, baseline1, array)
    engine = pyuvsim.UVEngine(task1)

    visibility1 = engine.make_visibility()

    task2 = pyuvsim.UVTask(source, time, freq, baseline2, array)
    engine = pyuvsim.UVEngine(task2)

    visibility2 = engine.make_visibility()

    nt.assert_true(np.allclose(visibility1, visibility2))


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
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos, _ = hera_uv.get_ENU_antpos()

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    # create_mock_catalog uses azimuth of 90
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)
    source = source_arr[0]

    src_alt_az = source.alt_az_calc(time, array_location)
    nt.assert_true(np.isclose(src_alt_az[0], src_alt.rad))
    nt.assert_true(np.isclose(src_alt_az[1], src_az.rad))

    src_lmn = source.pos_lmn(time, array_location)
    nt.assert_true(np.isclose(src_lmn[0], src_l))
    nt.assert_true(np.isclose(src_lmn[1], src_m))
    nt.assert_true(np.isclose(src_lmn[2], src_n))

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

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

    nt.assert_true(np.isclose(beam_za, beam_za2))
    nt.assert_true(np.isclose(beam_az, beam_az2))

    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array([beam_az]),
                                                         za_array=np.array([beam_za]),
                                                         freq_array=np.array([freq.to('Hz').value]))
    jones = np.zeros((2, 2), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    beam_jones = antenna1.get_beam_jones(array, src_alt_az, freq)
    print(beam_jones)
    print(jones)
    print(beam_jones - jones)
    nt.assert_true(np.allclose(beam_jones, jones))

    uvw_wavelength_array = hera_uv.uvw_array * units.m / const.c * freq.to('1/s')

    vis_analytic = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(2j * np.pi * (uvw_wavelength_array[0, 0] * src_l + uvw_wavelength_array[0, 1] * src_m + uvw_wavelength_array[0, 2] * src_n))
    vis_analytic = np.array([vis_analytic[0, 0], vis_analytic[1, 1], vis_analytic[0, 1], vis_analytic[1, 0]])

    nt.assert_true(np.allclose(baseline.uvw.to('m').value, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4))

    print(vis_analytic)
    print(visibility)
    print(vis_analytic - visibility)
    nt.assert_true(np.allclose(visibility, vis_analytic, atol=1e-4))


def test_offzenith_source_multibl_uvfits():
    """Test single off-zenith source using test uvdata file.
        Calculate visibilities for a baseline triangle.
    """
    hera_uv = UVData()
    hera_uv.read_uvfits(longbl_uvfits_file, ant_str='cross')   # consists of a right triangle of baselines with w term
    hera_uv.unphase_to_drift(use_ant_pos=True)

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos, ants = uvtest.checkWarnings(hera_uv.get_ENU_antpos,
                                        category=[DeprecationWarning],
                                        message=['The default for the `center` keyword has changed',
                                                 'The xyz array in ENU_from_ECEF is being interpreted as (Npts, 3)'],
                                        nwarnings=2)
    antenna1 = pyuvsim.Antenna('ant1', ants[0], np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', ants[1], np.array(antpos[1, :]), 0)
    antenna3 = pyuvsim.Antenna('ant3', ants[2], np.array(antpos[2, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    # create_mock_catalog uses azimuth of 90
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', alt=src_alt.deg)
    source = source_arr[0]

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
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array([src_az.rad]), za_array=np.array([src_za.rad]), freq_array=np.array([freq.to('Hz').value]))
    jones = np.zeros((2, 2), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    uvw_wavelength_array = hera_uv.uvw_array[0:hera_uv.Nbls] * units.m / const.c * freq.to('1/s')

    visibilities_analytic = []
    for u, v, w in uvw_wavelength_array:
        vis = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(2j * np.pi * (u * src_l + v * src_m + w * src_n))
        visibilities_analytic.append(np.array([vis[0, 0], vis[1, 1], vis[1, 0], vis[0, 1]]))

    # the file used different phasing code than the test uses -- increase the tolerance
    nt.assert_true(np.allclose(uvws, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4))

    # the file used different phasing code than the test uses -- increase the tolerance
    nt.assert_true(np.allclose(visibilities, visibilities_analytic, atol=1e-4))


def test_file_to_tasks():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    beam = simtest.make_cst_beams()
    beam_list = [beam]

    Nblts = hera_uv.Nblts
    Nbls = hera_uv.Nbls
    Ntimes = hera_uv.Ntimes
    Nfreqs = hera_uv.Nfreqs
    Nsrcs = len(sources)

    Ntasks = Nblts * Nfreqs * Nsrcs
    beam_dict = None

    taskiter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict)
    uvtask_list = uvtest.checkWarnings(list, [taskiter],
                                       message=['The default for the `center` keyword has changed'],
                                       category=DeprecationWarning)

    tlist = copy.deepcopy(uvtask_list)

    # Test task comparisons
    tlist.sort()
    nt.assert_true(np.all([tlist[i + 1] > tlist[i] for i in range(Ntasks - 1)]))
    nt.assert_true(np.all([tlist[i + 1] >= tlist[i] for i in range(Ntasks - 1)]))
    task0 = copy.deepcopy(tlist[0])
    task1 = copy.deepcopy(tlist[1])
    uvind0 = copy.copy(task0.uvdata_index)
    uvind1 = copy.copy(task1.uvdata_index)
    task0.baseline = task1.baseline
    task0.uvdata_index = (0, 0, 0)
    task1.uvdata_index = (1, 0, 0)
    nt.assert_true(task1 > task0)
    nt.assert_true(task1 >= task0)

    task0.uvdata_index = (0, 0, 0)
    task1.uvdata_index = (0, 0, 1)
    nt.assert_true(task1 > task0)
    nt.assert_true(task1 >= task0)
    nt.assert_true(task0 <= task1)

    tel_loc = EarthLocation.from_geocentric(*hera_uv.telescope_location, unit='m')

    telescope = pyuvsim.Telescope(hera_uv.telescope_name, tel_loc, beam_list)

    ant_pos = hera_uv.antenna_positions + hera_uv.telescope_location
    ant_pos_enu = uvutils.ENU_from_ECEF(ant_pos,
                                        *hera_uv.telescope_location_lat_lon_alt)

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

    for idx, antenna1 in enumerate(antennas1):
        antenna2 = antennas2[idx]
        baseline = pyuvsim.Baseline(antenna1, antenna2)
        task = pyuvsim.UVTask(sources[0], time.jd, hera_uv.freq_array[0, 0], baseline, telescope)
        task.uvdata_index = (idx, 0, 0)
        expected_task_list.append(task)
    expected_task_list.sort()
    for idx, task in enumerate(uvtask_list):
        exp_task = expected_task_list[idx]
        nt.assert_equal(task, exp_task)


def test_uvdata_init():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_10time10chan)
    hera_uv.unphase_to_drift(use_ant_pos=True)
    uvdata_out = pyuvsim.init_uvdata_out(hera_uv, 'zenith_source',
                                         obs_param_file='', telescope_config_file='', antenna_location_file='')

    hera_uv.data_array = np.zeros_like(hera_uv.data_array, dtype=np.complex)
    hera_uv.flag_array = np.zeros_like(hera_uv.data_array, dtype=bool)
    hera_uv.nsample_array = np.ones_like(hera_uv.data_array)
    obs_param_file = ''
    telescope_config_file = ''
    antenna_location_file = ''
    hera_uv.history = (pyuvsim.get_version_string()
                       + 'Sources from source list: zenith_source. '
                       + ' Based on config files: ' + obs_param_file + ', '
                       + telescope_config_file + ', ' + antenna_location_file + ' Npus = 1.'
                       + hera_uv.pyuvdata_version_str)
    hera_uv.instrument = hera_uv.telescope_name
    nt.assert_true(np.allclose(hera_uv.antenna_positions, uvdata_out.antenna_positions))
    nt.assert_true(uvdata_out.__eq__(hera_uv, check_extra=False))


def test_uvdata_init_errors():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 1.0)
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str')
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', obs_param_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', telescope_config_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', antenna_location_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', obs_param_file='string')
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', obs_param_file='string', telescope_config_file='string')
    nt.assert_raises(ValueError, pyuvsim.init_uvdata_out, hera_uv, 'source_list_str', obs_param_file='string', telescope_config_file=1.0)


def test_gather():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

    Nblts = hera_uv.Nblts
    Nbls = hera_uv.Nbls
    Ntimes = hera_uv.Ntimes
    Nfreqs = hera_uv.Nfreqs
    Nsrcs = len(sources)

    Ntasks = Nblts * Nfreqs * Nsrcs
    beam_dict = dict(zip(hera_uv.antenna_names, [0] * hera_uv.Nants_data))

    taskiter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict)
    uvtask_list = uvtest.checkWarnings(list, [taskiter],
                                       message=['The default for the `center` keyword has changed'],
                                       category=DeprecationWarning)

    uv_out = pyuvsim.init_uvdata_out(hera_uv, 'zenith_source',
                                     obs_param_file='', telescope_config_file='', antenna_location_file='')
    for task in uvtask_list:
        engine = pyuvsim.UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uv_out = pyuvsim.serial_gather(uvtask_list, uv_out)

    nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3))


def test_local_task_gen():
    # Confirm I get the same results looping over the task list as I do with the generator function.
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_10time10chan)
    hera_uv.select(times=np.unique(hera_uv.time_array)[0:3], freq_chans=range(3))
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, kwds = pyuvsim.create_mock_catalog(time, arrangement='random', Nsrcs=5)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

    Nblts = hera_uv.Nblts
    Nbls = hera_uv.Nbls
    Ntimes = hera_uv.Ntimes
    Nfreqs = hera_uv.Nfreqs
    Nsrcs = len(sources)
    Ntasks = Nblts * Nfreqs * Nsrcs
    beam_dict = None

    # Copy sources and beams so we don't accidentally reuse quantities.
    taskiter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, sources, beam_list, beam_dict)
    uvtask_list = uvtest.checkWarnings(list, [taskiter],
                                       message=['The default for the `center` keyword has changed'],
                                       category=DeprecationWarning)
    uvtask_iter = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, copy.deepcopy(sources), copy.deepcopy(beam_list), beam_dict)

    # Check error conditions
    uv_iter0 = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), 'not_uvdata', sources, beam_list, beam_dict)
    nt.assert_raises(TypeError, next, uv_iter0)
    uv_iter1 = pyuvsim.uvdata_to_task_iter(np.arange(Ntasks), hera_uv, 'not_ndarray', beam_list, beam_dict)
    nt.assert_raises(TypeError, next, uv_iter1)

    engine0 = pyuvsim.UVEngine(reuse_spline=False)
    for tki, task0 in enumerate(uvtask_iter):
        task1 = uvtask_list[tki]
        engine1 = pyuvsim.UVEngine(task1, reuse_spline=True)
        engine0.set_task(task0)
        nt.assert_true(np.allclose(engine1.make_visibility(), engine0.make_visibility()))
