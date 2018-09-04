# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import os
import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy import units
from pyuvdata import UVBeam, UVData
import pyuvdata.utils as uvutils
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
import astropy.constants as const
from memory_profiler import profile
import yaml
import sys
from scipy.special import spherical_jn as jn
import pickle

cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
triangle_uvfits_file = os.path.join(SIM_DATA_PATH, '28m_triangle_10time_10chan.uvfits')
longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')
longbl_yaml_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.yaml')
laptop_size_sim = os.path.join(SIM_DATA_PATH, 'laptop_size_sim.yaml')


@profile
def create_zenith_source(time, name):
    """Create pyuvsim Source object at zenith.

    Input: Astropy Time object
        sample: Time('2018-03-01 00:00:00', scale='utc')
    Returns: Pyuvsim Source object
    """
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    del(icrs_coord)
    return pyuvsim.Source(name, ra, dec, freq, [1, 0, 0, 0])


def create_offzenith_source(time, name, az, alt):
    """Create pyuvsim Source object off zenith at az/alt.

    Inputs: Astropy Time object
        sample: Time('2018-03-01 00:00:00', scale='utc')
            az/alt as Angles
    Returns: Pyuvsim Source object
    """
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=alt, az=az,
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    return pyuvsim.Source(name, ra, dec, freq, [1.0, 0, 0, 0])


def test_source_zenith_icrs():
    """Test single source position at zenith constructed using icrs."""
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time = Time('2018-03-01 00:00:00', scale='utc', location=array_location)

    freq = (150e6 * units.Hz)

    lst = time.sidereal_time('apparent')

    tee_ra = lst
    cirs_ra = pyuvsim.tee_to_cirs_ra(tee_ra, time)

    cirs_source_coord = SkyCoord(ra=cirs_ra, dec=array_location.lat,
                                 obstime=time, frame='cirs', location=array_location)

    icrs_coord = cirs_source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    zenith_source = pyuvsim.Source('icrs_zen', ra, dec, freq, [1, 0, 0, 0])

    zenith_source_lmn = zenith_source.pos_lmn(time, array_location)

    nt.assert_true(np.allclose(zenith_source_lmn, np.array([0, 0, 1]), atol=1e-5))


def test_source_zenith():
    """Test single source position at zenith."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    zenith_source = create_zenith_source(time, 'zensrc')
    zenith_source_lmn = zenith_source.pos_lmn(time, array_location)
    print('Zenith Source lmn')
    print(zenith_source_lmn)

    nt.assert_true(np.allclose(zenith_source_lmn, np.array([0, 0, 1])))


@profile
def test_single_zenith_source():
    """Test single zenith source."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location

    freq = (150e6 * units.Hz)
    source = create_zenith_source(time, 'zensrc')

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([.5, .5, 0, 0]), atol=5e-3))


def test_source_below_horizon():
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location

    freq = (150e6 * units.Hz)

    source = create_offzenith_source(time, 'src_down', az=Angle('0d'), alt=Angle('-40d'))

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([0, 0, 0, 0])))

    # redo with RA/Dec defined source
    time = Time(2458098.27471265, format='jd')
    time.location = array_location

    source_coord = SkyCoord(ra=Angle('13h20m'), dec=Angle('-30d43m17.5s'),
                            obstime=time, frame='icrs', location=array_location)

    source = pyuvsim.Source('src_down', source_coord.ra, source_coord.dec, freq,
                            [1.0, 0, 0, 0])

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([0, 0, 0, 0])))


def test_single_zenith_source_uvdata():
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
    antpos = uvutils.ENU_from_ECEF(antpos.T, *hera_uv.telescope_location_lat_lon_alt).T

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source at zenith
    time.location = array_location
    source = create_zenith_source(time, 'zensrc')

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

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

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

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
    source = create_offzenith_source(time, 'offzensrc', az=src_az, alt=src_alt)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    baseline1 = pyuvsim.Baseline(antenna1, antenna2)
    baseline2 = pyuvsim.Baseline(antenna3, antenna4)

    print('Baseline1 uvw: ', baseline1.uvw)
    print('Baseline2 uvw: ', baseline2.uvw)
    print('Baseline1 enus: ', baseline1.antenna1.pos_enu, baseline1.antenna2.pos_enu)
    print('Baseline2 enus: ', baseline2.antenna1.pos_enu, baseline2.antenna2.pos_enu)

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
#    antpos = hera_uv.antenna_positions[0:2, :] + hera_uv.telescope_location
#    antpos = uvutils.ENU_from_ECEF(antpos.T, *hera_uv.telescope_location_lat_lon_alt).T
    antpos, _ = hera_uv.get_ENU_antpos()

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    source = create_offzenith_source(time, 'offzensrc', az=src_az, alt=src_alt)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    # analytically calculate visibility

    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array([src_az.rad]), za_array=np.array([src_za.rad]), freq_array=np.array([freq.to('Hz').value]))
    jones = np.zeros((2, 2), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    uvw_wavelength_array = hera_uv.uvw_array * units.m / const.c * freq.to('1/s')

    vis_analytic = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(2j * np.pi * (uvw_wavelength_array[0, 0] * src_l + uvw_wavelength_array[0, 1] * src_m + uvw_wavelength_array[0, 2] * src_n))
    vis_analytic = np.array([vis_analytic[0, 0], vis_analytic[1, 1], vis_analytic[0, 1], vis_analytic[1, 0]])

    print('Analytic visibility', vis_analytic)
    print('Calculated visibility', visibility)
    print baseline.uvw.to('m').value, hera_uv.uvw_array[0:hera_uv.Nbls]
    nt.assert_true(np.allclose(baseline.uvw.to('m').value, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4))
    nt.assert_true(np.allclose(visibility, vis_analytic, atol=1e-4))


def test_uniform_beam():
    beam = pyuvsim.AnalyticBeam('uniform')
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    time = Time('2018-03-01 00:00:00', scale='utc')
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    source_list, mock_keywords = pyuvsim.create_mock_catalog(time, 'hera_text', array_location=array_location)

    nsrcs = len(source_list)
    az_vals = []
    za_vals = []
    freq_vals = []
    for src in source_list:
        src_coord = SkyCoord(ra=src.ra, dec=src.dec, frame='icrs', obstime=time,
                             location=array_location)
        src_coord_altaz = src_coord.transform_to('altaz')
        az_vals.append(src_coord_altaz.az.to('rad').value)
        za_vals.append((Angle('90d') - src_coord_altaz.alt).to('rad').value)
        if len(freq_vals) > 0:
            if src.freq.to('Hz').value != freq_vals[0]:
                freq_vals.append(src.freq.to('Hz').value)
        else:
            freq_vals.append(src.freq.to('Hz').value)

    n_freqs = len(freq_vals)
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array(az_vals),
                                                         za_array=np.array(za_vals),
                                                         freq_array=np.array(freq_vals))
    expected_data = np.zeros((2, 1, 2, n_freqs, nsrcs), dtype=np.float)
    expected_data[1, 0, 0, :, :] = 1
    expected_data[0, 0, 1, :, :] = 1
    expected_data[1, 0, 1, :, :] = 1
    expected_data[0, 0, 0, :, :] = 1
    nt.assert_true(np.allclose(interpolated_beam, expected_data))


def test_airy_beam():
    diameter_m = 14.
    beam = pyuvsim.AnalyticBeam('airy', diameter=diameter_m)
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    time = Time('2018-03-01 00:00:00', scale='utc')
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    source_list, mock_keywords = pyuvsim.create_mock_catalog(time, 'hera_text', array_location=array_location)

    nsrcs = len(source_list)
    az_vals = []
    za_vals = []
    freq_vals = []
    for src in source_list:
        src_coord = SkyCoord(ra=src.ra, dec=src.dec, frame='icrs', obstime=time,
                             location=array_location)
        src_coord_altaz = src_coord.transform_to('altaz')
        az_vals.append(src_coord_altaz.az.to('rad').value)
        za_vals.append((Angle('90d') - src_coord_altaz.alt).to('rad').value)
        if len(freq_vals) > 0:
            if src.freq.to('Hz').value != freq_vals[0]:
                freq_vals.append(src.freq.to('Hz').value)
        else:
            freq_vals.append(src.freq.to('Hz').value)

    n_freqs = len(freq_vals)
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array(az_vals),
                                                         za_array=np.array(za_vals),
                                                         freq_array=np.array(freq_vals))

    expected_data = np.zeros((2, 1, 2, n_freqs, nsrcs), dtype=np.float)
    interp_zas = np.zeros((n_freqs, nsrcs), dtype=np.float)
    for f_ind in range(n_freqs):
        interp_zas[f_ind, :] = np.array(za_vals)
    za_grid, f_grid = np.meshgrid(interp_zas, freq_vals)
    xvals = diameter_m / 2. * np.sin(za_grid) * 2. * np.pi * f_grid / 3e8
    airy_vals = np.zeros_like(xvals)
    airy_vals[xvals > 0.] = 2. * jn(1, xvals[xvals > 0.]) / xvals[xvals > 0.]
    airy_vals[xvals == 0.] = 1.

    expected_data[1, 0, 0, :, :] = airy_vals
    expected_data[0, 0, 1, :, :] = airy_vals
    expected_data[1, 0, 1, :, :] = airy_vals
    expected_data[0, 0, 0, :, :] = airy_vals
    nt.assert_true(np.allclose(interpolated_beam, expected_data))


def test_gaussian_beam():
    sigma_rad = Angle('5d').to('rad').value
    beam = pyuvsim.AnalyticBeam('gaussian', sigma=sigma_rad)
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    time = Time('2018-03-01 00:00:00', scale='utc')
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    source_list, mock_keywords = pyuvsim.create_mock_catalog(time, 'hera_text', array_location=array_location)

    nsrcs = len(source_list)
    az_vals = []
    za_vals = []
    freq_vals = []
    for src in source_list:
        src_coord = SkyCoord(ra=src.ra, dec=src.dec, frame='icrs', obstime=time,
                             location=array_location)
        src_coord_altaz = src_coord.transform_to('altaz')
        az_vals.append(src_coord_altaz.az.to('rad').value)
        za_vals.append((Angle('90d') - src_coord_altaz.alt).to('rad').value)
        if len(freq_vals) > 0:
            if src.freq.to('Hz').value != freq_vals[0]:
                freq_vals.append(src.freq.to('Hz').value)
        else:
            freq_vals.append(src.freq.to('Hz').value)

    n_freqs = len(freq_vals)
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array(az_vals),
                                                         za_array=np.array(za_vals),
                                                         freq_array=np.array(freq_vals))

    expected_data = np.zeros((2, 1, 2, n_freqs, nsrcs), dtype=np.float)
    interp_zas = np.zeros((n_freqs, nsrcs), dtype=np.float)
    for f_ind in range(n_freqs):
        interp_zas[f_ind, :] = np.array(za_vals)
    gaussian_vals = np.exp(-(interp_zas**2) / (2 * sigma_rad**2))

    expected_data[1, 0, 0, :, :] = gaussian_vals
    expected_data[0, 0, 1, :, :] = gaussian_vals
    expected_data[1, 0, 1, :, :] = gaussian_vals
    expected_data[0, 0, 0, :, :] = gaussian_vals

    nt.assert_true(np.allclose(interpolated_beam, expected_data))


def test_offzenith_source_multibl_uvfits():
    """Test single off-zenith source using test uvdata file.
        Calculate visibilities for a baseline triangle.
    """
    hera_uv = UVData()
    hera_uv.read_uvfits(longbl_uvfits_file, ant_str='cross')   # consists of a right triangle of baselines with w term
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
    antpos, ants = hera_uv.get_ENU_antpos()
    antenna1 = pyuvsim.Antenna('ant1', ants[0], np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', ants[1], np.array(antpos[1, :]), 0)
    antenna3 = pyuvsim.Antenna('ant3', ants[2], np.array(antpos[2, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    source = create_offzenith_source(time, 'offzensrc', az=src_az, alt=src_alt)

    # beam = UVBeam()
    # beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
    #                    telescope_name='HERA',
    #                    feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
    #                    model_name='E-field pattern - Rigging height 4.9m',
    #                    model_version='1.0')

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

    print('pyuvsim uvws: ', np.around(uvws))
    print('file uvws: ', np.around(hera_uv.uvw_array[0:hera_uv.Nbls]))
    print('Difference: ', uvws - hera_uv.uvw_array[0:hera_uv.Nbls])

    # the file used different phasing code than the test uses -- increase the tolerance
    nt.assert_true(np.allclose(uvws, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4))

    print('Analytic visibility', visibilities_analytic)
    print('Calculated visibility', visibilities)
    print(np.array(visibilities) - np.array(visibilities_analytic))

    # the file used different phasing code than the test uses -- increase the tolerance
    nt.assert_true(np.allclose(visibilities, visibilities_analytic, atol=1e-4))


# This test is supposed to see if miriad works for conjugation, but right now there are too
# many unknowns with the HERA test file to understand why it doesn't pass.
@nt.nottest
def test_single_offzenith_source_miriad():
    """Test single off-zenith source using test uvdata file."""
    miriad_uv = UVData()
    miriad_uv.read_miriad(os.path.join(DATA_PATH, 'hera_testfile'), ant_str='9_10')
    miriad_uv.select(times=miriad_uv.time_array[0])

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    time = Time(miriad_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(miriad_uv.telescope_location[0],
                                                   miriad_uv.telescope_location[1],
                                                   miriad_uv.telescope_location[2],
                                                   unit='m')
    freq = miriad_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos = miriad_uv.antenna_positions[0:2, :] + miriad_uv.telescope_location
    antpos = uvutils.ENU_from_ECEF(antpos.T, *miriad_uv.telescope_location_lat_lon_alt).T

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source off zenith
    time.location = array_location
    source = create_offzenith_source(time, 'offzensrc', az=src_az, alt=src_alt)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    # analytically calculate visibility

    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'
    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array([src_az.rad]), za_array=np.array([src_za.rad]), freq_array=np.array([freq.to('Hz').value]))
    jones = np.zeros((2, 2), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    uvw_wavelength_array = hera_uv.uvw_array * units.m / const.c * freq.to('1/s')

    vis_analytic = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(-2j * np.pi * (uvw_wavelength_array[0, 0] * src_l + uvw_wavelength_array[0, 1] * src_m + uvw_wavelength_array[0, 2] * src_n))
    vis_analytic = np.array([vis_analytic[0, 0], vis_analytic[1, 1], vis_analytic[1, 0], vis_analytic[0, 1]])

    print('Analytic visibility', vis_analytic)
    print('Calculated visibility', visibility)

    nt.assert_true(np.allclose(visibility, vis_analytic, atol=5e-3))


@profile
def test_yaml_to_tasks():
    #    params = yaml.safe_load(open(longbl_yaml_file))
    params = yaml.safe_load(open(laptop_size_sim))
    params['config_path'] = SIM_DATA_PATH
    output = pyuvsim.simsetup.initialize_uvdata_from_params(params)
    input_uv, beam_list, beam_dict, beam_ids = \
        pyuvsim.simsetup.initialize_uvdata_from_params(params)
    time = Time('2018-03-01 00:00:00', scale='utc')
    HERA_location = EarthLocation(lat='-30d43m17.5s',
                                  lon='21d25m41.9s',
                                  height=1073.)
    catalog, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='zenith',
                                                         Nsrcs=10, array_location=HERA_location)
    print("Size of catalog:", sys.getsizeof(catalog), " bytes with ",
          len(catalog), "entries")
    uvtask_list = pyuvsim.uvdata_to_task_list(input_uv, catalog, beam_list)
    uvtask_pickle = pickle.dumps(uvtask_list)
    print("getsizeof(uvtask_pickle) = ", sys.getsizeof(uvtask_pickle), " bytes with ",
          len(uvtask_list), "entries")
    print("len(uvtask_pickle) = ", len(uvtask_pickle))


@profile
def test_file_to_tasks():

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time, 'zensrc')])

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources, beam_list)

    tel_loc = EarthLocation.from_geocentric(*hera_uv.telescope_location, unit='m')
    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]
    telescope = pyuvsim.Telescope(hera_uv.telescope_name, tel_loc, beam_list)

    ant_pos = hera_uv.antenna_positions + hera_uv.telescope_location
    ant_pos_enu = uvutils.ENU_from_ECEF(ant_pos.T,
                                        *hera_uv.telescope_location_lat_lon_alt).T

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

    for idx, task in enumerate(uvtask_list):
        exp_task = expected_task_list[idx]
        nt.assert_equal(task, exp_task)


def test_uvdata_init():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    hera_uv.unphase_to_drift(use_ant_pos=True)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time, 'zensrc')])

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources, beam_list)
    # for task in uvtask_list:
    #    task.time = Time(task.time, format='jd')
    #    task.freq = task.freq * units.Hz

    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 1.0)
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str', uvdata_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str', uvdata_file='testfile', telescope_config_file='tconfig')
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str')
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str', obs_param_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str', telescope_config_file=1.0)
    nt.assert_raises(ValueError, pyuvsim.initialize_uvdata, uvtask_list, 'source_list_str', antenna_location_file=1.0)

    uvdata_out = pyuvsim.initialize_uvdata(uvtask_list, 'zenith_source',
                                           uvdata_file=EW_uvfits_file)

    hera_uv.data_array = np.zeros_like(hera_uv.data_array, dtype=np.complex)
    hera_uv.flag_array = np.zeros_like(hera_uv.data_array, dtype=bool)
    hera_uv.nsample_array = np.ones_like(hera_uv.data_array)
    hera_uv.history = (pyuvsim.get_version_string()
                       + 'Sources from source list: zenith_source. '
                       'Based on UVData file: ' + EW_uvfits_file + '. Npus = 1.'
                       + hera_uv.pyuvdata_version_str)
    hera_uv.instrument = hera_uv.telescope_name
    nt.assert_equal(hera_uv._antenna_positions, uvdata_out._antenna_positions)
    nt.assert_true(uvdata_out.__eq__(hera_uv, check_extra=False))


def test_gather():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time, 'src1')
                        # create_zenith_source(time, 'src2'),
                        # create_zenith_source(time, 'src3')
                        ])

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]

    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources, beam_list)
    uv_out = pyuvsim.initialize_uvdata(uvtask_list, 'zenith_source',
                                       uvdata_file=EW_uvfits_file)

    for task in uvtask_list:
        engine = pyuvsim.UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uv_out = pyuvsim.serial_gather(uvtask_list, uv_out)

    nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3))


def test_sources_equal():
    time = Time('2018-03-01 00:00:00', scale='utc')
    src1 = create_zenith_source(time, 'src')
    src2 = create_zenith_source(time, 'src')
    nt.assert_equal(src1, src2)


def test_mock_catalog():

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file, ant_str='cross')

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')

    test_source = create_offzenith_source(time, 'src0', az=src_az, alt=src_alt)
    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', zen_ang=src_za.deg)
    cat_source = cat[0]
    for k in cat_source.__dict__:
        print 'Cat: ', k, cat_source.__dict__[k]
        print 'Test: ', k, test_source.__dict__[k]
        print '\n'

    nt.assert_equal(cat_source, test_source)


def test_gaussbeam_values():
    """
        Make the long-line point sources up to 10 degrees from zenith.
        Obtain visibilities
        Confirm that the values match the expected beam values at those zenith angles.
    """
    sigma = 0.05
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
    freq = hera_uv.freq_array[0, 0] * units.Hz

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')

    catalog, mock_keywords = pyuvsim.create_mock_catalog(time=time, arrangement='long-line', Nsrcs=41,
                                                         max_za=10., array_location=array_location)

    beam = pyuvsim.AnalyticBeam('gaussian', sigma=sigma)
    array = pyuvsim.Telescope('telescope_name', array_location, [beam])

    # Need a dummy baseline for this test.
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    coherencies = []
    zenith_angles = []
    for src in catalog:
        task = pyuvsim.UVTask(src, time, freq, baseline, array)
        engine = pyuvsim.UVEngine(task)
        engine.apply_beam()
#        task.source.az_za_calc(time, array_location)
        zenith_angles.append(task.source.az_za[1])   # In radians.
        coherencies.append(np.real(engine.apparent_coherency[0, 0]).astype(float))  # All four components should be identical

    coherencies = np.array(coherencies)
    zenith_angles = np.array(zenith_angles)

    # Confirm the coherency values (ie., brightnesses) match the beam values.

    beam_values = np.exp(-(zenith_angles)**2 / (2 * beam.sigma**2))
    nt.assert_true(np.all(beam_values**2 == coherencies))


def test_tee_ra_loop():
    time = Time(2457458.1739, scale='utc', format='jd')
    tee_ra = Angle(np.pi / 4., unit='rad')  # rad
    cirs_ra = pyuvsim.utils.tee_to_cirs_ra(tee_ra, time)
    new_tee_ra = pyuvsim.utils.cirs_to_tee_ra(cirs_ra, time)
    nt.assert_equal(new_tee_ra, tee_ra)
