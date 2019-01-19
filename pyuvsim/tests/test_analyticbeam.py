# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import numpy as np
import nose.tools as nt
from scipy.special import j1
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy import units

from pyuvdata import UVData

import pyuvsim
import pyuvsim.utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')


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

        beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_coord_altaz.alt.to('rad').value,
                                                                 src_coord_altaz.az.to('rad').value)
        az_vals.append(beam_az)
        za_vals.append(beam_za)

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

        beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_coord_altaz.alt.to('rad').value,
                                                                 src_coord_altaz.az.to('rad').value)
        az_vals.append(beam_az)
        za_vals.append(beam_za)

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
    nz = xvals != 0.
    ze = xvals == 0.
    airy_vals[nz] = (2. * j1(xvals[nz]) / xvals[nz])**2
    airy_vals[ze] = 1.

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

        beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_coord_altaz.alt.to('rad').value,
                                                                 src_coord_altaz.az.to('rad').value)
        az_vals.append(beam_az)
        za_vals.append(beam_za)

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
                                                         min_alt=80., array_location=array_location)

    beam = pyuvsim.AnalyticBeam('gaussian', sigma=sigma)
    array = pyuvsim.Telescope('telescope_name', array_location, [beam])

    # Need a dummy baseline for this test.
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    coherencies = []
    altitudes = []
    for src in catalog:
        task = pyuvsim.UVTask(src, time, freq, baseline, array)
        engine = pyuvsim.UVEngine(task)
        engine.apply_beam()
        altitudes.append(task.source.alt_az[0])   # In radians.
        coherencies.append(np.real(engine.apparent_coherency[0, 0]).astype(float))  # All four components should be identical

    coherencies = np.array(coherencies)
    zenith_angles, _ = simutils.altaz_to_zenithangle_azimuth(altitudes,
                                                             np.zeros_like(np.array(altitudes)))

    # Confirm the coherency values (ie., brightnesses) match the beam values.

    beam_values = np.exp(-(zenith_angles)**2 / (2 * beam.sigma**2))
    nt.assert_true(np.all(beam_values**2 == coherencies))


def test_beamerrs():
    """
    Error cases.
    """
    nt.assert_raises(ValueError, pyuvsim.AnalyticBeam, 'unsupported_type')
    beam = pyuvsim.AnalyticBeam('gaussian')
    az, za = np.random.uniform(0.0, np.pi, (2, 5))
    freq_arr = np.linspace(1e8, 1.5e8, 10)
    nt.assert_raises(ValueError, beam.interp, az, za, freq_arr)
    beam.type = 'airy'
    nt.assert_raises(ValueError, beam.interp, az, za, freq_arr)
    beam.type = 'noninterpretable'
    nt.assert_raises(ValueError, beam.interp, az, za, freq_arr)
