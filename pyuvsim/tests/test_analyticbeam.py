# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from astropy import units
from astropy.coordinates import Angle, EarthLocation
from astropy.time import Time
from pyuvdata import UVData
from scipy.special import j1

import pyuvsim
import pyuvsim.tests as simtest
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
    sources, mock_keywords = pyuvsim.create_mock_catalog(
        time, 'hera_text', array_location=array_location
    )

    nsrcs = sources.Ncomponents

    sources.update_positions(time, array_location)
    za_vals = np.pi / 2. - sources.alt_az[1]  # rad
    az_vals = sources.alt_az[1]
    freq_vals = sources.freq

    n_freqs = len(freq_vals)
    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=az_vals, za_array=za_vals, freq_array=freq_vals
    )
    expected_data = np.zeros((2, 1, 2, n_freqs, nsrcs), dtype=np.float)
    expected_data[1, 0, 0, :, :] = 1
    expected_data[0, 0, 1, :, :] = 1
    # expected_data[1, 0, 1, :, :] = 1
    # expected_data[0, 0, 0, :, :] = 1
    assert np.allclose(interpolated_beam, expected_data)


def test_airy_beam_values():
    diameter_m = 14.
    beam = pyuvsim.AnalyticBeam('airy', diameter=diameter_m)
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    time = Time('2018-03-01 00:00:00', scale='utc')
    array_location = EarthLocation(
        lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.
    )
    sources, mock_keywords = pyuvsim.create_mock_catalog(
        time, 'hera_text', array_location=array_location
    )

    sources.update_positions(time, array_location)
    za_vals = np.pi / 2. - sources.alt_az[1]  # rad
    az_vals = sources.alt_az[1]
    freq_vals = sources.freq.to("Hz").value

    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=az_vals, za_array=za_vals, freq_array=freq_vals
    )

    expected_data = np.zeros((2, 1, 2, freq_vals.size, az_vals.size), dtype=np.float)
    za_grid, f_grid = np.meshgrid(za_vals, freq_vals)
    xvals = diameter_m / 2. * np.sin(za_grid) * 2. * np.pi * f_grid / 3e8
    airy_values = np.zeros_like(xvals)
    nz = xvals != 0.
    ze = xvals == 0.
    airy_values[nz] = 2. * j1(xvals[nz]) / xvals[nz]
    airy_values[ze] = 1.
    expected_data[1, 0, 0, :, :] = airy_values
    expected_data[0, 0, 1, :, :] = airy_values

    assert np.allclose(interpolated_beam, expected_data)


def test_uv_beam_widths():
    diameter_m = 25.0
    beam = pyuvsim.AnalyticBeam('airy', diameter=diameter_m)
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    Nfreqs = 20
    freq_vals = np.linspace(100e6, 130e6, Nfreqs)
    lams = 3e8 / freq_vals

    N = 250
    Npix = 500
    zmax = np.radians(90)  # Degrees
    arr = np.arange(-N, N)
    x, y = np.meshgrid(arr, arr)
    r = np.sqrt(x ** 2 + y ** 2) / float(N)
    zas = r * zmax
    azs = np.arctan2(y, x)
    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array(azs), za_array=np.array(zas), freq_array=np.array(freq_vals)
    )

    ebeam = interpolated_beam[0, 0, 1, :, :]
    ebeam = ebeam.reshape(Nfreqs, Npix, Npix)
    beam_kern = np.fft.fft2(ebeam, axes=(1, 2))
    beam_kern = np.fft.fftshift(beam_kern, axes=(1, 2))
    for i, bk in enumerate(beam_kern):
        # Cutoff at half a % of the maximum value in Fourier space.
        thresh = np.max(np.abs(bk)) * 0.005
        points = np.sum(np.abs(bk) >= thresh)
        upix = 1 / (2 * np.sin(zmax))  # 2*sin(zmax) = fov extent projected onto the xy plane
        area = np.sum(points) * upix ** 2
        kern_radius = np.sqrt(area / np.pi)
        assert np.isclose(diameter_m / lams[i], kern_radius, rtol=0.5)


def test_achromatic_gaussian_beam():
    sigma_rad = Angle('5d').to('rad').value
    beam = pyuvsim.AnalyticBeam('gaussian', sigma=sigma_rad)
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'

    time = Time('2018-03-01 00:00:00', scale='utc')
    array_location = EarthLocation(
        lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.
    )
    sources, mock_keywords = pyuvsim.create_mock_catalog(
        time, 'hera_text', array_location=array_location
    )

    nsrcs = sources.Ncomponents

    sources.update_positions(time, array_location)
    za_vals = np.pi / 2. - sources.alt_az[1]  # rad
    az_vals = sources.alt_az[1]
    freq_vals = sources.freq.to("Hz").value

    n_freqs = len(freq_vals)
    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array(az_vals), za_array=np.array(za_vals), freq_array=np.array(freq_vals)
    )

    expected_data = np.zeros((2, 1, 2, n_freqs, nsrcs), dtype=np.float)
    interp_zas = np.zeros((n_freqs, nsrcs), dtype=np.float)
    for f_ind in range(n_freqs):
        interp_zas[f_ind, :] = np.array(za_vals)
    gaussian_vals = np.exp(-(interp_zas ** 2) / (2 * sigma_rad ** 2))

    expected_data[1, 0, 0, :, :] = gaussian_vals
    expected_data[0, 0, 1, :, :] = gaussian_vals

    assert np.allclose(interpolated_beam, expected_data)


def test_gaussbeam_values():
    """
    Make the long-line point sources up to 10 degrees from zenith.
    Obtain visibilities
    Confirm that the values match the expected beam values at those zenith angles.
    """
    sigma = 0.05
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    array_location = EarthLocation.from_geocentric(
        *hera_uv.telescope_location, unit='m'
    )
    freq = hera_uv.freq_array[0, 0] * units.Hz

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')

    catalog, mock_keywords = pyuvsim.create_mock_catalog(
        time=time, arrangement='long-line', Nsrcs=41, min_alt=80., array_location=array_location
    )

    catalog.update_positions(time, array_location)
    beam = pyuvsim.AnalyticBeam('gaussian', sigma=sigma)
    array = pyuvsim.Telescope('telescope_name', array_location, [beam])

    # Need a dummy baseline for this test.
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    task = pyuvsim.UVTask(catalog, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)
    engine.apply_beam()
    altitudes = task.sources.alt_az[0]  # In radians.
    # All four components should be identical
    coherencies = np.real(
        engine.apparent_coherency[0, 0] + engine.apparent_coherency[1, 1]
    ).astype(float)

    zenith_angles, _ = simutils.altaz_to_zenithangle_azimuth(
        altitudes, np.zeros_like(np.array(altitudes))
    )

    # Confirm the coherency values (ie., brightnesses) match the beam values.

    beam_values = np.exp(-(zenith_angles) ** 2 / (2 * beam.sigma ** 2))
    assert np.all(beam_values ** 2 == coherencies)


def test_chromatic_gaussian():
    """
    test_chromatic_gaussian
    Defining a gaussian beam with a spectral index and reference frequency.
    Check that beam width follows prescribed power law.
    """
    freqs = np.arange(120e6, 160e6, 4e6)
    Nfreqs = len(freqs)
    Npix = 1000
    alpha = -1.5
    sigma = np.radians(15.0)

    az = np.zeros(Npix)
    za = np.linspace(0, np.pi / 2., Npix)

    # Error if trying to define chromatic beam without a reference frequency

    simtest.assert_raises_message(
        ValueError, 'ref_freq must be set for nonzero gaussian beam spectral index',
        pyuvsim.AnalyticBeam, 'gaussian', sigma=sigma, spectral_index=alpha
    )
    A = pyuvsim.AnalyticBeam('gaussian', sigma=sigma, ref_freq=freqs[0], spectral_index=alpha)

    # Get the widths at each frequency.

    vals, _ = A.interp(az, za, freqs)

    vals = vals[0, 0, 1]

    for fi in range(Nfreqs):
        hwhm = za[np.argmin(np.abs(vals[fi] - 0.5))]
        sig_f = sigma * (freqs[fi] / freqs[0]) ** alpha
        assert np.isclose(sig_f, 2 * hwhm / 2.355, atol=1e-3)


def test_power_analytic_beam():
    freqs = np.arange(120e6, 160e6, 4e6)
    Npix = 1000
    diam = 14.0

    az = np.zeros(Npix)
    za = np.linspace(0, np.pi / 2., Npix)

    eb = pyuvsim.AnalyticBeam('gaussian', diameter=diam)
    pb = pyuvsim.AnalyticBeam('gaussian', diameter=diam)
    pb.efield_to_power()
    evals = eb.interp(az, za, freqs)[0][0, 0, 1]
    pvals = pb.interp(az, za, freqs)[0][0, 0, 0]
    assert np.allclose(evals ** 2, pvals)

    # Ensure uniform beam works
    pb = pyuvsim.AnalyticBeam('uniform')
    pb.efield_to_power()
    pb.interp(az, za, freqs)


def test_comparison():
    """
    Beam __eq__ method
    """
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam2.type = 'undefined'

    not_beam = UVData()
    assert beam1 != not_beam
    assert beam2 != beam1


def test_beamerrs():
    """
    Error cases.
    """
    simtest.assert_raises_message(
        ValueError, 'type not recognized', pyuvsim.AnalyticBeam, 'unsupported_type'
    )
    beam = pyuvsim.AnalyticBeam('gaussian')
    az, za = np.random.uniform(0.0, np.pi, (2, 5))
    freq_arr = np.linspace(1e8, 1.5e8, 10)
    simtest.assert_raises_message(
        ValueError, 'Dish diameter needed for gaussian beam -- units: meters',
        beam.interp, az, za, freq_arr
    )
    beam.type = 'airy'
    simtest.assert_raises_message(
        ValueError, 'Dish diameter needed for airy beam -- units: meters',
        beam.interp, az, za, freq_arr
    )
    beam.type = 'noninterpolable'
    simtest.assert_raises_message(
        ValueError, 'no interp for this type: noninterpolable', beam.interp, az, za, freq_arr
    )


def test_diameter_to_sigma():
    diameter_m = 25.0
    abm = pyuvsim.AnalyticBeam('airy', diameter=diameter_m)
    gbm = pyuvsim.AnalyticBeam('gaussian', diameter=diameter_m)

    Nfreqs = 20
    freq_vals = np.linspace(100e6, 130e6, Nfreqs)
    lams = 3e8 / freq_vals

    N = 250
    Npix = 501
    zmax = np.radians(40)  # Degrees

    zas = np.linspace(-zmax, zmax, Npix)
    azs = np.array([0.0] * (N + 1) + [np.pi] * N)

    shape = (2, 1, 2, Nfreqs,) + azs.shape
    airy_vals, interp_basis_vector = abm.interp(
        az_array=azs.flatten(), za_array=zas.flatten(), freq_array=freq_vals
    )

    gauss_vals, interp_basis_vector = gbm.interp(
        az_array=azs.flatten(), za_array=zas.flatten(), freq_array=freq_vals
    )

    airy_vals = airy_vals.reshape(shape)
    gauss_vals = gauss_vals.reshape(shape)

    airy_vals = airy_vals[0, 0, 0] * airy_vals[0, 0, 1]
    gauss_vals = gauss_vals[0, 0, 0] * gauss_vals[
        0, 0, 1]  # Remove pol/spw/feed axes. Make power beam.

    for fi in range(Nfreqs):
        null = 1.22 * lams[fi] / diameter_m
        inds = np.where(np.abs(zas) < null)

        # Assert integral of power beams within the first Airy null are close
        assert np.isclose(np.sum(airy_vals[fi, inds]), np.sum(gauss_vals[fi, inds]), rtol=1e-2)
