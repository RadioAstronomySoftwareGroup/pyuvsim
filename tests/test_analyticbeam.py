# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import re

import numpy as np
import pytest
from astropy import units
from astropy.coordinates import Angle, EarthLocation
from astropy.time import Time
from pyuvdata import UVData
from scipy.special import j1

import pyuvsim
import pyuvsim.utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

EW_uvfits_file = os.path.join(SIM_DATA_PATH, "28mEWbl_1time_1chan.uvfits")

c_ms = pyuvsim.analyticbeam.c_ms


@pytest.fixture
def heratext_posfreq():
    time = Time("2018-03-01 00:00:00", scale="utc")
    array_location = EarthLocation(lat="-30d43m17.5s", lon="21d25m41.9s", height=1073.0)
    sources, _ = pyuvsim.create_mock_catalog(
        time, "hera_text", array_location=array_location
    )

    sources.update_positions(time, array_location)
    za_vals = np.pi / 2.0 - sources.alt_az[1]  # rad
    az_vals = sources.alt_az[1]

    freq_vals = np.array([10**8])

    return az_vals, za_vals, freq_vals


def test_uniform_beam(heratext_posfreq):
    beam = pyuvsim.AnalyticBeam("uniform")
    beam.peak_normalize()

    az_vals, za_vals, freqs = heratext_posfreq

    nsrcs = az_vals.size
    n_freqs = freqs.size

    interpolated_beam, _ = beam.interp(
        az_array=az_vals, za_array=za_vals, freq_array=freqs
    )
    expected_data = np.zeros((2, 2, n_freqs, nsrcs), dtype=float)
    expected_data[1, 0, :, :] = 1
    expected_data[0, 1, :, :] = 1
    np.testing.assert_allclose(interpolated_beam, expected_data)


def test_airy_beam_values(heratext_posfreq):
    diameter_m = 14.0
    beam = pyuvsim.AnalyticBeam("airy", diameter=diameter_m)
    beam.peak_normalize()

    az_vals, za_vals, freq_vals = heratext_posfreq

    interpolated_beam, _ = beam.interp(
        az_array=az_vals, za_array=za_vals, freq_array=freq_vals
    )

    expected_data = np.zeros((2, 2, 1, az_vals.size), dtype=float)
    za_grid, f_grid = np.meshgrid(za_vals, freq_vals)
    xvals = diameter_m / 2.0 * np.sin(za_grid) * 2.0 * np.pi * f_grid / c_ms
    airy_values = np.zeros_like(xvals)
    nz = xvals != 0.0
    ze = xvals == 0.0
    airy_values[nz] = 2.0 * j1(xvals[nz]) / xvals[nz]
    airy_values[ze] = 1.0
    expected_data[1, 0, :, :] = airy_values
    expected_data[0, 1, :, :] = airy_values

    np.testing.assert_allclose(interpolated_beam, expected_data)


def test_interp_errors(heratext_posfreq):
    diameter_m = 14.0
    beam = pyuvsim.AnalyticBeam("airy", diameter=diameter_m)
    beam.peak_normalize()

    az_vals, za_vals, freq_vals = heratext_posfreq

    az_mesh, za_mesh = np.meshgrid(az_vals, za_vals)

    with pytest.raises(
        ValueError,
        match="az_array, za_array and freq_array must all be one dimensional.",
    ):
        beam.interp(az_array=az_mesh, za_array=za_mesh, freq_array=freq_vals)

    with pytest.raises(
        ValueError, match="az_array and za_array must have the same shape."
    ):
        beam.interp(az_array=az_vals, za_array=za_vals[0:-1], freq_array=freq_vals)


def test_uv_beam_widths():
    # Check that the width of the Airy disk beam in UV space corresponds with the dish diameter.
    diameter_m = 25.0
    beam = pyuvsim.AnalyticBeam("airy", diameter=diameter_m)
    beam.peak_normalize()

    Nfreqs = 20
    freq_vals = np.linspace(100e6, 130e6, Nfreqs)
    lams = c_ms / freq_vals

    N = 250
    Npix = 500
    zmax = np.radians(90)  # Degrees
    arr = np.arange(-N, N)
    x_arr, y_arr = np.meshgrid(arr, arr)
    x_arr = x_arr.flatten()
    y_arr = y_arr.flatten()
    radius = np.sqrt(x_arr**2 + y_arr**2) / float(N)
    zas = radius * zmax
    azs = np.arctan2(y_arr, x_arr)
    interpolated_beam, _ = beam.interp(
        az_array=np.asarray(azs),
        za_array=np.asarray(zas),
        freq_array=np.array(freq_vals),
    )

    ebeam = interpolated_beam[0, 1, :, :]
    ebeam = ebeam.reshape(Nfreqs, Npix, Npix)
    beam_kern = np.fft.fft2(ebeam, axes=(1, 2))
    beam_kern = np.fft.fftshift(beam_kern, axes=(1, 2))
    for i, bk in enumerate(beam_kern):
        # Cutoff at half a % of the maximum value in Fourier space.
        thresh = np.max(np.abs(bk)) * 0.005
        points = np.sum(np.abs(bk) >= thresh)
        upix = 1 / (
            2 * np.sin(zmax)
        )  # 2*sin(zmax) = fov extent projected onto the xy plane
        area = np.sum(points) * upix**2
        kern_radius = np.sqrt(area / np.pi)
        assert np.isclose(diameter_m / lams[i], kern_radius, rtol=0.5)


def test_achromatic_gaussian_beam(heratext_posfreq):
    sigma_rad = Angle("5d").to_value("rad")
    beam = pyuvsim.AnalyticBeam("gaussian", sigma=sigma_rad)
    beam.peak_normalize()

    az_vals, za_vals, freq_vals = heratext_posfreq
    nsrcs = az_vals.size
    n_freqs = freq_vals.size

    interpolated_beam, interp_basis_vector = beam.interp(
        az_array=np.array(az_vals),
        za_array=np.array(za_vals),
        freq_array=np.array(freq_vals),
    )

    expected_data = np.zeros((2, 2, n_freqs, nsrcs), dtype=float)
    interp_zas = np.zeros((n_freqs, nsrcs), dtype=float)
    for f_ind in range(n_freqs):
        interp_zas[f_ind, :] = np.array(za_vals)
    gaussian_vals = np.exp(-(interp_zas**2) / (2 * sigma_rad**2))

    expected_data[1, 0, :, :] = gaussian_vals
    expected_data[0, 1, :, :] = gaussian_vals

    np.testing.assert_allclose(interpolated_beam, expected_data)


@pytest.mark.filterwarnings("ignore:UVW orientation appears to be flipped")
@pytest.mark.filterwarnings("ignore:The shapes of several attributes will be changing")
def test_gaussbeam_values():
    """
    Make the long-line point sources up to 10 degrees from zenith.
    Confirm that the coherencies match the expected beam values at those zenith angles.
    """
    sigma = 0.05
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    if hasattr(hera_uv, "use_current_array_shapes"):
        hera_uv.use_future_array_shapes()

    if hasattr(hera_uv, "telescope"):
        array_location = hera_uv.telescope.location
    else:
        # this can be removed when we require pyuvdata >= 3.0
        array_location = EarthLocation.from_geocentric(
            *hera_uv.telescope_location, unit="m"
        )
    freq = hera_uv.freq_array[0] * units.Hz

    time = Time(hera_uv.time_array[0], scale="utc", format="jd")

    catalog, _ = pyuvsim.create_mock_catalog(
        time=time,
        arrangement="long-line",
        Nsrcs=41,
        min_alt=80.0,
        array_location=array_location,
    )

    catalog.update_positions(time, array_location)
    beam = pyuvsim.AnalyticBeam("gaussian", sigma=sigma)
    array = pyuvsim.Telescope("telescope_name", array_location, [beam])

    # Need a dummy baseline for this test.
    antenna1 = pyuvsim.Antenna("ant1", 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna("ant2", 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    task = pyuvsim.UVTask(catalog, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)
    engine.apply_beam()
    altitudes = task.sources.alt_az[0]  # In radians.
    # All four components should be identical
    if isinstance(engine.apparent_coherency, units.Quantity):
        coherency_use = engine.apparent_coherency.to_value("Jy")
    else:
        coherency_use = engine.apparent_coherency

    coherencies = np.real(coherency_use[0, 0] + coherency_use[1, 1]).astype(float)

    zenith_angles, _ = simutils.altaz_to_zenithangle_azimuth(
        altitudes, np.zeros_like(np.array(altitudes))
    )

    # Confirm the coherency values (ie., brightnesses) match the beam values.
    beam_values = np.exp(-((zenith_angles) ** 2) / (2 * beam.sigma**2))
    assert np.all(beam_values**2 == coherencies)


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
    za = np.linspace(0, np.pi / 2.0, Npix)

    gauss = pyuvsim.AnalyticBeam(
        "gaussian", sigma=sigma, ref_freq=freqs[0], spectral_index=alpha
    )

    # Get the widths at each frequency.

    vals, _ = gauss.interp(az, za, freqs)

    vals = vals[0, 1]

    for fi in range(Nfreqs):
        hwhm = za[np.argmin(np.abs(vals[fi] - 0.5))]
        sig_f = sigma * (freqs[fi] / freqs[0]) ** alpha
        assert np.isclose(sig_f, 2 * hwhm / 2.355, atol=1e-3)


def test_chromatic_gaussian_error():
    # Error if trying to define chromatic beam without a reference frequency
    with pytest.raises(
        ValueError,
        match="ref_freq must be set for nonzero gaussian beam spectral index",
    ):
        pyuvsim.AnalyticBeam("gaussian", sigma=np.radians(15.0), spectral_index=-1.5)


def test_power_analytic_beam():
    # Check that power beam evaluation matches electric field amp**2 for analytic beams.
    freqs = np.arange(120e6, 160e6, 4e6)
    Npix = 1000
    diam = 14.0

    az = np.zeros(Npix)
    za = np.linspace(0, np.pi / 2.0, Npix)

    for b in ["gaussian", "uniform", "airy"]:
        eb = pyuvsim.AnalyticBeam(b, diameter=diam)
        pb = pyuvsim.AnalyticBeam(b, diameter=diam)
        pb.efield_to_power()
        evals = eb.interp(az, za, freqs)[0][0, 1]
        pvals = pb.interp(az, za, freqs)[0][0, 0]
        np.testing.assert_allclose(evals**2, pvals)


def test_comparison():
    """
    Beam __eq__ method
    """
    beam1 = pyuvsim.AnalyticBeam("uniform")
    beam2 = pyuvsim.AnalyticBeam("gaussian", sigma=0.02)
    beam2.type = "undefined"

    not_beam = UVData()
    assert beam1 != not_beam
    assert beam2 != beam1


def test_beam_init_errs():
    """
    Error cases.
    """
    with pytest.raises(ValueError, match="type not recognized"):
        pyuvsim.AnalyticBeam("unsupported_type")


@pytest.mark.parametrize(
    ("beam_type", "error_msg"),
    [
        (
            "gaussian",
            re.escape(
                "Antenna diameter (meters) or sigma (radians) needed for gaussian beams."
            ),
        ),
        ("airy", "Antenna diameter needed for airy beam"),
        ("noninterpolable", "no interp for this type: noninterpolable"),
    ],
)
def test_beam_interp_errs(beam_type, error_msg):
    if beam_type == "noninterpolable":
        beam = pyuvsim.AnalyticBeam("uniform")
        beam.type = "noninterpolable"
    else:
        beam = pyuvsim.AnalyticBeam(beam_type)
    az, za = np.random.uniform(0.0, np.pi, (2, 5))
    freq_arr = np.linspace(1e8, 1.5e8, 10)
    with pytest.raises(ValueError, match=error_msg):
        beam.interp(az, za, freq_arr)


def test_diameter_to_sigma():
    # The integrals of an Airy power beam and a Gaussian power beam, within
    # the first Airy null, should be close if the Gaussian width is set to the Airy width.
    diameter_m = 25.0
    abm = pyuvsim.AnalyticBeam("airy", diameter=diameter_m)
    gbm = pyuvsim.AnalyticBeam("gaussian", diameter=diameter_m)

    Nfreqs = 20
    freq_vals = np.linspace(100e6, 130e6, Nfreqs)
    lams = c_ms / freq_vals

    N = 250
    Npix = 501
    zmax = np.radians(40)  # Degrees

    zas = np.linspace(-zmax, zmax, Npix)
    azs = np.array([0.0] * (N + 1) + [np.pi] * N)

    shape = (2, 1, 2, Nfreqs) + azs.shape
    airy_vals, interp_basis_vector = abm.interp(
        az_array=azs.flatten(), za_array=zas.flatten(), freq_array=freq_vals
    )

    gauss_vals, interp_basis_vector = gbm.interp(
        az_array=azs.flatten(), za_array=zas.flatten(), freq_array=freq_vals
    )

    airy_vals = airy_vals.reshape(shape)
    gauss_vals = gauss_vals.reshape(shape)

    airy_vals = airy_vals[0, 0, 0] * airy_vals[0, 0, 1]
    gauss_vals = (
        gauss_vals[0, 0, 0] * gauss_vals[0, 0, 1]
    )  # Remove pol/spw/feed axes. Make power beam.

    for fi in range(Nfreqs):
        null = 1.22 * lams[fi] / diameter_m
        inds = np.where(np.abs(zas) < null)

        # Assert integral of power beams within the first Airy null are close
        assert np.isclose(
            np.sum(airy_vals[fi, inds]), np.sum(gauss_vals[fi, inds]), rtol=1e-2
        )
