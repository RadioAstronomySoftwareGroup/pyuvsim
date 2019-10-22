# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
#import healpy as hp
import os
import astropy.constants as const
import astropy_healpix
from astropy_healpix import HEALPix
from astropy import units
from astropy.coordinates import SkyCoord, EarthLocation, Angle, AltAz
from astropy.time import Time

import pyuvsim
import pyuvsim.tests as simtest
from pyuvdata import UVData
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim import utils as simutils

EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
single_offzenith_healpix_hdf5_file = os.path.join(SIM_DATA_PATH, 'single_off_zenith_healpix.hdf5')


EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
single_offzenith_healpix_hdf5_file = os.path.join(SIM_DATA_PATH, 'single_off_zenith_healpix.hdf5')


def test_source_zenith_from_icrs():
    """Test single source position at zenith constructed using icrs."""
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time = Time('2018-03-01 00:00:00', scale='utc', location=array_location)

    lst = time.sidereal_time('apparent')

    tee_ra = lst
    cirs_ra = pyuvsim.tee_to_cirs_ra(tee_ra, time)

    cirs_source_coord = SkyCoord(ra=cirs_ra, dec=array_location.lat,
                                 obstime=time, frame='cirs', location=array_location)

    icrs_coord = cirs_source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    # Check error cases
    simtest.assert_raises_message(
        ValueError, 'ra must be an astropy Angle object. value was: 3.14',
        pyuvsim.SkyModel, 'icrs_zen', ra.rad, dec.rad, [1, 0, 0, 0]
    )
    simtest.assert_raises_message(
        ValueError, 'dec must be an astropy Angle object. value was: -0.53',
        pyuvsim.SkyModel, 'icrs_zen', ra, dec.rad, [1, 0, 0, 0]
    )
    zenith_source = pyuvsim.SkyModel('icrs_zen', ra, dec, [1, 0, 0, 0])

    zenith_source.update_positions(time, array_location)
    zenith_source_lmn = zenith_source.pos_lmn.squeeze()
    assert np.allclose(zenith_source_lmn, np.array([0, 0, 1]), atol=1e-5)


def test_source_zenith():
    """Test single source position at zenith."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    source_arr, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith',
                                                array_location=array_location)
    zenith_source = source_arr

    zenith_source.update_positions(time, array_location)
    zenith_source_lmn = zenith_source.pos_lmn.squeeze()
    assert np.allclose(zenith_source_lmn, np.array([0, 0, 1]))


def test_sources_equal():
    time = Time('2018-03-01 00:00:00', scale='utc')
    src1, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    src2, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    assert src1 == src2


def test_pol_coherency_calc():
    Ncomp = 100
    ra = Angle(np.linspace(0.0, 2 * np.pi, Ncomp), unit='rad')
    dec = Angle(np.linspace(-np.pi / 2., np.pi / 2., Ncomp), unit='rad')
    names = np.arange(Ncomp).astype('str')
    freq = 1e8 * units.Hz

    stokes = np.zeros((4, 1, Ncomp))
    stokes[0, :, :] = 1.0
    stokes[1, :, :Ncomp // 2] = 2.5

    sky = pyuvsim.source.SkyModel(names, ra, dec, stokes, freq_array=freq)

    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='0.0d', lon='21d25m41.9s',
                                   height=1073.)
    sky.update_positions(time, array_location)
    coh_loc = sky.coherency_calc(array_location)

    inds = np.ones(Ncomp).astype(bool)
    inds[Ncomp // 2:] = False

    unpol_srcs_up = (sky.alt_az[0] > 0) * (~inds)

    assert np.allclose(coh_loc[0, 0, 0, unpol_srcs_up], 0.5)
    assert np.allclose(coh_loc[1, 1, 0, unpol_srcs_up], 0.5)
    assert np.allclose(coh_loc[1, 0, 0, unpol_srcs_up], 0.0)
    assert np.allclose(coh_loc[0, 1, 0, unpol_srcs_up], 0.0)


def test_calc_basis_rotation_matrix():
    """
    This tests whether the 3-D rotation matrix from RA/Dec to Alt/Az is
    actually a rotation matrix (R R^T = R^T R = I)
    """

    time = Time('2018-01-01 00:00')
    telescope_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)

    source = pyuvsim.source.SkyModel('Test', Angle(12. * units.hr),
                                     Angle(-30. * units.deg),
                                     150 * units.MHz, [1., 0., 0., 0.])
    source.update_positions(time, telescope_location)

    basis_rot_matrix = source._calc_average_rotation_matrix(telescope_location)

    assert np.allclose(np.matmul(basis_rot_matrix, basis_rot_matrix.T), np.eye(3))
    assert np.allclose(np.matmul(basis_rot_matrix.T, basis_rot_matrix), np.eye(3))


def test_calc_vector_rotation():
    """
    This checks that the 2-D coherency rotation matrix is unit determinant.
    I suppose we could also have checked (R R^T = R^T R = I)
    """

    time = Time('2018-01-01 00:00')
    telescope_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)

    source = pyuvsim.source.SkyModel('Test', Angle(12. * units.hr),
                                     Angle(-30. * units.deg),
                                     150 * units.MHz, [1., 0., 0., 0.])
    source.update_positions(time, telescope_location)

    coherency_rotation = np.squeeze(source._calc_coherency_rotation(telescope_location))

    assert(np.isclose(np.linalg.det(coherency_rotation), 1))


def analytic_beam_jones(za, az, sigma=0.3):
    """
    Analytic beam with sensible polarization response.

    Required for testing polarized sources.
    """
    # B = np.exp(-np.tan(za/2.)**2. / 2. / sigma**2.)
    B = 1
    J = np.array([[np.cos(za) * np.sin(az), np.cos(az)],
                  [np.cos(az) * np.cos(za), -np.sin(az)]])
    return B * J


def test_polarized_source_visibilities():
    """Test that visibilities of a polarized source match prior calculations."""
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)
    time0 = Time('2018-03-01 18:00:00', scale='utc', location=array_location)

    ha_off = 1 / 6.
    ha_delta = 0.1
    time_offsets = np.arange(-ha_off, ha_off + ha_delta, ha_delta)
    zero_indx = np.argmin(np.abs(time_offsets))
    # make sure we get a true zenith time
    time_offsets[zero_indx] = 0.
    times = time0 + time_offsets * units.hr
    ntimes = times.size

    zenith = SkyCoord(alt=90. * units.deg, az=0 * units.deg, frame='altaz',
                      obstime=time0, location=array_location)
    zenith_icrs = zenith.transform_to('icrs')

    src_astropy = SkyCoord(ra=zenith_icrs.ra, dec=zenith_icrs.dec,
                           obstime=times, location=array_location)
    src_astropy_altaz = src_astropy.transform_to('altaz')
    assert np.isclose(src_astropy_altaz.alt.rad[zero_indx], np.pi / 2)

    freq = (150e6 * units.Hz)
    stokes_radec = [1, -0.2, 0.3, 0.1]

    decoff = 0.0 * units.arcmin  # -0.17 * units.arcsec
    raoff = 0.0 * units.arcsec

    source = pyuvsim.SkyModel('icrs_zen', zenith_icrs.ra + raoff,
                              zenith_icrs.dec + decoff, freq, stokes_radec)

    coherency_matrix_local = np.zeros([2, 2, ntimes], dtype='complex128')
    alts = np.zeros(ntimes)
    azs = np.zeros(ntimes)
    for ti, time in enumerate(times):
        source.update_positions(time, telescope_location=array_location)
        alt, az = source.alt_az
        assert alt == src_astropy_altaz[ti].alt.radian
        assert az == src_astropy_altaz[ti].az.radian
        alts[ti] = alt
        azs[ti] = az

        coherency_tmp = source.coherency_calc(array_location).squeeze()
        coherency_matrix_local[:, :, ti] = coherency_tmp

    zas = np.pi / 2. - alts
    Jbeam = analytic_beam_jones(zas, azs)
    coherency_instr_local = np.einsum('ab...,bc...,dc...->ad...', Jbeam,
                                      coherency_matrix_local, np.conj(Jbeam))

    expected_instr_local = np.array(
        [[[0.60632557 + 0.00000000e+00j, 0.6031185 - 2.71050543e-20j,
           0.60059597 + 0.00000000e+00j, 0.59464231 + 5.42101086e-20j,
           0.58939657 + 0.00000000e+00j],
          [0.14486082 + 4.99646382e-02j, 0.14776209 + 4.99943414e-02j,
           0.14960097 + 5.00000000e-02j, 0.15302905 + 4.99773672e-02j,
           0.15536376 + 4.99307015e-02j]],
         [[0.14486082 - 4.99646382e-02j, 0.14776209 - 4.99943414e-02j,
           0.14960097 - 5.00000000e-02j, 0.15302905 - 4.99773672e-02j,
           0.15536376 - 4.99307015e-02j],
          [0.39282051 + 0.00000000e+00j, 0.39674527 + 0.00000000e+00j,
           0.39940403 + 0.00000000e+00j, 0.40481652 + 5.42101086e-20j,
           0.40895287 + 0.00000000e+00j]]])

    assert np.allclose(coherency_instr_local, expected_instr_local)


def test_polarized_source_smooth_visibilities():
    """Test that visibilities change smoothly as a polarized source transits."""
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)
    time0 = Time('2018-03-01 18:00:00', scale='utc', location=array_location)

    ha_off = 1
    ha_delta = 0.01
    time_offsets = np.arange(-ha_off, ha_off + ha_delta, ha_delta)
    zero_indx = np.argmin(np.abs(time_offsets))
    # make sure we get a true zenith time
    time_offsets[zero_indx] = 0.
    times = time0 + time_offsets * units.hr
    ntimes = times.size

    zenith = SkyCoord(alt=90. * units.deg, az=0 * units.deg, frame='altaz',
                      obstime=time0, location=array_location)
    zenith_icrs = zenith.transform_to('icrs')

    src_astropy = SkyCoord(ra=zenith_icrs.ra, dec=zenith_icrs.dec,
                           obstime=times, location=array_location)
    src_astropy_altaz = src_astropy.transform_to('altaz')
    assert np.isclose(src_astropy_altaz.alt.rad[zero_indx], np.pi / 2)

    freq = (150e6 * units.Hz)
    stokes_radec = [1, -0.2, 0.3, 0.1]

    source = pyuvsim.SkyModel('icrs_zen', zenith_icrs.ra,
                              zenith_icrs.dec, freq, stokes_radec)

    coherency_matrix_local = np.zeros([2, 2, ntimes], dtype='complex128')
    alts = np.zeros(ntimes)
    azs = np.zeros(ntimes)
    for ti, time in enumerate(times):
        source.update_positions(time, telescope_location=array_location)
        alt, az = source.alt_az
        assert alt == src_astropy_altaz[ti].alt.radian
        assert az == src_astropy_altaz[ti].az.radian
        alts[ti] = alt
        azs[ti] = az

        coherency_tmp = source.coherency_calc(array_location).squeeze()
        coherency_matrix_local[:, :, ti] = coherency_tmp

    zas = np.pi / 2. - alts
    Jbeam = analytic_beam_jones(zas, azs)
    coherency_instr_local = np.einsum('ab...,bc...,dc...->ad...', Jbeam,
                                      coherency_matrix_local, np.conj(Jbeam))

    # test that all the instrumental coherencies are smooth
    t_diff_sec = np.diff(times.jd) * 24 * 3600
    for pol_i in [0, 1]:
        for pol_j in [0, 1]:
            real_coherency = coherency_instr_local[pol_i, pol_j, :].real
            real_derivative = np.diff(real_coherency) / t_diff_sec
            real_derivative_diff = np.diff(real_derivative)
            assert np.max(np.abs(real_derivative_diff)) < 1e-6
            imag_coherency = coherency_instr_local[pol_i, pol_j, :].imag
            imag_derivative = np.diff(imag_coherency) / t_diff_sec
            imag_derivative_diff = np.diff(imag_derivative)
            assert np.max(np.abs(imag_derivative_diff)) < 1e-6

    # test that the stokes coherencies are smooth
    stokes_instr_local = simutils.coherency_to_stokes(coherency_instr_local)
    for pol_i in range(4):
        real_stokes = stokes_instr_local[pol_i, :].real
        real_derivative = np.diff(real_stokes) / t_diff_sec
        real_derivative_diff = np.diff(real_derivative)
        assert np.max(np.abs(real_derivative_diff)) < 1e-6
        imag_stokes = stokes_instr_local[pol_i, :].imag
        assert np.all(imag_stokes == 0)


def test_read_healpix_hdf5():
    Nside = 32
    # hp_obj = HEALPix(nside=Nside)
    Npix = astropy_healpix.nside_to_npix(Nside)
    # ipix_disc = hp_obj.cone_search_lonlat((np.pi / 2) * units.rad,
    # (np.pi * 3 / 4) * units.rad, radius = 10 * units.rad)
    # Npix = hp.nside2npix(Nside)
    # vec = astropy_healpix.healpy.ang2vec(np.pi / 2, np.pi * 3 / 4)
    # vec = hp.ang2vec(np.pi / 2, np.pi * 3 / 4)
    # ipix_disc = hp.query_disc(nside=32, vec=vec, radius=np.radians(10))
    m = np.arange(Npix)
    ipix_disc = [5103, 5104, 5231, 5232, 5233, 5358, 5359, 5360, 5361, 5486, 5487, 5488, 5489, 5490,
                 5613, 5614, 5615, 5616, 5617, 5618, 5741, 5742, 5743, 5744, 5745, 5746, 5747, 5869,
                 5870, 5871, 5872, 5873, 5874, 5997, 5998, 5999, 6000, 6001, 6002, 6003, 6124, 6125,
                 6126, 6127, 6128, 6129, 6130, 6131, 6253, 6254, 6255, 6256, 6257, 6258, 6259, 6381,
                 6382, 6383, 6384, 6385, 6386, 6509, 6510, 6511, 6512, 6513, 6514, 6515, 6637, 6638,
                 6639, 6640, 6641, 6642, 6766, 6767, 6768, 6769, 6770, 6894, 6895, 6896, 6897, 7023,
                 7024, 7025, 7151, 7152]
    m[ipix_disc] = m.max()

    indices = np.arange(Npix)

    frequencies = np.linspace(100, 110, 10)

    hpmap, inds, freqs = pyuvsim.source.read_healpix_hdf5(
        os.path.join(SIM_DATA_PATH, 'test_file.hdf5')
    )

    assert np.allclose(hpmap[0, :], m)
    assert np.allclose(inds, indices)
    assert np.allclose(freqs, frequencies)


def test_healpix_to_sky():
    Nside = 32
    Npix = astropy_healpix.nside_to_npix(Nside)
    #vec = hp.ang2vec(np.pi / 2, np.pi * 3 / 4)
    #ipix_disc = hp.query_disc(nside=32, vec=vec, radius=np.radians(10))
    m = np.arange(Npix)
    ipix_disc = [5103, 5104, 5231, 5232, 5233, 5358, 5359, 5360, 5361, 5486, 5487, 5488, 5489, 5490,
                 5613, 5614, 5615, 5616, 5617, 5618, 5741, 5742, 5743, 5744, 5745, 5746, 5747, 5869,
                 5870, 5871, 5872, 5873, 5874, 5997, 5998, 5999, 6000, 6001, 6002, 6003, 6124, 6125,
                 6126, 6127, 6128, 6129, 6130, 6131, 6253, 6254, 6255, 6256, 6257, 6258, 6259, 6381,
                 6382, 6383, 6384, 6385, 6386, 6509, 6510, 6511, 6512, 6513, 6514, 6515, 6637, 6638,
                 6639, 6640, 6641, 6642, 6766, 6767, 6768, 6769, 6770, 6894, 6895, 6896, 6897, 7023,
                 7024, 7025, 7151, 7152]
    m[ipix_disc] = m.max()

    m = np.repeat(m[None, :], 10, axis=0)
    hpmap, inds, freqs = pyuvsim.source.read_healpix_hdf5(
        os.path.join(SIM_DATA_PATH, 'test_file.hdf5')
    )
    m = (m.T / simutils.jy2Tsr(freqs, bm=astropy_healpix.nside_to_pixel_area(Nside), mK=False)).T
    sky = pyuvsim.source.healpix_to_sky(hpmap, inds, freqs)
    print(m)
    print(sky.stokes[0])
    assert np.allclose(sky.stokes[0], m.value)


def test_units_healpix_to_sky():
    Nside = 32
    beam_area = astropy_healpix.nside_to_pixel_area(Nside)  # * units.sr
    #beam_area = hp.pixelfunc.nside2pixarea(Nside) * units.sr
    hpmap, inds, freqs = pyuvsim.source.read_healpix_hdf5(
        os.path.join(SIM_DATA_PATH, 'test_file.hdf5')
    )
    freqs = freqs * units.Hz
    stokes = (hpmap.T * units.K).to(units.Jy, units.brightness_temperature(beam_area, freqs)).T
    sky = pyuvsim.source.healpix_to_sky(hpmap, inds, freqs)

    assert np.allclose(sky.stokes[0, 0], stokes.value[0])


def test_single_offzenith_source_uvfits():
    """Test single off-zenith source in a healpix map using test uvdata file."""
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file, ant_str='cross')
    hera_uv.unphase_to_drift()

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

    Nside = 8
    Npix = astropy_healpix.nside_to_npix(Nside)
    freqs = np.arange(100, 100.5, 0.1)
    Nfreqs = len(freqs)
    m = np.zeros((Npix, Nfreqs))
    ipix = 357
    m[ipix, :] = [3254807413.337762, 3248307549.9303513,
                  3241827137.4792957, 3235366098.452164, 3228924355.702452]

    Nskies = 1
    inds = np.arange(Npix)
    dataset = np.zeros((Nskies, Nfreqs, len(m)))
    for j in range(0, len(m)):
        dataset[0, :, j] = freqs
    for i in range(0, Nfreqs):
        dataset[0, i, :] = m[:, i]

    ra, dec = astropy_healpix.healpix_to_lonlat(ipix, Nside)
    skycoord_use = SkyCoord(ra, dec, frame='icrs')
    source_altaz = skycoord_use.transform_to(AltAz(obstime=time, location=array_location))
    alt_az = np.array([source_altaz.alt.value, source_altaz.az.value])

    src_az = Angle(alt_az[1], unit='deg')
    src_alt = Angle(alt_az[0], unit='deg')
    src_za = Angle('90.d') - src_alt

    src_l = np.sin(src_az.rad) * np.sin(src_za.rad)
    src_m = np.cos(src_az.rad) * np.sin(src_za.rad)
    src_n = np.cos(src_za.rad)

    hpmap_hpx, indices_hpx, freqs_hpx = pyuvsim.source.read_healpix_hdf5(
        single_offzenith_healpix_hdf5_file)
    sky_hpx = pyuvsim.source.healpix_to_sky(hpmap_hpx, indices_hpx, freqs_hpx)

    time.location = array_location

    sky_hpx.update_positions(time, array_location)
    src_alt_az = sky_hpx.alt_az
    assert np.isclose(src_alt_az[0][ipix], src_alt.rad)
    assert np.isclose(src_alt_az[1][ipix], src_az.rad)

    src_lmn = sky_hpx.pos_lmn
    assert np.isclose(src_lmn[0][ipix], src_l)
    assert np.isclose(src_lmn[1][ipix], src_m)
    assert np.isclose(src_lmn[2][ipix], src_n)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beam_list = [beam]

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    task = pyuvsim.UVTask(sky_hpx, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    # analytically calculate visibility
    beam.peak_normalize()
    beam.interpolation_function = 'az_za_simple'
    beam_za, beam_az = simutils.altaz_to_zenithangle_azimuth(src_alt.rad, src_az.rad)
    beam_za2, beam_az2 = simutils.altaz_to_zenithangle_azimuth(
        src_alt_az[0][ipix], src_alt_az[1][ipix])

    assert np.isclose(beam_za, beam_za2)
    assert np.isclose(beam_az, beam_az2)

    print(beam_az)

    interpolated_beam, interp_basis_vector = beam.interp(az_array=np.array([beam_az]),
                                                         za_array=np.array([beam_za]),
                                                         freq_array=np.array([freq.to('Hz').value]))
    jones = np.zeros((2, 2, 1), dtype=np.complex64)
    jones[0, 0] = interpolated_beam[1, 0, 0, 0, 0]
    jones[1, 1] = interpolated_beam[0, 0, 1, 0, 0]
    jones[1, 0] = interpolated_beam[1, 0, 1, 0, 0]
    jones[0, 1] = interpolated_beam[0, 0, 0, 0, 0]

    beam_jones = antenna1.get_beam_jones(array, src_alt_az, freq)
    print(beam_jones[:, :, ipix])
    print(jones)
    # print(beam_jones - jones)
    assert np.allclose(beam_jones[:, :, ipix], jones.squeeze())

    uvw_wavelength_array = hera_uv.uvw_array * units.m / const.c * freq.to('1/s')
    jones_T = np.swapaxes(jones, 0, 1)
    # Remove source axis from jones matrix
    jones = jones.squeeze()
    vis_analytic = 0.5 * np.dot(jones, np.conj(jones).T) * np.exp(2j * np.pi * (
        uvw_wavelength_array[0, 0] * src_l + uvw_wavelength_array[0, 1] * src_m + uvw_wavelength_array[0, 2] * src_n))
    vis_analytic = np.array([vis_analytic[0, 0], vis_analytic[1, 1],
                             vis_analytic[0, 1], vis_analytic[1, 0]])

    assert np.allclose(baseline.uvw.to('m').value, hera_uv.uvw_array[0:hera_uv.Nbls], atol=1e-4)

    print(vis_analytic)
    print(visibility)
    print(vis_analytic - visibility)
    assert np.allclose(visibility, vis_analytic, atol=1e-4)
