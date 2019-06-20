# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import healpy as hp
import h5py
import os
from astropy import units
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy.time import Time

import pyuvsim
import pyuvsim.tests as simtest
<<<<<<< HEAD
import pyuvsim.utils as simutils

=======
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
>>>>>>> 5b1f657... Added tests for healpix_to_sky and read_hdf5 functions

def test_source_zenith_from_icrs():
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
    # Check error cases
    simtest.assert_raises_message(
        ValueError, 'ra must be an astropy Angle object. value was: 3.14',
        pyuvsim.SkyModel, 'icrs_zen', ra.rad, dec.rad, freq.value, [1, 0, 0, 0]
    )
    simtest.assert_raises_message(
        ValueError, 'dec must be an astropy Angle object. value was: -0.53',
        pyuvsim.SkyModel, 'icrs_zen', ra, dec.rad, freq.value, [1, 0, 0, 0]
    )
    simtest.assert_raises_message(
        ValueError, 'freq must be an astropy Quantity object. value was: 150000000.0',
        pyuvsim.SkyModel, 'icrs_zen', ra, dec, freq.value, [1, 0, 0, 0]
    )
    zenith_source = pyuvsim.SkyModel('icrs_zen', ra, dec, freq, [1, 0, 0, 0])

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
    freq = np.ones(Ncomp) * 1e8 * units.Hz

    stokes = np.zeros((4, Ncomp))
    stokes[0, :] = 1.0
    stokes[1, :Ncomp // 2] = 2.5

    sky = pyuvsim.source.SkyModel(names, ra, dec, freq, stokes)

    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='0.0d', lon='21d25m41.9s',
                                   height=1073.)
    sky.update_positions(time, array_location)
    coh_loc = sky.coherency_calc(array_location)

    inds = np.ones(Ncomp).astype(bool)
    inds[Ncomp // 2:] = False

    unpol_srcs_up = (sky.alt_az[0] > 0) * (~inds)

    assert np.allclose(coh_loc[0, 0, unpol_srcs_up], 0.5)
    assert np.allclose(coh_loc[1, 1, unpol_srcs_up], 0.5)
    assert np.allclose(coh_loc[1, 0, unpol_srcs_up], 0.0)
    assert np.allclose(coh_loc[0, 1, unpol_srcs_up], 0.0)


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


def test_read_hdf5():
    Nside = 32
    Npix = hp.nside2npix(Nside)
    vec = hp.ang2vec(np.pi / 2, np.pi * 3 / 4)
    ipix_disc = hp.query_disc(nside=32, vec=vec, radius=np.radians(10))
    m = np.arange(Npix)
    m[ipix_disc] = m.max()

    indices = np.arange(Npix)

    frequencies = 100000000*np.ones(len(indices))

    hpmap, inds, freqs = pyuvsim.source.read_hdf5(os.path.join(SIM_DATA_PATH,'test_file.hdf5'), 0)

    nt.assert_true(np.allclose(hpmap, m))
    nt.assert_true(np.allclose(inds, indices))
    nt.assert_true(np.allclose(freqs, frequencies))
    
    
def test_healpix_to_sky():
    Nside = 32
    Npix = hp.nside2npix(Nside)
    vec = hp.ang2vec(np.pi / 2, np.pi * 3 / 4)
    ipix_disc = hp.query_disc(nside=32, vec=vec, radius=np.radians(10))
    m = np.arange(Npix)
    m[ipix_disc] = m.max()
    
    hpmap, inds, freqs = pyuvsim.source.read_hdf5(os.path.join(SIM_DATA_PATH,'test_file.hdf5'), 0)
    sky = pyuvsim.source.healpix_to_sky(hpmap, inds, freqs)
    
    nt.assert_true(np.allclose(sky.stokes[0], m))
