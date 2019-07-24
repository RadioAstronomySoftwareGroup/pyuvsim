# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import shutil

import numpy as np
import pytest
import pyuvdata.tests as uvtest
from astropy.coordinates import Angle
from astropy.time import Time
from pyuvdata import UVData

import pyuvsim.tests as simtest
from pyuvsim import utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

triangle_uvfits_file = os.path.join(SIM_DATA_PATH, '28m_triangle_10time_10chan.uvfits')


def test_tee_ra_loop():
    time = Time(2457458.1739, scale='utc', format='jd')
    tee_ra = Angle(np.pi / 4., unit='rad')  # rad
    cirs_ra = simutils.tee_to_cirs_ra(tee_ra, time)
    new_tee_ra = simutils.cirs_to_tee_ra(cirs_ra, time)
    assert new_tee_ra == tee_ra


def test_altaz_to_za_az():
    alts = [90, 75, 30, 0, -10, -45, -90]
    # 0=North, 90=East
    azs = [0, 45, 90, 135, 180, 270, 350]

    zas = [0, 15, 60, 90, 100, 135, 180]
    # 0=East, 90=North
    beam_azs = [90, 45, 0, 315, 270, 180, 100]

    calc_za, calc_az = simutils.altaz_to_zenithangle_azimuth(
        np.deg2rad(alts), np.deg2rad(azs)
    )
    assert np.allclose(calc_za, np.deg2rad(zas))
    assert np.allclose(calc_az, np.deg2rad(beam_azs))


def test_single_altaz_to_za_az():
    alts = 90
    # 0=North, 90=East
    azs = 135

    zas = 0
    # 0=East, 90=North
    beam_azs = 315

    calc_za, calc_az = simutils.altaz_to_zenithangle_azimuth(
        np.deg2rad(alts), np.deg2rad(azs)
    )
    assert np.isclose(calc_za, np.deg2rad(zas))
    assert np.isclose(calc_az, np.deg2rad(beam_azs))


def test_za_az_to_altaz():
    # 5 degrees off zenith in the East direction
    zas = [0, 5, 45, 90, 120, 150, 180]
    # 0=East, 90=North
    azs = [0, 45, 90, 135, 180, 270, 350]

    alts = [90, 85, 45, 0, -30, -60, -90]
    # 0=North, 90=East
    astropy_azs = [90, 45, 0, 315, 270, 180, 100]

    calc_alt, calc_az = simutils.zenithangle_azimuth_to_altaz(
        np.deg2rad(zas), np.deg2rad(azs)
    )
    assert np.allclose(calc_alt, np.deg2rad(alts))
    assert np.allclose(calc_az, np.deg2rad(astropy_azs))


def test_single_za_az_to_altaz():
    # 5 degrees off zenith in the East direction
    zas = 5
    # 0=East, 90=North
    azs = 180

    alts = 85
    # 0=North, 90=East
    astropy_azs = 270

    calc_alt, calc_az = simutils.zenithangle_azimuth_to_altaz(
        np.deg2rad(zas), np.deg2rad(azs)
    )
    assert np.isclose(calc_alt, np.deg2rad(alts))
    assert np.isclose(calc_az, np.deg2rad(astropy_azs))


def test_altaz_za_az_errors():
    simtest.assert_raises_message(
        ValueError, 'number of altitude and azimuth values must match.',
        simutils.altaz_to_zenithangle_azimuth, 0, [0, np.pi / 2]
    )
    simtest.assert_raises_message(
        ValueError, 'number of zenith_angle and azimuth values must match.',
        simutils.zenithangle_azimuth_to_altaz, 0, [0, np.pi / 2]
    )


def test_file_namer():
    """
    File name incrementer utility
    """
    fnames = []
    for i in range(111):
        fname = os.path.join(simtest.TESTDATA_PATH, 'file_' + str(i))
        with open(fname, 'w') as f:
            f.write(' ')
        fnames.append(fname)
    existing_file = fnames[0]
    new_filepath = simutils.check_file_exists_and_increment(existing_file)
    for fn in fnames:
        os.remove(fn)
    assert new_filepath.endswith("_111")


def test_file_namer_extensions():
    """
    File name incrementer with specified extension
    """
    os.mkdir(os.path.join(simtest.TESTDATA_PATH, 'tempfiles'))
    fnames = []
    for i in range(111):
        fname = os.path.join(simtest.TESTDATA_PATH, 'tempfiles', 'file_' + str(i) + '.ext')
        with open(fname, 'w') as f:
            f.write(' ')
        fnames.append(fname)
    existing_file = fnames[0]
    new_filepath = simutils.check_file_exists_and_increment(existing_file, 'ext')
    shutil.rmtree(os.path.join(simtest.TESTDATA_PATH, 'tempfiles'))
    assert new_filepath.endswith("_111.ext")


@pytest.mark.parametrize("save_format", [None, 'uvfits', 'miriad', 'uvh5'])
def test_write_uvdata(save_format):
    """ Test function that defines filenames from parameter dict """
    if save_format == 'uvh5':
        pytest.importorskip('h5py')

    uv = UVData()
    uvtest.checkWarnings(
        uv.read_uvfits, [triangle_uvfits_file],
        message='Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.'
    )

    ofname = os.path.join(simtest.TESTDATA_PATH, 'test_file')
    filing_dict = {'outfile_name': ofname}
    expected_ofname = simutils.write_uvdata(uv, filing_dict,
                                            return_filename=True,
                                            out_format=save_format)
    ofname = os.path.join('.', ofname)

    if save_format == 'uvfits' or save_format is None:
        assert ofname + '.uvfits' == expected_ofname
        os.remove(ofname + '.uvfits')
    elif save_format == 'uvh5':
        assert ofname + '.uvh5' == expected_ofname
        os.remove(ofname + '.uvh5')
    else:
        assert ofname == expected_ofname
        shutil.rmtree(ofname)


def test_write_error_with_no_format():
    """Test write_uvdata will error if no format is given."""
    uv = UVData()
    uvtest.checkWarnings(
        uv.read_uvfits, [triangle_uvfits_file],
        message='Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.'
    )

    ofname = os.path.join(simtest.TESTDATA_PATH, 'test_file')
    filing_dict = {'outfile_name': ofname}
    simtest.assert_raises_message(
        ValueError, 'Invalid output format. Options are " uvfits", "uvh5", or "miriad"',
        simutils.write_uvdata, uv, filing_dict, return_filename=True, out_format=''
    )


def test_file_format_in_filing_dict():
    """Test file is written out when output_format is set in filing dict."""
    uv = UVData()
    uvtest.checkWarnings(
        uv.read_uvfits, [triangle_uvfits_file],
        message='Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.'
    )

    ofname = os.path.join(simtest.TESTDATA_PATH, 'test_file')
    filing_dict = {'outfile_name': ofname}
    filing_dict['output_format'] = 'uvfits'
    expected_ofname = simutils.write_uvdata(uv, filing_dict, return_filename=True)
    assert ofname + '.uvfits' == expected_ofname

    # Cleanup
    os.remove(ofname + '.uvfits')
