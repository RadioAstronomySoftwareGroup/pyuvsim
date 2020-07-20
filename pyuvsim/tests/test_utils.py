# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os

import numpy as np
import pytest
from pyuvdata import UVData

from pyuvsim import utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

triangle_uvfits_file = os.path.join(SIM_DATA_PATH, '28m_triangle_10time_10chan.uvfits')


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
    with pytest.raises(ValueError, match='number of altitude and azimuth values must match.'):
        simutils.altaz_to_zenithangle_azimuth(0, [0, np.pi / 2])
    with pytest.raises(ValueError, match='number of zenith_angle and azimuth values must match.'):
        simutils.zenithangle_azimuth_to_altaz(0, [0, np.pi / 2])


@pytest.mark.parametrize('ext', ['.ext', '.uvfits', '.uvh5', '.yaml', ''])
def test_file_namer(tmpdir, ext):
    """
    File name incrementer utility, with extensions.
    """
    fnames = []
    for i in range(111):
        fname = str(tmpdir.join(f"file_{i}{ext}"))
        with open(fname, 'w') as f:
            f.write(' ')
        fnames.append(fname)
    existing_file = fnames[0]
    if ext == '.ext':
        new_filepath = simutils.check_file_exists_and_increment(existing_file, 'ext')
    else:
        new_filepath = simutils.check_file_exists_and_increment(existing_file)
    assert new_filepath.endswith(f"_111{ext}")


@pytest.mark.parametrize("save_format", [None, 'uvfits', 'miriad', 'uvh5'])
def test_write_uvdata(save_format, tmpdir):
    """ Test function that defines filenames from parameter dict """
    uv = UVData()
    uv.read_uvfits(triangle_uvfits_file)

    ofname = str(tmpdir.join('test_file'))
    filing_dict = {'outfile_name': ofname}
    expected_ofname = simutils.write_uvdata(uv, filing_dict,
                                            return_filename=True,
                                            out_format=save_format)
    ofname = os.path.join('.', ofname)

    if save_format == 'uvfits' or save_format is None:
        assert ofname + '.uvfits' == expected_ofname
    elif save_format == 'uvh5':
        assert ofname + '.uvh5' == expected_ofname
    else:
        assert ofname == expected_ofname


def test_write_error_with_no_format(tmpdir):
    """Test write_uvdata will error if no format is given."""
    uv = UVData()
    uv.read_uvfits(triangle_uvfits_file)

    ofname = str(tmpdir.join('test_file'))
    filing_dict = {'outfile_name': ofname}
    with pytest.raises(ValueError,
                       match='Invalid output format. Options are'):
        simutils.write_uvdata(uv, filing_dict, return_filename=True, out_format='')


def test_file_format_in_filing_dict(tmpdir):
    """Test file is written out when output_format is set in filing dict."""
    uv = UVData()
    uv.read_uvfits(triangle_uvfits_file)

    ofname = str(tmpdir.join('test_file'))
    filing_dict = {'outfile_name': ofname}
    filing_dict['output_format'] = 'uvfits'
    expected_ofname = simutils.write_uvdata(uv, filing_dict, return_filename=True)
    assert ofname + '.uvfits' == expected_ofname

    # Cleanup
    os.remove(ofname + '.uvfits')
