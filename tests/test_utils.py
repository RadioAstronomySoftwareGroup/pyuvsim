# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import sys
import warnings

import numpy as np
import pytest
from pyuvdata import UVData
from pyuvdata.testing import check_warnings

from pyuvsim import utils as simutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

triangle_uvfits_file = os.path.join(SIM_DATA_PATH, "28m_triangle_10time_10chan.uvfits")


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
    np.testing.assert_allclose(calc_za, np.deg2rad(zas))
    np.testing.assert_allclose(calc_az, np.deg2rad(beam_azs))


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
    np.testing.assert_allclose(calc_alt, np.deg2rad(alts))
    np.testing.assert_allclose(calc_az, np.deg2rad(astropy_azs))


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
    with pytest.raises(
        ValueError, match="number of altitude and azimuth values must match."
    ):
        simutils.altaz_to_zenithangle_azimuth(0, [0, np.pi / 2])
    with pytest.raises(
        ValueError, match="number of zenith_angle and azimuth values must match."
    ):
        simutils.zenithangle_azimuth_to_altaz(0, [0, np.pi / 2])


@pytest.mark.parametrize("ext", [".ext", ".uvfits", ".uvh5", ".yaml", ""])
def test_file_namer(tmpdir, ext):
    """
    File name incrementer utility, with extensions.
    """
    fnames = []
    for i in range(111):
        fname = str(tmpdir.join(f"file_{i}{ext}"))
        with open(fname, "w") as f:
            f.write(" ")
        fnames.append(fname)
    existing_file = fnames[0]
    new_filepath = simutils.check_file_exists_and_increment(existing_file)
    assert new_filepath.endswith(f"_111{ext}")


@pytest.mark.parametrize("save_format", [None, "uvfits", "miriad", "uvh5", "ms"])
def test_write_uvdata(save_format, tmpdir):
    """Test function that defines filenames from parameter dict"""
    if save_format == "ms":
        pytest.importorskip("casacore")
    if save_format == "miriad" and sys.platform.startswith("win"):
        pytest.skip("miriad is not supported on Windows")

    uv = UVData.from_file(triangle_uvfits_file)

    ofname = str(tmpdir.join("test_file"))
    filing_dict = {"outfile_name": ofname}
    if save_format is None:
        warn_str = [
            "No out format specified for uvdata file. Defaulting to uvh5 (note "
            "this is a defaulting change, it used to default to uvfits)"
        ]
        warn_type = [UserWarning]
    elif save_format == "miriad":
        warn_str = [
            "writing default values for restfreq, vsource, veldop, jyperk, and systemp"
        ]
        warn_type = [UserWarning]
    else:
        warn_type = None
        warn_str = ""
    with check_warnings(warn_type, match=warn_str):
        expected_ofname = simutils.write_uvdata(
            uv, filing_dict, return_filename=True, out_format=save_format
        )

    ofname = os.path.join(".", ofname)

    if save_format == "uvfits":
        assert ofname + ".uvfits" == expected_ofname
    elif save_format == "uvh5" or save_format is None:
        assert ofname + ".uvh5" == expected_ofname
    elif save_format == "ms":
        assert ofname + ".ms" == expected_ofname
    else:
        assert ofname == expected_ofname


@pytest.mark.filterwarnings("ignore:writing default values for restfreq, vsource")
@pytest.mark.parametrize("save_format", [None, "uvfits", "miriad", "uvh5", "ms"])
def test_write_uvdata_clobber(save_format, tmpdir):
    """Test overwriting a uvdata object yields the expected results."""
    if save_format == "ms":
        pytest.importorskip("casacore")

    if save_format == "miriad" and sys.platform.startswith("win"):
        pytest.skip("miriad is not supported on windows")

    uv = UVData.from_file(triangle_uvfits_file)

    uv.set_lsts_from_time_array()
    filing_dict = {
        "outdir": tmpdir.join("test_dir"),
        "outfile_prefix": "test",
        "outfile_suffix": "file",
    }
    if save_format is None:
        warn_str = [
            "No out format specified for uvdata file. Defaulting to uvh5 (note "
            "this is a defaulting change, it used to default to uvfits)"
        ]
        warn_type = [UserWarning]
    elif save_format == "miriad":
        warn_str = [
            "writing default values for restfreq, vsource, veldop, jyperk, and systemp"
        ]
        warn_type = [UserWarning]
    else:
        warn_type = None
        warn_str = ""

    with check_warnings(warn_type, match=warn_str), warnings.catch_warnings():
        warnings.filterwarnings("ignore", "`np.int` is a deprecated alias")
        warnings.filterwarnings("ignore", "`np.bool` is a deprecated alias")
        expected_ofname = simutils.write_uvdata(
            uv, filing_dict, return_filename=True, out_format=save_format
        )

    assert os.path.exists(expected_ofname)

    uv2 = UVData.from_file(expected_ofname)

    if save_format == "ms":
        # MS adds some stuff to history & extra keywords
        uv2.history = uv.history
        uv2.extra_keywords = uv.extra_keywords

        # when we write to MS we set some things that weren't set before.
        # We intend to set these on UVData objects generically in the future...
        uv.dut1 = uv2.dut1
        uv.earth_omega = uv2.earth_omega
        uv.gst0 = uv2.gst0
        uv.rdate = uv2.rdate
        uv.timesys = uv2.timesys

        # for some reason, the vis_units also change. This is more problematic...
        uv2.vis_units = uv.vis_units
    elif save_format == "miriad":
        # ordering gets changes in `write_miriad`
        uv.reorder_blts()
        uv2.reorder_blts()
    elif save_format == "uvfits":
        uv.rdate = uv2.rdate

    uv2._consolidate_phase_center_catalogs(other=uv, ignore_name=True)

    assert uv == uv2

    uv.data_array += 1

    filing_dict["clobber"] = True
    with check_warnings(warn_type, match=warn_str):
        simutils.write_uvdata(uv, filing_dict, out_format=save_format)

    uv2.read(expected_ofname)

    if save_format == "ms":
        # MS adds some stuff to history & extra keywords
        uv2.history = uv.history
        uv2.extra_keywords = uv.extra_keywords
        # for some reason, the vis_units also change. This is more problematic...
        uv2.vis_units = uv.vis_units

    uv2._consolidate_phase_center_catalogs(other=uv, ignore_name=True)
    assert uv2 == uv


def test_write_fix_autos(tmpdir):
    uv = UVData.from_file(triangle_uvfits_file)

    uv.set_lsts_from_time_array()

    auto_screen = uv.ant_1_array == uv.ant_2_array
    assert np.all(np.abs(uv.data_array[auto_screen]) == 0)
    uv.data_array[auto_screen] = 1
    uv.data_array[auto_screen] += 1e-11 * complex(0, 1)

    ofname = str(tmpdir.join("test_file"))
    filing_dict = {"outfile_name": ofname}

    with check_warnings(
        UserWarning, match="Fixing auto-correlations to be be real-only"
    ):
        simutils.write_uvdata(uv, filing_dict, return_filename=True, out_format="uvh5")


def test_write_error_with_no_format(tmpdir):
    """Test write_uvdata will error if no format is given."""
    uv = UVData.from_file(triangle_uvfits_file)

    ofname = str(tmpdir.join("test_file"))
    filing_dict = {"outfile_name": ofname}
    with pytest.raises(ValueError, match="Invalid output format. Options are"):
        simutils.write_uvdata(uv, filing_dict, return_filename=True, out_format="")


def test_file_format_in_filing_dict(tmpdir):
    """Test file is written out when output_format is set in filing dict."""
    uv = UVData.from_file(triangle_uvfits_file)

    ofname = str(tmpdir.join("test_file"))
    filing_dict = {"outfile_name": ofname}
    filing_dict["output_format"] = "uvfits"
    expected_ofname = simutils.write_uvdata(uv, filing_dict, return_filename=True)
    assert ofname + ".uvfits" == expected_ofname

    # Cleanup
    os.remove(ofname + ".uvfits")


def test_progsteps_error():
    with pytest.raises(ValueError, match="Maximum value is needed."):
        simutils.progsteps()
