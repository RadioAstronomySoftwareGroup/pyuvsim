# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import ast
import copy
import os
import re
import shutil
import subprocess  # nosec
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy import units
from astropy.coordinates import Angle, EarthLocation, Latitude, Longitude, SkyCoord
from astropy.time import Time
from astropy.utils.data import download_file, get_cached_urls
from pyradiosky import SkyModel
from pyradiosky.data import DATA_PATH as SKY_DATA_PATH
from pyuvdata import (
    AiryBeam,
    GaussianBeam,
    ShortDipoleBeam,
    UniformBeam,
    UVBeam,
    UVData,
    utils as uvutils,
)
from pyuvdata.analytic_beam import AnalyticBeam
from pyuvdata.testing import check_warnings

try:
    from pyuvsim import mpi
except ImportError:
    mpi = None

import pyuvsim
import pyuvsim.utils as simutils
from pyuvsim import simsetup
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

from . import compare_dictionaries

herabeam_default = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")

pytestmark = pytest.mark.filterwarnings(
    "ignore:antenna_diameters are not set or are being overwritten."
)

# Five different test configs
param_filenames = [
    os.path.join(SIM_DATA_PATH, "test_config", f"param_10time_10chan_{x}.yaml")
    for x in range(6)
]

longbl_uvfits_file = os.path.join(SIM_DATA_PATH, "5km_triangle_1time_1chan.uvfits")
triangle_uvfits_file = os.path.join(SIM_DATA_PATH, "28m_triangle_10time_10chan.uvfits")
manytimes_config = os.path.join(
    SIM_DATA_PATH, "test_config", "param_100times_1.5days_triangle.yaml"
)
gleam_param_file = os.path.join(
    SIM_DATA_PATH, "test_config", "param_1time_testgleam.yaml"
)


@pytest.fixture(scope="module")
def times_and_freqs():
    freqs, fstep = np.linspace(100, 200, 1024, retstep=True)
    times, tstep = np.linspace(2458570, 2458570 + 0.5, 239, retstep=True)
    # yield the time and frequency arrays to the tests
    # then delete after
    yield times, tstep, freqs, fstep

    del (times, freqs)


@pytest.fixture
def time_dict_base():
    tdict_base = {
        "Ntimes": 25,
        "duration_hours": 25.0,
        "start_time": 2457458.0,
        "end_time": 2457459.0,
        "integration_time": 3600,
    }

    time_array = np.linspace(
        tdict_base["start_time"],
        tdict_base["end_time"],
        tdict_base["Ntimes"],
        endpoint=True,
    )

    tdict_base["time_array"] = time_array

    return tdict_base


@pytest.fixture(scope="module")
def cat_with_some_pols():
    # Mock catalog with a couple sources polarized.
    Nsrcs = 30
    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    pol_inds = range(12, 15)
    stokes = np.zeros((4, Nfreqs, Nsrcs))
    stokes[0, :, :] = 1.0

    stokes[1, :, pol_inds] = 0.2
    stokes[2, :, pol_inds] = 1.2
    stokes[3, :, pol_inds] = 0.3

    ra = Longitude(np.linspace(0, 2 * np.pi, Nsrcs), "rad")
    dec = Latitude(np.linspace(-np.pi / 2, np.pi / 3, Nsrcs), "rad")

    sky = SkyModel(
        name=np.arange(Nsrcs).astype(str),
        ra=ra,
        dec=dec,
        frame="icrs",
        stokes=stokes * units.Jy,
        spectral_type="full",
        freq_array=freqs,
    )

    return sky


@pytest.fixture
def uvdata_keyword_dict():
    init_dict = {
        "antenna_layout_filepath": os.path.join(
            SIM_DATA_PATH, "test_config/triangle_bl_layout.csv"
        ),
        "telescope_location": (
            -30.72152777777791,
            21.428305555555557,
            1073.0000000093132,
        ),
        "telescope_name": "HERA",
        "Nfreqs": 10,
        "start_freq": 1e8,
        "bandwidth": 1e8,
        "Ntimes": 60,
        "integration_time": 100.0,
        "start_time": 2458101.0,
        "polarization_array": ["xx"],
        "no_autos": True,
        "conjugation_convention": "ant1<ant2",
        "blt_order": ["time", "baseline"],
        "write_files": False,
        "run_check": True,
        "feed_array": ["x", "y"],
        "feed_angle": [np.pi / 2, 0],
        "mount_type": "fixed",
    }

    return init_dict


@pytest.fixture(scope="session")
def mwa_beam_path():
    try:
        from pyuvdata.datasets import fetch_data

        path = fetch_data("mwa_full_EE")
    except ImportError:
        # This can be removed once we require pyuvdata > 3.2.3
        from pyuvdata.data import DATA_PATH as UV_DATA_PATH

        path = os.path.join(UV_DATA_PATH, "mwa_full_EE_test.h5")

    yield path


@pytest.mark.parametrize("return_data", [True, False])
def test_mock_catalog_zenith_source(hera_loc, return_data):
    time = Time(2457458.65410, scale="utc", format="jd")

    array_location = hera_loc

    source_coord = SkyCoord(
        alt=Angle(90 * units.deg),
        az=Angle(0 * units.deg),
        obstime=time,
        frame="altaz",
        location=array_location,
    )
    icrs_coord = source_coord.transform_to("icrs")

    test_source = SkyModel(
        name="src0",
        skycoord=icrs_coord,
        stokes=units.Quantity([1, 0, 0, 0], "Jy"),
        spectral_type="flat",
    )

    cat, _ = simsetup.create_mock_catalog(
        time, arrangement="zenith", return_data=return_data
    )

    if return_data:
        test_data = simsetup.SkyModelData(test_source)

        for attr in cat.__dict__:
            assert getattr(cat, attr) == getattr(test_data, attr)

        new_cat = cat.get_skymodel()

        assert new_cat == test_source

    else:
        assert cat == test_source


def test_shared_mpierr():
    time = Time(2457458.65410, scale="utc", format="jd")
    cat, _ = simsetup.create_mock_catalog(time, arrangement="zenith")
    cat_data = simsetup.SkyModelData(cat)

    if mpi is None:
        with pytest.raises(ImportError, match="You need mpi4py to use this method."):
            cat_data.share(root=0)


def test_mock_catalog_off_zenith_source(hera_loc):
    src_az = Angle("90.0d")
    src_alt = Angle("85.0d")

    time = Time(2457458.65410, scale="utc", format="jd")

    array_location = hera_loc

    source_coord = SkyCoord(
        alt=src_alt, az=src_az, obstime=time, frame="altaz", location=array_location
    )
    icrs_coord = source_coord.transform_to("icrs")

    test_source = SkyModel(
        name="src0",
        skycoord=icrs_coord,
        stokes=units.Quantity([1.0, 0, 0, 0], "Jy"),
        spectral_type="flat",
    )

    cat, _ = simsetup.create_mock_catalog(
        time, arrangement="off-zenith", alt=src_alt.deg
    )

    assert cat == test_source


def test_mock_diffuse_maps_errors():
    analytic_diffuse = simsetup.analytic_diffuse
    astropy_healpix = simsetup.astropy_healpix
    if (analytic_diffuse is not None) and (astropy_healpix is not None):
        # Error cases:
        with pytest.raises(ValueError, match="Diffuse arrangement selected"):
            simsetup.create_mock_catalog(Time.now(), arrangement="diffuse")

        with check_warnings(UserWarning, match="No nside chosen"):
            simsetup.create_mock_catalog(
                Time.now(), arrangement="diffuse", diffuse_model="monopole"
            )

    else:
        with pytest.raises(ValueError, match="analytic_diffuse and astropy_healpix"):
            simsetup.create_mock_catalog(Time.now(), arrangement="diffuse")


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize(
    ("modname", "modkwargs"),
    [("monopole", {}), ("gauss", {"a": 0.05}), ("polydome", {"n": 4})],
)
@pytest.mark.parametrize("location", ["earth", "moon"])
def test_mock_diffuse_maps(modname, modkwargs, hera_loc, apollo_loc, location):
    analytic_diffuse = pytest.importorskip("analytic_diffuse")
    pytest.importorskip("astropy_healpix")
    if location == "earth":
        loc = hera_loc
    else:
        pytest.importorskip("lunarsky")
        loc = apollo_loc
    map_nside = 128
    t0 = Time.now()
    cat, _ = simsetup.create_mock_catalog(
        t0,
        arrangement="diffuse",
        array_location=loc,
        diffuse_model=modname,
        map_nside=map_nside,
        diffuse_params=modkwargs,
    )

    cat.update_positions(t0, loc)

    modfunc = analytic_diffuse.get_model(modname)
    alt, az = cat.alt_az
    za = np.pi / 2 - alt

    vals = modfunc(az, za, **modkwargs)

    assert cat.nside == map_nside
    np.testing.assert_allclose(cat.stokes[0, 0].to_value("K"), vals, rtol=0, atol=1e-12)


@pytest.mark.parametrize(
    ("horizon_buffer", "pass_time", "pass_array_loc", "pass_uv", "return_catname"),
    [
        (True, True, True, False, True),
        (False, False, False, True, False),
        (True, True, True, True, False),
    ],
)
def test_initialize_catalog_from_params(
    horizon_buffer, pass_time, pass_array_loc, pass_uv, return_catname, hera_loc
):
    # Pass in parameter dictionary as dict
    uv_in = UVData.from_file(triangle_uvfits_file)

    source_dict = {"catalog": "mock", "mock_arrangement": "zenith", "Nsrcs": 5}
    if horizon_buffer:
        source_dict["horizon_buffer"] = 0.04364

    if pass_array_loc:
        source_dict["array_location"] = ",".join(
            [str(coord) for coord in uv_in.telescope.location_lat_lon_alt_degrees]
        )

    warn_type = []
    warn_str = []
    if pass_time:
        source_dict["time"] = uv_in.time_array[0]
    else:
        warn_type += [UserWarning]
        warn_str += [
            "No julian date given for mock catalog. Defaulting to first time step."
        ]

    if pass_uv:
        uv_use = uv_in
    else:
        uv_use = None
        if not pass_array_loc:
            warn_type += [UserWarning]
            warn_str += ["No array_location specified. Defaulting to the HERA site."]

    if len(warn_type) == 0:
        warn_type = None
        warn_str = ""

    with check_warnings(warn_type, match=warn_str):
        catalog_uv = simsetup.initialize_catalog_from_params(
            {"sources": source_dict}, input_uv=uv_use, return_catname=return_catname
        )

    if return_catname:
        catalog_uv, cat_name = catalog_uv
        assert cat_name.startswith("mock")

    exp_cat, _ = simsetup.create_mock_catalog(
        uv_in.time_array[0], arrangement="zenith", array_location=hera_loc, Nsrcs=5
    )

    assert exp_cat == catalog_uv


@pytest.mark.parametrize(
    ("source_dict", "input_uv", "err_type", "err_msg", "warn_msg"),
    [
        ({}, None, KeyError, "No catalog defined.", []),
        (
            {"catalog": "mock", "mock_arrangement": "zenith", "Nsrcs": 5},
            "foo",
            TypeError,
            "input_uv must be UVData object",
            [],
        ),
        (
            {"catalog": "mock", "mock_arrangement": "zenith", "Nsrcs": 5},
            None,
            ValueError,
            "input_uv must be supplied if using mock catalog without specified julian date",
            ["No array_location specified. Defaulting to the HERA site."],
        ),
    ],
)
def test_initialize_catalog_from_params_errors(
    source_dict, input_uv, err_type, err_msg, warn_msg
):
    if len(warn_msg) == 0:
        warn_type = None
    else:
        warn_type = UserWarning
    with (
        check_warnings(warn_type, match=warn_msg),
        pytest.raises(err_type, match=err_msg),
    ):
        simsetup.initialize_catalog_from_params(
            {"sources": source_dict}, input_uv=input_uv
        )


@pytest.mark.parametrize("use_filetype", [True, False])
def test_vot_catalog(use_filetype):
    filetype = None

    vot_param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testvot.yaml"
    )
    if use_filetype:
        filetype = "vot"
    vot_catalog = simsetup.initialize_catalog_from_params(
        vot_param_filename, filetype=filetype
    )

    if use_filetype:
        filetype = "text"
    txt_param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml"
    )
    txt_catalog = simsetup.initialize_catalog_from_params(
        txt_param_filename, filetype=filetype
    )

    assert vot_catalog == txt_catalog


def test_vot_catalog_errors():
    vot_param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testvot.yaml"
    )
    with pytest.raises(ValueError, match="Invalid filetype. Filetype options are:"):
        simsetup.initialize_catalog_from_params(vot_param_filename, filetype="foo")


@pytest.mark.filterwarnings("ignore:Some stokes I values are negative.")
@pytest.mark.parametrize(
    ("filetype", "flux_cut", "nan_cut", "neg_cut"),
    [
        ("gleam", True, True, True),
        (None, True, False, False),
        (None, False, True, True),
        (None, False, False, True),
        (None, False, False, False),
    ],
)
def test_gleam_catalog(filetype, flux_cut, nan_cut, neg_cut):
    params_use = gleam_param_file
    if flux_cut:
        expected_ncomp = 9
    elif neg_cut:
        expected_ncomp = 32
    else:
        expected_ncomp = 50

    if not flux_cut or not nan_cut or not neg_cut:
        with open(gleam_param_file) as pfile:
            param_dict = yaml.safe_load(pfile)
        param_dict["config_path"] = os.path.dirname(gleam_param_file)
        if not flux_cut:
            param_dict["sources"].pop("min_flux")
            param_dict["sources"].pop("max_flux")
        if not nan_cut:
            param_dict["sources"].pop("non_nan")
        if not neg_cut:
            param_dict["sources"].pop("non_negative")
        params_use = param_dict

    gleam_catalog = simsetup.initialize_catalog_from_params(
        params_use, filetype=filetype
    )

    assert gleam_catalog.Ncomponents == expected_ncomp


@pytest.mark.parametrize("use_filetype", [True, False])
@pytest.mark.parametrize("yaml_filetype", [True, False])
def test_skyh5_catalog(use_filetype, yaml_filetype, tmp_path):
    filetype = None
    gleam_filename = os.path.join(SIM_DATA_PATH, "gleam_50srcs.vot")
    skyobj = SkyModel.from_gleam_catalog(gleam_filename, run_check=False)
    skyobj.select(non_negative=True)
    assert skyobj.Ncomponents == 32

    skyh5_file = os.path.join(tmp_path, "gleam.skyh5")
    skyobj.write_skyh5(skyh5_file, clobber=True)

    starting_param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml"
    )

    param_filename = os.path.join(tmp_path, "param_test_skyh5_gleam.yaml")

    with open(starting_param_filename) as yf:
        param_dict = yaml.safe_load(yf)

    param_dict["sources"]["catalog"] = skyh5_file
    if yaml_filetype:
        param_dict["sources"]["filetype"] = "skyh5"
    with open(param_filename, "w") as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)

    if use_filetype:
        filetype = "skyh5"
    skyh5_catalog = simsetup.initialize_catalog_from_params(
        param_filename, filetype=filetype
    )

    assert skyh5_catalog == skyobj


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
def test_healpix_catalog():
    pytest.importorskip("astropy_healpix")
    path = os.path.join(SKY_DATA_PATH, "healpix_disk.skyh5")
    sky = SkyModel.from_file(path)

    params = {"sources": {"catalog": path}}
    hpx_sky = simsetup.initialize_catalog_from_params(params)
    assert hpx_sky == sky


@pytest.mark.parametrize("spectral_type", ["flat", "subband", "spectral_index"])
def test_gleam_catalog_spectral_type(spectral_type):
    with open(gleam_param_file) as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict["config_path"] = os.path.dirname(gleam_param_file)
    param_dict["sources"].pop("min_flux")
    param_dict["sources"].pop("max_flux")
    param_dict["sources"]["spectral_type"] = spectral_type

    gleam_catalog = simsetup.initialize_catalog_from_params(param_dict)
    assert gleam_catalog.spectral_type == spectral_type
    if spectral_type == "flat":
        assert gleam_catalog.Ncomponents == 50
    else:
        assert gleam_catalog.Ncomponents == 32


@pytest.mark.parametrize("telparam_in_obsparam", [True, False])
def test_param_reader(telparam_in_obsparam, tmpdir):
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_10time_10chan_0.yaml"
    )
    with open(param_filename) as fhandle:
        param_dict = yaml.safe_load(fhandle)

    if telparam_in_obsparam:
        new_param_file = os.path.join(tmpdir, "new_obsparam.yaml")
        new_telconfig_file = os.path.join(tmpdir, "new_tel_config.yaml")
        new_layout_file = os.path.join(tmpdir, "triangle_bl_layout.csv")

        orig_tel_config = os.path.join(
            SIM_DATA_PATH,
            "test_config",
            param_dict["telescope"]["telescope_config_name"],
        )
        # add an instrument here to make sure it gets picked up properly
        param_dict["telescope"]["instrument"] = "foo"
        param_dict["telescope"]["telescope_config_name"] = new_telconfig_file

        with open(orig_tel_config) as fhandle:
            tel_config = yaml.safe_load(fhandle)
            tel_config["beam_paths"][0].filename = os.path.join(
                SIM_DATA_PATH, tel_config["beam_paths"][0].filename[0]
            )
            # add a mount type
            for beam in tel_config["beam_paths"].values():
                if issubclass(type(beam), AnalyticBeam):
                    beam.mount_type = "fixed"
            assert tel_config["beam_paths"][1].mount_type is not None
        with open(new_telconfig_file, "w") as fhandle:
            yaml.safe_dump(tel_config, fhandle)
        with open(new_param_file, "w") as fhandle:
            yaml.safe_dump(param_dict, fhandle)

        shutil.copyfile(
            os.path.join(SIM_DATA_PATH, "test_config", "triangle_bl_layout.csv"),
            new_layout_file,
        )

    else:
        new_param_file = param_filename

    uv_in = UVData.from_file(triangle_uvfits_file)
    uv_in.unproject_phase()

    beam0 = UVBeam.from_file(herabeam_default, mount_type="fixed")
    beam1 = UniformBeam()
    beam2 = GaussianBeam(sigma=0.02)
    beam3 = AiryBeam(diameter=14.6)

    beam_list = pyuvsim.BeamList([beam0, beam1, beam2, beam3])

    beam_dict = {"ANT1": 0, "ANT2": 1, "ANT3": 2, "ANT4": 3}

    warn_list = [
        "The reorder_blt_kw parameter is deprecated in favor of setting "
        "obs_param['ordering']['blt_order']. This will become an error in "
        "version 1.5"
    ]

    # Check default configuration
    with check_warnings(DeprecationWarning, match=warn_list):
        uv_obj, new_beam_list, new_beam_dict = simsetup.initialize_uvdata_from_params(
            new_param_file,
            reorder_blt_kw={"order": "time", "minor_order": "baseline"},
            check_kw={"run_check_acceptability": True},
            return_beams=True,
        )

    for bi in new_beam_list:
        assert bi.beam.mount_type is not None

    assert np.all(uv_obj.telescope.feed_angle[:, 0] == np.pi / 2)
    assert np.all(uv_obj.telescope.feed_angle[:, 1] == 0)

    assert uv_obj.telescope.mount_type is not None
    assert uv_obj.telescope.get_x_orientation_from_feeds() == "east"

    if telparam_in_obsparam:
        assert uv_obj.telescope.instrument == "foo"
        # reset it for equality check
        uv_obj.telescope.instrument = "Triangle"

    simsetup._complete_uvdata(uv_obj, inplace=True)

    with check_warnings(
        UserWarning,
        match="No out format specified for uvdata file. Defaulting to uvh5 "
        "(note this is a defaulting change, it used to default to uvfits).",
    ):
        expected_ofilepath = simutils.write_uvdata(
            uv_obj, param_dict, return_filename=True, dryrun=True
        )
    assert Path(expected_ofilepath).name == "sim_results.uvh5"

    # Spoof attributes that won't match.
    uv_obj.history = uv_in.history

    uvfits_required_extra = ["_gst0", "_rdate", "_earth_omega", "_dut1", "_timesys"]
    for attr in uvfits_required_extra:
        param = getattr(uv_obj, attr)
        if param.value is None:
            param.value = param.spoof_val
            setattr(uv_obj, attr, param)

    assert new_beam_dict == beam_dict

    for b_ind, bi in enumerate(new_beam_list):
        assert bi.beam == beam_list[b_ind].beam
        assert bi == beam_list[b_ind]
    assert new_beam_list == beam_list

    # renumber/rename the phase centers so the equality check will pass.
    uv_obj._consolidate_phase_center_catalogs(other=uv_in, ignore_name=True)

    # file is missing extra_keywords info
    assert uv_obj._extra_keywords != uv_in._extra_keywords
    if telparam_in_obsparam:
        uv_in.extra_keywords = {
            "obsparam": "new_obsparam.yaml",
            "telecfg": new_telconfig_file,
            "layout": "triangle_bl_layout.csv",
        }
    else:
        uv_in.extra_keywords = {
            "obsparam": "param_10time_10chan_0.yaml",
            "telecfg": "28m_triangle_10time_10chan.yaml",
            "layout": "triangle_bl_layout.csv",
        }

    assert uv_obj.telescope.mount_type is not None
    uv_in.telescope.mount_type = uv_obj.telescope.mount_type
    assert uv_obj == uv_in


@pytest.mark.filterwarnings("ignore:Entries in 'beam_paths' should be specified")
@pytest.mark.parametrize(
    ("subdict", "error", "msg"),
    [
        (
            {
                "config_path": os.path.join(
                    SIM_DATA_PATH, "nonexistent_directory", "nonexistent_file"
                )
            },
            ValueError,
            "nonexistent_directory is not a directory",
        ),
        (
            {
                "config_path": os.path.join(SIM_DATA_PATH, "test_config"),
                "telescope": {"array_layout": "nonexistent_file"},
            },
            ValueError,
            "nonexistent_file from yaml does not exist",
        ),
        (
            {
                "config_path": os.path.join(SIM_DATA_PATH, "test_config"),
                "telescope": {"telescope_config_name": "nonexistent_file"},
            },
            ValueError,
            "telescope_config_name file from yaml does not exist",
        ),
        (
            {
                "telescope": {
                    "telescope_config_name": os.path.join(
                        SIM_DATA_PATH,
                        "test_config",
                        "28m_triangle_10time_10chan_nofile.yaml",
                    )
                }
            },
            ValueError,
            "Unrecognized beam file or model",
        ),
        (
            {"polarization_array": [-1]},
            ValueError,
            re.escape(
                "Specified polarization array [-1] requires feeds {'r'} but the "
                "beams have feeds ['x' 'y']."
            ),
        ),
    ],
)
def test_param_reader_errors(subdict, error, msg):
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_10time_10chan_0.yaml"
    )

    # Error conditions:
    params_bad = simsetup._config_str_to_dict(param_filename)

    for key, value in subdict.items():
        if isinstance(value, dict):
            for key2, value2 in value.items():
                params_bad[key][key2] = value2
        else:
            params_bad[key] = value

    with pytest.raises(error, match=msg):
        simsetup.initialize_uvdata_from_params(params_bad, return_beams=True)


@pytest.mark.parametrize(
    ("world", "selenoid"),
    [
        (None, None),
        ("earth", None),
        ("moon", "SPHERE"),
        ("moon", "GSFC"),
        ("moon", "GRAIL23"),
        ("moon", "CE-1-LAM-GEO"),
        ("moon", None),
    ],
)
def test_tele_parser(world, selenoid):
    """
    Test minimal dict passed (not covered by param reader tests)
    """
    tdict = {
        "array_layout": os.path.join(SIM_DATA_PATH, "test_layout_6ant.csv"),
        "telescope_location": "(-30.72152777777791, 21.428305555555557, 1073.0000000093132)",
        "telescope_name": "foo",
    }
    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    if world is not None:
        if world != "earth":
            pytest.importorskip("lunarsky")
        tdict["world"] = world
        if selenoid is not None:
            tdict["ellipsoid"] = selenoid

    warn_type = UserWarning
    warn_str = [
        "No beam information, so cannot determine telescope mount_type, "
        "feed_array or feed_angle. Specify a telescope config file in the "
        "obs param file to get beam information."
    ]

    with check_warnings(warn_type, match=warn_str):
        tpars, blist, _ = simsetup.parse_telescope_params(tdict, freq_array=freqs)

    assert tpars["Nants_data"] == 6
    assert len(blist) == 0
    if world is not None:
        assert tpars["world"] == world
        if selenoid is not None:
            assert tpars["ellipsoid"] == selenoid
        elif world != "earth":
            assert tpars["ellipsoid"] == "SPHERE"
        else:
            assert tpars["ellipsoid"] is None


@pytest.mark.parametrize(
    ("tele_dict", "err_type", "err_msg"),
    [
        (
            {"array_layout": os.path.join(SIM_DATA_PATH, "test_layout_6ant.csv")},
            KeyError,
            "If telescope_config_name not provided in `telescope` obsparam section, "
            "you must provide telescope_location",
        ),
        (
            {
                "array_layout": os.path.join(SIM_DATA_PATH, "test_layout_6ant.csv"),
                "telescope_location": (
                    "(-30.72152777777791, 21.428305555555557, 1073.0000000093132)"
                ),
            },
            KeyError,
            "If telescope_config_name not provided in `telescope` obsparam section, "
            "you must provide telescope_name",
        ),
        (
            {
                "telescope_name": "foo",
                "telescope_location": (
                    "(-30.72152777777791, 21.428305555555557, 1073.0000000093132)"
                ),
            },
            KeyError,
            "array_layout must be provided.",
        ),
        (
            {
                "array_layout": 5,
                "telescope_name": "foo",
                "telescope_location": (
                    "(-30.72152777777791, 21.428305555555557, 1073.0000000093132)"
                ),
            },
            ValueError,
            "array_layout must be a string or have options that parse as a dict.",
        ),
        (
            {
                "array_layout": os.path.join(SIM_DATA_PATH, "test_layout_6ant.csv"),
                "telescope_name": "foo",
                "telescope_location": (
                    "(-30.72152777777791, 21.428305555555557, 1073.0000000093132)"
                ),
                "world": "tatooine",
            },
            ValueError,
            "Invalid world tatooine",
        ),
    ],
)
def test_tele_parser_errors(tele_dict, err_type, err_msg):
    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    with pytest.raises(err_type, match=err_msg):
        simsetup.parse_telescope_params(tele_dict, freq_array=freqs)


@pytest.mark.parametrize(
    "bpass_kwds",
    [("start_freq", "end_freq"), ("channel_width", "Nfreqs"), ("bandwidth",)],
)
@pytest.mark.parametrize("chwid_kwds", [("Nfreqs",), ("channel_width",)])
@pytest.mark.parametrize("ref_freq_kwds", [("start_freq",), ("end_freq",)])
def test_freq_parser(bpass_kwds, chwid_kwds, ref_freq_kwds):
    """
    Check all valid input parameter cases for frequencies.
    """

    fdict_base = {
        "Nfreqs": 10,
        "channel_width": 0.5,
        "start_freq": 0.0,
        "end_freq": 4.5,
        "bandwidth": 5.0,
    }

    freq_array = np.linspace(
        fdict_base["start_freq"],
        fdict_base["start_freq"]
        + fdict_base["bandwidth"]
        - fdict_base["channel_width"],
        fdict_base["Nfreqs"],
        endpoint=True,
    )

    fdict_base["freq_array"] = freq_array

    # As long as one tuple from each set is represented,
    # the param parser should work.
    keys = tuple(set(bpass_kwds + chwid_kwds + (ref_freq_kwds)))  # Get unique keys
    subdict = {key: fdict_base[key] for key in keys}
    test = simsetup.parse_frequency_params(subdict)
    np.testing.assert_allclose(test["freq_array"], freq_array)

    cw_array = np.full((10,), fdict_base["channel_width"], dtype=float)
    np.testing.assert_allclose(test["channel_width"], cw_array)


def test_freq_parse_linspace_bug():
    fdict = {
        "Nfreqs": 500,
        "end_freq": 150000000.0,
        "start_freq": 100000000.0,
        "bandwidth": 50000000.0,
    }

    freq_array = np.linspace(100000000, 150000000, num=500)

    with check_warnings(
        UserWarning,
        match="The bandwidth is not consistent with the start_freq, end_freq, "
        "Nfreqs specified. Using the values calculated from the start_freq, "
        "end_freq, Nfreqs parameters.",
    ):
        new_fdict = simsetup.parse_frequency_params(fdict)
    np.testing.assert_allclose(new_fdict["freq_array"], freq_array)


@pytest.mark.parametrize(
    "fdict_add", [{"channel_width": 50000000.0}, {"bandwidth": 50000000.0}, {}]
)
def test_freq_parser_single(fdict_add):
    fdict = {"Nfreqs": 1, "end_freq": 150000000.0, "start_freq": 100000000.0}
    fdict = fdict | fdict_add

    if len(fdict_add) == 0:
        with pytest.raises(
            ValueError,
            match="If Nfreqs is 1 then either channel_width or bandwidth must be specified.",
        ):
            new_fdict = simsetup.parse_frequency_params(fdict)
        return
    with check_warnings(
        UserWarning,
        match="The end_freq is not consistent with the start_freq, Nfreqs, "
        f"{', '.join(fdict_add.keys())} specified.",
    ):
        new_fdict = simsetup.parse_frequency_params(fdict)

    freq_array = np.asarray([100000000.0])
    channel_width = np.asarray([50000000.0])
    np.testing.assert_allclose(new_fdict["freq_array"], freq_array)
    np.testing.assert_allclose(new_fdict["channel_width"], channel_width)


@pytest.mark.parametrize(
    ("freq_array", "channel_width"),
    [
        (np.linspace(0.0, 4.5, 10), 0.5),
        (np.asarray([0.0, 0.5, 2, 4]), 0.5),
        (np.linspace(0.0, 4.5, 10), np.full((10,), 0.5, dtype=float)),
        (np.asarray([0.0, 0.5, 2, 4]), np.full((4,), 0.5, dtype=float)),
        (np.asarray([0.0, 0.5, 2, 4]), np.asarray([0.5, 0.5, 1, 2])),
    ],
)
def test_freq_parser_freq_array(freq_array, channel_width):
    """
    Check parsing works with a vector of channel widths.
    """
    subdict = {"freq_array": freq_array, "channel_width": channel_width}
    test = simsetup.parse_frequency_params(subdict)
    np.testing.assert_allclose(test["freq_array"], freq_array)
    if not isinstance(channel_width, np.ndarray):
        np.testing.assert_allclose(
            test["channel_width"], np.ones_like(freq_array) * 0.5
        )
    else:
        np.testing.assert_allclose(test["channel_width"], channel_width)


@pytest.mark.parametrize(
    ("freq_dict", "msg"),
    [
        (
            {"bandwidth": 5.0},
            "Either start_freq or end_freq must be specified. The parameters "
            "that were specified were: bandwidth",
        ),
        (
            {"start_freq": 0.0, "Nfreqs": 10},
            "If only one of start_freq and end_freq is specified, two more of Nfreqs, "
            "channel_width and bandwidth must be specified as well. The "
            "parameters that were specified were: Nfreqs, start_freq",
        ),
        (
            {"start_freq": 0.0, "channel_width": 0.5},
            "If only one of start_freq and end_freq is specified, two more of "
            "Nfreqs, channel_width and bandwidth must be specified as well. The "
            "parameters that were specified were: channel_width, start_freq",
        ),
        (
            {"start_freq": 0.0, "end_freq": 4.5},
            "If both start_freq and end_freq are specified, either Nfreqs or "
            "channel_width must be specified as well. The parameters "
            "that were specified were: end_freq, start_freq",
        ),
        (
            {"freq_array": [0.0]},
            "channel_width must be specified if freq_array has length 1",
        ),
        (
            {"freq_array": [0.0, 0.5, 2, 4]},
            "channel_width must be specified if spacing in freq_array is uneven.",
        ),
        (
            {"freq_array": [0.0, 0.5, 2, 4], "channel_width": [0.5, 0.5, 0.5]},
            "If channel_width has multiple elements, the channel_width must be "
            "the same length as freq_array",
        ),
        (
            {"channel_width": 3.14, "start_freq": 1.0, "end_freq": 8.3},
            "end_freq - start_freq must be evenly divisible by channel_width",
        ),
        (
            {
                "start_freq": 0.0,
                "Nfreqs": 10,
                "channel_width": np.full((10,), 0.5, dtype=float),
            },
            "channel_width must be a scalar if freq_array is not specified",
        ),
    ],
)
def test_freq_parser_errors(freq_dict, msg):
    # Now check error cases
    with pytest.raises(ValueError, match=msg):
        simsetup.parse_frequency_params(freq_dict)


@pytest.mark.parametrize(
    ("freq_dict", "msg"),
    [
        (
            {
                "freq_array": [100.0, 101.0, 102.0],
                "start_freq": 100.5,
                "end_freq": 101.5,
                "Nfreqs": 5,
                "bandwidth": 10.0,
            },
            "The start_freq, end_freq, Nfreqs, bandwidth is not consistent with "
            "the freq_array specified. Using the values calculated from the "
            "freq_array parameters.",
        ),
        (
            {
                "start_freq": 100,
                "end_freq": 105,
                "Nfreqs": 5,
                "bandwidth": 10.0,
                "channel_width": 0.5,
            },
            "The bandwidth, channel_width is not consistent with the start_freq, "
            "end_freq, Nfreqs specified. Using the values calculated from the "
            "start_freq, end_freq, Nfreqs parameters.",
        ),
        (
            {"start_freq": 100, "Nfreqs": 5, "bandwidth": 10.0, "channel_width": 0.5},
            "The bandwidth is not consistent with the Nfreqs, channel_width, "
            "start_freq specified. Using the values calculated from the Nfreqs, "
            "channel_width, start_freq parameters. Input values were: [10.0], "
            "calculated values were: [2.5].",
        ),
    ],
)
def test_freq_parser_inconsistency_warnings(freq_dict, msg):
    # Now check error cases
    with check_warnings(UserWarning, match=msg):
        simsetup.parse_frequency_params(freq_dict)


@pytest.mark.parametrize(
    "time_keys",
    [
        [
            "start_time",
            "end_time",
            "integration_time",
            "duration_hours",
            "Ntimes",
            "time_array",
        ],
        ["start_time", "end_time", "duration_hours", "Ntimes"],
        ["start_time", "end_time", "integration_time", "Ntimes"],
        ["start_time", "integration_time", "duration_hours"],
        ["end_time", "integration_time", "duration_hours"],
        ["start_time", "integration_time", "Ntimes"],
        ["end_time", "integration_time", "Ntimes"],
        ["start_time", "duration_hours", "Ntimes"],
        ["end_time", "duration_hours", "Ntimes"],
        ["start_time", "duration_hours", "integration_time"],
        ["end_time", "duration_hours", "integration_time"],
        ["time_array"],
    ],
)
def test_time_parser(time_keys, time_dict_base):
    """
    Check a variety of cases for the time parser.
    """
    dayspersec = 1 / (24 * 3600.0)

    subdict = {key: time_dict_base[key] for key in time_keys}
    print(subdict)
    test = simsetup.parse_time_params(subdict)
    print(test)
    np.testing.assert_allclose(
        test["time_array"], time_dict_base["time_array"], atol=dayspersec
    )


@pytest.mark.parametrize(
    ("time_keys", "err_msg"),
    [
        (
            ("duration_hours",),
            "Either start_time or end_time must be specified. The parameters "
            "that were specified were: duration",
        ),
        (
            ("start_time", "Ntimes"),
            "If only one of start_time and end_time is specified, two more of "
            "Ntimes, integration_time and duration must be specified as well. "
            "The parameters that were specified were: Ntimes, start_time",
        ),
        (
            ("start_time", "integration_time"),
            "If only one of start_time and end_time is specified, two more of "
            "Ntimes, integration_time and duration must be specified as well. "
            "The parameters that were specified were: integration_time, start_time",
        ),
        (
            ("start_time", "end_time"),
            "If both start_time and end_time are specified, either Ntimes or "
            "integration_time must be specified as well. The parameters that were "
            "specified were: end_time, start_time",
        ),
    ],
)
def test_time_parser_errors(time_keys, err_msg, time_dict_base):
    # Now check error cases
    subdict = {key: time_dict_base[key] for key in time_keys}
    with pytest.raises(ValueError, match=err_msg):
        simsetup.parse_time_params(subdict)


@pytest.mark.parametrize(
    ("time_dict", "add_base", "rm_keys", "msg"),
    [
        ({}, True, [], None),
        (
            {
                "integration_time": 3.14,
                "start_time": 10000.0,
                "end_time": 80000.3,
                "Ntimes": 30,
            },
            False,
            [],
            "The integration_time is not consistent with the start_time, "
            "end_time, Ntimes specified. Using the values calculated from the "
            "start_time, end_time, Ntimes parameters.",
        ),
        (
            {"Ntimes": 7},
            True,
            ["time_array"],
            "The duration_hours, integration_time is not consistent with the "
            "start_time, end_time, Ntimes specified. Using the values calculated "
            "from the start_time, end_time, Ntimes parameters.",
        ),
        (
            {"time_offset": 2457457.5},
            True,
            [],
            [
                "time_offset is present, but start_time is larger than time_offset, "
                "so not adding time_offset to start_time.",
                "time_offset is present, but end_time is larger than time_offset, "
                "so not adding time_offset to end_time.",
                "time_offset is present, but time_array is larger than time_offset, "
                "so not adding time_offset to time_array.",
            ],
        ),
        (
            {"duration_days": 25.0 / 24.0},
            True,
            [],
            "Both duration_hours and duration_days are specified, using duration_days.",
        ),
    ],
)
def test_time_parser_inconsistent(time_dict, add_base, rm_keys, msg, time_dict_base):
    if add_base:
        time_dict = {**time_dict_base, **time_dict}
        for key in rm_keys:
            del time_dict[key]

    warn_type = UserWarning
    if msg is None:
        warn_type = None
    with check_warnings(warn_type, match=msg):
        simsetup.parse_time_params(time_dict)


def test_single_input_time():
    time_dict = simsetup.time_array_to_params(np.asarray([1.0]), np.asarray([1.0]))
    assert time_dict["integration_time"] == [1.0]


def test_freq_time_params_match(times_and_freqs):
    times, tstep, freqs, fstep = times_and_freqs
    daysperhour = 1 / 24.0
    hourspersec = 1 / 60.0**2
    dayspersec = daysperhour * hourspersec
    int_time_days = np.full_like(times, tstep)
    int_times_sec = int_time_days / dayspersec
    channel_width = np.full_like(freqs, fstep)
    time_dict = simsetup.time_array_to_params(times, int_times_sec)
    assert "start_time" in time_dict
    assert "end_time" in time_dict
    assert "Ntimes" in time_dict
    assert "time_array" not in time_dict

    freq_dict = simsetup.freq_array_to_params(freqs, channel_width)
    assert "start_freq" in freq_dict
    assert "end_freq" in freq_dict
    assert "Nfreqs" in freq_dict
    assert "freq_array" not in freq_dict

    ftest = simsetup.parse_frequency_params(freq_dict)
    ttest = simsetup.parse_time_params(time_dict)
    np.testing.assert_allclose(ftest["freq_array"], freqs)
    np.testing.assert_allclose(ftest["channel_width"], channel_width)
    np.testing.assert_allclose(ttest["time_array"], times)
    np.testing.assert_allclose(ttest["integration_time"], int_times_sec)


def test_uneven_time_array_to_params():
    times = np.linspace(2458570, 2458570 + 0.5, 239)
    # Check that this works for unevenly-spaced times
    times = np.random.choice(times, 150, replace=False)
    times.sort()
    daysperhour = 1 / 24.0
    hourspersec = 1 / 60.0**2
    dayspersec = daysperhour * hourspersec
    time_diff = np.diff(times)
    int_time_days = np.concatenate((time_diff, [time_diff[-1]]))
    int_times_sec = int_time_days / dayspersec
    time_dict = simsetup.time_array_to_params(times, int_times_sec)
    assert "time_array" in time_dict
    assert "integration_time" in time_dict

    ttest = simsetup.parse_time_params(time_dict)
    np.testing.assert_allclose(ttest["time_array"], times)
    np.testing.assert_allclose(ttest["integration_time"], int_times_sec)


def test_single_time_array_to_params():
    # check Ntimes = 1 case
    times = np.linspace(2458570, 2458570.5, 1)
    inttimes = np.asarray([0.5])
    tdict = simsetup.time_array_to_params(times, inttimes)
    assert tdict["time_offset"] == 2458569.5
    assert tdict["time_array"] == times - tdict["time_offset"]
    assert tdict["integration_time"] == inttimes


def test_single_freq_array_to_params():
    freqs = np.linspace(100, 200, 1)
    channel_width = np.asarray([2])
    fdict = simsetup.freq_array_to_params(freqs, channel_width)
    assert fdict["freq_array"] == freqs
    assert fdict["channel_width"] == channel_width


def test_param_select_cross():
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_mwa_nocore.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)

    uv_obj_full = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)

    # test only keeping cross pols
    param_dict["select"] = {"ant_str": "cross"}
    uv_obj_cross = simsetup.initialize_uvdata_from_params(
        param_dict, return_beams=False
    )
    uv_obj_cross2 = uv_obj_full.select(ant_str="cross", inplace=False)

    # histories are different because of time stamp from UVData.new() method
    uv_obj_cross2.history = uv_obj_cross.history

    assert uv_obj_cross == uv_obj_cross2


def test_param_select_bls():
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_mwa_nocore.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)
    uv_obj_full = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)

    # test only keeping certain baselines
    param_dict["select"] = {"bls": "[(40, 41), (42, 43), (44, 45)]"}  # Test as string
    uv_obj_bls = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)

    uv_obj_bls2 = uv_obj_full.select(bls=[(40, 41), (42, 43), (44, 45)], inplace=False)
    uv_obj_bls.history, uv_obj_bls2.history = "", ""
    assert uv_obj_bls == uv_obj_bls2


def test_param_select_redundant():
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_hex37_14.6m.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)
    uv_obj_full = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)

    # test only keeping one baseline per redundant group
    param_dict["select"] = {"redundant_threshold": 0.1}
    uv_obj_red = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)
    uv_obj_red2 = uv_obj_full.compress_by_redundancy(
        tol=0.1, inplace=False, use_grid_alg=True
    )
    uv_obj_red.history, uv_obj_red2.history = "", ""

    assert uv_obj_red == uv_obj_red2
    assert uv_obj_red.Nbls < uv_obj_full.Nbls


@pytest.mark.parametrize(
    "order_dict",
    [
        None,
        {"conjugation_convention": "v>0", "blt_order": ["ant1", "ant2"]},
        {"conjugation_convention": "ant2<ant1", "blt_order": ["baseline"]},
        {"conjugation_convention": "ant2<ant1", "blt_order": "baseline"},
    ],
)
def test_param_ordering(order_dict):
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_mwa_nocore.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)
    uv_obj_orig = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)

    param_dict2 = copy.deepcopy(param_dict)
    if order_dict is None:
        param_dict2.pop("ordering")
        warn_str = [
            "The default baseline conjugation convention has changed. In the past "
            "it was 'ant2<ant1', it now defaults to 'ant1<ant2'. You can specify "
            "the baseline conjugation convention in `obs_param` by setting the "
            "obs_param['ordering']['conjugation_convention'] field. This warning "
            "will go away in version 1.5."
        ]
        warn_type = UserWarning
    else:
        warn_type = None
        warn_str = []
        param_dict2["ordering"] = order_dict

    with check_warnings(warn_type, match=warn_str):
        uv_obj2 = simsetup.initialize_uvdata_from_params(
            param_dict2, return_beams=False
        )

    if order_dict is not None:
        if order_dict["blt_order"] == "baseline" or order_dict["blt_order"] == [
            "baseline"
        ]:
            expected_blt_order = ("baseline", "time")
        else:
            expected_blt_order = tuple(order_dict["blt_order"])
        assert uv_obj2.blt_order == expected_blt_order

        temp_uv = uv_obj_orig.copy()
        temp_uv.conjugate_bls(order_dict["conjugation_convention"])
        temp_uv.reorder_blts(
            order=expected_blt_order[0], minor_order=expected_blt_order[1]
        )
        temp_uv.history = uv_obj2.history
        assert uv_obj2 == temp_uv

        uv_obj2.reorder_blts(
            order="time", minor_order="baseline", conj_convention="ant1<ant2"
        )

    uv_obj2.history = uv_obj_orig.history
    assert uv_obj2 == uv_obj_orig


@pytest.mark.parametrize("key", ["cat_name", "object_name"])
def test_param_set_cat_name(key):
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_mwa_nocore.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)

    param_dict[key] = "foo"
    uv_obj = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)
    assert uv_obj.phase_center_catalog[0]["cat_name"] == "foo"


@pytest.mark.parametrize("pol_array", [[-5, -6], None])
@pytest.mark.parametrize("set_x_orient", [True, False])
def test_param_polarization_array(pol_array, set_x_orient):
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_mwa_nocore.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)

    if pol_array is None:
        pol_array = [-5, -6, -7, -8]
    else:
        param_dict["polarization_array"] = pol_array

    if set_x_orient:
        param_dict["telescope"]["telescope_config_name"] = (
            "../mwa88_nocore_config_dipole.yaml"
        )

    uv_obj, _, _ = simsetup.initialize_uvdata_from_params(param_dict, return_beams=True)
    assert uv_obj.polarization_array.tolist() == pol_array


def test_uvdata_keyword_init(uvdata_keyword_dict):
    # check it runs through

    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    lla_deg = uvd.telescope.location_lat_lon_alt_degrees
    tel_name = uvd.telescope.name

    np.testing.assert_allclose(uvdata_keyword_dict["telescope_location"], lla_deg)
    np.testing.assert_allclose(
        uvdata_keyword_dict["integration_time"],
        uvd.integration_time,
        rtol=uvd._integration_time.tols[0],
        atol=uvd._integration_time.tols[1],
    )
    assert uvdata_keyword_dict["telescope_name"] == tel_name
    assert uvdata_keyword_dict["start_freq"] == uvd.freq_array[0]
    assert uvdata_keyword_dict["start_time"] == uvd.time_array[0]
    assert uvdata_keyword_dict["Ntimes"] == uvd.Ntimes
    assert uvdata_keyword_dict["Nfreqs"] == uvd.Nfreqs
    assert uvdata_keyword_dict["polarization_array"] == uvutils.polnum2str(
        uvd.polarization_array
    )
    assert not np.any(uvd.ant_1_array == uvd.ant_2_array)


def test_uvdata_keyword_init_select_bls(uvdata_keyword_dict):
    # check bls and antenna_nums selections work
    bls = [(0, 1), (0, 2), (0, 3)]
    uvdata_keyword_dict["bls"] = bls
    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    antpairs = uvd.get_antpairs()
    assert antpairs == bls


@pytest.mark.parametrize("pols", [["xx", "yy"], ["xx", "yy", "xy", "yx"], ["yy"], None])
def test_uvdata_keyword_init_polarization_array(uvdata_keyword_dict, pols):
    if pols is None:
        del uvdata_keyword_dict["polarization_array"]
        pols = ["xx", "yy", "xy", "yx"]
    else:
        uvdata_keyword_dict["polarization_array"] = pols

    for feed_angle in [[np.pi / 2, 0], [0, np.pi / 2], [np.pi / 4, np.pi * 3 / 4]]:
        uvdata_keyword_dict["feed_angle"] = feed_angle
        uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
        exp_pols = copy.deepcopy(pols)
        if feed_angle[0] == np.pi / 2:
            exp_pols = [pol.replace("x", "e") for pol in exp_pols]
            exp_pols = [pol.replace("y", "n") for pol in exp_pols]
        elif feed_angle[0] == 0:
            exp_pols = [pol.replace("x", "n") for pol in exp_pols]
            exp_pols = [pol.replace("y", "e") for pol in exp_pols]

        assert uvd.get_pols() == exp_pols


def test_uvdata_keyword_init_select_antnum_str(uvdata_keyword_dict):
    # check that '1' gets converted to [1]
    uvdata_keyword_dict["no_autos"] = False
    uvdata_keyword_dict["antenna_nums"] = "1"
    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)

    assert uvd.Nbls == 1
    assert uvd.polarization_array.tolist() == [-5]


@pytest.mark.filterwarnings("ignore:The Nfreqs, bandwidth is not consistent")
@pytest.mark.filterwarnings("ignore:The Ntimes is not consistent")
def test_uvdata_keyword_init_time_freq_override(uvdata_keyword_dict):
    # check time and freq array definitions supersede other parameters
    freq_array = np.linspace(100, 200, 11) * 1e6
    time_array = np.linspace(2458101, 2458102, 21)
    uvdata_keyword_dict["freq_array"] = freq_array
    uvdata_keyword_dict["time_array"] = time_array
    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)

    np.testing.assert_allclose(uvd.time_array[:: uvd.Nbls], time_array)
    np.testing.assert_allclose(uvd.freq_array, freq_array)


def test_uvdata_keyword_init_layout_dict(uvdata_keyword_dict, tmpdir):
    # test feeding array layout as dictionary
    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    antpos = uvd.telescope.get_enu_antpos()
    ants = uvd.telescope.antenna_numbers
    antpos_d = dict(zip(ants, antpos, strict=False))
    layout_fname = "temp_layout.csv"
    obsparam_fname = "temp_obsparam.yaml"

    uvdata_keyword_dict["output_layout_filename"] = layout_fname
    uvdata_keyword_dict["output_yaml_filename"] = obsparam_fname
    uvdata_keyword_dict["array_layout"] = antpos_d
    uvdata_keyword_dict["path_out"] = str(tmpdir)
    uvdata_keyword_dict["write_files"] = True

    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    layout_path = str(tmpdir.join(layout_fname))
    obsparam_path = str(tmpdir.join(obsparam_fname))
    assert os.path.exists(layout_path)
    assert os.path.exists(obsparam_path)

    assert uvd.Nbls == 6
    assert uvd.Nants_data == 4
    antpos2 = uvd.telescope.get_enu_antpos()
    ants2 = uvd.telescope.antenna_numbers
    antpos_d2 = dict(zip(ants2, antpos2, strict=False))
    assert np.all([np.isclose(antpos_d[ant], antpos_d2[ant]) for ant in ants])


@pytest.mark.parametrize("pass_layout", [True, False])
def test_uvdata_keyword_init_write(pass_layout, uvdata_keyword_dict, tmpdir):
    # Check defaults when writing to file.

    new_layout_file = os.path.join(tmpdir, "triangle_bl_layout.csv")
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "triangle_bl_layout.csv"),
        new_layout_file,
    )

    uvd = simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    antpos = uvd.telescope.get_enu_antpos()
    ants = uvd.telescope.antenna_numbers
    antpos_d = dict(zip(ants, antpos, strict=False))
    # Checking -- Default to a copy of the original layout, if layout is provided.
    obsparam_fname = os.path.join(tmpdir, "obsparam.yaml")

    if pass_layout:
        layout_fname = os.path.join(tmpdir, "triangle_bl_layout.csv")
    else:
        layout_fname = os.path.join(tmpdir, "antenna_layout.csv")
        uvdata_keyword_dict.pop("antenna_layout_filepath")
        uvdata_keyword_dict["array_layout"] = antpos_d
        uvdata_keyword_dict["complete"] = True

    uvdata_keyword_dict["write_files"] = True

    cwd = os.getcwd()
    os.chdir(tmpdir)
    simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)
    assert os.path.exists(layout_fname)
    assert os.path.exists(obsparam_fname)

    os.chdir(cwd)


def test_initialize_uvdata_from_keywords_errors(uvdata_keyword_dict):
    del uvdata_keyword_dict["antenna_layout_filepath"]

    with pytest.raises(
        ValueError,
        match="Either array_layout or antenna_layout_filepath must be passed.",
    ):
        simsetup.initialize_uvdata_from_keywords(**uvdata_keyword_dict)


@pytest.mark.parametrize("int_time_varies", [False, "single", "single_time"])
@pytest.mark.parametrize("use_cli", [True, False])
def test_uvfits_to_config(tmp_path, int_time_varies, use_cli):
    """
    Loopback test of reading parameters from uvfits file, generating uvfits file, and reading
    in again.
    """
    opath = os.path.join(tmp_path, "uvfits_yaml_temp")
    param_filename = "obsparam.yaml"
    second_param_filename = "test2_config.yaml"
    if not os.path.exists(opath):
        os.makedirs(opath)  # Directory will be deleted when test completed.

    # Read uvfits file to params.
    uv0 = UVData.from_file(triangle_uvfits_file)

    if use_cli:
        path = opath
        layout_fname = os.path.join(opath, "test_layout.csv")
        telescope_config = os.path.join(opath, "test_tel_config.yaml")
        subprocess.check_output(  # nosec
            [
                "uvdata_to_telescope_config",
                triangle_uvfits_file,
                "-b",
                herabeam_default,
                "-l",
                layout_fname,
                "-t",
                telescope_config,
            ]
        )
    else:
        path, telescope_config, layout_fname = simsetup.uvdata_to_telescope_config(
            uv0, herabeam_default, path_out=opath, return_names=True
        )
    tel_config_full = os.path.join(path, telescope_config)
    with open(tel_config_full) as yf:
        telconfig = yaml.safe_load(yf)

    assert 0 in telconfig["beam_paths"]

    warn_type = None
    msg = ""
    if int_time_varies == "single":
        # Test case of a single non-uniform integration times: no effect
        uv0.integration_time[-1] += 2
        warn_type = UserWarning
        msg = (
            "integration time varies for unique times, using the shortest "
            "integration time for each unique time."
        )
    elif int_time_varies == "single_time":
        # Test case of a non-uniform integration time for last time
        assert uv0.Ntimes > 1
        time_array, _, unique_inverse = np.unique(
            uv0.time_array, return_index=True, return_inverse=True
        )
        assert time_array.size < uv0.time_array.size
        uv0.integration_time[np.nonzero(unique_inverse == unique_inverse[-1])] += 2
        assert not uvutils.tools._test_array_constant(
            uv0.integration_time, tols=UVData()._integration_time.tols
        )
    with check_warnings(warn_type, match=msg):
        simsetup.uvdata_to_config_file(
            uv0,
            telescope_config_name=os.path.join(path, telescope_config),
            layout_csv_name=os.path.join(path, layout_fname),
            path_out=opath,
        )

    # From parameters, generate a uvdata object.
    param_dict = simsetup._config_str_to_dict(os.path.join(opath, param_filename))
    # make a copy because the parameter dictionary gets modified in the function below.
    orig_param_dict = copy.deepcopy(param_dict)
    uv1 = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)
    # Generate parameters from new uvfits and compare with old.
    path, telescope_config, layout_fname = simsetup.uvdata_to_telescope_config(
        uv1,
        herabeam_default,
        telescope_config_name=telescope_config,
        layout_csv_name=layout_fname,
        path_out=opath,
        return_names=True,
    )

    if use_cli:
        # add the data-like arrays so it can be written to disk
        shape = (uv1.Nblts, uv1.Nfreqs, uv1.Npols)
        uv1.data_array = np.zeros(shape, dtype=complex)
        uv1.flag_array = np.zeros(shape, dtype=bool)
        uv1.nsample_array = np.ones(shape, dtype=float)

        uvd_file = os.path.join(tmp_path, "temp.uvh5")
        uv1.write_uvh5(uvd_file)

        subprocess.check_output(  # nosec
            [
                "uvdata_to_config",
                uvd_file,
                "-p",
                second_param_filename,
                "-t",
                os.path.join(path, telescope_config),
                "-l",
                os.path.join(path, layout_fname),
                "--outpath",
                opath,
            ]
        )
    else:
        simsetup.uvdata_to_config_file(
            uv1,
            param_filename=second_param_filename,
            telescope_config_name=os.path.join(path, telescope_config),
            layout_csv_name=os.path.join(path, layout_fname),
            path_out=opath,
        )

    del param_dict

    param_dict = simsetup._config_str_to_dict(
        os.path.join(opath, second_param_filename)
    )
    shutil.rmtree(opath)
    assert param_dict["obs_param_file"] == second_param_filename
    assert orig_param_dict["obs_param_file"] == param_filename
    orig_param_dict["obs_param_file"] = second_param_filename
    assert compare_dictionaries(param_dict, orig_param_dict)


@pytest.mark.parametrize(
    ("arrangement", "text_cat"),
    [
        ("cross", "mock_cross_2458098.27471.txt"),
        ("hera_text", "mock_hera_text_2458098.27471.txt"),
        ("long-line", "mock_long-line_2458098.27471.txt"),
        ("off-zenith", "mock_off-zenith_2458098.27471.txt"),
        ("triangle", "mock_triangle_2458098.27471.txt"),
        ("random", "mock_random_2458098.27471.txt"),
        ("zenith", "mock_zenith_2458098.27471.txt"),
    ],
)
def test_mock_catalogs(arrangement, text_cat):
    time = Time(2458098.27471265, scale="utc", format="jd")
    cat, mock_kwds = simsetup.create_mock_catalog(time, arrangement, rseed=2458098)

    # For each mock catalog, verify the Ra/Dec source positions against a saved text catalog.
    radec_catalog = SkyModel.from_file(
        os.path.join(SIM_DATA_PATH, "test_catalogs", text_cat)
    )
    assert np.all(radec_catalog == cat)


@pytest.mark.filterwarnings("ignore:UVW orientation appears to be flipped")
def test_mock_catalog_gaussbeam_values():
    """
    Make the long-line point sources up to 10 degrees from zenith.
    Confirm that the coherencies match the expected beam values at those zenith angles.
    """
    sigma = 0.05
    hera_uv = UVData()
    EW_uvfits_file = os.path.join(SIM_DATA_PATH, "28mEWbl_1time_1chan.uvfits")
    hera_uv.read_uvfits(EW_uvfits_file)

    array_location = hera_uv.telescope.location

    freq = hera_uv.freq_array[0] * units.Hz

    time = Time(hera_uv.time_array[0], scale="utc", format="jd")

    catalog, _ = simsetup.create_mock_catalog(
        time=time,
        arrangement="long-line",
        Nsrcs=41,
        min_alt=80.0,
        array_location=array_location,
    )

    catalog.update_positions(time, array_location)
    beam = GaussianBeam(sigma=sigma)
    array = pyuvsim.Telescope(
        "telescope_name", array_location, pyuvsim.BeamList([beam])
    )

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
    np.testing.assert_allclose(beam_values**2, coherencies)


def test_saved_mock_catalog(tmpdir):
    time = Time(2458098.27471265, scale="utc", format="jd")

    cwd = os.getcwd()
    os.chdir(tmpdir)

    cat, mock_kwds = simsetup.create_mock_catalog(time, "random", Nsrcs=100, save=True)
    loc = ast.literal_eval(mock_kwds["array_location"])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    fname = "mock_catalog_random.npz"
    alts_reload = np.load(fname)["alts"]
    cat.update_positions(time, loc)
    alt, _ = cat.alt_az

    np.testing.assert_allclose(alts_reload, np.degrees(alt))

    os.chdir(cwd)


@pytest.mark.parametrize("min_alt", [-20, 0, None, 50])
def test_randsource_minalt(min_alt):
    time = Time(2458098.27471265, scale="utc", format="jd")
    cat, mock_kwds = simsetup.create_mock_catalog(
        time, "random", Nsrcs=100, min_alt=min_alt
    )
    loc = ast.literal_eval(mock_kwds["array_location"])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    cat.update_positions(time, loc)
    alt, _ = cat.alt_az
    if min_alt is None:
        min_alt = 30  # Checking default
    assert np.all(alt >= np.radians(min_alt))


def test_randsource_distribution():
    # Check that random sources are uniformly distributed per solid angle

    astropy_healpix = pytest.importorskip("astropy_healpix")
    Nsrcs = 40000
    time = Time(2458098.27471265, scale="utc", format="jd")
    cat, mock_kwds = simsetup.create_mock_catalog(
        time, "random", Nsrcs=Nsrcs, min_alt=-90, rseed=2458098
    )
    loc = ast.literal_eval(mock_kwds["array_location"])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    cat.update_positions(time, loc)
    alt, az = cat.alt_az

    # Bin into low-res HEALPix pixels.
    nside = 32
    npix = 12 * nside**2
    hp = astropy_healpix.HEALPix(nside)

    inds = hp.lonlat_to_healpix(Longitude(az, unit="rad"), Latitude(alt, unit="rad"))
    un, counts = np.unique(inds, return_counts=True)

    # counts should be Poisson-distributed with rate lambda = nsrcs / npix
    # variance and mean should be close to lambda
    lam = Nsrcs / npix
    assert np.isclose(np.mean(counts), lam, atol=1.0)
    assert np.isclose(np.var(counts), lam, atol=1.0)


def test_mock_catalog_error():
    time = Time(2458098.27471265, scale="utc", format="jd")
    with pytest.raises(
        KeyError, match="Invalid mock catalog arrangement: invalid_catalog_name"
    ):
        simsetup.create_mock_catalog(time, "invalid_catalog_name")


def test_keyword_param_loop(tmpdir):
    # Check that yaml/csv files made by intialize_uvdata_from_keywords will work
    # on their own.
    layout_fname = "temp_layout_kwdloop.csv"
    obsparam_fname = "temp_obsparam_kwdloop.yaml"
    path_out = str(tmpdir)
    # add some jiggle so you get non-zero uvws
    antpos_enu = (np.ones(30) + np.random.uniform(-10, 10, 30)).reshape((10, 3))
    antnums = np.arange(10)
    antpos_d = dict(zip(antnums, antpos_enu, strict=False))

    uvd = simsetup.initialize_uvdata_from_keywords(
        array_layout=antpos_d,
        telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
        telescope_name="HERA",
        Nfreqs=10,
        start_freq=1e8,
        bandwidth=1e8,
        Ntimes=60,
        integration_time=100.0,
        start_time=2458101.0,
        no_autos=True,
        conjugation_convention="ant1<ant2",
        path_out=path_out,
        antenna_layout_filepath=layout_fname,
        output_yaml_filename=obsparam_fname,
        feed_array=["x", "y"],
        feed_angle=[np.pi / 2, 0],
        mount_type="fixed",
    )

    uv2 = simsetup.initialize_uvdata_from_params(
        os.path.join(path_out, obsparam_fname), return_beams=False
    )

    uv2.extra_keywords = {}
    uvd.extra_keywords = {}  # These will not match

    uv2.history = uvd.history

    assert uv2 == uvd


def test_multi_analytic_beams(tmpdir):
    # Test inline definitions of beam attributes.
    # eg. (in beam configuration file):
    #
    # beam_paths:
    #   0 : airy, diameter=14
    #   1 : airy, diameter=20
    #   2 : gaussian, sigma=0.5
    par_fname = str(tmpdir.join("test_teleconfig.yaml"))
    layout_fname = str(tmpdir.join("test_layout_5ant.csv"))

    telescope_location = (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
    telescope_name = "SKA"
    mt_kwargs = {"mount_type": None}
    beam_specs = {
        0: AiryBeam(diameter=14, **mt_kwargs),
        1: AiryBeam(diameter=20),
        2: GaussianBeam(sigma=0.5),
    }
    expected_beams = [
        AiryBeam(diameter=14, **mt_kwargs),
        AiryBeam(diameter=20),
        GaussianBeam(sigma=0.5),
    ]

    Nants = 5
    antenna_numbers = np.arange(Nants)
    antpos = np.zeros((Nants, 3))
    antpos[:, 0] = np.arange(Nants)
    names = antenna_numbers.astype(str)
    beam_ids = [0, 1, 2, 2, 0]
    simsetup._write_layout_csv(layout_fname, antpos, names, antenna_numbers, beam_ids)

    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    # Write tele config to file.
    pdict = {
        "telescope_location": str(telescope_location),
        "telescope_name": telescope_name,
        "beam_paths": beam_specs,
    }
    with open(par_fname, "w") as yfile:
        yaml.safe_dump(pdict, yfile, default_flow_style=False)

    param_dict = {"telescope_config_name": par_fname, "array_layout": layout_fname}

    pdict, beam_list, beam_dict = simsetup.parse_telescope_params(
        param_dict, config_path=str(tmpdir), freq_array=freqs
    )

    exp_mount_type = []
    for i, nm in enumerate(names):
        bid = beam_ids[i]
        assert beam_dict[nm] == bid
        assert beam_list[bid].beam == expected_beams[bid]
        if expected_beams[bid].mount_type is None:
            exp_mount_type.append("other")
        else:
            exp_mount_type.append(expected_beams[bid].mount_type)

    assert np.all(pdict["mount_type"] == exp_mount_type)


def test_direct_fname(tmpdir):
    new_tel_config = os.path.join(tmpdir, "28m_triangle_10time_10chan.yaml")
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "28m_triangle_10time_10chan.yaml"),
        new_tel_config,
    )
    new_param_file = os.path.join(tmpdir, "param_100times_1.5days_triangle.yaml")
    shutil.copyfile(
        os.path.join(
            SIM_DATA_PATH, "test_config", "param_100times_1.5days_triangle.yaml"
        ),
        new_param_file,
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "triangle_bl_layout.csv"),
        os.path.join(tmpdir, "triangle_bl_layout.csv"),
    )

    cwd = os.getcwd()
    os.chdir(tmpdir)

    simsetup.initialize_uvdata_from_params(
        "param_100times_1.5days_triangle.yaml", return_beams=False
    )

    os.chdir(cwd)


@pytest.mark.parametrize(
    ("input_yaml", "beam_id", "err_msg"),
    [
        (
            """
            beam_paths:
                0: 1.35
            spline_interp_opts:
                    kx: 4
                    ky: 4
            freq_interp_kind: 'cubic'
            telescope_location: (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
            telescope_name: BLLITE
            """,
            0,
            "Beam model is not properly specified",
        ),
        (
            """
            beam_paths:
                0:
                    diameter: 12
            spline_interp_opts:
                    kx: 4
                    ky: 4
            freq_interp_kind: 'cubic'
            telescope_location: (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
            telescope_name: BLLITE
            """,
            0,
            "Beam model must have either a 'filename' field for UVBeam files "
            "or a 'type' field for analytic beams.",
        ),
        (
            """
            beam_paths:
                0:
                    type: unsupported_type
                    diameter: 12
            spline_interp_opts:
                    kx: 4
                    ky: 4
            freq_interp_kind: 'cubic'
            telescope_location: (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
            telescope_name: BLLITE
            """,
            0,
            "Undefined beam model type: unsupported_type",
        ),
        (
            """
            beam_paths:
                0: !AnalyticBeam
                    class: AiryBeam
                    diameter: 12
            spline_interp_opts:
                    kx: 4
                    ky: 4
            freq_interp_kind: 'cubic'
            telescope_location: (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
            telescope_name: BLLITE
            """,
            1,
            "beam_id 1 is not listed in the telconfig.",
        ),
        (
            """
            beam_paths:
                0: !AnalyticBeam
                    class: AiryBeam
                    diameter: 12
            diameter: 12
            spline_interp_opts:
                    kx: 4
                    ky: 4
            freq_interp_kind: 'cubic'
            telescope_location: (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
            telescope_name: BLLITE
            """,
            0,
            "Beam shape options diameter and sigma should be specified per "
            "beamID in the 'beam_paths' section not as globals. For examples "
            "see the parameter_files documentation.",
        ),
    ],
)
def test_beamlist_init_errors(input_yaml, beam_id, err_msg):
    telconfig = yaml.safe_load(input_yaml)

    warn_msg = []
    if ("!AnalyticBeam" not in input_yaml) and not isinstance(
        telconfig["beam_paths"][0], float
    ):
        warn_msg += [
            "Entries in 'beam_paths' should be specified using either the UVBeam "
            "or AnalyticBeam constructors or using a dict syntax for UVBeams. "
            "For examples see the parameter_files documentation. Specifying "
            "analytic beam without the AnalyticBeam constructors will cause an "
            "error in version 1.6"
        ]
        warn_type = DeprecationWarning
    if len(warn_msg) == 0:
        warn_type = None

    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    with (
        check_warnings(warn_type, match=warn_msg),
        pytest.raises(ValueError, match=err_msg),
    ):
        simsetup._construct_beam_list([beam_id], telconfig, freq_array=freqs)


def test_beamlist_init_spline():
    telescope_config_name = os.path.join(SIM_DATA_PATH, "bl_lite_mixed.yaml")
    with open(telescope_config_name) as yf:
        telconfig = yaml.safe_load(yf)

    # The path for beam 0 is invalid, and it's not needed for this test.
    del telconfig["beam_paths"][0]
    beam_ids = np.arange(1, 6)

    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    # Check that spline_interp_opts is passed along correctly to BeamList
    telconfig["spline_interp_opts"] = {"kx": 2, "ky": 2}
    with check_warnings(
        DeprecationWarning,
        match=[
            "Entries in 'beam_paths' should be specified using either the UVBeam "
            "or AnalyticBeam constructors or using a dict syntax for UVBeams. "
            "For examples see the parameter_files documentation. Specifying "
            "analytic beam without the AnalyticBeam constructors will cause an "
            "error in version 1.6"
        ]
        * 5,
    ):
        beam_list = simsetup._construct_beam_list(beam_ids, telconfig, freq_array=freqs)
    assert isinstance(beam_list, pyuvsim.BeamList)
    assert beam_list.spline_interp_opts == {"kx": 2, "ky": 2}


@pytest.mark.parametrize(
    ("rename_beamfits", "pass_beam_type"), [(False, True), (True, True), (False, False)]
)
def test_beamlist_init(rename_beamfits, pass_beam_type, tmp_path, mwa_beam_path):
    telescope_config_name = os.path.join(SIM_DATA_PATH, "bl_lite_mixed.yaml")
    with open(telescope_config_name) as yf:
        telconfig = yaml.safe_load(yf)

    beamfits_file = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")

    if rename_beamfits:
        new_beam_file = os.path.join(tmp_path, "HERA_NicCST.uvbeam")
        shutil.copyfile(beamfits_file, new_beam_file)

        telconfig["beam_paths"][0] = {"filename": new_beam_file}
    else:
        telconfig["beam_paths"][0] = {"filename": "HERA_NicCST.beamfits"}

    if pass_beam_type:
        beam_file = beamfits_file
        telconfig["beam_paths"][0] = {"filename": beam_file, "file_type": "beamfits"}

    telconfig["beam_paths"][0]["mount_type"] = "fixed"

    telconfig["beam_paths"][6]["filename"] = mwa_beam_path

    nbeams = len(telconfig["beam_paths"])
    entries_warnings = 6

    warn_list = [
        "Entries in 'beam_paths' should be specified using either the UVBeam "
        "or AnalyticBeam constructors or using a dict syntax for UVBeams. "
        "For examples see the parameter_files documentation. Specifying "
        "analytic beam without the AnalyticBeam constructors will cause an "
        "error in version 1.6"
    ] * entries_warnings

    warn_types = DeprecationWarning

    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    with check_warnings(warn_types, match=warn_list):
        beam_list = simsetup._construct_beam_list(
            np.arange(nbeams), telconfig, freq_array=freqs
        )

    # How the beam attributes should turn out for this file:
    assert isinstance(beam_list[0].beam, UVBeam)
    assert isinstance(beam_list[1].beam, AiryBeam)
    assert beam_list[1].beam.diameter == 16
    assert isinstance(beam_list[2].beam, GaussianBeam)
    assert beam_list[2].beam.sigma == 0.03
    assert isinstance(beam_list[3].beam, AiryBeam)
    assert beam_list[3].beam.diameter == 12
    assert isinstance(beam_list[4].beam, GaussianBeam)
    assert beam_list[4].beam.diameter == 14
    assert isinstance(beam_list[5].beam, GaussianBeam)
    assert beam_list[5].beam.diameter == 12
    assert isinstance(beam_list[6].beam, UVBeam)
    assert isinstance(beam_list[7].beam, ShortDipoleBeam)


@pytest.mark.parametrize("sel_type", ["freq_range", "freq_range_yaml", "freq_buff"])
@pytest.mark.parametrize(
    "tel_config", ["bl_lite_mixed.yaml", "bl_lite_mixed_constructors.yaml"]
)
def test_beamlist_init_freqrange(sel_type, tel_config, mwa_beam_path):
    telescope_config_name = os.path.join(SIM_DATA_PATH, tel_config)
    hera_beam_path = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")

    with open(telescope_config_name) as yf:
        lines = yf.readlines()
    lines[2] = lines[2].replace("hera.beamfits", hera_beam_path)
    lines[20] = lines[2].replace("mwa_full_EE_test.h5", mwa_beam_path)

    telconfig = yaml.safe_load("\n".join(lines))

    Nfreqs = 10
    freqs = np.linspace(120, 145, Nfreqs) * 1e6

    freq_range = [117e6, 148e6]
    if sel_type == "freq_range":
        freq_range_pass = freq_range
    elif sel_type == "freq_range_yaml":
        telconfig["select"] = {"freq_range": freq_range}
        freq_range_pass = None
    else:
        telconfig["select"] = {"freq_buffer": 0}
        freq_range_pass = None

    warn_list = []
    if "constructors" not in tel_config:
        warn_list += [
            "Entries in 'beam_paths' should be specified using either the UVBeam "
            "or AnalyticBeam constructors or using a dict syntax for UVBeams. "
            "For examples see the parameter_files documentation. Specifying "
            "analytic beam without the AnalyticBeam constructors will cause an "
            "error in version 1.6"
        ] * 5 + ["The mount_type parameter must be set for UVBeam objects"]

    if len(warn_list) == 0:
        warn_use = None
    else:
        warn_use = DeprecationWarning
    with check_warnings(warn_use, match=warn_list):
        beam_list = simsetup._construct_beam_list(
            np.arange(6), telconfig, freq_range=freq_range_pass, freq_array=freqs
        )

    # How the beam attributes should turn out for this file:
    assert isinstance(beam_list[0].beam, UVBeam)
    assert len(beam_list[0].beam.freq_array) == 2

    with pytest.raises(ValueError, match="If passed, freq_range have 2 elements"):
        beam_list = simsetup._construct_beam_list(
            np.arange(6), telconfig, freq_range=[117e6], freq_array=freqs
        )


def test_moon_lsts():
    # Check that setting lsts for a Moon simulation works as expected.
    pytest.importorskip("lunarsky")
    from lunarsky import Time as LTime

    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_tranquility_hex.yaml"
    )
    param_dict = simsetup._config_str_to_dict(param_filename)
    uv_obj = simsetup.initialize_uvdata_from_params(param_dict, return_beams=False)
    assert "world" in uv_obj.extra_keywords
    assert uv_obj.extra_keywords["world"] == "moon"

    # Check ordering on lsts -- unique lsts should correspond with unique times.
    lsts, lst_inds = np.unique(uv_obj.lst_array, return_inverse=True)
    times, time_inds = np.unique(uv_obj.time_array, return_inverse=True)
    assert np.all(lst_inds == time_inds)

    # Confirm that the LSTs match the zenith RAs for the telescope_location
    loc = uv_obj.telescope.location
    times_quant = LTime(times, format="jd", location=loc)
    skycrds = SkyCoord(
        alt=["90d"] * times.size,
        az=["0d"] * times.size,
        frame="lunartopo",
        obstime=times_quant,
        location=loc,
    ).transform_to("icrs")

    np.testing.assert_allclose(skycrds.ra.rad, lsts, atol=1e-7, rtol=0)

    # Unset the lst array and confirm that the call from _complete_uvdata returns the same.
    backup_lst_array = uv_obj.lst_array.copy()
    uv_obj.lst_array = None

    new_obj = simsetup._complete_uvdata(uv_obj)

    np.testing.assert_allclose(new_obj.lst_array, backup_lst_array)
    assert new_obj.check()


def test_mock_catalog_moon():
    # A mock catalog made with a MoonLocation.
    pytest.importorskip("lunarsky")
    from lunarsky import MoonLocation, Time as LTime

    time = LTime.now()
    loc = MoonLocation.from_selenodetic(24.433333333, 0.687500000)
    mmock, mkwds = simsetup.create_mock_catalog(time, "hera_text", array_location=loc)
    eloc = EarthLocation.from_geodetic(24.433, 0.6875)
    emock, ekwds = simsetup.create_mock_catalog(time, "hera_text", array_location=eloc)

    assert mkwds["world"] == "moon"
    assert ekwds["world"] == "earth"

    # Simple check that the given lat/lon were interpreted differently in each call.
    # use not == rather than != to avoid a pyradiosky bug, fixed in v0.1.4
    assert mmock != emock


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize("unit", ["Jy", "K", "Jy/sr", "K sr"])
def test_skymodeldata_with_quantity_stokes(unit, cat_with_some_pols):
    # Support for upcoming pyradiosky change setting SkyModel.stokes
    # to an astropy Quantity.
    if unit in ["Jy", "K sr"]:
        sky = cat_with_some_pols
        if unit != "Jy":
            sky.jansky_to_kelvin()
    else:
        pytest.importorskip("analytic_diffuse")
        pytest.importorskip("astropy_healpix")
        sky, _ = simsetup.create_mock_catalog(
            Time.now(), arrangement="diffuse", diffuse_model="monopole", map_nside=16
        )
        sky.spectral_type = "spectral_index"
        sky.reference_frequency = np.full(sky.Ncomponents, 1e8) * units.Hz
        sky.spectral_index = np.full(sky.Ncomponents, -0.8)
        sky.check()

        if unit != "K":
            sky.kelvin_to_jansky()

    if not isinstance(sky.stokes, units.Quantity):
        sky.stokes *= units.Unit(unit)

    smd = simsetup.SkyModelData(sky)
    assert np.all(sky.stokes.to_value(unit)[0] == smd.stokes_I)
    assert units.Unit(smd.flux_unit) == units.Unit(unit)

    sky2 = smd.get_skymodel()
    assert sky2 == sky


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize("component_type", ["point", "healpix"])
@pytest.mark.parametrize("select", [True, False])
def test_skymodeldata(component_type, select, cat_with_some_pols):
    # Test that SkyModelData class can properly recreate a SkyModel and subselect.
    if component_type == "point":
        sky = cat_with_some_pols
        filename_use = "mock_with_pol"
    else:
        pytest.importorskip("astropy_healpix")
        path = os.path.join(SKY_DATA_PATH, "healpix_disk.skyh5")
        sky = SkyModel.from_file(path)
        filename_use = ["healpix_disk"]

    smd = simsetup.SkyModelData(sky, filename=filename_use)
    if isinstance(filename_use, str):
        assert smd.filename == [filename_use]
    else:
        assert smd.filename == filename_use

    sky_ra, sky_dec = sky.get_lon_lat()
    assert (smd.ra == sky_ra.deg).all()
    assert (smd.dec == sky_dec.deg).all()

    if isinstance(sky.stokes, units.Quantity):
        smd.stokes_I *= units.Unit(smd.flux_unit)
        if smd.polarized is not None:
            smd.stokes_Q *= units.Unit(smd.flux_unit)
            smd.stokes_U *= units.Unit(smd.flux_unit)
            smd.stokes_V *= units.Unit(smd.flux_unit)

    assert (smd.stokes_I == sky.stokes[0]).all()

    if smd.polarized is not None:
        assert (smd.stokes_Q == sky.stokes[..., smd.polarized][1]).all()
        assert (smd.stokes_U == sky.stokes[..., smd.polarized][2]).all()
        assert (smd.stokes_V == sky.stokes[..., smd.polarized][3]).all()

    if select:
        sky1_sub = smd.get_skymodel(range(8, 13))

        assert sky1_sub.Ncomponents == 5
        if smd.polarized is not None:
            assert sky1_sub._n_polarized == 1

        sky_sub = sky.select(component_inds=list(range(8, 13)), inplace=False)
        sky1_sub.history = sky_sub.history
        assert sky_sub == sky1_sub
    else:
        # Make skymodel from SkyModelData.
        sky1 = smd.get_skymodel()
        # history is not copied into SkyModelData.
        sky1.history = sky.history

        assert sky1 == sky


@pytest.mark.parametrize("inds", [range(30), range(5), np.arange(9, 14)])
def test_skymodeldata_pol_select(inds, cat_with_some_pols):
    # When running SkyModelData.subselect, confirm that the
    # polarization array and Q, U, V are properly selected.

    smd = simsetup.SkyModelData(cat_with_some_pols)
    sub_smd = smd.subselect(inds)

    test_q = np.zeros((smd.Nfreqs, smd.Ncomponents))
    temp = np.zeros((sub_smd.Nfreqs, sub_smd.Ncomponents))
    temp[..., sub_smd.polarized] = sub_smd.stokes_Q
    test_q[..., inds] = temp[()]

    full_q = np.zeros_like(test_q)
    full_q[..., smd.polarized] = smd.stokes_Q

    assert np.all(full_q[..., inds] == test_q[..., inds])


def test_skymodeldata_non_icrs(cat_with_some_pols):
    ra, dec = cat_with_some_pols.get_lon_lat()
    gcrs_coord = SkyCoord(ra, dec, frame="gcrs")
    icrs_coord = gcrs_coord.transform_to("icrs")

    sky = SkyModel(
        name=cat_with_some_pols.name,
        ra=ra,
        dec=dec,
        frame="gcrs",
        stokes=cat_with_some_pols.stokes,
        spectral_type=cat_with_some_pols.spectral_type,
        freq_array=cat_with_some_pols.freq_array,
    )
    smd = simsetup.SkyModelData(sky)

    np.testing.assert_allclose(icrs_coord.ra.deg, smd.ra)
    np.testing.assert_allclose(icrs_coord.dec.deg, smd.dec)


@pytest.mark.parametrize("spectral_type", ["spectral_index", "subband"])
def test_skymodel_data_spectype_select(spectral_type):
    sky = SkyModel.from_file(
        os.path.join(SIM_DATA_PATH, "gleam_50srcs.vot"),
        spectral_type=spectral_type,
        run_check=False,
    )
    sky.select(non_nan="all", non_negative=True)

    smd = simsetup.SkyModelData(sky)

    sky1_sub = smd.get_skymodel(range(8, 13))

    assert sky1_sub.Ncomponents == 5

    sky_sub = sky.select(component_inds=list(range(8, 13)), inplace=False)
    sky1_sub.history = sky_sub.history
    assert sky_sub == sky1_sub


@pytest.mark.parametrize("inds", [range(30), range(5)])
def test_skymodeldata_attr_bases(inds, cat_with_some_pols):
    # Check that downselecting doesn't copy length-Ncomponent arrays.

    smd = simsetup.SkyModelData(cat_with_some_pols)
    smd_copy = smd.subselect(inds)
    assert smd_copy.ra.base is smd.ra.base
    assert smd_copy.dec.base is smd.dec.base
    assert smd_copy.stokes_I.base is smd.stokes_I.base


def test_simsetup_with_obsparam_freq_buffer():
    fl = os.path.join(SIM_DATA_PATH, "test_config", "obsparam_diffuse_sky_freqbuf.yaml")

    warn_list = [
        "Beam selections should be specified in the telescope "
        "configuration, not in the obsparam. This will become an error in "
        "version 1.5"
    ]

    with check_warnings(DeprecationWarning, match=warn_list):
        _, beams, _ = simsetup.initialize_uvdata_from_params(fl, return_beams=True)

    assert beams[0].beam.freq_array.max() < 101e6


def test_simsetup_with_freq_buffer():
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_diffuse_sky_freqbuf_tel.yaml"
    )

    params = simsetup._config_str_to_dict(param_filename)
    _, beams, _ = simsetup.initialize_uvdata_from_params(params, return_beams=True)

    assert beams[0].beam.freq_array.max() < 101e6


@pytest.mark.filterwarnings("ignore:Some Stokes I values are negative")
def test_simsetup_with_cached_catalog():
    # url of catalog file to download and cache, and its filetype
    url = "https://repository.library.brown.edu/storage/bdr:eafzyycj/content/"
    filetype = "skyh5"

    # If file already cached then `cache=True` will not redownload it.
    # Returns path to downloaded file content.
    filename = download_file(url, cache=True, pkgname="pyuvsim")

    # check that file exists in cache
    assert url in get_cached_urls("pyuvsim")
    # use filename to prevent warning
    assert os.path.basename(filename) == "contents"

    # create a source dictionary
    source_dict = {"catalog": url, "filetype": filetype}

    # try to initialize catalog using the cached catalog url
    catalog_uv, catname = simsetup.initialize_catalog_from_params(
        {"sources": source_dict}, return_catname=True
    )

    # first assert should be true for all astropy downloaded files
    # second and third assert specific to the downloaded healpix map
    assert catname == "contents"
    assert catalog_uv.Ncomponents == 196608
    assert catalog_uv.freq_array == [100 * units.MHz]
