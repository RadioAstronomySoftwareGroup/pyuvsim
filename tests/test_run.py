# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import copy
import os
import shutil

import numpy as np
import pyradiosky
import pytest
import pyuvdata
import pyuvdata.utils as uvutils
import yaml
from astropy import units
from astropy.coordinates import Latitude, Longitude
from astropy.time import Time
from pyradiosky.utils import jy_to_ksr, stokes_to_coherency
from pyuvdata import UVData

import pyuvsim
from pyuvsim.analyticbeam import c_ms
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamList

try:
    from pyuvdata.testing import check_warnings
except ImportError:
    # this can be removed once we require pyuvdata >= v3.0
    from pyuvdata.tests import check_warnings

pytest.importorskip("mpi4py")  # noqa
# everything in this file requires 2 PUs
pytestmark = pytest.mark.parallel(2)

future_shapes_options = [True]
if hasattr(UVData(), "use_current_array_shapes"):
    future_shapes_options += [False]


@pytest.fixture
def goto_tempdir(tmpdir):
    # Run test within temporary directory.
    newpath = str(tmpdir)
    cwd = os.getcwd()
    os.chdir(newpath)

    yield newpath

    os.chdir(cwd)


@pytest.mark.filterwarnings("ignore:antenna_diameters are not set")
@pytest.mark.filterwarnings("ignore:Telescope Triangle is not in known_telescopes.")
@pytest.mark.filterwarnings("ignore:Fixing auto-correlations to be be real-only")
@pytest.mark.parametrize(
    "paramfile", ["param_1time_1src_testcat.yaml", "param_1time_1src_testvot.yaml"]
)
def test_run_paramfile_uvsim(goto_tempdir, paramfile):
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    uv_ref = UVData()
    uv_ref.read(os.path.join(SIM_DATA_PATH, "testfile_singlesource.uvh5"))
    if hasattr(uv_ref, "use_current_array_shapes"):
        uv_ref.use_future_array_shapes()

    param_filename = os.path.join(SIM_DATA_PATH, "test_config", paramfile)
    # This test obsparam file has "single_source.txt" as its catalog.
    pyuvsim.uvsim.run_uvsim(param_filename)

    # Loading the file and comparing is only done on rank 0.
    if pyuvsim.mpi.rank != 0:
        return

    path = goto_tempdir
    ofilepath = os.path.join(path, "tempfile.uvfits")

    uv_new = UVData.from_file(ofilepath)
    if hasattr(uv_new, "use_current_array_shapes"):
        uv_new.use_future_array_shapes()

    uv_new.unproject_phase(use_ant_pos=True)
    uv_new._consolidate_phase_center_catalogs(other=uv_ref)

    if pyuvsim.mpi.rank == 0:
        assert uvutils._check_history_version(uv_new.history, pyradiosky.__version__)
        assert uvutils._check_history_version(uv_new.history, pyuvdata.__version__)
        assert uvutils._check_history_version(uv_new.history, pyuvsim.__version__)
        assert uvutils._check_history_version(uv_new.history, paramfile)
        assert uvutils._check_history_version(uv_new.history, "triangle_bl_layout.csv")
        assert uvutils._check_history_version(
            uv_new.history, "28m_triangle_10time_10chan.yaml"
        )
        assert uvutils._check_history_version(uv_new.history, "Npus =")

        # Reset parts that will deviate
        uv_new.history = uv_ref.history
        uv_new.object_name = uv_ref.object_name
        uv_ref.dut1 = uv_new.dut1
        uv_ref.gst0 = uv_new.gst0
        uv_ref.rdate = uv_new.rdate

        uv_ref.conjugate_bls()
        uv_ref.reorder_blts()
        uv_ref.integration_time = np.full_like(uv_ref.integration_time, 11.0)

        # remove filename attribute to ensure equality
        uv_new.filename = None
        uv_ref.filename = None
        assert uv_new == uv_ref


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
# Set the tolerances as low as we can achieve currently. Ideally these tolerances
# would be lower, but it's complicated.
# See Lanman, Murray and Jacobs, 2022, DOI: 10.3847/1538-4365/ac45fd
@pytest.mark.parametrize(
    ("model", "tol"),
    [
        ("monopole", 3e-4),
        ("cosza", 2e-4),
        ("quaddome", 8e-5),
        ("monopole-nonflat", 3e-4),
    ],
)
@pytest.mark.parametrize("backend", ["rma", "send_recv"])
@pytest.mark.parametrize("progbar", ["progsteps", "tqdm"])
@pytest.mark.parallel(2, timeout=90)
def test_analytic_diffuse(model, tol, tmpdir, backend, progbar):
    # Generate the given model and simulate for a few baselines.
    # Import from analytic_diffuse  (consider moving to rasg_affiliates?)
    analytic_diffuse = pytest.importorskip("analytic_diffuse")
    pytest.importorskip("astropy_healpix")
    if progbar == "tqdm":
        pytest.importorskip("tqdm")

    modname = model
    use_w = False
    params = {}
    if model == "quaddome":
        modname = "polydome"
        params["n"] = 2
    elif model == "monopole-nonflat":
        modname = "monopole"
        use_w = True
        params["order"] = 50  # Expansion order for the non-flat monopole solution.

    # Making configuration files for this simulation.
    template_path = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_diffuse_sky.yaml"
    )
    obspar_path = str(tmpdir.join("obsparam_diffuse_sky.yaml"))
    layout_path = str(tmpdir.join("threeant_layout.csv"))
    herauniform_path = str(tmpdir.join("hera_uniform.yaml"))

    teleconfig = {
        "beam_paths": {0: {"type": "uniform"}},
        "telescope_location": "(-30.72153, 21.42830, 1073.0)",
        "telescope_name": "HERA",
    }
    if not use_w:
        antpos_enu = np.array([[0, 0, 0], [0, 3, 0], [5, 0, 0]], dtype=float)
    else:
        antpos_enu = np.array([[0, 0, 0], [0, 3, 0], [0, 3, 5]], dtype=float)

    pyuvsim.simsetup._write_layout_csv(
        layout_path, antpos_enu, np.arange(3).astype(str), np.arange(3)
    )
    with open(herauniform_path, "w") as ofile:
        yaml.dump(teleconfig, ofile, default_flow_style=False)

    with open(template_path, "r") as yfile:
        obspar = yaml.safe_load(yfile)
    obspar["telescope"]["array_layout"] = layout_path
    obspar["telescope"]["telescope_config_name"] = herauniform_path
    obspar["sources"]["diffuse_model"] = modname
    obspar["sources"].update(params)
    if model == "monopole":
        # use a higher nside for monopole to improve the accuracy
        obspar["sources"]["map_nside"] = 256
    obspar["filing"]["outfile_name"] = "diffuse_sim.uvh5"
    obspar["filing"]["output_format"] = "uvh5"
    obspar["filing"]["outdir"] = str(tmpdir)

    with open(obspar_path, "w") as ofile:
        yaml.dump(obspar, ofile, default_flow_style=False)

    uv_out = pyuvsim.run_uvsim(
        obspar_path, return_uv=True, backend=backend, progbar=progbar
    )
    if pyuvsim.mpi.rank == 0:
        # Convert from Jy to K sr
        dat = uv_out.data_array[:, 0, 0] * jy_to_ksr(uv_out.freq_array[0]).value
        # Evaluate the solution and compare to visibilities.
        soln = analytic_diffuse.get_solution(modname)
        uvw_lam = uv_out.uvw_array * uv_out.freq_array[0] / c_ms
        ana = soln(uvw_lam, **params)
        np.testing.assert_allclose(ana / 2, dat, atol=tol, rtol=0)


@pytest.mark.filterwarnings("ignore:Fixing auto polarization power beams")
def test_powerbeam_sim(cst_beam):
    new_cst = copy.deepcopy(cst_beam)
    new_cst.efield_to_power()
    beams = BeamList([new_cst] * 4)
    cfg = os.path.join(SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml")
    input_uv = pyuvsim.simsetup.initialize_uvdata_from_params(cfg, return_beams=False)
    sky_model = pyuvsim.simsetup.initialize_catalog_from_params(
        cfg, return_catname=False
    )
    sky_model = pyuvsim.simsetup.SkyModelData(sky_model)
    beam_dict = {n: 0 for n in range(4)}

    with pytest.raises(ValueError, match="Beam type must be efield!"):
        pyuvsim.run_uvdata_uvsim(input_uv, beams, beam_dict, catalog=sky_model)


@pytest.mark.parametrize("rename_beamfits", [True, False])
def test_run_paramdict_uvsim(rename_beamfits, tmp_path):
    # Running a simulation from parameter dictionary.
    param_file = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml"
    )

    if rename_beamfits:
        os.makedirs(os.path.join(tmp_path, "test_config"))
        new_param_file = os.path.join(
            tmp_path, "test_config", "param_1time_1src_testcat.yaml"
        )
        shutil.copyfile(param_file, new_param_file)

        telescope_param_file = os.path.join(
            SIM_DATA_PATH, "test_config", "28m_triangle_10time_10chan.yaml"
        )
        new_telescope_param_file = os.path.join(
            tmp_path, "test_config", "28m_triangle_10time_10chan.yaml"
        )
        shutil.copyfile(telescope_param_file, new_telescope_param_file)

        telescope_layout_file = os.path.join(
            SIM_DATA_PATH, "test_config", "triangle_bl_layout.csv"
        )
        new_telescope_layout_file = os.path.join(
            tmp_path, "test_config", "triangle_bl_layout.csv"
        )
        shutil.copyfile(telescope_layout_file, new_telescope_layout_file)

        source_file = os.path.join(SIM_DATA_PATH, "single_source.txt")
        new_source_file = os.path.join(tmp_path, "single_source.txt")
        shutil.copyfile(source_file, new_source_file)

        beamfits_file = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")
        new_beam_file = os.path.join(tmp_path, "test_config", "HERA_NicCST.uvbeam")
        shutil.copyfile(beamfits_file, new_beam_file)

        # change the beam file name to .uvbeam
        with open(new_telescope_param_file, "r") as pfile:
            tele_param_dict = yaml.safe_load(pfile)
            tele_param_dict["beam_paths"][0] = {"filename": new_beam_file}

        with open(new_telescope_param_file, "w") as yfile:
            yaml.dump(tele_param_dict, yfile, default_flow_style=False)

        n_beam_warnings = 3
        warn_type = [DeprecationWarning] * n_beam_warnings
        msg = [pyuvsim.telescope.weird_beamfits_extension_warning] * n_beam_warnings

        params = pyuvsim.simsetup._config_str_to_dict(new_param_file)
    else:
        warn_type = None
        msg = ""
        params = pyuvsim.simsetup._config_str_to_dict(param_file)

    if pyuvsim.mpi.rank > 0:
        with check_warnings(warn_type, match=msg):
            pyuvsim.run_uvsim(params, return_uv=True)
    else:
        pyuvsim.run_uvsim(params, return_uv=True)


@pytest.mark.filterwarnings("ignore:Telescope Triangle is not in known_telescopes.")
@pytest.mark.filterwarnings("ignore:The shapes of several attributes will be changing")
@pytest.mark.parametrize("spectral_type", ["flat", "subband", "spectral_index"])
def test_run_gleam_uvsim(spectral_type):
    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, "test_config", "param_1time_testgleam.yaml")
    )
    params["sources"]["spectral_type"] = spectral_type
    params["sources"].pop("min_flux")
    params["sources"].pop("max_flux")

    uv_out = pyuvsim.run_uvsim(params, return_uv=True)

    if pyuvsim.mpi.rank == 0:
        assert uv_out.telescope_name == "Triangle"

        file_name = f"gleam_triangle_{spectral_type}.uvh5"
        uv_in = UVData.from_file(os.path.join(SIM_DATA_PATH, file_name))
        # This can be removed when we require pyuvdata >= 3.0
        if hasattr(uv_in, "use_current_array_shapes"):
            uv_in.use_future_array_shapes()
        uv_in.conjugate_bls()
        uv_in.reorder_blts()
        uv_in.integration_time = np.full_like(uv_in.integration_time, 11.0)
        # This just tests that we get the same answer as an earlier run, not that
        # the data are correct (that's covered in other tests)
        uv_out.history = uv_in.history
        if hasattr(uv_out, "telescope"):
            assert uv_in.telescope._location == uv_out.telescope._location
        else:
            # this can be removed when we require pyuvdata >= 3.0
            assert uv_in._telescope_location == uv_out._telescope_location
        assert uv_in == uv_out


@pytest.mark.filterwarnings("ignore:The reference_frequency is aliased as `frequency`")
@pytest.mark.parametrize("spectral_type", ["subband", "spectral_index"])
@pytest.mark.parametrize("backend", ["rma", "send_recv"])
@pytest.mark.parametrize("progbar", ["progsteps", "tqdm"])
def test_zenith_spectral_sim(spectral_type, tmpdir, backend, progbar):
    # Make a power law source at zenith in three ways.
    # Confirm that simulated visibilities match expectation.
    if progbar == "tqdm":
        pytest.importorskip("tqdm")

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml")
    )

    alpha = -0.5
    ref_freq = 111e6
    Nfreqs = 20
    freqs = np.linspace(110e6, 115e6, Nfreqs)
    freq_params = pyuvsim.simsetup.freq_array_to_params(freqs)
    freqs = pyuvsim.simsetup.parse_frequency_params(freq_params)["freq_array"]
    freqs *= units.Hz
    spectrum = (freqs.value / ref_freq) ** alpha

    source, kwds = pyuvsim.create_mock_catalog(
        Time.now(), arrangement="zenith", Nsrcs=1
    )
    source.spectral_type = spectral_type
    if spectral_type == "spectral_index":
        source.reference_frequency = np.array([ref_freq]) * units.Hz
        source.spectral_index = np.array([alpha])
    else:
        freq_lower = freqs - freq_params["channel_width"] * units.Hz / 2.0
        freq_upper = freqs + freq_params["channel_width"] * units.Hz / 2.0
        source.freq_edge_array = np.concatenate(
            (freq_lower[np.newaxis, :], freq_upper[np.newaxis, :]), axis=0
        )
        source.Nfreqs = Nfreqs
        source.freq_array = freqs
        source.stokes = np.repeat(source.stokes, Nfreqs, axis=1)
        source.stokes[0, :, 0] *= spectrum
        source.coherency_radec = stokes_to_coherency(source.stokes)

    catpath = str(tmpdir.join("spectral_test_catalog.skyh5"))
    source.write_skyh5(catpath)
    params["sources"] = {"catalog": catpath}
    params["filing"]["outdir"] = str(tmpdir)
    params["freq"] = freq_params
    params["time"]["start_time"] = kwds["time"]
    params["select"] = {"antenna_nums": [1, 2]}

    uv_out = pyuvsim.run_uvsim(params, return_uv=True, backend=backend, progbar=progbar)
    if pyuvsim.mpi.rank == 0:
        for ii in range(uv_out.Nbls):
            assert np.allclose(uv_out.data_array[ii, :, 0], spectrum / 2)


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    hera_uv = UVData()

    hera_uv.polarizations = ["xx"]

    with pytest.raises(ValueError, match="input_uv must have XX,YY,XY,YX polarization"):
        pyuvsim.run_uvdata_uvsim(
            hera_uv, ["beamlist"], {}, catalog=pyuvsim.SkyModelData()
        )


def test_input_uv_error():
    with pytest.raises(TypeError, match="input_uv must be UVData object"):
        pyuvsim.run_uvdata_uvsim(None, None, {}, catalog=pyuvsim.SkyModelData())


# several of these filters should be removed once we require pyuvdata>=3.0
@pytest.mark.filterwarnings("ignore:Setting the location attribute post initialization")
@pytest.mark.filterwarnings("ignore:This method will be removed in version 3.0")
@pytest.mark.filterwarnings("ignore:The lst_array is not self-consistent")
@pytest.mark.filterwarnings("ignore:Telescope apollo11 is not in known_telescopes.")
@pytest.mark.filterwarnings("ignore:The shapes of several attributes will be changing")
@pytest.mark.parametrize("future_shapes", future_shapes_options)
@pytest.mark.parametrize("selenoid", ["SPHERE", "GSFC", "GRAIL23", "CE-1-LAM-GEO"])
def test_sim_on_moon(future_shapes, goto_tempdir, selenoid):
    pytest.importorskip("lunarsky")
    from lunarsky import MoonLocation

    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_tranquility_hex.yaml"
    )

    tmpdir = goto_tempdir
    # copy the config files and modify the lunar ellipsoid
    os.makedirs(os.path.join(tmpdir, "test_config"))
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "obsparam_tranquility_hex.yaml"),
        os.path.join(tmpdir, "test_config", "obsparam_tranquility_hex.yaml"),
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "tranquility_config.yaml"),
        os.path.join(tmpdir, "tranquility_config.yaml"),
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "layout_hex37_14.6m.csv"),
        os.path.join(tmpdir, "test_config", "layout_hex37_14.6m.csv"),
    )

    with open(os.path.join(tmpdir, "tranquility_config.yaml"), "r") as pfile:
        tele_param_dict = yaml.safe_load(pfile)
    tele_param_dict["ellipsoid"] = selenoid
    with open(os.path.join(tmpdir, "tranquility_config.yaml"), "w") as pfile:
        yaml.dump(tele_param_dict, pfile, default_flow_style=False)

    param_filename = os.path.join(
        tmpdir, "test_config", "obsparam_tranquility_hex.yaml"
    )

    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    param_dict["select"] = {"redundant_threshold": 0.1}
    uv_obj, beam_list, beam_dict = pyuvsim.initialize_uvdata_from_params(
        param_dict, return_beams=True
    )
    if hasattr(uv_obj, "telescope"):
        assert uv_obj.telescope.location.ellipsoid == selenoid
    else:
        # this can be removed when we require pyuvdata >= 3.0
        assert uv_obj._telescope_location.ellipsoid == selenoid

    # set the filename to make sure it ends up in the history,
    # remove the parameter file info from extra_keywords
    uv_obj.filename = ["moon_sim"]
    uv_obj._filename.form = (1,)
    uv_obj.extra_keywords.pop("obsparam")
    uv_obj.extra_keywords.pop("telecfg")
    uv_obj.extra_keywords.pop("layout")
    uv_obj.check()

    uv_obj.select(times=uv_obj.time_array[0])

    if hasattr(uv_obj, "telescope"):
        tranquility_base = uv_obj.telescope.location
    else:
        # this can be removed when we require pyuvdata >= 3.0
        tranquility_base = MoonLocation.from_selenocentric(
            *uv_obj.telescope_location, "meter"
        )
        tranquility_base.ellipsoid = selenoid

    time = Time(uv_obj.time_array[0], format="jd", scale="utc")
    sources, _ = pyuvsim.create_mock_catalog(
        time,
        array_location=tranquility_base,
        arrangement="zenith",
        Nsrcs=30,
        return_data=True,
    )
    # Run simulation.
    if not future_shapes:
        uv_obj.use_current_array_shapes()
    uv_out = pyuvsim.uvsim.run_uvdata_uvsim(
        uv_obj, beam_list, beam_dict, catalog=sources, quiet=True
    )
    if pyuvsim.mpi.rank == 0:
        assert uvutils._check_history_version(uv_out.history, pyradiosky.__version__)
        assert uvutils._check_history_version(uv_out.history, pyuvdata.__version__)
        assert uvutils._check_history_version(uv_out.history, pyuvsim.__version__)
        assert uvutils._check_history_version(uv_out.history, uv_obj.filename[0])
        assert uvutils._check_history_version(uv_out.history, "Npus =")

        assert uv_out.extra_keywords["world"] == "moon"
        assert uv_out._telescope_location.ellipsoid == selenoid
        assert np.allclose(uv_out.data_array[:, :, 0], 0.5)

        # Lunar Frame Roundtripping
        param_dict["filing"]["outdir"] = str(tmpdir)

        uv_filename = pyuvsim.utils.write_uvdata(
            uv_out, param_dict, return_filename=True, quiet=True
        )
        uv_compare = UVData()
        uv_compare.read(uv_filename, use_future_array_shapes=future_shapes)
        assert np.allclose(uv_out.telescope_location, uv_compare.telescope_location)
        assert uv_out._telescope_location.frame == uv_compare._telescope_location.frame
        assert uv_compare._telescope_location.ellipsoid == selenoid

        # Cleanup
        os.remove(uv_filename)


@pytest.mark.filterwarnings("ignore:This method will be removed in version 3.0")
@pytest.mark.filterwarnings("ignore:The lst_array is not self-consistent")
@pytest.mark.parametrize("selenoid", ["SPHERE", "GSFC", "GRAIL23", "CE-1-LAM-GEO"])
def test_lunar_gauss(goto_tempdir, selenoid):
    pytest.importorskip("lunarsky")
    from lunarsky import MoonLocation
    from spiceypy.utils.exceptions import SpiceUNKNOWNFRAME

    # Make a gaussian source that passes through zenith
    # Confirm that simulated visibilities match expectation.

    tmpdir = goto_tempdir
    # copy the config files and modify the lunar ellipsoid
    os.makedirs(os.path.join(tmpdir, "test_config"))
    os.makedirs(os.path.join(tmpdir, "test_catalogs"))
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "obsparam_lunar_gauss.yaml"),
        os.path.join(tmpdir, "test_config", "obsparam_lunar_gauss.yaml"),
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "bl_single_gauss.yaml"),
        os.path.join(tmpdir, "test_config", "bl_single_gauss.yaml"),
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "baseline_moon.csv"),
        os.path.join(tmpdir, "test_config", "baseline_moon.csv"),
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_catalogs", "one_distant_point_2458178.5.txt"),
        os.path.join(tmpdir, "test_catalogs", "one_distant_point_2458178.5.txt"),
    )

    with open(
        os.path.join(tmpdir, "test_config", "bl_single_gauss.yaml"), "r"
    ) as pfile:
        tele_param_dict = yaml.safe_load(pfile)
    tele_param_dict["ellipsoid"] = selenoid
    with open(
        os.path.join(tmpdir, "test_config", "bl_single_gauss.yaml"), "w"
    ) as pfile:
        yaml.dump(tele_param_dict, pfile, default_flow_style=False)

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(tmpdir, "test_config", "obsparam_lunar_gauss.yaml")
    )

    params["filing"]["outdir"] = str(tmpdir)

    try:
        uv_out = pyuvsim.run_uvsim(params, return_uv=True, quiet=True)
    except SpiceUNKNOWNFRAME as err:
        pytest.skip("SpiceUNKNOWNFRAME error: " + str(err))

    if pyuvsim.mpi.rank == 0:
        assert uv_out._telescope_location.ellipsoid == selenoid
        # Skymodel and update positions

        # Init sky model
        sm = pyradiosky.SkyModel(
            name="source0",
            ra=Longitude(308.32686, unit="deg"),
            dec=Latitude(-21, unit="deg"),
            stokes=units.Quantity([1, 0, 0, 0], unit="Jy"),
            spectral_type="flat",
            frame="fk5",
        )

        if hasattr(uv_out, "telescope"):
            pos = uv_out.telescope.location_lat_lon_alt_degrees
        else:
            # this can be removed when we require pyuvdata >= 3.0
            pos = uv_out.telescope_location_lat_lon_alt_degrees

        # Creating the analytical gaussian
        Alt = np.zeros(uv_out.Ntimes)
        Az = np.zeros(uv_out.Ntimes)
        refTimes = uv_out.get_times(0, 1)
        telescope_location_obj = MoonLocation(
            Longitude(pos[1], unit="deg"),
            Latitude(pos[0], unit="deg"),
            units.Quantity(pos[2], unit="m"),
            ellipsoid=selenoid,
        )
        for t in range(uv_out.Ntimes):
            sm.update_positions(Time(refTimes[t], format="jd"), telescope_location_obj)
            Alt[t] = sm.alt_az[0, 0]
            Az[t] = sm.alt_az[1, 0]

        sigma = 0.5
        Vis = 0.5 * np.exp(-np.power((Alt - np.pi / 2) / sigma, 2))

        # Check that the analytical visibility agrees with the simulation
        assert np.allclose(
            Vis, np.abs(uv_out.get_data(0, 1)[:, 0, 0]), rtol=1e-04, atol=1e-04
        )
