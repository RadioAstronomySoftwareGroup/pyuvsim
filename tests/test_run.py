# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import re
import shutil
import subprocess  # nosec

import numpy as np
import pyradiosky
import pytest
import pyuvdata
import pyuvdata.utils.history as history_utils
import yaml
from astropy import units
from astropy.constants import c as speed_of_light
from astropy.coordinates import EarthLocation, Latitude, Longitude
from astropy.time import Time
from pyradiosky.utils import jy_to_ksr
from pyuvdata import ShortDipoleBeam, UniformBeam, UVData

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamList

c_ms = speed_of_light.to("m/s").value

pytest.importorskip("mpi4py")  # noqa


@pytest.mark.filterwarnings("ignore:Fixing auto-correlations to be be real-only")
@pytest.mark.parametrize(
    "paramfile", ["param_1time_1src_testcat.yaml", "param_1time_1src_testvot.yaml"]
)
@pytest.mark.parametrize("use_uvdata", [True, False])
@pytest.mark.skipif(
    not pytest.pyuvsim_can_parallel,
    reason="mpi-pytest is not installed. Cannot run parallel tests.",
)
@pytest.mark.parallel(2)
def test_run_one_source(goto_tempdir, paramfile, use_uvdata):
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    ref_file = os.path.join(SIM_DATA_PATH, "testfile_singlesource.uvh5")
    uv_ref = UVData.from_file(ref_file)

    param_filename = os.path.join(SIM_DATA_PATH, "test_config", paramfile)
    path = goto_tempdir
    ofilepath = os.path.join(path, "tempfile.uvh5")

    if use_uvdata:
        # only run most of the set up on rank 0
        if pyuvsim.mpi.rank == 0:
            uvd, beam_list, beam_dict = pyuvsim.simsetup.initialize_uvdata_from_params(
                param_filename, return_beams=True
            )
        else:
            uvd = None
            beam_dict = None
            beam_list = None

        # initialize the catalog on all nodes. For non-zero ranks, this skymodel
        # is made special so that it will break if line 877 in uvsim.py does not exist
        catalog = pyuvsim.simsetup.initialize_catalog_from_params(param_filename)
        # This test obsparam file has "single_source.txt" as its catalog.

        uv_out = pyuvsim.run_uvdata_uvsim(
            input_uv=uvd, beam_list=beam_list, beam_dict=beam_dict, catalog=catalog
        )
        if pyuvsim.mpi.rank == 0:
            # write the file out so we can make the rest of the
            # test work together nicely
            uv_out.write_uvh5(ofilepath, clobber=True)
    else:
        pyuvsim.uvsim.run_uvsim(param_filename)

    # Loading the file and comparing is only done on rank 0.
    if pyuvsim.mpi.rank != 0:
        return

    uv_new = UVData.from_file(ofilepath)

    # harmonize the phase centers
    uv_new._consolidate_phase_center_catalogs(other=uv_ref, ignore_name=True)

    assert history_utils._check_history_version(uv_new.history, pyradiosky.__version__)
    assert history_utils._check_history_version(uv_new.history, pyuvdata.__version__)
    assert history_utils._check_history_version(uv_new.history, pyuvsim.__version__)
    assert history_utils._check_history_version(uv_new.history, paramfile)
    assert history_utils._check_history_version(
        uv_new.history, "triangle_bl_layout.csv"
    )
    assert history_utils._check_history_version(
        uv_new.history, "28m_triangle_10time_10chan.yaml"
    )
    assert history_utils._check_history_version(uv_new.history, "Npus =")

    # Reset history and extra_keywords because they will deviate
    uv_new.history = uv_ref.history
    uv_new.extra_keywords = uv_ref.extra_keywords
    uv_ref.telescope.mount_type = uv_new.telescope.mount_type
    uv_new.filename = uv_ref.filename

    assert uv_new == uv_ref


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.filterwarnings("ignore:antenna_diameters are not set")
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
def test_analytic_diffuse(model, tol, tmpdir):
    # Generate the given model and simulate for a few baselines.
    # Import from analytic_diffuse  (consider moving to rasg_affiliates?)
    pytest.importorskip("analytic_diffuse")
    pytest.importorskip("astropy_healpix")
    import analytic_diffuse

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
        "beam_paths": {0: UniformBeam()},
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

    with open(template_path) as yfile:
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

    uv_out = pyuvsim.run_uvsim(obspar_path, return_uv=True)
    # Convert from Jy to K sr
    dat = uv_out.data_array[:, 0, 0] * jy_to_ksr(uv_out.freq_array[0]).value
    # Evaluate the solution and compare to visibilities.
    soln = analytic_diffuse.get_solution(modname)
    uvw_lam = uv_out.uvw_array * uv_out.freq_array[0] / c_ms
    ana = soln(uvw_lam, **params)
    np.testing.assert_allclose(ana / 2, dat, atol=tol, rtol=0)


@pytest.mark.filterwarnings("ignore:antenna_diameters are not set")
def test_diffuse_units(tmpdir):
    pytest.importorskip("analytic_diffuse")
    pytest.importorskip("astropy_healpix")

    sky_k, _, _, _, _ = pyuvsim.simsetup._create_catalog_diffuse(
        map_nside=128,
        diffuse_model="polydome",
        diffuse_params={"n": 2},
        time=Time(2457458.1738949567, format="jd"),
        array_location=EarthLocation.from_geodetic(
            lat=-30.72153, lon=21.42830, height=1073.0
        ),
    )

    sky_k.spectral_type = "spectral_index"
    sky_k.reference_frequency = np.full(sky_k.Ncomponents, 1e8) * units.Hz
    sky_k.spectral_index = np.full(sky_k.Ncomponents, -0.8)
    sky_k.check()

    sky_jy = sky_k.copy()
    sky_jy.kelvin_to_jansky()
    assert sky_jy.stokes.unit == units.Unit("Jy/sr")

    catpath_k = str(tmpdir.join("diffuse_units_k.skyh5"))
    sky_k.write_skyh5(catpath_k)

    catpath_jy = str(tmpdir.join("diffuse_units_jy.skyh5"))
    sky_jy.write_skyh5(catpath_jy)

    # Making configuration files for this simulation.
    template_path = os.path.join(
        SIM_DATA_PATH, "test_config", "obsparam_diffuse_sky.yaml"
    )
    obspar_path = str(tmpdir.join("obsparam_diffuse_sky.yaml"))
    layout_path = str(tmpdir.join("threeant_layout.csv"))
    herauniform_path = str(tmpdir.join("hera_uniform.yaml"))

    teleconfig = {
        "beam_paths": {0: UniformBeam()},
        "telescope_location": "(-30.72153, 21.42830, 1073.0)",
        "telescope_name": "HERA",
    }
    antpos_enu = np.array([[0, 0, 0], [0, 3, 0], [5, 0, 0]], dtype=float)

    pyuvsim.simsetup._write_layout_csv(
        layout_path, antpos_enu, np.arange(3).astype(str), np.arange(3)
    )
    with open(herauniform_path, "w") as ofile:
        yaml.dump(teleconfig, ofile, default_flow_style=False)

    with open(template_path) as yfile:
        obspar = yaml.safe_load(yfile)
    obspar["telescope"]["array_layout"] = layout_path
    obspar["telescope"]["telescope_config_name"] = herauniform_path
    obspar["sources"] = {"catalog": catpath_k}

    obspar["filing"]["outfile_name"] = "diffuse_units_k_sim.uvh5"
    obspar["filing"]["output_format"] = "uvh5"
    obspar["filing"]["outdir"] = str(tmpdir)

    with open(obspar_path, "w") as ofile:
        yaml.dump(obspar, ofile, default_flow_style=False)

    uv_out_k = pyuvsim.run_uvsim(obspar_path, return_uv=True)

    obspar["sources"] = {"catalog": catpath_jy}
    obspar["filing"]["outfile_name"] = "diffuse_units_jy_sim.uvh5"
    with open(obspar_path, "w") as ofile:
        yaml.dump(obspar, ofile, default_flow_style=False)

    uv_out_jy = pyuvsim.run_uvsim(obspar_path, return_uv=True)

    # make the histories the same for comparison
    uv_out_jy.history = uv_out_k.history
    # harmonize phase center catalog info
    uv_out_jy._consolidate_phase_center_catalogs(other=uv_out_k, ignore_name=True)

    assert uv_out_k == uv_out_jy


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
        new_beam_file = os.path.join(tmp_path, "test_config", "HERA_NicCST.beamfits")
        shutil.copyfile(beamfits_file, new_beam_file)

        params = pyuvsim.simsetup._config_str_to_dict(new_param_file)
    else:
        params = pyuvsim.simsetup._config_str_to_dict(param_file)

    pyuvsim.run_uvsim(params, return_uv=True)


@pytest.mark.parametrize("spectral_type", ["flat", "subband", "spectral_index"])
@pytest.mark.parametrize("nfeeds", [1, 2])
def test_run_gleam_uvsim(spectral_type, nfeeds):
    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, "test_config", "param_1time_testgleam.yaml")
    )
    params["sources"]["spectral_type"] = spectral_type
    params["sources"].pop("min_flux")
    params["sources"].pop("max_flux")

    if nfeeds == 1:
        if spectral_type == "subband":
            tel_config = "28m_triangle_10time_10chan_yfeed.yaml"
            select_pol = ["yy"]
        else:
            tel_config = "28m_triangle_10time_10chan_xfeed.yaml"
            select_pol = ["xx"]
        params["telescope"]["telescope_config_name"] = tel_config

    input_uv, beam_list, _ = pyuvsim.simsetup.initialize_uvdata_from_params(
        params, return_beams=True
    )

    assert input_uv.Npols == nfeeds**2

    for beam_inter in beam_list:
        assert beam_inter.Nfeeds == nfeeds

    uv_out = pyuvsim.run_uvsim(params, return_uv=True)
    assert uv_out.telescope.name == "Triangle"

    file_name = f"gleam_triangle_{spectral_type}.uvh5"
    uv_in = UVData.from_file(os.path.join(SIM_DATA_PATH, file_name))
    uv_in.rename_phase_center(0, "unprojected")
    uv_in.conjugate_bls()
    uv_in.reorder_blts()
    uv_in.integration_time = np.full_like(uv_in.integration_time, 11.0)

    assert uv_out.Npols == nfeeds**2

    if nfeeds == 1:
        assert uv_out.Npols == 1
        select_feed = select_pol[0][0]
        uv_in.select(polarizations=select_pol)
        if hasattr(uv_in.telescope, "feed_array"):
            feed_index = np.nonzero(uv_in.telescope.feed_array[0] == select_feed)[0]
            uv_in.telescope._select_along_param_axis({"Nfeeds": feed_index})

    # This just tests that we get the same answer as an earlier run, not that
    # the data are correct (that's covered in other tests)
    uv_out.history = uv_in.history
    uv_in.extra_keywords = uv_out.extra_keywords
    uv_in.telescope.mount_type = uv_out.telescope.mount_type
    assert uv_in.telescope._location == uv_out.telescope._location
    assert uv_in == uv_out


@pytest.mark.parametrize(
    ["spectral_type", "use_cli", "use_uvdata"],
    (
        ["subband", False, False],
        ["subband", True, True],
        ["spectral_index", True, False],
        ["spectral_index", False, True],
    ),
)
def test_zenith_spectral_sim(spectral_type, tmpdir, use_cli, use_uvdata):
    # Make a power law source at zenith in three ways.
    # Confirm that simulated visibilities match expectation.

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml")
    )

    alpha = -0.5
    ref_freq = 111e6
    Nfreqs = 20
    freqs, ch_width = np.linspace(110e6, 115e6, Nfreqs, retstep=True)
    freq_params = pyuvsim.simsetup.freq_array_to_params(
        freq_array=freqs, channel_width=ch_width
    )
    freq_params["Nspws"] = 1  # added to testing some logic to handle this.

    parsed_freq_params = pyuvsim.simsetup.parse_frequency_params(freq_params)
    freqs = parsed_freq_params["freq_array"]
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
        freq_lower = freqs - parsed_freq_params["channel_width"] * units.Hz / 2.0
        freq_upper = freqs + parsed_freq_params["channel_width"] * units.Hz / 2.0
        source.freq_edge_array = np.concatenate(
            (freq_lower[np.newaxis, :], freq_upper[np.newaxis, :]), axis=0
        )
        source.Nfreqs = Nfreqs
        source.freq_array = freqs
        source.stokes = np.repeat(source.stokes, Nfreqs, axis=1)
        source.stokes[0, :, 0] *= spectrum

    catpath = str(tmpdir.join("spectral_test_catalog.skyh5"))
    source.write_skyh5(catpath)
    params["sources"] = {"catalog": catpath}
    params["filing"]["outdir"] = str(tmpdir)
    params["filing"]["outfile_name"] = "zenith_spectral.uvh5"
    params["freq"] = freq_params
    new_time_params = params["time"]
    new_time_params["start_time"] = kwds["time"]
    new_time_params = pyuvsim.simsetup.parse_time_params(new_time_params)
    params["time"] = pyuvsim.simsetup.time_array_to_params(
        time_array=new_time_params["time_array"],
        integration_time=new_time_params["integration_time"],
    )
    params["select"] = {"antenna_nums": [1, 2]}

    new_telescope_param_file = str(
        tmpdir / params["telescope"]["telescope_config_name"]
    )
    shutil.copyfile(
        os.path.join(
            SIM_DATA_PATH, "test_config", params["telescope"]["telescope_config_name"]
        ),
        new_telescope_param_file,
    )

    new_telescope_layout_file = str(tmpdir / params["telescope"]["array_layout"])
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", params["telescope"]["array_layout"]),
        new_telescope_layout_file,
    )

    param_file = str(tmpdir / "zenith_spectral.yaml")
    with open(param_file, "w") as pfile:
        yaml.dump(params, pfile, default_flow_style=False)

    if use_uvdata:
        uvd, beam_list, beam_dict = pyuvsim.simsetup.initialize_uvdata_from_params(
            param_file, return_beams=True
        )

        catalog = pyuvsim.simsetup.initialize_catalog_from_params(param_file)

    if use_cli:
        uvd_out = str(tmpdir / params["filing"]["outfile_name"])
        if use_uvdata:
            uvd_in = str(tmpdir / "zenith_spectral_input.uvh5")
            pyuvsim.simsetup._complete_uvdata(uvd, inplace=True)
            uvd.write_uvh5(uvd_in)
            cat_in = str(tmpdir / "zenith_spectral_input.skyh5")
            catalog.write_skyh5(cat_in)
            uvb_in = str(tmpdir / "uniform_beam.beamfits")
            abeam = UniformBeam()
            uvb = abeam.to_uvbeam(
                freq_array=np.linspace(100e6, 120e6, 4),
                axis1_array=np.linspace(0, 2 * np.pi, 36),
                axis2_array=np.linspace(0, np.pi / 2, 18),
            )
            uvb.write_beamfits(uvb_in)

            subprocess.check_output(  # nosec
                [
                    "run_pyuvsim",
                    "--uvdata",
                    uvd_in,
                    "--uvbeam",
                    uvb_in,
                    "--skymodel",
                    cat_in,
                    "--outfile",
                    uvd_out,
                ]
            )
        else:
            profile_file = str(tmpdir / "profile.out")
            profile_file_meta = str(tmpdir / "profile_meta.out")
            subprocess.check_output(  # nosec
                [
                    "run_param_pyuvsim",
                    param_file,
                    "--profile",
                    profile_file,
                    "--quiet",
                    "--keep_nonroot_stdout",
                    "--raw_profile",
                ]
            )
            assert os.path.exists(profile_file)
            assert os.path.exists(profile_file_meta)
        uv_out = UVData.from_file(uvd_out)
    else:
        if use_uvdata:
            uv_out = pyuvsim.run_uvdata_uvsim(
                input_uv=uvd, beam_list=beam_list, beam_dict=beam_dict, catalog=catalog
            )
        else:
            uv_out = pyuvsim.run_uvsim(param_file, return_uv=True)

    if use_cli and use_uvdata:
        # this is different because we are writing out the analytic beam to a
        # UVBeam, which then gets peak normalized.
        exp_spectrum = spectrum
    else:
        exp_spectrum = spectrum / 2
    for ii in range(uv_out.Nbls):
        np.testing.assert_allclose(uv_out.data_array[ii, :, 0], exp_spectrum)


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    hera_uv = UVData()

    hera_uv.polarization_array = np.asarray([-5])
    beam_list = BeamList([ShortDipoleBeam()])

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Input UVData object/simulation parameters and beams polarizations "
            "do not agree. Input beams have output polarizations: "
            "['xx', 'yy', 'xy', 'yx'], Simulation has expected polarizations ['xx']"
        ),
    ):
        pyuvsim.run_uvdata_uvsim(hera_uv, beam_list, {}, catalog=pyuvsim.SkyModelData())


@pytest.mark.filterwarnings("ignore:Setting the location attribute post initialization")
@pytest.mark.parametrize("selenoid", ["SPHERE", "GSFC", "GRAIL23", "CE-1-LAM-GEO"])
def test_sim_on_moon(goto_tempdir, selenoid):
    pytest.importorskip("lunarsky")
    from spiceypy.utils.exceptions import SpiceUNKNOWNFRAME

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

    with open(os.path.join(tmpdir, "tranquility_config.yaml")) as pfile:
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
    assert uv_obj.telescope.location.ellipsoid == selenoid

    # set the filename to make sure it ends up in the history,
    # remove the parameter file info from extra_keywords
    uv_obj.filename = ["moon_sim"]
    uv_obj._filename.form = (1,)
    uv_obj.extra_keywords.pop("obsparam")
    uv_obj.extra_keywords.pop("telecfg")
    uv_obj.extra_keywords.pop("layout")
    uv_obj.check()

    uv_obj.select(times=uv_obj.time_array[0])

    tranquility_base = uv_obj.telescope.location

    time = Time(uv_obj.time_array[0], format="jd", scale="utc")
    sources, _ = pyuvsim.create_mock_catalog(
        time,
        array_location=tranquility_base,
        arrangement="zenith",
        Nsrcs=30,
        return_data=True,
    )
    # Run simulation.
    try:
        uv_out = pyuvsim.uvsim.run_uvdata_uvsim(
            uv_obj, beam_list, beam_dict, catalog=sources, quiet=True
        )
    except SpiceUNKNOWNFRAME as err:
        pytest.skip("SpiceUNKNOWNFRAME error: " + str(err))

    assert history_utils._check_history_version(uv_out.history, pyradiosky.__version__)
    assert history_utils._check_history_version(uv_out.history, pyuvdata.__version__)
    assert history_utils._check_history_version(uv_out.history, pyuvsim.__version__)
    assert history_utils._check_history_version(uv_out.history, uv_obj.filename[0])
    assert history_utils._check_history_version(uv_out.history, "Npus =")

    assert uv_out.extra_keywords["world"] == "moon"
    assert uv_out.telescope.location.ellipsoid == selenoid
    assert np.allclose(uv_out.data_array[:, :, 0], 0.5)

    # Lunar Frame Roundtripping
    param_dict["filing"]["outdir"] = str(tmpdir)
    uv_filename = pyuvsim.utils.write_uvdata(
        uv_out, param_dict, return_filename=True, quiet=True
    )
    uv_compare = UVData()
    uv_compare.read(uv_filename)
    assert uv_out.telescope._location == uv_compare.telescope._location


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

    with open(os.path.join(tmpdir, "test_config", "bl_single_gauss.yaml")) as pfile:
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

    assert uv_out.telescope.location.ellipsoid == selenoid

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

    pos = uv_out.telescope.location_lat_lon_alt_degrees

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
