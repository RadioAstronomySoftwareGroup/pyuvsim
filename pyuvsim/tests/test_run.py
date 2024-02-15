# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import copy
import os
import resource
import shutil
import sys

import numpy as np
import pyradiosky
import pytest
import pyuvdata
import pyuvdata.tests as uvtest
import pyuvdata.utils as uvutils
import yaml
from astropy import units
from astropy.time import Time
from packaging import version  # packaging is installed with setuptools
from pyradiosky.utils import jy_to_ksr, stokes_to_coherency
from pyuvdata import UVData

import pyuvsim
from pyuvsim.analyticbeam import c_ms
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamList

pytest.importorskip('mpi4py')  # noqa


@pytest.fixture()
def goto_tempdir(tmpdir):
    # Run test within temporary directory.
    newpath = str(tmpdir)
    cwd = os.getcwd()
    os.chdir(newpath)

    yield newpath

    os.chdir(cwd)


@pytest.mark.filterwarnings("ignore:antenna_diameters are not set")
@pytest.mark.parametrize('paramfile', ['param_1time_1src_testcat.yaml',
                                       'param_1time_1src_testvot.yaml'])
@pytest.mark.parallel(2)
def test_run_paramfile_uvsim(goto_tempdir, paramfile):
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    uv_ref = UVData()
    if version.parse(pyuvdata.__version__) > version.parse("2.2.12"):
        # if we can, use the new file that has many things fixed
        uv_ref.read(os.path.join(SIM_DATA_PATH, 'testfile_singlesource.uvh5'))
    else:
        uv_ref.read(os.path.join(SIM_DATA_PATH, 'testfile_singlesource.uvfits'))
        uv_ref.unphase_to_drift(use_ant_pos=True)
        uv_ref.reorder_blts()
        # This is an old file with the bug that added one to the
        # antenna numbers for uvfits files. Fix them (if pyuvdata is recent)
        if np.min(np.union1d(uv_ref.ant_1_array, uv_ref.ant_2_array)) > 0:
            uv_ref.ant_1_array = uv_ref.ant_1_array - 1
            uv_ref.ant_2_array = uv_ref.ant_2_array - 1
            uv_ref.antenna_numbers = uv_ref.antenna_numbers - 1
            uv_ref.baseline_array = uv_ref.antnums_to_baseline(
                uv_ref.ant_1_array, uv_ref.ant_2_array
            )

        # set the x_orientation
        uv_ref.x_orientation = "east"

        # fix the channel width, which doesn't match the channel width in the parameter file
        uv_ref.channel_width = 1000000.0

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', paramfile)
    # This test obsparam file has "single_source.txt" as its catalog.
    pyuvsim.uvsim.run_uvsim(param_filename)

    # Loading the file and comparing is only done on rank 0.
    if pyuvsim.mpi.rank != 0:
        return

    path = goto_tempdir
    ofilepath = os.path.join(path, 'tempfile.uvfits')

    uv_new = UVData()
    uv_new.read(ofilepath)

    if version.parse(pyuvdata.__version__) > version.parse("2.2.12"):
        uv_new.unproject_phase(use_ant_pos=True)
        uv_new._consolidate_phase_center_catalogs(other=uv_ref)
    else:
        uv_new.unphase_to_drift(use_ant_pos=True)

    assert uvutils._check_history_version(uv_new.history, pyradiosky.__version__)
    assert uvutils._check_history_version(uv_new.history, pyuvdata.__version__)
    assert uvutils._check_history_version(uv_new.history, pyuvsim.__version__)
    assert uvutils._check_history_version(uv_new.history, paramfile)
    assert uvutils._check_history_version(uv_new.history, 'triangle_bl_layout.csv')
    assert uvutils._check_history_version(uv_new.history, '28m_triangle_10time_10chan.yaml')
    assert uvutils._check_history_version(uv_new.history, "Npus =")

    # Reset parts that will deviate
    uv_new.history = uv_ref.history
    uv_new.object_name = uv_ref.object_name
    uv_ref.dut1 = uv_new.dut1
    uv_ref.gst0 = uv_new.gst0
    uv_ref.rdate = uv_new.rdate

    assert uv_new == uv_ref


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
@pytest.mark.parametrize('model', ['monopole', 'cosza', 'quaddome', 'monopole-nonflat'])
def test_analytic_diffuse(model, tmpdir):
    # Generate the given model and simulate for a few baselines.
    # Import from analytic_diffuse  (consider moving to rasg_affiliates?)
    pytest.importorskip('analytic_diffuse')
    pytest.importorskip('astropy_healpix')
    import analytic_diffuse

    modname = model
    use_w = False
    params = {}
    if model == 'quaddome':
        modname = 'polydome'
        params['n'] = 2
    elif model == 'monopole-nonflat':
        modname = 'monopole'
        use_w = True
        params['order'] = 30    # Expansion order for the non-flat monopole solution.

    # Making configuration files for this simulation.
    template_path = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_diffuse_sky.yaml')
    obspar_path = str(tmpdir.join('obsparam_diffuse_sky.yaml'))
    layout_path = str(tmpdir.join('threeant_layout.csv'))
    herauniform_path = str(tmpdir.join('hera_uniform.yaml'))

    teleconfig = {
        'beam_paths': {0: {"type": "uniform"}},
        'telescope_location': "(-30.72153, 21.42830, 1073.0)",
        'telescope_name': 'HERA'
    }
    if not use_w:
        antpos_enu = np.array([[0, 0, 0], [0, 3, 0], [5, 0, 0]], dtype=float)
    else:
        antpos_enu = np.array([[0, 0, 0], [0, 3, 0], [0, 3, 5]], dtype=float)

    pyuvsim.simsetup._write_layout_csv(
        layout_path, antpos_enu, np.arange(3).astype(str), np.arange(3)
    )
    with open(herauniform_path, 'w') as ofile:
        yaml.dump(teleconfig, ofile, default_flow_style=False)

    with open(template_path, 'r') as yfile:
        obspar = yaml.safe_load(yfile)
    obspar['telescope']['array_layout'] = layout_path
    obspar['telescope']['telescope_config_name'] = herauniform_path
    obspar['sources']['diffuse_model'] = modname
    obspar['sources'].update(params)
    obspar['filing']['outfile_name'] = 'diffuse_sim.uvh5'
    obspar['filing']['output_format'] = 'uvh5'
    obspar['filing']['outdir'] = str(tmpdir)

    with open(obspar_path, 'w') as ofile:
        yaml.dump(obspar, ofile, default_flow_style=False)

    uv_out = pyuvsim.run_uvsim(obspar_path, return_uv=True)
    # Convert from Jy to K sr
    dat = uv_out.data_array[:, 0, 0] * jy_to_ksr(uv_out.freq_array[0]).value
    # Evaluate the solution and compare to visibilities.
    soln = analytic_diffuse.get_solution(modname)
    uvw_lam = uv_out.uvw_array * uv_out.freq_array[0] / c_ms
    ana = soln(uvw_lam, **params)
    assert np.allclose(ana / 2, dat, atol=1e-2)


@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.filterwarnings("ignore:Fixing auto polarization power beams")
def test_powerbeam_sim(cst_beam):
    new_cst = copy.deepcopy(cst_beam)
    if hasattr(new_cst, "_freq_interp_kind"):
        new_cst.freq_interp_kind = 'nearest'  # otherwise we get an error about freq interpolation
    new_cst.efield_to_power()
    beams = BeamList([new_cst] * 4)
    cfg = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    input_uv = pyuvsim.simsetup.initialize_uvdata_from_params(cfg, return_beams=False)
    sky_model = pyuvsim.simsetup.initialize_catalog_from_params(
        cfg, return_catname=False
    )
    sky_model = pyuvsim.simsetup.SkyModelData(sky_model)

    with pytest.raises(ValueError, match="Beam type must be efield!"):
        pyuvsim.run_uvdata_uvsim(input_uv, beams, catalog=sky_model)


@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.parametrize("rename_beamfits", [True, False])
def test_run_paramdict_uvsim(rename_beamfits, tmp_path):
    # Running a simulation from parameter dictionary.
    param_file = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')

    msg = ["Cannot check consistency of a string-mode BeamList"]
    warn_type = [UserWarning]
    if rename_beamfits:
        os.makedirs(os.path.join(tmp_path, 'test_config'))
        new_param_file = os.path.join(tmp_path, 'test_config', 'param_1time_1src_testcat.yaml')
        shutil.copyfile(param_file, new_param_file)

        telescope_param_file = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan.yaml'
        )
        new_telescope_param_file = os.path.join(
            tmp_path, 'test_config', '28m_triangle_10time_10chan.yaml'
        )
        shutil.copyfile(telescope_param_file, new_telescope_param_file)

        telescope_layout_file = os.path.join(
            SIM_DATA_PATH, 'test_config', 'triangle_bl_layout.csv'
        )
        new_telescope_layout_file = os.path.join(tmp_path, 'test_config', 'triangle_bl_layout.csv')
        shutil.copyfile(telescope_layout_file, new_telescope_layout_file)

        source_file = os.path.join(SIM_DATA_PATH, 'single_source.txt')
        new_source_file = os.path.join(tmp_path, 'single_source.txt')
        shutil.copyfile(source_file, new_source_file)

        beamfits_file = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.beamfits')
        new_beam_file = os.path.join(tmp_path, 'test_config', 'HERA_NicCST.uvbeam')
        shutil.copyfile(beamfits_file, new_beam_file)

        # change the beam file name to .uvbeam
        with open(new_telescope_param_file, 'r') as pfile:
            tele_param_dict = yaml.safe_load(pfile)
            tele_param_dict["beam_paths"][0] = {"filename": new_beam_file}

        with open(new_telescope_param_file, 'w') as yfile:
            yaml.dump(tele_param_dict, yfile, default_flow_style=False)

        n_beam_warnings = 3
        warn_type += [DeprecationWarning] * n_beam_warnings
        msg += [pyuvsim.telescope.weird_beamfits_extension_warning] * n_beam_warnings

        params = pyuvsim.simsetup._config_str_to_dict(new_param_file)
    else:
        params = pyuvsim.simsetup._config_str_to_dict(param_file)

    with uvtest.check_warnings(warn_type, match=msg):
        pyuvsim.run_uvsim(params, return_uv=True)


@pytest.mark.filterwarnings("ignore:No julian date given for mock catalog")
@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
def test_run_nsky_parts(capsys):
    # there parameters were hand picked and fine-tuned to create nsky_parts = 2
    #  this test feels very wonky just to ensure the nsky_parts is printed
    scale = 1.0
    if 'linux' in sys.platform:
        scale = 2**10
    # Running a simulation from parameter dictionary.
    os.environ["SLURM_MEM_PER_NODE"] = str(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 1e6 + 5
    )  # Only 5MB of memory more than used right now
    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    )
    # just piggy backing on the setup but we need more sources for this to work
    params['sources']["catalog"] = "mock"
    params["sources"]["mock_arrangement"] = "random"
    params["sources"]["Nsrcs"] = 5000

    pyuvsim.run_uvsim(params, return_uv=True)

    captured = capsys.readouterr()
    assert "The source list has been split into Nsky_parts" in captured.out
    del os.environ["SLURM_MEM_PER_NODE"]


@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.parametrize(
    "spectral_type",
    ["flat", "subband", "spectral_index"])
def test_run_gleam_uvsim(spectral_type):
    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testgleam.yaml')
    )
    params["sources"]["spectral_type"] = spectral_type
    params["sources"].pop("min_flux")
    params["sources"].pop("max_flux")

    pyuvsim.run_uvsim(params, return_uv=True)


@pytest.mark.filterwarnings("ignore:The reference_frequency is aliased as `frequency`")
@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.parametrize(
    "spectral_type",
    ["subband", "spectral_index"])
def test_zenith_spectral_sim(spectral_type, tmpdir):
    # Make a power law source at zenith in three ways.
    # Confirm that simulated visibilities match expectation.

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    )

    alpha = -0.5
    ref_freq = 111e6
    Nfreqs = 20
    freqs = np.linspace(110e6, 115e6, Nfreqs)
    freq_params = pyuvsim.simsetup.freq_array_to_params(freqs)
    freqs = pyuvsim.simsetup.parse_frequency_params(freq_params)['freq_array']
    freqs *= units.Hz
    spectrum = (freqs.value / ref_freq)**alpha

    source, kwds = pyuvsim.create_mock_catalog(Time.now(), arrangement='zenith', Nsrcs=1)
    source.spectral_type = spectral_type
    if spectral_type == 'spectral_index':
        source.reference_frequency = np.array([ref_freq]) * units.Hz
        source.spectral_index = np.array([alpha])
    else:
        freq_lower = freqs - freq_params["channel_width"] * units.Hz / 2.
        freq_upper = freqs + freq_params["channel_width"] * units.Hz / 2.
        source.freq_edge_array = np.concatenate(
            (freq_lower[np.newaxis, :], freq_upper[np.newaxis, :]), axis=0
        )
        source.Nfreqs = Nfreqs
        source.freq_array = freqs
        source.stokes = np.repeat(source.stokes, Nfreqs, axis=1)
        source.stokes[0, :, 0] *= spectrum
        source.coherency_radec = stokes_to_coherency(source.stokes)

    catpath = str(tmpdir.join('spectral_test_catalog.skyh5'))
    source.write_skyh5(catpath)
    params['sources'] = {"catalog" : catpath}
    params['filing']['outdir'] = str(tmpdir)
    params['freq'] = freq_params
    params['time']['start_time'] = kwds['time']
    params['select'] = {'antenna_nums' : [1, 2]}

    uv_out = pyuvsim.run_uvsim(params, return_uv=True)

    for ii in range(uv_out.Nbls):
        assert np.allclose(uv_out.data_array[ii, :, 0], spectrum / 2)


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    hera_uv = UVData()

    hera_uv.polarizations = ['xx']

    with pytest.raises(ValueError, match='input_uv must have XX,YY,XY,YX polarization'):
        pyuvsim.run_uvdata_uvsim(hera_uv, ['beamlist'])


def test_input_uv_error():
    with pytest.raises(TypeError, match="input_uv must be UVData object"):
        pyuvsim.run_uvdata_uvsim(None, None)


@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
@pytest.mark.filterwarnings("ignore:This method will be removed in version 3.0")
@pytest.mark.filterwarnings("ignore:The lst_array is not self-consistent")
@pytest.mark.parametrize("future_shapes", [True, False])
def test_sim_on_moon(future_shapes):
    pytest.importorskip("lunarsky")
    from lunarsky import MoonLocation
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_tranquility_hex.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    param_dict['select'] = {'redundant_threshold': 0.1}
    uv_obj, beam_list, beam_dict = pyuvsim.initialize_uvdata_from_params(
        param_dict, return_beams=True
    )

    # set the filename to make sure it ends up in the history,
    # remove the parameter file info from extra_keywords
    uv_obj.filename = ["moon_sim"]
    uv_obj._filename.form = (1,)
    uv_obj.extra_keywords.pop('obsparam')
    uv_obj.extra_keywords.pop('telecfg')
    uv_obj.extra_keywords.pop('layout')
    uv_obj.check()

    uv_obj.select(times=uv_obj.time_array[0])
    tranquility_base = MoonLocation.from_selenocentric(*uv_obj.telescope_location, 'meter')

    time = Time(uv_obj.time_array[0], format='jd', scale='utc')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, array_location=tranquility_base, arrangement='zenith', Nsrcs=30, return_data=True
    )
    # Run simulation.
    if not future_shapes:
        uv_obj.use_current_array_shapes()
    uv_out = pyuvsim.uvsim.run_uvdata_uvsim(
        uv_obj, beam_list, beam_dict, catalog=sources, quiet=True
    )
    assert uvutils._check_history_version(uv_out.history, pyradiosky.__version__)
    assert uvutils._check_history_version(uv_out.history, pyuvdata.__version__)
    assert uvutils._check_history_version(uv_out.history, pyuvsim.__version__)
    assert uvutils._check_history_version(uv_out.history, uv_obj.filename[0])
    assert uvutils._check_history_version(uv_out.history, "Npus =")

    assert np.allclose(uv_out.data_array[:, :, 0], 0.5)
    assert uv_out.extra_keywords['world'] == 'moon'
