# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
import os
import pytest
import yaml
from astropy import units

from pyuvdata import UVData
from pyradiosky.utils import jy_to_ksr, stokes_to_coherency

import pyuvsim
from pyuvsim.astropy_interface import Time
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.analyticbeam import c_ms


pytest.importorskip('mpi4py')  # noqa


@pytest.mark.filterwarnings("ignore:The frequency field is included in the recarray")
def test_run_paramfile_uvsim():
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.

    uv_ref = UVData()
    uv_ref.read_uvfits(os.path.join(SIM_DATA_PATH, 'testfile_singlesource.uvfits'))
    uv_ref.unphase_to_drift(use_ant_pos=True)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    with open(param_filename) as pfile:
        params_dict = yaml.safe_load(pfile)
    tempfilename = params_dict['filing']['outfile_name']

    # This test obsparam file has "single_source.txt" as its catalog.
    pyuvsim.uvsim.run_uvsim(param_filename)

    uv_new_txt = UVData()
    with pytest.warns(UserWarning, match='antenna_diameters is not set'):
        uv_new_txt.read_uvfits(tempfilename)

    uv_new_txt.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')
    pyuvsim.uvsim.run_uvsim(param_filename)

    uv_new_vot = UVData()
    with pytest.warns(UserWarning, match='antenna_diameters is not set'):
        uv_new_vot.read_uvfits(tempfilename)

    uv_new_vot.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)
    uv_new_txt.history = uv_ref.history  # History includes irrelevant info for comparison
    uv_new_vot.history = uv_ref.history
    uv_new_txt.object_name = uv_ref.object_name
    uv_new_vot.object_name = uv_ref.object_name
    assert uv_new_txt == uv_ref
    assert uv_new_vot == uv_ref


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
        'beam_paths': {0: 'uniform'},
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
    dat = uv_out.data_array[:, 0, 0, 0] * jy_to_ksr(uv_out.freq_array[0, 0]).value
    # Evaluate the solution and compare to visibilities.
    soln = analytic_diffuse.get_solution(modname)
    uvw_lam = uv_out.uvw_array * uv_out.freq_array[0, 0] / c_ms
    ana = soln(uvw_lam, **params)
    assert np.allclose(ana / 2, dat, atol=1e-2)


@pytest.mark.filterwarnings("ignore:The frequency field is included in the recarray")
def test_run_paramdict_uvsim():
    # Running a simulation from parameter dictionary.

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    )

    pyuvsim.run_uvsim(params, return_uv=True)


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
    freqs = pyuvsim.simsetup.parse_frequency_params(freq_params)['freq_array'][0, :]
    freqs *= units.Hz
    spectrum = (freqs.value / ref_freq)**alpha

    source, kwds = pyuvsim.create_mock_catalog(Time.now(), arrangement='zenith', Nsrcs=1)
    source.spectral_type = spectral_type
    if spectral_type == 'spectral_index':
        source.reference_frequency = np.array([ref_freq]) * units.Hz
        source.spectral_index = np.array([alpha])
    else:
        source.Nfreqs = Nfreqs
        source.freq_array = freqs
        source.stokes = np.repeat(source.stokes, Nfreqs, axis=1)
        source.stokes[0, :, 0] *= spectrum
        source.coherency_radec = stokes_to_coherency(source.stokes)

    catpath = str(tmpdir.join('spectral_test_catalog.txt'))
    source.write_text_catalog(catpath)
    params['sources'] = {"catalog" : catpath}
    params['filing']['outdir'] = str(tmpdir)
    params['freq'] = freq_params
    params['time']['start_time'] = kwds['time']
    params['select'] = {'antenna_nums' : [1, 2]}

    uv_out = pyuvsim.run_uvsim(params, return_uv=True)

    for ii in range(uv_out.Nbls):
        assert np.allclose(uv_out.data_array[ii, 0, :, 0], spectrum / 2)


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    hera_uv = UVData()

    hera_uv.polarizations = ['xx']

    with pytest.raises(ValueError, match='input_uv must have XX,YY,XY,YX polarization'):
        pyuvsim.run_uvdata_uvsim(hera_uv, ['beamlist'])


def test_input_uv_error():
    with pytest.raises(TypeError, match="input_uv must be a UVData object"):
        pyuvsim.run_uvdata_uvsim(None, None)


@pytest.mark.skipif('not pyuvsim.astropy_interface.hasmoon')
def test_sim_on_moon():
    from pyuvsim.astropy_interface import MoonLocation
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_tranquility_hex.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    param_dict['select'] = {'redundant_threshold': 0.1}
    uv_obj, beam_list, beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    uv_obj.select(times=uv_obj.time_array[0])
    tranquility_base = MoonLocation.from_selenocentric(*uv_obj.telescope_location, 'meter')

    time = Time(uv_obj.time_array[0], format='jd', scale='utc')
    sources, kwds = pyuvsim.create_mock_catalog(
        time, array_location=tranquility_base, arrangement='zenith', Nsrcs=30, return_data=True
    )
    # Run simulation.
    uv_out = pyuvsim.uvsim.run_uvdata_uvsim(
        uv_obj, beam_list, beam_dict, catalog=sources, quiet=True
    )
    assert np.allclose(uv_out.data_array[:, 0, :, 0], 0.5)
    assert uv_out.extra_keywords['world'] == 'moon'
