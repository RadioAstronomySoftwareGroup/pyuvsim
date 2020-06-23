# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
import os
import pytest
import yaml

from pyuvdata import UVData

import pyuvsim
from pyuvsim.astropy_interface import Time
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


@pytest.mark.filterwarnings("ignore:The frequency field is included in the recarray")
def test_run_paramfile_uvsim():
    # Test vot and txt catalogs for parameter simulation

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


def test_pol_error():
    # Check that running with a uvdata object without the proper polarizations will fail.
    pytest.importorskip('mpi4py')

    hera_uv = UVData()

    hera_uv.polarizations = ['xx']

    with pytest.raises(ValueError, match='input_uv must have XX,YY,XY,YX polarization'):
        pyuvsim.run_uvdata_uvsim(hera_uv, ['beamlist'])


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
