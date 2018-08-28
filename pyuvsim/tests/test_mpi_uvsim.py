# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import numpy as np
import nose.tools as nt
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
import yaml
from mpi4py import MPI

from pyuvsim import mpi


cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
param_filenames = [os.path.join(SIM_DATA_PATH, 'test_config', 'param_10time_10chan_{}.yaml'.format(x)) for x in range(4)]   # Five different test configs
singlesource_vot = os.path.join(SIM_DATA_PATH, 'single_source.vot')
singlesource_txt = os.path.join(SIM_DATA_PATH, 'single_source.txt')


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]
    mock_keywords = {"Nsrcs": 3}
    uv_out = pyuvsim.run_uvsim(hera_uv, beam_list, catalog_file=None, mock_keywords=mock_keywords,
                               uvdata_file=EW_uvfits_file)
    rank = mpi.get_rank()
    if rank == 0:
        nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3))


def test_run_param_uvsim():
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    with open(param_filename, 'r') as pfile:
        params_dict = yaml.safe_load(pfile)
    uv_in, beam_list, beam_dict, beam_ids = pyuvsim.simsetup.initialize_uvdata_from_params(param_filename)
    beam_list[0] = pyuvsim.analyticbeam.AnalyticBeam('uniform')  # Replace the one that's a HERA beam
    # This test obsparam file has "single_source.txt" as its catalog.
    catalog = os.path.join(SIM_DATA_PATH, params_dict['sources']['catalog'])
    uv_out = pyuvsim.run_uvsim(uv_in, beam_list, catalog_file=catalog, beam_dict=beam_dict)
    pyuvsim.simsetup.write_uvfits(uv_out, params_dict)
    tempfilename = params_dict['filing']['outfile_name']

    uv_new = UVData()
    uv_new.read_uvfits(tempfilename)
    os.remove(tempfilename)
    uv_ref = UVData()
    uv_ref.read_uvfits(os.path.join(SIM_DATA_PATH, 'testfile_singlesource.uvfits'))
    uv_new.history = uv_ref.history  # History includes irrelevant info for comparison
    nt.assert_equal(uv_new, uv_ref)

    # Test votable catalog
    catalog = singlesource_vot
    del uv_out, uv_new
    uv_out = pyuvsim.run_uvsim(uv_in, beam_list, catalog_file=catalog, beam_dict=beam_dict)
    pyuvsim.simsetup.write_uvfits(uv_out, params_dict)
    uv_new = UVData()
    uv_new.read_uvfits(tempfilename)
    print(tempfilename)
    os.remove(tempfilename)
    uv_new.history = uv_ref.history  # History includes irrelevant info for comparison
    nt.assert_equal(uv_new, uv_ref)


def test_mpi_funcs():
    mpi.start_mpi()
    nt.assert_true(mpi.get_rank() == 0)
    nt.assert_true(mpi.get_Npus() == 1)
    nt.assert_true(isinstance(mpi.get_comm(), MPI.Intracomm))
