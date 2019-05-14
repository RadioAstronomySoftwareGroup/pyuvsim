# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import numpy as np
import yaml
import nose.tools as nt
from mpi4py import MPI
import astropy

from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
import pyuvdata.tests as uvtest

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
from pyuvsim import mpi
import pyuvsim.tests as simtest


cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')
param_filenames = [os.path.join(SIM_DATA_PATH, 'test_config', 'param_10time_10chan_{}.yaml'.format(x)) for x in range(4)]   # Five different test configs
singlesource_vot = os.path.join(SIM_DATA_PATH, 'single_source.vot')
singlesource_txt = os.path.join(SIM_DATA_PATH, 'single_source.txt')


def test_mpi_version():
    nt.assert_true(MPI.VERSION == 3)


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beamfile = os.path.join(simtest.TESTDATA_PATH, "temp.uvbeam")
    beam.write_beamfits(beamfile)

    beam_list = [beamfile]
    mock_keywords = {"Nsrcs": 3}
    nt.assert_raises(TypeError, pyuvsim.run_uvdata_uvsim, 'not_uvdata', beam_list)
    mpi.start_mpi()
    catalog, mock_kwds = pyuvsim.simsetup.create_mock_catalog(hera_uv.time_array[0], **mock_keywords)
    uv_out = pyuvsim.run_uvdata_uvsim(hera_uv, beam_list, catalog=catalog, source_list_name='mock', obs_param_file='', telescope_config_file='', antenna_location_file='')
    rank = mpi.get_rank()
    if rank == 0:
        nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3))
    os.remove(beamfile)


def test_run_param_uvsim():
    # Test vot and txt catalogs for parameter simulation

    uv_ref = UVData()
    uvtest.checkWarnings(uv_ref.read_uvfits, [os.path.join(SIM_DATA_PATH, 'testfile_singlesource.uvfits')],
                         nwarnings=1, message='antenna_diameters is not set')
    uv_ref.unphase_to_drift(use_ant_pos=True)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')

    with open(param_filename) as pfile:
        params_dict = yaml.safe_load(pfile)

    tempfilename = params_dict['filing']['outfile_name']
    # This test obsparam file has "single_source.txt" as its catalog.

    uvtest.checkWarnings(pyuvsim.uvsim.run_uvsim, [param_filename], nwarnings=1,
                         message=['The default for the `center` keyword has changed'],
                         category=DeprecationWarning)
    uv_new_txt = UVData()
    uvtest.checkWarnings(uv_new_txt.read_uvfits, [tempfilename], message='antenna_diameters is not set')
    uv_new_txt.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')

    uvtest.checkWarnings(pyuvsim.uvsim.run_uvsim, [param_filename],
                         nwarnings=11,
                         message=([SIM_DATA_PATH] * 10
                                  + ['The default for the `center` keyword has changed']),
                         category=([astropy.io.votable.exceptions.W50] * 10
                                   + [DeprecationWarning]))

    uv_new_vot = UVData()
    uvtest.checkWarnings(uv_new_vot.read_uvfits, [tempfilename], message='antenna_diameters is not set')
    uv_new_vot.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)
    uv_new_txt.history = uv_ref.history  # History includes irrelevant info for comparison
    uv_new_vot.history = uv_ref.history
    uv_new_txt.object_name = uv_ref.object_name
    nt.assert_equal(uv_new_txt, uv_ref)
    nt.assert_equal(uv_new_vot, uv_ref)


def test_mpi_funcs():
    mpi.start_mpi()
    nt.assert_true(mpi.get_rank() == 0)
    nt.assert_true(mpi.get_Npus() == 1)
    nt.assert_true(isinstance(mpi.get_comm(), MPI.Intracomm))
