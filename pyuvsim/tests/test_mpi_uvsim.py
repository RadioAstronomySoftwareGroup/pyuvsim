# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os

import mpi4py
import numpy as np
import yaml
import resource
import time
import six

mpi4py.rc.initialize = False  # noqa
import pytest
from mpi4py import MPI

from pyuvdata import UVData
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
param_filenames = [
    os.path.join(SIM_DATA_PATH, 'test_config', 'param_10time_10chan_{}.yaml'.format(x))
    for x in range(4)
]  # Five different test configs
singlesource_vot = os.path.join(SIM_DATA_PATH, 'single_source.vot')
singlesource_txt = os.path.join(SIM_DATA_PATH, 'single_source.txt')


def test_mpi_version():
    assert MPI.VERSION == 3


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beamfile = os.path.join(simtest.TESTDATA_PATH, "temp.uvbeam")
    beam.write_beamfits(beamfile)
    beam_list = [beamfile]
    mock_keywords = {"Nsrcs": 3}
    simtest.assert_raises_message(
        TypeError, 'input_uv must be UVData object',
        pyuvsim.run_uvdata_uvsim, 'not_uvdata', beam_list
    )
    mpi.start_mpi()
    catalog, mock_kwds = pyuvsim.simsetup.create_mock_catalog(
        hera_uv.time_array[0], return_table=True, **mock_keywords
    )
    uv_out = pyuvsim.run_uvdata_uvsim(hera_uv, beam_list, catalog=catalog)
    rank = mpi.get_rank()
    if rank == 0:
        assert np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3)
    os.remove(beamfile)


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
    uvtest.checkWarnings(pyuvsim.uvsim.run_uvsim, [param_filename], nwarnings=1,
                         message=['The default for the `center` keyword'],
                         category=[DeprecationWarning])
    uv_new_txt = UVData()
    uvtest.checkWarnings(uv_new_txt.read_uvfits, [tempfilename],
                         message='antenna_diameters is not set')
    uv_new_txt.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')
    if six.PY2:
        pyuvsim.uvsim.run_uvsim(param_filename)
    else:
        uvtest.checkWarnings(
            pyuvsim.uvsim.run_uvsim, [param_filename],
            nwarnings=1,
            message=['The default for the `center` keyword has changed'],
            category=[DeprecationWarning]
        )

    uv_new_vot = UVData()
    uvtest.checkWarnings(uv_new_vot.read_uvfits, [tempfilename], nwarnings=1,
                         category=[UserWarning], message=['antenna_diameters is not set'])
    uv_new_vot.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)
    uv_new_txt.history = uv_ref.history  # History includes irrelevant info for comparison
    uv_new_vot.history = uv_ref.history
    uv_new_txt.object_name = uv_ref.object_name
    uv_new_vot.object_name = uv_ref.object_name
    assert uv_new_txt == uv_ref
    assert uv_new_vot == uv_ref


def test_run_paramdict_uvsim():
    # Running a simulation from parameter dictionary.

    params = pyuvsim.simsetup._config_str_to_dict(
        os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    )

    pyuvsim.run_uvsim(params, return_uv=True)


def test_mpi_funcs():
    mpi.start_mpi()
    assert mpi.get_rank() == 0
    assert mpi.get_Npus() == 1
    assert isinstance(mpi.get_comm(), MPI.Intracomm)
    assert isinstance(mpi.get_comm(), MPI.Intracomm)
    assert isinstance(mpi.get_node_comm(), MPI.Intracomm)


def test_shared_mem():
    mpi.start_mpi()
    N = 200
    shape = (20, 10)
    A = np.arange(N, dtype=float).reshape(shape)

    sA = mpi.shared_mem_bcast(A)

    # Equivalent to original
    assert np.all(sA == A)

    # Not the original object:
    assert hex(id(sA)) != hex(id(A))

    # Shared array should be read-only
    pytest.raises(ValueError, sA.itemset, 0, 3.0)


def test_mem_usage():
    # Check that the mpi-enabled memory check is consistent
    # with a local memory check.

    # Also check that making a variable of a given size
    # increases memory usage by the expected amount.

    mpi.start_mpi()

    memory_usage_GiB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 2**10 / 2**30
    assert memory_usage_GiB == mpi.get_max_node_rss()
    incsize = 50 * 2**20    # 50 MiB
    arr = bytearray(incsize)
    time.sleep(1)
    change = mpi.get_max_node_rss() - memory_usage_GiB
    assert np.isclose(change, incsize / 2**30, atol=5e-2)
    del arr


@pytest.mark.skip
def test_mpi_counter():
    # This test should be run in parallel to check likely bug sources.
    mpi.start_mpi()

    count = mpi.Counter()
    N = 20
    for i in range(N):
        count.next()
    assert count.current_value() == N
    count.free()
