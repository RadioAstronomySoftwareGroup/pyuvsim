# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os

import numpy as np
import yaml
import resource
import time
import pytest
import sys

pytest.importorskip('mpi4py')  # noqa
import mpi4py
mpi4py.rc.initialize = False  # noqa
from mpi4py import MPI
from astropy import units
from astropy.coordinates import Latitude, Longitude

from pyuvdata import UVData
from pyuvdata.data import DATA_PATH

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
from pyuvsim import mpi
from pyuvsim.astropy_interface import Time
import pyuvsim.tests as simtest
import pyradiosky


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


@pytest.fixture
def fake_tasks():
    sky = pyradiosky.SkyModel(
        'src', Longitude('10d'), Latitude('5d'), [1, 0, 0, 0], 'spectral_index',
        reference_frequency=np.array([100e6]) * units.Hz, spectral_index=np.array([-0.74])
    )
    n_tasks = 30
    t0 = Time.now()
    freq = 50 * units.Hz
    objs = [pyuvsim.UVTask(sky, t0, freq, None, None, freq_i=ii) for ii in range(n_tasks)]
    for ti, task in enumerate(objs):
        task.visibility_vector = np.random.uniform(0, 3) + np.random.uniform(0, 3) * 1j
        task.uvdata_index = (ti, 0, 0)
        task.sources = 1     # Replace with something easier to compare later.
        objs[ti] = task

    return objs


def test_mpi_version():
    assert MPI.VERSION == 3


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    beam = simtest.make_cst_beams(freqs=[100e6, 123e6])
    beamfile = os.path.join(simtest.TESTDATA_PATH, "temp.uvbeam")
    beam.write_beamfits(beamfile)
    beam_list = pyuvsim.BeamList([beamfile])
    mock_keywords = {"Nsrcs": 3}

    with pytest.raises(TypeError, match='input_uv must be UVData object'):
        pyuvsim.run_uvdata_uvsim('not_uvdata', beam_list, quiet=True)

    mpi.start_mpi()
    catalog, mock_kwds = pyuvsim.simsetup.create_mock_catalog(
        hera_uv.time_array[0], return_data=True, **mock_keywords
    )
    uv_out = pyuvsim.run_uvdata_uvsim(hera_uv, beam_list, catalog=catalog)
    rank = mpi.get_rank()
    if rank == 0:
        assert np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3)
    os.remove(beamfile)


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
    with pytest.warns(UserWarning) as antdiam:
        uv_new_txt.read_uvfits(tempfilename)
    assert str(antdiam.pop().message).startswith('antenna_diameters is not set')
    uv_new_txt.unphase_to_drift(use_ant_pos=True)
    os.remove(tempfilename)

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')
    pyuvsim.uvsim.run_uvsim(param_filename)

    uv_new_vot = UVData()
    with pytest.warns(UserWarning) as antdiam:
        uv_new_vot.read_uvfits(tempfilename)
    assert str(antdiam.pop().message).startswith('antenna_diameters is not set')
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

    scale = 1.0
    if 'linux' in sys.platform:
        scale = 2**10

    memory_usage_GiB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**30
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


@pytest.mark.parametrize('MAX_BYTES', [mpi.INT_MAX, 100])
def test_big_gather(MAX_BYTES, fake_tasks):
    mpi.start_mpi()

    objs = fake_tasks
    n_tasks = len(objs)
    result, split_info = mpi.big_gather(
        mpi.world_comm, objs, root=0, return_split_info=True, MAX_BYTES=MAX_BYTES
    )
    assert all(objs[ii].freq_i == ii for ii in range(n_tasks))
    assert all(result[0][ii].freq == objs[ii].freq for ii in range(n_tasks))
    assert all(result[0][ii].time == objs[ii].time for ii in range(n_tasks))

    # Compare with normal gather:
    result2 = mpi.world_comm.gather(objs, root=0)

    assert result2 == result

    assert split_info['MAX_BYTES'] == MAX_BYTES
    if MAX_BYTES < 200:
        assert len(split_info['ranges']) > 1


@pytest.mark.parametrize('MAX_BYTES', [mpi.INT_MAX, 100])
def test_big_bcast(MAX_BYTES, fake_tasks):
    mpi.start_mpi()

    objs = fake_tasks
    n_tasks = len(objs)

    result, split_info = mpi.big_bcast(
        mpi.world_comm, objs, root=0, return_split_info=True, MAX_BYTES=MAX_BYTES
    )
    assert all(result[ii].freq_i == ii for ii in range(n_tasks))
    assert all(result[ii].freq == objs[ii].freq for ii in range(n_tasks))
    assert all(result[ii].time == objs[ii].time for ii in range(n_tasks))

    # Compare with normal gather:
    result2 = mpi.world_comm.bcast(objs, root=0)

    assert result2 == result

    assert split_info['MAX_BYTES'] == MAX_BYTES
    if MAX_BYTES < 200:
        assert len(split_info['ranges']) > 1


def test_big_bcast_gather_loop(fake_tasks):
    mpi.start_mpi()

    objs = fake_tasks

    broadcast = mpi.big_bcast(mpi.world_comm, objs, root=0, MAX_BYTES=35)
    gathered = mpi.big_gather(mpi.world_comm, broadcast, root=0, MAX_BYTES=27)

    assert broadcast == gathered[0]


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
    print(uv_obj.extra_keywords['world'])
    # Run simulation.
    uv_out = pyuvsim.uvsim.run_uvdata_uvsim(
        uv_obj, beam_list, beam_dict, catalog=sources, quiet=True
    )
    assert np.allclose(uv_out.data_array[:, 0, :, 0], 0.5)
    assert uv_out.extra_keywords['world'] == 'moon'


def test_sharedmem_bcast_with_quantities():
    # Use mpi.quantity_shared_bcast and check returned objects.

    mpi.start_mpi()
    lats = Latitude(np.linspace(-np.pi / 2, np.pi / 2, 10), 'rad')
    freqs = np.linspace(100e6, 130e6, 15) * units.Hz
    lat_return = mpi.quantity_shared_bcast(lats)
    freq_return = mpi.quantity_shared_bcast(freqs)

    assert np.all(lat_return == lats)
    assert np.all(freq_return == freqs)

    assert np.all(freq_return.to("MHz") == freqs.to("MHz"))
