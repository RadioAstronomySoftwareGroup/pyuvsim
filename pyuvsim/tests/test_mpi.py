# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import numpy as np
import resource
import time
import sys

import pytest
from astropy import units  # noqa
from astropy.coordinates import Latitude, Longitude
import pyradiosky

pytest.importorskip('mpi4py')
import mpi4py   # noqa
mpi4py.rc.initialize = False  # noqa
from mpi4py import MPI  # noqa

import pyuvsim  # noqa
from pyuvsim import mpi  # noqa
from pyuvsim.astropy_interface import Time  # noqa


@pytest.fixture(scope='module', autouse=True)
def start_mpi():
    mpi.start_mpi(False)


@pytest.fixture(scope='module')
def fake_tasks():
    sky = pyradiosky.SkyModel(
        'src', Longitude('10d'), Latitude('5d'), np.array([1, 0, 0, 0]) * units.Jy,
        'spectral_index', reference_frequency=np.array([100e6]) * units.Hz,
        spectral_index=np.array([-0.74])
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
    assert MPI.VERSION >= 3


@pytest.mark.parallel(2)
def test_mpi_funcs():
    assert mpi.get_rank() == MPI.COMM_WORLD.rank
    assert mpi.get_Npus() == MPI.COMM_WORLD.size
    assert isinstance(mpi.get_comm(), MPI.Intracomm)
    assert isinstance(mpi.get_comm(), MPI.Intracomm)
    assert isinstance(mpi.get_node_comm(), MPI.Intracomm)


def test_shared_mem():
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

    scale = 1.0
    if 'linux' in sys.platform:
        scale = 2**10

    memory_usage_GiB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**30
    assert np.isclose(memory_usage_GiB, mpi.get_max_node_rss())
    incsize = 50 * 2**20    # 50 MiB
    arr = bytearray(incsize)
    time.sleep(1)
    change = mpi.get_max_node_rss() - memory_usage_GiB
    del arr
    assert np.isclose(change, incsize / 2**30, atol=5e-2)


@pytest.mark.parametrize("count_rank", [0, 2])
@pytest.mark.parallel(4)
def test_mpi_counter(count_rank):
    # Warning -- This test has been flaky in the past.
    mpi.start_mpi()
    count = mpi.Counter(count_rank=count_rank)
    N = 20
    for i in range(N):
        count.next()
    mpi.world_comm.Barrier()
    if mpi.world_comm.rank == count_rank:
        assert count.current_value() == N * mpi.world_comm.size
    count.free()


@pytest.mark.parametrize('max_bytes', [mpi.INT_MAX, 100])
def test_big_gather(max_bytes, fake_tasks):

    objs = fake_tasks
    n_tasks = len(objs)
    result, split_info = mpi.big_gather(
        mpi.world_comm, objs, root=0, return_split_info=True, MAX_BYTES=max_bytes
    )
    if mpi.rank == 0:
        assert all(objs[ii].freq_i == ii for ii in range(n_tasks))
        assert all(result[0][ii].freq == objs[ii].freq for ii in range(n_tasks))
        assert all(result[0][ii].time == objs[ii].time for ii in range(n_tasks))

    # Compare with normal gather:
    result2 = mpi.world_comm.gather(objs, root=0)

    if mpi.rank == 0:
        assert result2 == result

        assert split_info['MAX_BYTES'] == max_bytes
        if max_bytes < 200:
            assert len(split_info['ranges']) > 1


@pytest.mark.parallel(3)
@pytest.mark.parametrize('max_bytes', [mpi.INT_MAX, 100])
def test_big_bcast(max_bytes, fake_tasks):

    objs = fake_tasks
    n_tasks = len(objs)

    result, split_info = mpi.big_bcast(
        mpi.world_comm, objs, root=0, return_split_info=True, MAX_BYTES=max_bytes
    )

    if mpi.rank == 0:
        assert all(result[ii].freq_i == ii for ii in range(n_tasks))
        assert all(result[ii].freq == objs[ii].freq for ii in range(n_tasks))
        assert all(result[ii].time == objs[ii].time for ii in range(n_tasks))

    # Compare with normal gather:
    result2 = mpi.world_comm.bcast(objs, root=0)

    if mpi.rank == 0:
        assert result2 == result
        assert split_info['MAX_BYTES'] == max_bytes
        if max_bytes < 200:
            assert len(split_info['ranges']) > 1


@pytest.mark.parallel(3)
def test_big_bcast_gather_loop(fake_tasks):

    objs = fake_tasks

    broadcast = mpi.big_bcast(mpi.world_comm, objs, root=0, MAX_BYTES=35)
    gathered = mpi.big_gather(mpi.world_comm, broadcast, root=0, MAX_BYTES=27)

    if mpi.rank == 0:
        assert broadcast == gathered[0]


@pytest.mark.parallel(3)
def test_sharedmem_bcast_with_quantities():
    # Use mpi.quantity_shared_bcast and check returned objects.

    lats = Latitude(np.linspace(-np.pi / 2, np.pi / 2, 10), 'rad')
    freqs = np.linspace(100e6, 130e6, 15) * units.Hz
    lat_return = mpi.quantity_shared_bcast(lats)
    freq_return = mpi.quantity_shared_bcast(freqs)

    if mpi.rank == 0:
        assert np.all(lat_return == lats)
        assert np.all(freq_return == freqs)
        assert np.all(freq_return.to("MHz") == freqs.to("MHz"))


@pytest.mark.parallel(3)
def test_skymodeldata_share():
    # Test the SkyModelData share method.
    sky = pyradiosky.SkyModel(
        'src', Longitude('10d'), Latitude('5d'), np.array([1, 0, 0, 0]) * units.Jy,
        'spectral_index', reference_frequency=np.array([100e6]) * units.Hz,
        spectral_index=np.array([-0.74])
    )

    smd = pyuvsim.simsetup.SkyModelData()
    if mpi.rank == 0:
        smd = pyuvsim.simsetup.SkyModelData(sky)

    smd.share()  # Shares among processes.

    sky2 = smd.get_skymodel()

    assert sky2 == sky
