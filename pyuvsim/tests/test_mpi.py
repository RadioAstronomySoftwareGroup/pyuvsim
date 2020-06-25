# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np
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

import pyuvsim
from pyuvsim import mpi
from pyuvsim.astropy_interface import Time
import pyradiosky


@pytest.fixture
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
    assert MPI.VERSION == 3


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


def test_mpi_counter(x):
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
