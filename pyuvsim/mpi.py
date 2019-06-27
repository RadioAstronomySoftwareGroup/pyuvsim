# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import sys
from threading import Thread
from array import array
import struct as _struct
import mpi4py
mpi4py.rc.initialize = False    # noqa
from mpi4py import MPI

rank = 0    # COMM_WORLD rank
Npus = 1
Npus_node = 1
world_comm = None
node_comm = None
rank_comm = None


def set_mpi_excepthook(mpi_comm):
    """Kill the whole job on an uncaught python exception"""

    def mpi_excepthook(exctype, value, traceback):      # pragma: no cover
        sys.__excepthook__(exctype, value, traceback)
        mpi_comm.Abort(1)

    sys.excepthook = mpi_excepthook


def start_mpi():
    """
    Check if MPI has already been initialized. If so, just set the communicators,
    Npus, and rank variables.
    """
    global world_comm, node_comm, rank_comm, rank, Npus
    if not MPI.Is_initialized():
        # Enable threading for the Counter class
        MPI.Init_thread(MPI.THREAD_MULTIPLE)
    world_comm = MPI.COMM_WORLD
    node_comm = world_comm.Split_type(MPI.COMM_TYPE_SHARED)
    rank_comm = world_comm.Split(color=node_comm.rank)

    Npus = world_comm.Get_size()
    Npus_node = node_comm.Get_size()
    rank = world_comm.Get_rank()
    set_mpi_excepthook(world_comm)

    world_comm.Barrier()

    if not rank == 0:       # pragma: no cover
        # For non-root ranks, do not print to stdout.
        # (Uncovered until we have multi-rank tests)
        global stdout
        sys.stdout = open('/dev/null', 'w')


def shared_mem_bcast(arr, root=0):
    """
    Allocate shared memory on each node and place contents of arr in it.

    Must be called from all PUs, but only the root process
    should pass in an array. Every other process should pass in None.
    """

    nbytes = 0
    itemsize = 0
    dtype = None
    Nitems = 0
    shape = tuple()

    if node_comm.rank == root:
        # Data cannot be shared across nodes.
        # Need to broadcast to the root PU on each node.
        arr = rank_comm.bcast(arr, root=root)
        Nitems = arr.size
        shape = arr.shape
        itemsize = sys.getsizeof(arr.flatten()[0])
        dtype = arr.dtype
        nbytes = itemsize * Nitems

    itemsize = node_comm.bcast(itemsize, root=root)
    dtype = node_comm.bcast(dtype, root=root)
    Nitems = node_comm.bcast(Nitems, root=root)
    shape = node_comm.bcast(shape, root=root)

    # Allocate a window if the node_comm rank is 0
    # Otherwise, make a handle to the window.
    # This will allocate nbytes on each node.

    win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=node_comm)
    buf, itemsize = win.Shared_query(0)
    sh_arr = np.ndarray(buffer=buf, dtype=dtype, shape=shape)

    if node_comm.rank == root:
        # Now fill the window on each node with the data.
        sh_arr[:] = arr

    sh_arr.flags['WRITEABLE'] = False   # Do not want ranks overwriting the data.

    world_comm.Barrier()
    return sh_arr


class Counter(object):
    """
    A basic parallelized counter class.

    Adapted from the mpi4py nxtval-threads.py demo.
    https://github.com/mpi4py/mpi4py/blob/master/demo/nxtval/nxtval-threads.py
    """

    def __init__(self, comm=None):
        """
        Create a new counter, and initialize to 0.
        """
        # duplicate communicator
        if comm is None:
            comm = world_comm
        assert not comm.Is_inter()
        self.comm = comm.Dup()
        # start counter thread
        self.thread = None
        rank = self.comm.Get_rank()
        if rank == 0:
            self.thread = Thread(target=self._counter_thread)
            self.thread.daemon = True
            self.thread.start()

    def _counter_thread(self):
        incr = array('i', [0])
        ival = array('i', [0])
        status = MPI.Status()
        while True:  # server loop
            self.comm.Recv([incr, MPI.INT],
                           MPI.ANY_SOURCE, MPI.ANY_TAG,
                           status)
            if status.Get_tag() == 1:
                return
            self.comm.Send([ival, MPI.INT],
                           status.Get_source(), 0)
            ival[0] += incr[0]

    def free(self):
        self.comm.Barrier()
        # stop counter thread
        rank = self.comm.Get_rank()
        if rank == 0:
            self.comm.Send([None, MPI.INT], 0, 1)
            self.thread.join()
        self.comm.Free()

    def next(self):
        """
        Increment counter.

        Returns:
            next value (integer)
        """
        incr = array('i', [1])
        ival = array('i', [0])
        self.comm.Send([incr, MPI.INT], 0, 0)
        self.comm.Recv([ival, MPI.INT], 0, 0)
        nxtval = ival[0]
        return nxtval

    def current_value(self):
        """
        Returns:
            current value of counter (integer)
        """
        incr = array('i', [0])
        ival = array('i', [0])
        self.comm.Send([incr, MPI.INT], 0, 0)
        self.comm.Recv([ival, MPI.INT], 0, 0)
        val = ival[0]
        return val


def get_rank():
    """
    Current rank on COMM_WORLD

    """
    return rank


def get_Npus():
    """
    Number of MPI processes.
    """
    return Npus


def get_comm():
    """
    world_comm, the communicator for all PUs
    """
    return world_comm


def get_node_comm():
    """
    node_comm : Communicator for all PUs on current node.
    """
    return node_comm
