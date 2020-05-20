# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import sys
from array import array
from threading import Thread
import resource
import atexit
from pickle import loads, dumps

import mpi4py
import numpy as np

mpi4py.rc.initialize = False  # noqa
from mpi4py import MPI

rank = 0  # COMM_WORLD rank
Npus = 1
world_comm = None
node_comm = None
rank_comm = None

# Maximum size of objects that can be handled by a single MPI operation.
INT_MAX = 2**32 - 1


def set_mpi_excepthook(mpi_comm):
    """Kill the whole job on an uncaught python exception"""

    def mpi_excepthook(exctype, value, traceback):  # pragma: no cover
        sys.__excepthook__(exctype, value, traceback)
        sys.stderr.flush()
        mpi_comm.Abort(1)

    sys.excepthook = mpi_excepthook


def start_mpi(block_nonroot_stdout=True):
    """
    Check if MPI has already been initialized. If so, just set the communicators,
    Npus, and rank variables.

    Parameters
    ----------

    block_nonroot_stdout : bool (True)
        Redirect stdout on nonzero ranks to /dev/null, for cleaner output.

    """
    global world_comm, node_comm, rank_comm, rank, Npus
    if not MPI.Is_initialized():
        # Enable threading for the Counter class
        MPI.Init_thread(MPI.THREAD_MULTIPLE)
        atexit.register(MPI.Finalize)
    world_comm = MPI.COMM_WORLD
    node_comm = world_comm.Split_type(MPI.COMM_TYPE_SHARED)
    rank_comm = world_comm.Split(color=node_comm.rank)

    Npus = world_comm.Get_size()
    rank = world_comm.Get_rank()
    set_mpi_excepthook(world_comm)

    world_comm.Barrier()

    if (not rank == 0) and block_nonroot_stdout:  # pragma: no cover
        # For non-root ranks, do not print to stdout.
        # (Uncovered until we have multi-rank tests)
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
    shape = tuple()  # noqa

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
        sh_arr[:] = arr[()]

    sh_arr.flags['WRITEABLE'] = False  # Do not want ranks overwriting the data.

    world_comm.Barrier()
    return sh_arr


def big_gather(comm, objs, root=0, return_split_info=False):
    """
    Gather operation that can exceed the MPI limit of ~4 GiB.

    MPI stores the total size of data to be gathered in a 32 bit integer,
    so that Gather operations which will gather more than 2**32 bytes will fail.

    This function replicates the behavior of mpi4py's `gather` method, while
    avoiding this size limit. The lowercase `gather` function first pickles
    Python objects before using the vectorized Gatherv. In this function,
    the pickled data are gathered in stages such that each step is smaller
    than the hard limit.

    Parameters
    ----------
    comm: mpi4py.MPI.Intracomm
        MPI communicator to use.
    objs: objects
        Data to gather from all processes.
    root: int
        Rank of process to receive the data.
    return_split_info: bool
        On root process, also a return a dictionary describing
        how the data were split.

    Returns
    -------
    list of objects:
        Length Npus list, such that the n'th entry is the data gathered from
        the n'th process.
        This is only filled on the root process. Other processes get None.
    dict:
        If return_split_info, the root process also gets a dictionary containing:
        - ranges: A list of tuples, giving the start and end byte of each chunk.
        - MAX_BYTES: The size limit that was used.
    """

    # The limit is on the integer describing the number of bytes gathered.
    MAX_BYTES = INT_MAX
    sbuf = dumps(objs)
    bytesize = len(sbuf)

    # Sizes of send buffers to be sent from each rank.
    counts = np.array(comm.allgather(bytesize))
    totsize = sum(counts)

    rbuf = None
    displ = None
    if comm.rank == 0:
        rbuf = np.empty(sum(counts), dtype=bytes)
        displ = np.array([sum(counts[1:p]) for p in range(comm.size)])

    # Position in the output buffer for the current send buffer.
    start_loc = sum(counts[:comm.rank])

    # Ranges of output bytes for each chunk.
    start = 0
    end = 0
    ranges = []
    while end < totsize:
        end = min(start + MAX_BYTES, totsize)
        ranges.append((start, end))
        start += MAX_BYTES + 1

    for start, end in ranges:
        # start/end indices of the local data to send, for this chunk.
        start_ind = min(max((start - start_loc), 0), bytesize)
        end_ind = min(max((end - start_loc), 0), bytesize)
        cur_sbuf = sbuf[start_ind:end_ind + 1]
        cur_counts = np.array(comm.gather(len(cur_sbuf), root=0))
        if len(cur_sbuf) > 0:
            loc_disp = max(start, start_loc)
        else:
            loc_disp = 0
        cur_displ = comm.gather(loc_disp, root=0)
        comm.Gatherv(sendbuf=cur_sbuf, recvbuf=(rbuf, cur_counts, cur_displ, MPI.BYTE), root=0)

    per_proc = None
    if comm.rank == 0:
        per_proc = [loads(rbuf[displ[ii]:displ[ii] + counts[ii]]) for ii in range(comm.size)]
    
    split_info_dict = None
    if comm.rank == 0:
        split_info_dict = {'MAX_BYTES': MAX_BYTES, 'ranges': ranges}
    
    if return_split_info:
        return per_proc, split_info_dict

    return per_proc

import time
class Counter(object):
    """
    A basic parallelized counter class.

    Adapted from the mpi4py nxtval-threads.py demo.
    https://github.com/mpi4py/mpi4py/blob/master/demo/nxtval/nxtval-threads.py
    """

    def __init__(self, comm=None, count_rank=0):
        """
        Create a new counter, and initialize to 0.
        """
        # duplicate communicator
        self.count_rank = count_rank
        if comm is None:
            comm = world_comm
        assert not comm.Is_inter()
        self.comm = comm.Dup()
        # start counter thread
        self.thread = None
        rank = self.comm.Get_rank()
        if rank == self.count_rank:
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
                           status.Get_source(), self.count_rank)
            ival[0] += incr[0]

    def free(self):
        self.comm.Barrier()
        # stop counter thread
        rank = self.comm.Get_rank()
        if rank == self.count_rank:
            self.comm.Send([None, MPI.INT], self.count_rank, 1)
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
        self.comm.Send([incr, MPI.INT], self.count_rank, 0)
        self.comm.Recv([ival, MPI.INT], self.count_rank, 0)
        nxtval = ival[0]
        return nxtval

    def current_value(self):
        """
        Returns:
            current value of counter (integer)
        """
        incr = array('i', [0])
        ival = array('i', [0])
        self.comm.Send([incr, MPI.INT], self.count_rank, 0)
        self.comm.Recv([ival, MPI.INT], self.count_rank, 0)
        val = ival[0]
        return val


def get_max_node_rss(return_per_node=False):
    """
    Find the maximum memory usage on any node in the job in GiB.

    Parameters
    ----------

    return_per_node : bool (optional)
        Return the total memory on the node to each rank on
        that node. (Default is False)

    Returns
    -------

    max_mem : float
        Maximum memory usage in GiB across the job.
        Only returns to the zero-th rank on the world_comm.
    """

    # On linux, getrusage returns in kiB
    # On Mac systems, getrusage returns in B
    scale = 1.0
    if 'linux' in sys.platform:
        scale = 2**10

    memory_usage_GiB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**30
    node_mem_tot = node_comm.allreduce(memory_usage_GiB, op=MPI.SUM)
    if return_per_node:
        return node_mem_tot

    max_mem = world_comm.reduce(node_mem_tot, op=MPI.MAX, root=0)
    return max_mem


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
