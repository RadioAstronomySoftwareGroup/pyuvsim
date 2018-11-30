# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import mpi4py
import sys
mpi4py.rc.initialize = False    # noqa
from mpi4py import MPI

rank = 0
Npus = 1
comm = None


def set_mpi_excepthook(mpi_comm):
    """Kill the whole job on an uncaught python exception"""

    def mpi_excepthook(exctype, value, traceback):      # pragma: no cover
        sys.__excepthook__(exctype, value, traceback)
        mpi_comm.Abort(1)

    sys.excepthook = mpi_excepthook


def start_mpi():
    """
    Check if MPI has already been initialized. If so, just set the comm,
    Npus, and rank variables.
    """
    global comm, rank, Npus
    if not MPI.Is_initialized():
        # Avoid accidentally doing MPI_INIT twice
        MPI.Init()
    comm = MPI.COMM_WORLD
    Npus = comm.Get_size()
    rank = comm.Get_rank()
    set_mpi_excepthook(comm)

    if not rank == 0:       # pragma: no cover
        # For non-root ranks, do not print to stdout.
        # (Uncovered until we have multi-rank tests)
        global stdout
        sys.stdout = open('/dev/null', 'w')


def get_rank():
    return rank


def get_Npus():
    return Npus


def get_comm():
    return comm
