# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

from __future__ import absolute_import, division, print_function

from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD
Npus = comm.Get_size()
rank = comm.Get_rank()

def get_mpi():
    """
    Initialize MPI, get the communicator, number of Processing Units (PUs)
    and the rank of this PU
    """

    return comm, rank, Npus


def set_mpi_excepthook(mpi_comm):
    """Kill the whole job on an uncaught python exception"""

    def mpi_excepthook(exctype, value, traceback):
        sys.__excepthook__(exctype, value, traceback)
        mpi_comm.Abort(1)

    sys.excepthook = mpi_excepthook
