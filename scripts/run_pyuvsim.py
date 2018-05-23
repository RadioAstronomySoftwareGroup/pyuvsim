#!/usr/bin/env python

# from pyuvsim import uvsim
import pyuvsim
import argparse
import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('--outdir', type=str, default='./')
parser.add_argument('--Nsrcs', type=int, default=3)
# parser.add_argument('--overwrite', action='store_true')


args = parser.parse_args()

for filename in args.file_in:
    uvdata_out = pyuvsim.uvsim.run_uvsim(filename, Nsrcs=args.Nsrcs)
    if rank == 0:
        outfile = os.path.join(args.outdir, 'sim_' + os.path.basename(filename))
        uvdata_out.write_uvfits(outfile, force_phase=True, spoof_nonessential=True)
