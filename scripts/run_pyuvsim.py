#!/usr/bin/env python

# from pyuvsim import uvsim
import pyuvsim
import argparse
import os

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation."))

parser.add_argument('file_in', metavar='<FILE>', type=str)
parser.add_argument('--outdir', type=str, nargs=1, default='./')
parser.add_argument('--Nsrcs', type=int, default=3)
parser.add_argument('--overwrite', action='store_true')


args = parser.parse_args()

for filename in args.file_in:
    uvdata_out = pyuvsim.uvsim.run_uvsim(filename, Nsrcs=args.Nsrcs)
    outfile = os.path.join(args.outdir, 'sim' + os.path.basename(filename))
    uvdata_out.write_uvfits(outfile, clobber=args.overwrite)
