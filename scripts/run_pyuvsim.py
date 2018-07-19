#!/usr/bin/env python

# from pyuvsim import uvsim
import pyuvsim
import argparse
import os
import numpy as np
from mpi4py import MPI
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import warnings

beam_file = os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if not rank == 0:
    warnings.simplefilter('ignore')

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('--outdir', type=str, default='./')
parser.add_argument('--Nsrcs', type=int, default=None)
parser.add_argument('--mock_arrangement', type=str, default='zenith')
parser.add_argument('--max_za', type=float, default=-1.0, help='Maximum zenith angle for mock arrangements')
parser.add_argument('--save_catalog', action='store_true', default=False, help='Save catalog')
# parser.add_argument('--overwrite', action='store_true')


args = parser.parse_args()

for filename in args.file_in:
    if rank == 0:
        print("Reading:", os.path.basename(filename))
    input_uv = UVData()
    input_uv.read_uvfits(filename, read_data=False)
    time0 = input_uv.time_array[0]
    input_uv.read_uvfits(filename, freq_chans=0, times=time0)
#    beam = UVBeam()
#    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[150e6],
#                       telescope_name='HERA',
#                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
#                       model_name='E-field pattern - Rigging height 4.9m',
#                       model_version='1.0')

    beam = pyuvsim.AnalyticBeam('gaussian', sigma=0.0222)
    beam_list = [beam]
    uvdata_out = pyuvsim.uvsim.run_uvsim(input_uv, beam_list=beam_list,
                                         mock_arrangement=args.mock_arrangement,
                                         max_za=args.max_za,
                                         save_catalog=args.save_catalog,
                                         Nsrcs=args.Nsrcs)
    if rank == 0:
        outfile = os.path.join(args.outdir, 'sim_' + os.path.basename(filename))
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        uvdata_out.write_uvfits(outfile, force_phase=True, spoof_nonessential=True)
