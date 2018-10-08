#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import argparse
import os
import numpy as np
import warnings
from mpi4py import MPI

from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

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
parser.add_argument('--min_alt', type=float, default=99.0, help='Maximum zenith angle for mock arrangements')
parser.add_argument('--save_catalog', action='store_true', default=False, help='Save catalog')


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
    extra_keywords = {'obs_param_file': 'uvfits_file=' + os.path.basename(filename),
                      'telescope_config_file': beam.type,
                      'antenna_location_file': os.path.basename(filename)}
    input_uv.extra_keywords = extra_keywords

    mock_keywords = {}
    mock_keywords['save'] = args.save_catalog
    mock_keywords['min_alt'] = args.save_catalog
    mock_keywords['arrangement'] = args.mock_arrangement

    uvdata_out = pyuvsim.uvsim.run_uvsim(input_uv, beam_list=beam_list,
                                         mock_keywords=mock_keywords)

    if rank == 0:
        outfile = os.path.join(args.outdir, 'sim_' + os.path.basename(filename))
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        uvdata_out.write_uvfits(outfile, force_phase=True, spoof_nonessential=True)
