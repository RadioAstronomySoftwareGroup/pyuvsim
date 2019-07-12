#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import argparse
import numpy as np
import yaml
import os
import sys
import resource
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim import mpi, profiling, simsetup, uvsim

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation "
                                              "for profiling purposes."))

paramsfile = os.path.join(SIM_DATA_PATH, 'profiling_params.yaml')
cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, f) for f in cst_files]

parser.add_argument('--Nsrcs', dest='Nsrcs', type=int, default=1)
parser.add_argument('--Ntimes', dest='Ntimes', type=int, default=1)
parser.add_argument('--Nfreqs', dest='Nfreqs', type=int, default=1)
parser.add_argument('--Nbls', dest='Nbls', type=int, default=1)
parser.add_argument('--beam', dest='beam', type=str, default='uniform')
parser.add_argument('--prof_out', dest='prof_out', type=str, default='time_profile.out')
parser.add_argument('--mem_out', dest='mem_out', type=str, default='memory_usage.out')

args = parser.parse_args()

with open(paramsfile, 'r') as pfile:
    params = yaml.safe_load(pfile)

params['config_path'] = os.path.dirname(paramsfile)

min_alt = 70  # Degrees

mpi.start_mpi()
rank = mpi.get_rank()

beam_list = None
beam_dict = None
input_uv = UVData()
mock_keywords = None
catalog = None
profiling.set_profiler(outfile_name=args.prof_out)
if rank == 0:

    params['freq']['Nfreqs'] = args.Nfreqs
    params['time']['Ntimes'] = args.Ntimes
    params['sources'] = {'catalog': 'mock'}

    input_uv, beam_list, beam_dict = simsetup.initialize_uvdata_from_params(params)

    # Baseline selection:
    input_uv.baseline_array = np.repeat(input_uv.baseline_array[:args.Nbls], args.Ntimes)
    input_uv.ant_1_array, input_uv.ant_2_array = input_uv.baseline_to_antnums(
        input_uv.baseline_array)
    ants_new = np.unique(input_uv.ant_1_array.tolist() + input_uv.ant_2_array.tolist())
    input_uv.antenna_numbers = ants_new
    # Antnames/numbers are going to be messed up by the baseline selection. Unimportant.
    input_uv.antenna_names = ants_new.astype(str)
    Nants = ants_new.size
    beam_dict = dict(zip(input_uv.antenna_names, np.zeros(Nants, dtype=int))
                     )  # For now, all use the same beam model
    input_uv.antenna_positions = input_uv.antenna_positions[:Nants, :]
    input_uv.Nants_data = Nants
    input_uv.Nants_telescope = Nants

    # Time selection:
    inds = np.array([np.arange(args.Nbls) + i
                     * input_uv.Nbls for i in range(args.Ntimes)]).flatten()
    input_uv.time_array = input_uv.time_array[inds]
    input_uv.lst_array = input_uv.lst_array[inds]
    input_uv.Nbls = args.Nbls
    input_uv.Nblts = args.Nbls * args.Ntimes

    # Beam selection:
    # Default is uniform
    if args.beam == 'hera':
        beam = UVBeam()
        beamfile = '/users/alanman/data/alanman/NickFagnoniBeams/HERA_NicCST_fullfreq.uvbeam'
        beam_list = [beamfile]

    mock_keywords = {'mock_arrangement': 'random', 'Nsrcs': args.Nsrcs,
                     'min_alt': min_alt, 'time': input_uv.time_array[0]}
    params['sources'].update(**mock_keywords)

    # Catalog setup
    catalog, _ = simsetup.initialize_catalog_from_params(params)

comm = mpi.world_comm
input_uv = comm.bcast(input_uv, root=0)
beam_list = comm.bcast(beam_list, root=0)
beam_dict = comm.bcast(beam_dict, root=0)
catalog = mpi.shared_mem_bcast(catalog, root=0)
if rank == 0:
    print("Starting simulation.")
    sys.stdout.flush()

uv_out = uvsim.run_uvdata_uvsim(input_uv, beam_list, beam_dict=beam_dict, catalog=catalog)

memory_usage_GB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
comm.Barrier()

memory_usage_GB = comm.gather(memory_usage_GB, root=0)

if rank == 0:
    memory_usage_GB = np.min(memory_usage_GB)
    print('Mem_usage: ' + str(memory_usage_GB))
    with open(args.mem_out, 'w') as memfile:
        memfile.write(str(memory_usage_GB))
