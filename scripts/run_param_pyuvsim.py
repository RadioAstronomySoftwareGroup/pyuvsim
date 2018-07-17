#!/usr/bin/env python

import pyuvsim
import argparse
import os
import numpy as np
import yaml
from mpi4py import MPI
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim import simsetup


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation from a parameter file."))

parser.add_argument('-p', '--paramsfile', dest='paramsfile', type=str, help='Parameter yaml file.')
args = vars(parser.parse_args())

if 'paramsfile' not in args:
    raise KeyError("Parameter file required")

with open(args['paramsfile'], 'r') as pfile:
    params = yaml.safe_load(pfile)

if params is None:
    params = {}

if rank == 0:

    if 'uvfile' in params:
        # simulate from a uvfits file if one is specified in the param file.

        filename = params['uvfile']
        print("Reading:", os.path.basename(filename))
        input_uv = UVData()
        input_uv.read_uvfits(filename)
        beam_ids = params['beam_files'].keys()
        beamfits_files = params['beam_files'].values()
        beam_list = []
        for bf in beamfits_files:
            uvb = UVBeam()
            uvb.read_beamfits()
            beam_list.append(bf)

        beam_list = (np.array(beam_list)[beam_ids]).tolist()
        outfile_name = os.path.join(params['outdir'], params['outfile_prefix'] + "_" + os.path.basename(filename))
        outfile_name = outfile_name + ".uvfits"

    else:
        # Not running off a uvfits.

        input_uv, beam_list, beam_ids = simsetup.initialize_uvdata_from_params(params)

        file_params = params['filing']
        params.update(file_params)
        del params['filing']
        if 'outfile_name' not in params or params['outfile_name'] == '':
            outfile_prefix = ""
            outfile_suffix = "_results"
            if 'outfile_prefix' in params:
                outfile_prefix = params['outfile_prefix'] + "_"
            if 'outfile_suffix' in params:
                outfile_suffix = "_" + params['outfile_suffix']
            outfile_name = os.path.join(params['outdir'], outfile_prefix
                                        + os.path.basename(args['paramsfile'])[:-5]
                                        + outfile_suffix)  # Strip .yaml extention
        else:
            outfile_name = params['outfile_name']

        outfile_name = outfile_name + ".uvfits"

        if 'clobber' not in params:
            outfile_name = simsetup.check_file_exists_and_increment(outfile_name)

        source_params = params['sources']
        if source_params['catalog'] == 'mock':
            params['mock_arrangement'] = source_params['mock_arrangement']
        if source_params['catalog'].endswith('txt'):
            catalog = simsetup.point_sources_from_params(source_params['catalog'])
        else:
            catalog = None



uvdata_out = pyuvsim.uvsim.run_uvsim(input_uv, beam_list=beam_list, mock_arrangement=params['mock_arrangement'], catalog=catalog)

if rank == 0:
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    uvdata_out.write_uvfits(outfile_name, force_phase=True, spoof_nonessential=True)
