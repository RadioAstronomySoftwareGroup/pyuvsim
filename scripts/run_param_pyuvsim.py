#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import pyuvsim
import argparse
import os
import numpy as np
import yaml
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim import simsetup
from astropy.coordinates import EarthLocation
from pyuvsim import mpi

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation from a parameter file."))

parser.add_argument('-p', '--paramsfile', dest='paramsfile', type=str, help='Parameter yaml file.')
args = vars(parser.parse_args())

if 'paramsfile' not in args:
    raise KeyError("Parameter file required")

with open(args['paramsfile'], 'r') as pfile:
    params = yaml.safe_load(pfile)

params['config_path'] = '.'

if params is None:
    params = {}

input_uv = UVData()

Nbeams = 0
beam_list = None
beam_dict = None
input_uv = UVData()
mock_keywords = None
catalog = None

mpi.start_mpi()
rank = mpi.get_rank()

if rank == 0:
    if 'uvfile' in params:
        # simulate from a uvfits file if one is specified in the param file.

        filename = params['uvfile']
        print("Reading:", os.path.basename(filename))
        input_uv.read_uvfits(filename)

        if 'beam_files' in params:
            beam_ids = params['beam_files'].keys()
            beamfits_files = params['beam_files'].values()
            beam_list = []
            for bf in beamfits_files:
                beam_list.append(bf)

        beam_list = (np.array(beam_list)[beam_ids]).tolist()
        Nbeams = len(beam_list)
        outfile_name = os.path.join(params['outdir'], params['outfile_prefix'] + "_" + os.path.basename(filename))
        outfile_name = outfile_name + ".uvfits"

    else:
        # Not running off a uvfits.
        print("Simulating from parameters")
        input_uv, beam_list, beam_dict, beam_ids = simsetup.initialize_uvdata_from_params(params)
        Nbeams = len(beam_list)
        print("Nfreqs: ", input_uv.Nfreqs)
        print("Ntimes: ", input_uv.Ntimes)
        source_params = params['sources']
        if source_params['catalog'] == 'mock':
            mock_keywords = {'time': input_uv.time_array[0], 'arrangement': source_params['mock_arrangement'],
                             'array_location': EarthLocation.from_geocentric(*input_uv.telescope_location, unit='m')}
            extra_mock_kwds = ['time', 'Nsrcs', 'zen_ang', 'save', 'max_za']
            for k in extra_mock_kwds:
                if k in source_params.keys():
                    mock_keywords[k] = source_params[k]
            catalog = 'mock'

        if 'catalog' in source_params:
            catalog = source_params['catalog']
        else:
            catalog = None
# mpi.start_mpi()
#beam_dict = mpi.get_comm().bcast(beam_dict, root=0)
uvdata_out = pyuvsim.uvsim.run_uvsim(input_uv, beam_list=beam_list, beam_dict=beam_dict, catalog_file=catalog, mock_keywords=mock_keywords)


if rank == 0:
    simsetup.write_uvfits(uvdata_out, params)
