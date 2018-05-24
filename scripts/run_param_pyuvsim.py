#!/usr/bin/env python

import pyuvsim
import argparse
import os, numpy as np
import yaml
from mpi4py import MPI
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


def fill_defaults(p):
    """
        If value is unset in the parameters dictionary p, set it to a default value.
        
    """
    assert( isinstance(p,dict) )
    sim_required_params = {'catalog' : 'mock', 'outdir' : '.', 'outfile_prefix' : 'sim', 'beam_file' : os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')}
    req_params_uvfile = { 'uvfile': os.path.join(SIM_DATA_PATH,'28mEWbl_1time_1chan.uvfits')}    #Required params for simulating off a uvfile
    req_params_nofile = { 'start_time' : None,
                            'end_time' : None,
                            'integration_time' : 1.0,
                            'Ntimes' : None,
                            'start_freq': None,
                            'end_freq' : None,
                            'Nfreqs' : None,
                            'channel_width': None,
                        }
    
    if  'uvfile' in p.keys() or len(p.keys()) == 0:
        sim_required_params.update(req_params_uvfile)
    else:
        sim_required_params.update(req_params_nofile)
    
    # Loop over the req_params_use and fill values in p that are missing.

    for k in sim_required_params.keys():
        if not k in p.keys():
            print("Missing parameter "+str(k)+", defaulting to "+str(sim_required_params[k]))
            p[k] = sim_required_params[k]

    # Certain parameters cannot be missing and have no default. Error if they're None.
    for k in p.keys():
        if p[k] is None: raise ValueError("Parameter "+str(k)+" must be defined.")

    if (p['catalog'] == 'mock') and ('mock_arrangement' not in p):
        print("Using mock catalog and arrangement unset. Defaulting to \'zenith\'")
        p['mock_arrangement'] = 'zenith'
        p['Nsrcs'] = 1

    return p

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to execute a pyuvsim simulation from a parameter file."))

parser.add_argument('-p','--paramsfile',dest='paramsfile', type=str, help='Parameter yaml file.')
args = vars(parser.parse_args())

if 'paramsfile' not in args:
    raise KeyError("Parameter file required")

with open(args['paramsfile'], 'r') as pfile:
    params = yaml.safe_load(pfile)

if params is None:
    params = {}

print params
params = fill_defaults(params)

if 'uvfile' in params:
    ### simulate from a uvfits file
    filename = params['uvfile']
    print("Reading:", os.path.basename(filename))
    input_uv = UVData()
    input_uv.read_uvfits(filename)
    beam = UVBeam()
    beam.read_cst_beam(params['beam_file'], beam_type='efield', frequency=150e6,
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')

    beam_list = [beam]
    uvdata_out = pyuvsim.uvsim.run_uvsim(input_uv, beam_list=beam_list, mock_arrangement=params['mock_arrangement'], Nsrcs=params['Nsrcs'])
    if rank == 0:
        outfile = os.path.join(params['outdir'], params['outfile_prefix'] + "_" + os.path.basename(filename))
        if not os.path.exists(params['outdir']): os.makedirs(params['outdir'])
        uvdata_out.write_uvfits(outfile, force_phase=True, spoof_nonessential=True)
