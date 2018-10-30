# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import os
import numpy as np
import psutil
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
import yaml
from pyuvsim import simsetup, create_mock_catalog, uvsim
from guppy import hpy

from mpi4py import MPI

# Initialize MPI variables
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print("rank = ", rank, "PID", os.getpid())
time = Time('2018-03-01 00:00:00', scale='utc')
array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
# pfile = '../data/5km_triangle_2time_1chan.yaml'
# pfile = '../data/5km_triangle_10time_1chan.yaml'
# pfile = '../data/5km_triangle_1time_1chan_6ant.yaml'
pfile = '../data/laptop_size_sim.yaml'
params = yaml.safe_load(open(pfile))
print(pfile)
input_uv, beam_list, beam_ids = simsetup.initialize_uvdata_from_params(params)
# for a range of Nsrcs measure ram usage
# Ns = np.array([1,10,100])
Ns = [10]
process = psutil.Process(os.getpid())
for Nsrcs in Ns:
    # make the sources
    catalog = create_mock_catalog(time, 'long-line', Nsrcs=Nsrcs)
    # run the thing
    print("Nsrcs,", Nsrcs, "rank,", rank, "pre uvdata:", process.get_memory_info()[0] / 1e6, "MB")
    uvdata_out = uvsim.run_uvsim(input_uv, beam_list=beam_list, catalog=catalog)
    print("Nsrcs,", Nsrcs, "rank,", rank, "post uvdata:", process.get_memory_info()[0] / 1e6, "MB")
