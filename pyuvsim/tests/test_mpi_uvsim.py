# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import os
import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from astropy import units
from pyuvdata import UVBeam, UVData
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
from mpi4py import MPI

# Initialize MPI variables
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
beam_files = [os.path.join(DATA_PATH, f) for f in cst_files]
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    beam = UVBeam()
    beam.read_cst_beam(beam_files, beam_type='efield', frequency=[100e6, 123e6],
                       telescope_name='HERA',
                       feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
                       model_name='E-field pattern - Rigging height 4.9m',
                       model_version='1.0')
    beam.write_beamfits("temp.uvbeam")
    beamfile = os.path.join(os.getcwd(), "temp.uvbeam")

    beam_list = [beamfile]
    mock_keywords = {"Nsrcs": 3}
    uv_out = pyuvsim.run_uvsim(hera_uv, beam_list, catalog_file=None, mock_keywords=mock_keywords,
                               uvdata_file=EW_uvfits_file)
    if rank == 0:
        nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array, atol=5e-3))
    os.remove(beamfile)
