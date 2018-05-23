import os
import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from astropy import units
from pyuvdata import UVBeam, UVData
import pyuvdata.utils as uvutils
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
from mpi4py import MPI

# Initialize MPI variables
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

beam_file = os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')


def test_run_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    uv_out = pyuvsim.run_uvsim(hera_uv, catalog=None, Nsrcs=3)
    if rank ==0:
        nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array))
