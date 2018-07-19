import os
import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from astropy import units
import pyuvdata
from pyuvdata import UVBeam, UVData
import pyuvdata.utils as uvutils
from pyuvdata.data import DATA_PATH
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
import astropy.constants as const
from memory_profiler import profile



time = Time('2018-03-01 00:00:00', scale='utc')
array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)

#for a range of Nsrcs measure ram usage
Ns = np.array([1,10,100])
for Nsrcs in Ns:
    #make the sources
    sources = pyuvsim.create_mock_catalog(time,'long-line',Nsrcs=Nsrcs)
    
