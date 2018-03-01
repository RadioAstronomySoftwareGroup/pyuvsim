import os
import numpy as np
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle
from astropy import units
from pyuvdata.uvbeam import UVBeam
from pyuvdata.data import DATA_PATH
import pyuvsim


beam_file = os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')


def test_single_source():
    """Test single zenith source"""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location
    lst = time.sidereal_time('mean')

    freq = (150e6 * units.Hz)
    source = pyuvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])

    ant_nums = range(2)

    antenna1 = pyuvsim.Antenna(np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna(np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]
    array = pyuvsim.Array(array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine([task])

    visibilities = engine.make_visibility()
