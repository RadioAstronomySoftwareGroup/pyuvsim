import os
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units
from pyuvdata.uvbeam import UVBeam
from pyuvdata.data import DATA_PATH
import uvsim


beam_file = os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')


def test_single_source():
    """Test single zenith source"""
    time = Time('2018-3-1 00:00:00')

    array_location = EarthLocation(latitude='-30d43m17.5s', longitude='21d25m41.9s',
                                   altitude=1073.)
    lst = time.sidereal_time('mean', longitude=array_location.lon)

    freq = (150e6 * units.Hz)
    source = uvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])

    ant_nums = range(2)

    antenna1 = uvsim.Antenna([0, 0, 0], 0)
    antenna2 = uvsim.Antenna([107, 0, 0], 0)

    baseline = uvsim.Baseline(antenna1, antenna2)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]
    array = uvsim.Array(array_location, beam_list)

    task = uvsim.UVTask(source, time, freq, baseline, array)

    engine = uvsim.UVEngine([task])

    visibilities = engine.make_visibility()
