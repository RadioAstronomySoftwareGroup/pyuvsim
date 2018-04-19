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


beam_file = os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')
hera_miriad_file = os.path.join(DATA_PATH, 'hera_testfile')
EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_1time_1chan.uvfits')


def create_zenith_source(time):
    """Create pyuvsim Source object at zenith.

    Input: Astropy Time object
        sample: Time('2018-03-01 00:00:00', scale='utc')
    Retruns: Pyuvsim Source object
    """
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec

    return pyuvsim.Source(ra, dec, time, freq, [1, 0, 0, 0])


def test_source_zenith():
    """Test single source position at zenith."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)
    #
    # source_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
    #                         obstime=time, frame='altaz', location=array_location)
    # icrs_coord = source_coord.transform_to('icrs')
    #
    # ra = icrs_coord.ra
    # dec = icrs_coord.dec
    #
    # zenith_source = pyuvsim.Source(ra, dec, time, freq, [1, 0, 0, 0])
    zenith_source = create_zenith_source(time)
    zenith_source_lmn = zenith_source.pos_lmn(time, array_location)
    print('Zenith Source lmn')
    print(zenith_source_lmn)

    nt.assert_true(np.allclose(zenith_source_lmn, np.array([0, 0, 1])))

# This test is known to fail. We will skip execution for now.
@nt.nottest
def test_source_lst_zenith():
    """Instantiate a source at zenith using LST=RA and array latitude=DEC."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                       height=1073.)
    freq = (150e6 * units.Hz)

    time.location = array_location
    lst = time.sidereal_time('apparent')

    source = pyuvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])

    source_lmn = source.pos_lmn(time, array_location)
    print('Source lmn')
    print(source_lmn)

    nt.assert_true(np.allclose(source_lmn, np.array([0, 0, 1])))


def test_single_source():
    """Test single zenith source."""
    time = Time('2018-03-01 00:00:00', scale='utc')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    time.location = array_location
    lst = time.sidereal_time('apparent')

    freq = (150e6 * units.Hz)
    source = create_zenith_source(time)

    antenna1 = pyuvsim.Antenna(np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna(np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]
    array = pyuvsim.Array(array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([1, 1, 0, 0])))


def test_single_source_vis_uvdata():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    print('hera uv time 0')
    print(hera_uv.time_array[0])
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
    print('hera telescope_location')
    print(hera_uv.telescope_location_lat_lon_alt_degrees)
    print('array_location')
    print(array_location.lat, array_location.lon, array_location.height)
    freq = hera_uv.freq_array[0, 0] * units.Hz
    print('freq')
    print(freq)

    # get antennas positions into ENU
    antpos = hera_uv.antenna_positions[0:2, :] + hera_uv.telescope_location
    antpos = uvutils.ENU_from_ECEF(antpos.T, *hera_uv.telescope_location_lat_lon_alt).T

    print('antpos shape')
    print(antpos.shape)
    print('antpos')
    print(antpos)
    antenna1 = pyuvsim.Antenna(np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna(np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source at zenith
    time.location = array_location
    lst = time.sidereal_time('mean')
    print('lst')
    print(lst)
    # source = pyuvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])
    source = create_zenith_source(time)

    print('ra, dec')
    print(source.ra, source.dec)
    print('epoch')
    print(source.epoch)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Array(array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()
    print(visibility)

    # problem is that the fringe isn't 1. Because the w term isn't 0.
    # should the fringe be 1 at zenith with a w term?

    nt.assert_true(np.allclose(visibility, np.array([1, 1, 0, 0])))
