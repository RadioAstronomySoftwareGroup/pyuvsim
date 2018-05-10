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

    antenna1 = pyuvsim.Antenna('ant1', np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]
    array = pyuvsim.Telescope(array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([1, 1, 0, 0])))


def test_single_source_vis_uvdata():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    array_location = EarthLocation.from_geocentric(hera_uv.telescope_location[0],
                                                   hera_uv.telescope_location[1],
                                                   hera_uv.telescope_location[2],
                                                   unit='m')
    freq = hera_uv.freq_array[0, 0] * units.Hz

    # get antennas positions into ENU
    antpos = hera_uv.antenna_positions[0:2, :] + hera_uv.telescope_location
    antpos = uvutils.ENU_from_ECEF(antpos.T, *hera_uv.telescope_location_lat_lon_alt).T

    antenna1 = pyuvsim.Antenna('ant1', np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source at zenith
    time.location = array_location
    lst = time.sidereal_time('mean')
    # source = pyuvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])
    source = create_zenith_source(time)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Telescope(array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([1, 1, 0, 0])))


def test_file_to_tasks():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time)])

    uvtask_list = pyuvsim.uvfile_to_task_list(EW_uvfits_file, sources)

    tel_loc = EarthLocation.from_geocentric(*hera_uv.telescope_location, unit='m')
    # beam list should be a list of UVBeam objects once we start using them
    beam_list = [0]
    telescope = pyuvsim.Telescope(tel_loc, beam_list)

    ant_pos = hera_uv.antenna_positions + hera_uv.telescope_location
    ant_pos_enu = uvutils.ENU_from_ECEF(ant_pos.T,
                                        *hera_uv.telescope_location_lat_lon_alt).T

    expected_task_list = []
    antenna_names = hera_uv.antenna_names
    antennas = []
    for num, antname in enumerate(antenna_names):
        beam_id = 0
        antennas.append(pyuvsim.Antenna(antname, ant_pos_enu[num], beam_id))

    antennas1 = []
    for antnum in hera_uv.ant_1_array:
        index = np.where(hera_uv.antenna_numbers == antnum)[0][0]
        antennas1.append(antennas[index])

    antennas2 = []
    for antnum in hera_uv.ant_2_array:
        index = np.where(hera_uv.antenna_numbers == antnum)[0][0]
        antennas2.append(antennas[index])

    for idx, antenna1 in enumerate(antennas1):
        antenna2 = antennas2[idx]
        baseline = pyuvsim.Baseline(antenna1, antenna2)
        task = pyuvsim.UVTask(sources[0], time, hera_uv.freq_array[0] * units.Hz, baseline, telescope)
        expected_task_list.append(task)

    for idx, task in enumerate(uvtask_list):
        exp_task = expected_task_list[idx]
        nt.assert_equal(task.baseline.antenna1.name, exp_task.baseline.antenna1.name)
        nt.assert_true(np.allclose(task.baseline.antenna1.pos_enu.to(units.m).value,
                                   exp_task.baseline.antenna1.pos_enu.to(units.m).value))
        nt.assert_equal(task.baseline.antenna1.beam_id, exp_task.baseline.antenna1.beam_id)
        nt.assert_equal(task.baseline.antenna2.name, exp_task.baseline.antenna2.name)
        nt.assert_true(np.allclose(task.baseline.antenna2.pos_enu.to(units.m).value,
                                   exp_task.baseline.antenna2.pos_enu.to(units.m).value))
        nt.assert_equal(task.baseline.antenna2.beam_id, exp_task.baseline.antenna2.beam_id)

        nt.assert_equal(task.time, exp_task.time)
        nt.assert_equal(task.freq, exp_task.freq)

        nt.assert_equal(task.source.ra, exp_task.source.ra)
        nt.assert_equal(task.source.dec, exp_task.source.dec)
        nt.assert_equal(task.source.epoch, exp_task.source.epoch)
        nt.assert_equal(task.source.stokes, exp_task.source.stokes)

        nt.assert_equal(task.telescope.telescope_location, exp_task.telescope.telescope_location)
        # this may break when we make beam_list contain UVBeam objects
        nt.assert_equal(task.telescope.beam_list, exp_task.telescope.beam_list)


# def test_loopback_file_to_sim():
