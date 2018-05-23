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


def create_zenith_source(time, name):
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
    return pyuvsim.Source(name, ra, dec, time, freq, [1, 0, 0, 0])


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
    zenith_source = create_zenith_source(time,'zensrc')
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

    source = pyuvsim.Source('testsrc',lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])

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
    source = create_zenith_source(time,'zensrc')

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 0, 0]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array([107, 0, 0]), 0)

    baseline = pyuvsim.Baseline(antenna1, antenna2)

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]
    array = pyuvsim.Telescope('telescope_name',array_location, beam_list)

    task = pyuvsim.UVTask(source, time, freq, baseline, array)

    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([.5, .5, 0, 0])))


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

    antenna1 = pyuvsim.Antenna('ant1', 1, np.array(antpos[0, :]), 0)
    antenna2 = pyuvsim.Antenna('ant2', 2, np.array(antpos[1, :]), 0)

    # setup the things that don't come from pyuvdata:
    # make a source at zenith
    time.location = array_location
    lst = time.sidereal_time('mean')
    # source = pyuvsim.Source(lst, Angle(array_location.lat), time, freq, [1, 0, 0, 0])
    source = create_zenith_source(time,'zensrc')

    # don't actually use a beam object for now because the thing you need to
    # calculate on the beam (the jones matrix at the source location) is bypassed for now
    beam_list = [0]

    baseline = pyuvsim.Baseline(antenna1, antenna2)
    array = pyuvsim.Telescope('telescope_name',array_location, beam_list)
    task = pyuvsim.UVTask(source, time, freq, baseline, array)
    engine = pyuvsim.UVEngine(task)

    visibility = engine.make_visibility()

    nt.assert_true(np.allclose(visibility, np.array([.5, .5, 0, 0])))


def test_file_to_tasks():

    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time,'zensrc')])

    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources)

    tel_loc = EarthLocation.from_geocentric(*hera_uv.telescope_location, unit='m')
    # beam list should be a list of UVBeam objects once we start using them
    beam_list = [0]
    telescope = pyuvsim.Telescope(hera_uv.telescope_name, tel_loc, beam_list)

    ant_pos = hera_uv.antenna_positions + hera_uv.telescope_location
    ant_pos_enu = uvutils.ENU_from_ECEF(ant_pos.T,
                                        *hera_uv.telescope_location_lat_lon_alt).T

    expected_task_list = []
    antenna_names = hera_uv.antenna_names
    antennas = []
    for num, antname in enumerate(antenna_names):
        beam_id = 0
        antennas.append(pyuvsim.Antenna(antname, num, ant_pos_enu[num], beam_id))

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
        task = pyuvsim.UVTask(sources[0], time, hera_uv.freq_array[0,0] * units.Hz, baseline, telescope)
        task.uvdata_index = (idx,0,0)
        expected_task_list.append(task)

    for idx, task in enumerate(uvtask_list):
        exp_task = expected_task_list[idx]
        nt.assert_equal(task, exp_task)


def test_uvdata_init():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    hera_uv.unphase_to_drift()
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time,'zensrc')])

    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources)

    uvdata_out = pyuvsim.initialize_uvdata(uvtask_list)

    hera_uv.data_array = np.zeros_like(hera_uv.data_array,dtype=np.complex)
    hera_uv.flag_array = np.zeros_like(hera_uv.data_array,dtype=bool)
    hera_uv.nsample_array = np.ones_like(hera_uv.data_array)
    hera_uv.history = 'UVSim'
    hera_uv.instrument = hera_uv.telescope_name
    hera_uv.integration_time = 1.

    ## FIX once pyuvdata gets rid of pyephem
    uvdata_out.zenith_ra = hera_uv.zenith_ra

    ## TODO Fix the test file's antenna_positions.
    nt.assert_equals(hera_uv._antenna_positions,uvdata_out._antenna_positions)
    nt.assert_true(uvdata_out.__eq__(hera_uv,check_extra=False))

def test_gather():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([
            create_zenith_source(time,'src1')
            # create_zenith_source(time,'src2'),
            # create_zenith_source(time,'src3')
            ])
    uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources)

    for task in uvtask_list:
        engine = pyuvsim.UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uv_out = pyuvsim.serial_gather(uvtask_list)

    nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array))

def test_run_serial_uvsim():
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)
    uv_out = pyuvsim.run_serial_uvsim(hera_uv, catalog=None, Nsrcs=1)

    nt.assert_true(np.allclose(uv_out.data_array, hera_uv.data_array))

def test_sources_equal():
    time = Time('2018-03-01 00:00:00', scale='utc')
    src1 = create_zenith_source(time,'src')
    src2 = create_zenith_source(time,'src')
    nt.assert_equal(src1,src2)
