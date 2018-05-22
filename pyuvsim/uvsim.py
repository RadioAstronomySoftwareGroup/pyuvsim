import numpy as np
import astropy.constants as const
import astropy.units as units
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from pyuvdata import UVData
import pyuvdata.utils as uvutils
# execution pseudo code

# if rank==0
    # load stuff
    # broadcast
    # gather
    # sum sources
    # save file
# else:
    # rankN
    # load array details
    # create object
    # recv
    # put details into object
    # object.calculate()
    # format results
    # send, return to recv

# TODO:
# loading array details, observation settings
# what does the MPI traffic look like?
#    who does the formatting to input mpi parameters into uvengine
#    ANS: each MPI process takes a list of inputs, creates a list of tasks, and
#      hands them to the UVEngine
# unittests

# FUTURE:
#   Baseline-dependent beams via Mueller matrix formalism

# read source catalog, generate Source objects,
#    set Source.calc to whatever
# 30 Nov 2017 - started source object, made _Blank_ test file
#    next steps: add basic parameter test for Source, flux calculation, ENU direction, coherency
#       -- HW Adam add some unittest boilerplate, next time testing Source attributes exist
#    after that, Baseline object, ECEF to ENU, empty beam hook
#    GOAL: 1 source on 1 baseline. -  first important unittest
# 7 Dec 2017 - Debated Jones vs Mueller and coordinate systems. Began defining Jones class and apply_beam method.
#    next steps: add in jones coordinate transformations and make unit tests for it.
#     refs: Smirnov Series 2011 papers. Nunhokee 2017
#    Future tie-in: UVBeam to get jones matrix option. Open github issue. Request support for coordinate transformations.
#    Zac --> Make us some diagrams of coord systems for different parts of the transformation.
# 18 Jan 2018
#    Zac -- Circulate what you have. Diagram coord systems.
#    Debate --- Do we rotate the beam or the sky? Calculate uvw in alt/az or ra/dec?
#    HW: Write down our requirements; Open up an issue to make first test (source at zenith).


class Source(object):
    freq = None
    stokes = None
    ra = None
    dec = None
    epoch = None
    coherency_radec = None

    def __init__(self, ra, dec, epoch, freq, stokes):
        """
        Initialize from source catalog

        Args:
            ra: astropy Angle object
                source RA at epoch
            dec: astropy Angle object
                source Dec at epoch
            epoch: astropy Time object
                epoch for RA/dec
            stokes:
                4 element vector giving the source [I, Q, U, V]
            freq: astropy quantity
                frequency of source catalog value
        """
        if not isinstance(ra, Angle):
            raise ValueError('ra must be an astropy Angle object. '
                             'value was: {ra}'.format(ra=ra))

        if not isinstance(dec, Angle):
            raise ValueError('dec must be an astropy Angle object. '
                             'value was: {dec}'.format(dec=dec))

        if not isinstance(freq, Quantity):
            raise ValueError('freq must be an astropy Quantity object. '
                             'value was: {f}'.format(f=freq))

        if not isinstance(epoch, Time):
            raise ValueError('epoch must be an astropy Time object. '
                             'value was: {t}'.format(t=epoch))

        self.freq = freq
        self.stokes = stokes
        self.ra = ra
        self.dec = dec
        self.epoch = epoch

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        self.coherency_radec = np.array([[self.stokes[0] + self.stokes[1],
                                          self.stokes[2] - 1j * self.stokes[3]],
                                         [self.stokes[2] + 1j * self.stokes[3],
                                          self.stokes[0] - self.stokes[1]]])

    def __eq__(self, other):
        return ((self.ra == other.ra) and
                (self.dec == other.dec) and
                (self.epoch == other.epoch) and
                (self.stokes == other.stokes) and
                (self.freq == other.freq))

    def coherency_calc(self, time, telescope_location):
        """
        Calculate the local coherency in az/za basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local az/za.

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            local coherency in az/za basis
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        if np.sum(np.abs(self.stokes[1:])) == 0:
            rotation_matrix = np.array([[1, 0], [0, 1]])
        else:
            # First need to calculate the sin & cos of the parallactic angle
            # See Meeus's astronomical algorithms eq 14.1
            # also see Astroplan.observer.parallactic_angle method
            time.location = telescope_location
            lst = time.sidereal_time('mean')
            hour_angle = (lst - self.ra).rad
            sinX = np.sin(hour_angle)
            cosX = np.tan(telescope_location.lat) * np.cos(self.dec) - np.sin(self.dec) * np.cos(hour_angle)

            rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]])

        print('rotation_matrix')
        print(rotation_matrix)
        coherency_local = np.einsum('ab,bc,cd->ad', rotation_matrix.T,
                                    self.coherency_radec, rotation_matrix)
        print('coherency_local')
        print(coherency_local)

        return coherency_local

    def az_za_calc(self, time, telescope_location):
        """
        calculate the azimuth & zenith angle for this source at a time & location

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            (azimuth, zenith_angle) in radians
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(telescope_location, EarthLocation):
            raise ValueError('telescope_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=telescope_location))

        source_coord = SkyCoord(self.ra, self.dec, frame='icrs', equinox=self.epoch)

        source_altaz = source_coord.transform_to(AltAz(obstime=time, location=telescope_location))

        az_za = (source_altaz.az.rad, source_altaz.zen.rad)
        return az_za

    def pos_lmn(self, time, telescope_location):
        """
        calculate the direction cosines of this source at a time & location

        Args:
            time: astropy Time object
            telescope_location: astropy EarthLocation object

        Returns:
            (l, m, n) direction cosine values
        """
        # calculate direction cosines of source at current time and array location
        az_za = self.az_za_calc(time, telescope_location)

        print('Source az_za')
        print(az_za)

        pos_l = np.cos(az_za[0]) * np.sin(az_za[1])
        pos_m = np.sin(az_za[0]) * np.sin(az_za[1])
        pos_n = np.cos(az_za[1])
        return (pos_l, pos_m, pos_n)

    # debate: will an object ever be _re_-used?
    # answer, that is not our initial intent. Remake objects for each t,f etc
    # new plan from 2/1: Source objects describe sources in celestial coords
    # but have methods to convert those coords to the local az/za frame


class Telescope(object):
    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.telescope_location = telescope_location
        self.telescope_name = telescope_name

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list

    def __eq__(self, other):
        return (self.telescope_location == other.telescope_location) and (self.beam_list == other.beam_list) and (self.telescope_name == other.telescope_name)


class Antenna(object):
    def __init__(self, name, number, enu_position, beam_id):
        self.name = name
        self.number = number
        # ENU position in meters relative to the telescope_location
        self.pos_enu = enu_position * units.m
        # index of beam for this antenna from array.beam_list
        self.beam_id = beam_id

    def get_beam_jones(self, array, source_lmn, frequency):
        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in Az/ZA
        # return array.beam_list[self.beam_id].get_direction_jones(source_lmn, frequency)
        return np.array([[1, 0], [0, 1]])

    def __eq__(self, other):
        return ((self.name == other.name) and
                np.all(self.pos_enu == other.pos_enu) and
                (self.beam_id == other.beam_id))

# Then there's the issue of how to pass this in the MPI context.
#   one option is for the UVEngine to get an array object and the tasks have
#   references to that object (i.e. antenna number which the array object can map to a beam)
#   Otherwise the task might need to carry the whole array, which would lead to
#   many copies of the array being passed to the UVEngine


class Baseline(object):
    def __init__(self, antenna1, antenna2):
        self.antenna1 = antenna1
        self.antenna2 = antenna2
        self.enu = antenna1.pos_enu - antenna2.pos_enu
        # we're using the local az/za frame so uvw is just enu
        self.uvw = self.enu

    def __eq__(self, other):
        return ((self.antenna1 == other.antenna1) and
                (self.antenna2 == other.antenna2) and
                np.all(self.enu == other.enu) and
                np.all(self.uvw == other.uvw))


class UVTask(object):
    # holds all the information necessary to calculate a single src, t, f, bl, array
    # need the array because we need an array location for mapping to locat az/za
    def __init__(self, source, time, freq, baseline, telescope):
        self.time = time
        self.freq = freq
        self.source = source
        self.baseline = baseline
        self.telescope = telescope
        self.visibility_vector = None

    def __eq__(self, other):
        return ((self.time == other.time) and
                (self.freq == other.freq) and
                (self.source == other.source) and
                (self.baseline == other.baseline) and
                (self.telescope == other.telescope))


class UVEngine(object):
    # inputs x,y,z,flux,baseline(u,v,w), time, freq
    # x,y,z in same coordinate system as uvws

    def __init__(self, task):   # task_array  = list of tuples (source,time,freq,uvw)
        # self.rank
        self.task = task
        # construct self based on MPI input

    # Debate --- Do we allow for baseline-defined beams, or stick with just antenna beams?
    #   This would necessitate the Mueller matrix formalism.
    #   As long as we stay modular, it should be simple to redefine things.

    def apply_beam(self):
        # Supply jones matrices and their coordinate. Knows current coords of visibilities, applies rotations.
        # for every antenna calculate the apparent jones flux
        # Beam --> Takes coherency matrix alt/az to ENU
        baseline = self.task.baseline
        source = self.task.source
        # coherency is a 2x2 matrix
        # (Ex^2 conj(Ex)Ey, conj(Ey)Ex Ey^2)
        # where x and y vectors along the local za/az coordinates.
        beam1_jones = baseline.antenna1.get_beam_jones(self.task.telescope, source.pos_lmn, self.task.freq)
        beam2_jones = baseline.antenna2.get_beam_jones(self.task.telescope, source.pos_lmn, self.task.freq)
        this_apparent_coherency = np.dot(beam1_jones,
                                         source.coherency_calc(self.task.time,
                                                               self.task.telescope.telescope_location))
        this_apparent_coherency = np.dot(this_apparent_coherency,
                                         (beam2_jones.conj().T))

        self.apparent_coherency = this_apparent_coherency

    def make_visibility(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.apply_beam()

        uvw_lambda = self.task.baseline.uvw * self.task.freq.to(1 / units.s) / (const.c)
        pos_lmn = self.task.source.pos_lmn(self.task.time, self.task.telescope.telescope_location)

        print('Pos lmn')
        print(pos_lmn)

        # This hard coding shouldn't be required.
        # pos_lmn = [0, 0, 1]

        fringe = np.exp(-2j * np.pi * np.dot(self.task.baseline.uvw, pos_lmn))
        pos_lmn = self.task.source.pos_lmn(self.task.time, self.task.telescope.telescope_location)

        print('fringe')
        print(fringe)

        vij = self.apparent_coherency * fringe

        # need to reshape to be [xx, yy, xy, yx]
        vis_vector = [vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]]

        return vis_vector

    def update_task(self):
        self.task.visibility_vector = self.make_visibility()


def uvfile_to_task_list(filename, sources, beam_dict=None):
    """Create task list from pyuvdata compatible input file.

    Returns: List of task parameters to be send to UVEngines
    List has task parameters defined in UVTask object
    This function extracts time, freq, Antenna1, Antenna2
    """

    if not isinstance(sources, np.ndarray):
        raise TypeError("sources must be a numpy array")

    input_uv = UVData()
    input_uv.read_uvfits(filename)

    freq = input_uv.freq_array[0, :] * units.Hz
    # TODO: Have this function make mega-list of [time,frequency, Ant1, Ant2]
    #  this may possibly be done with MPI-trickery later

    # beam_list should be a list of beam objects, once we have those.
    beam_list = [0]
    telescope = Telescope(input_uv.telescope_name,EarthLocation.from_geocentric(*input_uv.telescope_location, unit='m'), beam_list)

    times = Time(input_uv.time_array, format='jd', location=telescope.telescope_location)

    antpos_ECEF = input_uv.antenna_positions + input_uv.telescope_location
    antpos_ENU = uvutils.ENU_from_ECEF(antpos_ECEF.T,
                                       *input_uv.telescope_location_lat_lon_alt).T
    antenna_names = input_uv.antenna_names
    antennas = []
    for num, antname in enumerate(antenna_names):
        if beam_dict is None:
            beam_id = 0
        else:
            beam_id = beam_dict[antname]
        antennas.append(Antenna(antname, num, antpos_ENU[num], beam_id))

    antennas1 = []
    for antnum in input_uv.ant_1_array:
        index = np.where(input_uv.antenna_numbers == antnum)[0][0]
        antennas1.append(antennas[index])

    antennas2 = []
    for antnum in input_uv.ant_2_array:
        index = np.where(input_uv.antenna_numbers == antnum)[0][0]
        antennas2.append(antennas[index])

    baselines = []
    for index in range(len(antennas1)):
        baselines.append(Baseline(antennas1[index], antennas2[index]))
    baselines = np.array(baselines)

    blts_index = np.arange(input_uv.Nblts)
    frequency_index = np.arange(input_uv.Nfreqs)
    source_index = np.arange(len(sources))
    blts_ind, freq_ind, source_ind = np.meshgrid(blts_index, frequency_index, source_index)

    uvtask_list = []

    uvtask_params = zip(baselines[blts_ind].ravel(), freq[freq_ind].ravel(),
                        times[blts_ind].ravel(), sources[source_ind].ravel())
    for (bl, freq, t, source) in uvtask_params:
        task = UVTask(source, t, freq, bl, telescope)
        uvtask_list.append(task)

    return uvtask_list


def initialize_uvdata(uvtask_list):
    # writing 4 pol polarization files now

    # Assume all tasks have the same telescope.
    #   Enforce this generally!

    # Do after MPI_Reduce, summing over sources.

    task_freqs = []
    task_bls = []
    task_times = []
    task_antnames = []
    task_antnums = []
    task_antpos = []
    task_uvw = []
    ant_1_array = []
    ant_2_array = []
    telescope_name = uvtask_list[0].telescope.telescope_name
    telescope_location = uvtask_list[0].telescope.telescope_location.geocentric

    source_0 = uvtask_list[0].source
    freq_0 = uvtask_list[0].freq.to("Hz").value
    for task in uvtask_list:
        if not task.source == source_0: continue
        task_freqs.append(task.freq.to("Hz").value)

        if task.freq.to("Hz").value == freq_0:
            task_bls.append(task.baseline)
            task_times.append(task.time.jd)
            task_antnames.append(task.baseline.antenna1.name)
            task_antnames.append(task.baseline.antenna2.name)
            ant_1_array.append(task.baseline.antenna1.number)
            ant_2_array.append(task.baseline.antenna2.number)
            task_antnums.append(task.baseline.antenna1.number)
            task_antnums.append(task.baseline.antenna2.number)
            task_antpos.append(task.baseline.antenna1.pos_enu)
            task_antpos.append(task.baseline.antenna2.pos_enu)
            task_uvw.append(task.baseline.uvw)
    
    antnames, ant_indices = np.unique(task_antnames, return_index=True)
    task_antnums = np.array(task_antnums)
    task_antpos = np.array(task_antpos)
    antnums = task_antnums[ant_indices]
    antpos = task_antpos[ant_indices]

    freqs = np.unique(task_freqs)

    uv_obj = UVData()
    uv_obj.telescope_name = telescope_name
    uv_obj.telescope_location = np.array([tl.to('m').value for tl in telescope_location])
    uv_obj.instrument = telescope_name
    uv_obj.Nfreqs = freqs.size
    uv_obj.Ntimes = np.unique(task_times).size
    uv_obj.Nants_data = antnames.size
    uv_obj.Nants_telescope = uv_obj.Nants_data
    uv_obj.Nblts = len(ant_1_array)

    uv_obj.antenna_names = antnames.tolist()
    uv_obj.antenna_numbers = antnums
    antpos_ecef = uvutils.ECEF_from_ENU(antpos.T,*uv_obj.telescope_location_lat_lon_alt).T - uv_obj.telescope_location
    uv_obj.antenna_positions = antpos_ecef
    uv_obj.ant_1_array = np.array(ant_1_array)
    uv_obj.ant_2_array = np.array(ant_2_array)
    uv_obj.time_array = np.array(task_times)
    uv_obj.uvw_array = np.array(task_uvw)
    uv_obj.baseline_array = uv_obj.antnums_to_baseline(ant_1_array,ant_2_array)
    uv_obj.Nbls = np.unique(uv_obj.baseline_array).size
    if uv_obj.Nfreqs == 1:
        uv_obj.channel_width = 1.  #Hz
    else:
        uv_obj.channel_width = np.diff(freqs)[0]

    if uv_obj.Ntimes == 1:
        uv_obj.integration_time = 1.  #Second
    else:
        uv_obj.integration_time = np.diff(np.unique(task_times))[0]
    uv_obj.set_lsts_from_time_array()
    uv_obj.zenith_ra = uv_obj.lst_array
    uv_obj.zenith_dec = np.repeat(uv_obj.telescope_location_lat_lon_alt[0],uv_obj.Nblts)  #Latitude
    uv_obj.object_name = 'zenith'
    uv_obj.set_drift()
    uv_obj.vis_units = 'Jy'
    uv_obj.polarization_array = np.array([-5,-6,-7,-8])
    uv_obj.spw_array = np.array([0])
    uv_obj.freq_array = np.array([freqs])

    uv_obj.Nspws = uv_obj.spw_array.size
    uv_obj.Npols = uv_obj.polarization_array.size

    uv_obj.data_array = np.zeros((uv_obj.Nblts,uv_obj.Nspws,uv_obj.Nfreqs,uv_obj.Npols),dtype=np.complex)
    uv_obj.flag_array = np.zeros((uv_obj.Nblts,uv_obj.Nspws,uv_obj.Nfreqs,uv_obj.Npols),dtype=bool)
    uv_obj.nsample_array = np.ones_like(uv_obj.data_array,dtype=float)
    uv_obj.history = 'UVSim'

    uv_obj.check()

    return uv_obj
    

    


# TODO: make a gather function that puts the visibilities into a UVData object

# what a node does (pseudo code)
# def node:
#     if task_id > 0:
#         params_in = scatter()


# objects that are defined on input
# class Telescope(object):
#    telescope_location (maybe an astropy observer class)
#    contains antennas (locations), maps beams to antennas, cross-talk?
# class SkyModel(object):
#       """ """
#    can be derived from sources, a healpix map, or an image
#      (internally represented as lists of sources/components)
#      hands per freq, time, direction inputs to a broadcaster
# class Spectrum(object):
#
# class SkySelection(object):
#
# class Compute(object):
#
# class Observation(object):
#   defines user settable parameters
#   (start time, stop time, pointing, integration time)
# class Instrument(object):
#   limits of the instrument
#   available ranges of parameters
#   (frequencies, pointing ranges, integration times)
#   coordinate transforms
#    def Location(object):
#        #location of observatory, timekeeping?
# class Geometry(object):
#    # takes source positions and observation positions and calculates xyz
#    # could also calculate baselines
# class Beam(object)
