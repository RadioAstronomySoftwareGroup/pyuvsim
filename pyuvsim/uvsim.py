import numpy as np
import astropy.constants as const
import astropy.units as units
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from pyuvdata import UVData
import pyuvdata.utils as uvutils


class Source(object):
    name = None
    freq = None
    stokes = None
    ra = None
    dec = None
    epoch = None
    coherency_radec = None

    def __init__(self, name, ra, dec, epoch, freq, stokes):
        """
        Initialize from source catalog

        Args:
            name: Unique identifier for the source.
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

        self.name = name
        self.freq = freq
        self.stokes = stokes
        self.ra = ra
        self.dec = dec
        self.epoch = epoch

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        self.coherency_radec = .5*np.array([[self.stokes[0] + self.stokes[1],
                                             self.stokes[2] - 1j * self.stokes[3]],
                                           [self.stokes[2] + 1j * self.stokes[3],
                                            self.stokes[0] - self.stokes[1]]])

    def __eq__(self, other):
        return ((self.ra == other.ra) and
                (self.dec == other.dec) and
                (self.epoch == other.epoch) and
                (self.stokes == other.stokes) and
                (self.name == other.name) and
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

        coherency_local = np.einsum('ab,bc,cd->ad', rotation_matrix.T,
                                    self.coherency_radec, rotation_matrix)

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

        #print('Source az_za')
        #print(az_za)

        pos_l = np.cos(az_za[0]) * np.sin(az_za[1])
        pos_m = np.sin(az_za[0]) * np.sin(az_za[1])
        pos_n = np.cos(az_za[1])
        return (pos_l, pos_m, pos_n)


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
        self.enu = antenna2.pos_enu - antenna1.pos_enu
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
        self.uvdata_index = None        # Where to add the visibility in the uvdata object.

    def __eq__(self, other):
        return ((self.time == other.time) and
                (self.freq == other.freq) and
                (self.source == other.source) and
                (self.baseline == other.baseline) and
                (self.visibility_vector == other.visibility_vector) and
                (self.uvdata_index == other.uvdata_index) and
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

        #print('Pos lmn')
        #print(pos_lmn)

        # This hard coding shouldn't be required.
        # pos_lmn = [0, 0, 1]

        fringe = np.exp(-2j * np.pi * np.dot(self.task.baseline.uvw, pos_lmn))
        pos_lmn = self.task.source.pos_lmn(self.task.time, self.task.telescope.telescope_location)

        #print('fringe')
        #print(fringe)

        vij = self.apparent_coherency * fringe

        # need to reshape to be [xx, yy, xy, yx]
        vis_vector = [vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]]

        return np.array(vis_vector)

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
                        times[blts_ind].ravel(), sources[source_ind].ravel(),
                        blts_ind.ravel(), freq_ind.ravel())
    for (bl, freq, t, source, blti, fi) in uvtask_params:
        task = UVTask(source, t, freq, bl, telescope)
        task.uvdata_index = (blti,0,fi)  # 0 = spectral window index
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

def serial_gather(uvtask_list):
    """
        Initialize uvdata object, loop over uvtask list, acquire visibilities,
        and add to uvdata object.
    """
    uv_out = initialize_uvdata(uvtask_list)
    for task in uvtask_list:
        blt_ind,spw_ind,freq_ind = task.uvdata_index
        uv_out.data_array[blt_ind,spw_ind,freq_ind,:] += task.visibility_vector

    return uv_out

def create_mock_catalog(Nsrcs, time):
    """Create mock catalog with test sources."""
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec

    catalog = []
    # Divide totaly Stokes I intensity among all sources
    # Test file has Stokes I = 1 Jy
    for src_num in xrange(Nsrcs):
        catalog.append(Source('src'+str(src_num), ra, dec, time,
                                      freq, [1./Nsrcs, 0, 0, 0]))
    catalog = np.array(catalog)
    return catalog

def run_serial_uvsim(input_uvfits, catalog=None, Nsrcs=3):
    """Run uvsim."""
    input_uv = UVData()
    input_uv.read_uvfits(input_uvfits)

    time = Time(input_uv.time_array[0], scale='utc', format='jd')
    if catalog is None:
        catalog = create_mock_catalog(Nsrcs, time)

    uvtask_list = uvfile_to_task_list(input_uvfits, catalog)

    for task in uvtask_list:
        engine = UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uvdata_out = serial_gather(uvtask_list)

    return uvdata_out

def run_uvsim(input_uvfits, catalog=None, Nsrcs=3):
    """Run uvsim."""
    from mpi4py import MPI

    # Initialize MPI, get the communicator, number of Processing Units (PUs)
    # and the rank of this PU
    comm = MPI.COMM_WORLD
    Npus = comm.Get_size()
    rank = comm.Get_rank()

    # The Head node will initialize our simulation
    # Read input file and make uvtask list
    uvtask_list = []
    if rank==0:
        input_uv = UVData()
        input_uv.read_uvfits(input_uvfits)

        time = Time(input_uv.time_array[0], scale='utc', format='jd')
        if catalog is None:
            catalog = create_mock_catalog(Nsrcs, time)

        uvtask_list = uvfile_to_task_list(input_uvfits, catalog)
        # To split into PUs make a list of lists length NPUs
        uvtask_list = np.array_split(uvtask_list, Npus)
        uvtask_list = [ list(tl) for tl in uvtask_list]

    # Scatter the task list among all available PUs
    local_task_list = comm.scatter(uvtask_list, root=0)

    summed_task_dict = {}
    for task in local_task_list:
        engine = UVEngine(task)
        if task.uvdata_index not in summed_task_dict.keys():
            summed_task_dict[task.uvdata_index] = task
        if summed_task_dict[task.uvdata_index].visibility_vector is None:
            summed_task_dict[task.uvdata_index].visibility_vector = engine.make_visibility()
        else:
            summed_task_dict[task.uvdata_index].visibility_vector += engine.make_visibility()

    # All the sources in this summed list are foobar-ed
    # Source are summed over but only have 1 name
    # Some source may be correct
    summed_local_task_list = summed_task_dict.values()
    # gather all the finished local tasks into a list of list of len NPUs
    # gather is a blocking communication, have to wait for all PUs
    full_tasklist = comm.gather(summed_local_task_list, root=0)

    # Concatenate the list of lists into a flat list of tasks
    if rank == 0:
        uvtask_list = sum(full_tasklist,[])
        uvdata_out = serial_gather(uvtask_list)

        return uvdata_out
