import numpy as np
import astropy.constants as const
import astropy.units as units
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from pyuvdata import UVData
import pyuvdata.utils as uvutils
import os
import utils
from itertools import izip
from mpi4py import MPI
from astropy.io.votable import parse_single_table
import __builtin__

try:
    import progressbar
    progbar = True
except(ImportError):
    progbar = False

progsteps = False
try:
    if os.environ['PYUVSIM_BATCH_JOB'] == '1':
        progsteps = True
        progbar = False
except(KeyError):
    progbar = progbar

# Initialize MPI, get the communicator, number of Processing Units (PUs)
# and the rank of this PU

comm = MPI.COMM_WORLD
Npus = comm.Get_size()
rank = comm.Get_rank()


def blank(func):
    return func


try:
    __builtin__.profile
    if not rank == 0:
        __builtin__.profile = blank
except AttributeError:
    __builtin__.profile = blank


# The frame radio astronomers call the apparent or current epoch is the
# "true equator & equinox" frame, notated E_upsilon in the USNO circular
# astropy doesn't have this frame but it's pretty easy to adapt the CIRS frame
# by modifying the ra to reflect the difference between
# GAST (Grenwich Apparent Sidereal Time) and the earth rotation angle (theta)
def tee_to_cirs_ra(tee_ra, time):
    from astropy import _erfa as erfa
    from astropy.coordinates.builtin_frames.utils import get_jd12
    era = erfa.era00(*get_jd12(time, 'ut1'))
    theta_earth = Angle(era, unit='rad')

    assert(isinstance(time, Time))
    assert(isinstance(tee_ra, Angle))
    gast = time.sidereal_time('apparent', longitude=0)
    cirs_ra = tee_ra - (gast - theta_earth)
    return cirs_ra


def cirs_to_tee_ra(cirs_ra, time):
    from astropy import _erfa as erfa
    from astropy.coordinates.builtin_frames.utils import get_jd12
    era = erfa.era00(*get_jd12(time, 'ut1'))
    theta_earth = Angle(era, unit='rad')

    assert(isinstance(time, Time))
    assert(isinstance(cirs_ra, Angle))
    gast = time.sidereal_time('apparent', longitude=0)
    tee_ra = cirs_ra + (gast - theta_earth)
    return tee_ra


class Source(object):
    name = None
    freq = None
    stokes = None
    ra = None
    dec = None
    coherency_radec = None
    az_za = None

    @profile
    def __init__(self, name, ra, dec, freq, stokes):
        """
        Initialize from source catalog

        Args:
            name: Unique identifier for the source.
            ra: astropy Angle object
                source RA in J2000 (or ICRS) coordinates
            dec: astropy Angle object
                source Dec in J2000 (or ICRS) coordinates
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

        self.name = name
        self.freq = freq
        self.stokes = stokes
        self.ra = ra
        self.dec = dec

        self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        # The coherency is a 2x2 matrix giving electric field correlation in Jy.
        # Multiply by .5 to ensure that Trace sums to I not 2*I
        self.coherency_radec = .5 * np.array([[self.stokes[0] + self.stokes[1],
                                               self.stokes[2] - 1j * self.stokes[3]],
                                              [self.stokes[2] + 1j * self.stokes[3],
                                               self.stokes[0] - self.stokes[1]]])

    def __eq__(self, other):
        return ((self.ra == other.ra)
                and (self.dec == other.dec)
                and (self.stokes == other.stokes)
                and (self.name == other.name)
                and (self.freq == other.freq))

    @profile
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
            lst = time.sidereal_time('apparent')

            cirs_source_coord = self.skycoord.transform_to('cirs')
            tee_ra = cirs_to_tee_ra(cirs_source_coord.ra, time)

            hour_angle = (lst - tee_ra).rad
            sinX = np.sin(hour_angle)
            cosX = np.tan(telescope_location.lat) * np.cos(self.dec) - np.sin(self.dec) * np.cos(hour_angle)

            rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]])

        coherency_local = np.einsum('ab,bc,cd->ad', rotation_matrix.T,
                                    self.coherency_radec, rotation_matrix)

        return coherency_local

    @profile
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

        source_altaz = self.skycoord.transform_to(AltAz(obstime=time, location=telescope_location))

        az_za = (source_altaz.az.rad, source_altaz.zen.rad)
        self.az_za = az_za
        return az_za

    @profile
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

        pos_l = np.sin(az_za[0]) * np.sin(az_za[1])
        pos_m = np.cos(az_za[0]) * np.sin(az_za[1])
        pos_n = np.cos(az_za[1])
        return (pos_l, pos_m, pos_n)


class Telescope(object):
    @profile
    def __init__(self, telescope_name, telescope_location, beam_list):
        # telescope location (EarthLocation object)
        self.telescope_location = telescope_location
        self.telescope_name = telescope_name

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list

    def __eq__(self, other):
        return ((np.allclose(self.telescope_location.to('m').value, other.telescope_location.to("m").value, atol=1e-3))
                and (self.beam_list == other.beam_list)
                and (self.telescope_name == other.telescope_name))


class AnalyticBeam(object):
    supported_types = ['tophat', 'gaussian']

    def __init__(self, type, sigma=None):
        if type in self.supported_types:
            self.type = type
        else:
            raise ValueError('type not recognized')

        self.sigma = sigma
        self.data_normalization = 'peak'

    def peak_normalize(self):
        pass

    def interp(self, az_array, za_array, freq_array):
        # (Naxes_vec, Nspws, Nfeeds or Npols, freq_array.size, az_array.size)

        if self.type == 'tophat':
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            interp_data[1, 0, 0, :] = 1
            interp_data[0, 0, 1, :] = 1
            interp_basis_vector = None
        elif self.type == 'gaussian':
            if self.sigma is None:
                raise ValueError("Sigma needed for gaussian beam")
            interp_data = np.zeros((2, 1, 2, freq_array.size, az_array.size), dtype=np.float)
            values = np.exp(-(az_array**2) / (2 * self.sigma**2))
            interp_data[1, 0, 0, :] = values
            interp_data[0, 0, 1, :] = values
            interp_basis_vector = None
        else:
            raise ValueError('no interp for this type: ', self.type)

        return interp_data, interp_basis_vector


class Antenna(object):
    @profile
    def __init__(self, name, number, enu_position, beam_id):
        self.name = name
        self.number = number
        # ENU position in meters relative to the telescope_location
        self.pos_enu = enu_position * units.m
        # index of beam for this antenna from array.beam_list
        self.beam_id = beam_id

    @profile
    def get_beam_jones(self, array, source_az_za, frequency):
        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in Az/ZA
        # return array.beam_list[self.beam_id].get_direction_jones(source_lmn, frequency)

        source_az = np.array([source_az_za[0]])
        source_za = np.array([source_az_za[1]])
        freq = np.array([frequency.to('Hz').value])

        if array.beam_list[self.beam_id].data_normalization != 'peak':
            array.beam_list[self.beam_id].peak_normalize()
        array.beam_list[self.beam_id].interpolation_function = 'az_za_simple'

        interp_data, interp_basis_vector = \
            array.beam_list[self.beam_id].interp(az_array=source_az,
                                                 za_array=source_za,
                                                 freq_array=freq)

        # interp_data has shape: (Naxes_vec, Nspws, Nfeeds, 1 (freq), 1 (source position))
        jones_matrix = np.zeros((2, 2), dtype=np.complex)
        # first axis is feed, second axis is theta, phi (opposite order of beam!)
        jones_matrix[0, 0] = interp_data[1, 0, 0, 0, 0]
        jones_matrix[1, 1] = interp_data[0, 0, 1, 0, 0]
        jones_matrix[0, 1] = interp_data[0, 0, 0, 0, 0]
        jones_matrix[1, 0] = interp_data[1, 0, 1, 0, 0]

        return jones_matrix

    def __eq__(self, other):
        return ((self.name == other.name)
                and np.allclose(self.pos_enu.to('m').value, other.pos_enu.to('m').value, atol=1e-3)
                and (self.beam_id == other.beam_id))


class Baseline(object):
    @profile
    def __init__(self, antenna1, antenna2):
        self.antenna1 = antenna1
        self.antenna2 = antenna2
        self.enu = antenna2.pos_enu - antenna1.pos_enu
        # we're using the local az/za frame so uvw is just enu
        self.uvw = self.enu

    def __eq__(self, other):
        return ((self.antenna1 == other.antenna1)
                and (self.antenna2 == other.antenna2)
                and np.allclose(self.enu.to('m').value, other.enu.to('m').value, atol=1e-3)
                and np.allclose(self.uvw.to('m').value, other.uvw.to('m').value, atol=1e-3))


class UVTask(object):
    # holds all the information necessary to calculate a single src, t, f, bl, array
    # need the array because we need an array location for mapping to locat az/za

    @profile
    def __init__(self, source, time, freq, baseline, telescope):
        self.time = time
        self.freq = freq
        self.source = source
        self.baseline = baseline
        self.telescope = telescope
        self.visibility_vector = None
        self.uvdata_index = None        # Where to add the visibility in the uvdata object.

    def __eq__(self, other):
        return (np.isclose(self.time, other.time, atol=1e-4)
                and np.isclose(self.freq, other.freq, atol=1e-4)
                and (self.source == other.source)
                and (self.baseline == other.baseline)
                and (self.visibility_vector == other.visibility_vector)
                and (self.uvdata_index == other.uvdata_index)
                and (self.telescope == other.telescope))

    def __cmp__(self, other):
        # NB __cmp__ is not allowed in Python3.

        blti0, _, fi0 = self.uvdata_index
        blti1, _, fi1 = other.uvdata_index

        if blti0 > blti1:
            return 1
        elif blti0 == blti1:
            return 1 if fi0 > fi1 else -1
        else:
            return -1


class UVEngine(object):
    # inputs x,y,z,flux,baseline(u,v,w), time, freq
    # x,y,z in same coordinate system as uvws

    @profile
    def __init__(self, task):   # task_array  = list of tuples (source,time,freq,uvw)
        # self.rank
        self.task = task
        # Initialize task.time to a Time object.
        if isinstance(self.task.time, float):
            self.task.time = Time(self.task.time, format='jd')
        if isinstance(self.task.freq, float):
            self.task.freq = self.task.freq * units.Hz
        # construct self based on MPI input

    # Debate --- Do we allow for baseline-defined beams, or stick with just antenna beams?
    #   This would necessitate the Mueller matrix formalism.
    #   As long as we stay modular, it should be simple to redefine things.
    @profile
    def apply_beam(self):
        # Supply jones matrices and their coordinate. Knows current coords of visibilities, applies rotations.
        # for every antenna calculate the apparent jones flux
        # Beam --> Takes coherency matrix alt/az to ENU
        baseline = self.task.baseline
        source = self.task.source
        # coherency is a 2x2 matrix
        # (Ex^2 conj(Ex)Ey, conj(Ey)Ex Ey^2)
        # where x and y vectors along the local za/az coordinates.
        beam1_jones = baseline.antenna1.get_beam_jones(self.task.telescope,
                                                       source.az_za_calc(self.task.time,
                                                                         self.task.telescope.telescope_location),
                                                       self.task.freq)
        beam2_jones = baseline.antenna2.get_beam_jones(self.task.telescope,
                                                       source.az_za_calc(self.task.time,
                                                                         self.task.telescope.telescope_location),
                                                       self.task.freq)
        this_apparent_coherency = np.dot(beam1_jones,
                                         source.coherency_calc(self.task.time,
                                                               self.task.telescope.telescope_location))
        this_apparent_coherency = np.dot(this_apparent_coherency,
                                         (beam2_jones.conj().T))

        self.apparent_coherency = this_apparent_coherency

    @profile
    def make_visibility(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.apply_beam()

        pos_lmn = self.task.source.pos_lmn(self.task.time, self.task.telescope.telescope_location)

        # need to convert uvws from meters to wavelengths
        assert(isinstance(self.task.freq, Quantity))
        uvw_wavelength = self.task.baseline.uvw / const.c * self.task.freq.to('1/s')
        fringe = np.exp(2j * np.pi * np.dot(uvw_wavelength, pos_lmn))

        vij = self.apparent_coherency * fringe
        # need to reshape to be [xx, yy, xy, yx]
        vis_vector = [vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]]

        # Temporary -- write out task and other things to file.
        bl = str(self.task.baseline.antenna1.number) + "_" + str(self.task.baseline.antenna2.number)
        return np.array(vis_vector)

    @profile
    def update_task(self):
        self.task.visibility_vector = self.make_visibility()


@profile
def read_gleam_catalog(gleam_votable):
    """
    Creates a list of pyuvsim source objects from the GLEAM votable catalog.
    Despite the semi-standard votable format, there are enough differences that every catalog probably
    needs its own function.
    List of tested catalogs: GLEAM EGC catalog, version 2
    """

    table = parse_single_table(gleam_votable)
    data = table.array

    sourcelist = []
    for entry in data:
        source = Source(entry['GLEAM'], Angle(entry['RAJ2000'], unit=units.deg),
                        Angle(entry['DEJ2000'], unit=units.deg),
                        freq=(200e6 * units.Hz),
                        stokes=np.array([entry['Fintwide'], 0., 0., 0.]))
        sourcelist.append(source)
    return sourcelist


@profile
def uvdata_to_task_list(input_uv, sources, beam_list, beam_dict=None):
    """Create task list from pyuvdata compatible input file.

    Returns: List of task parameters to be send to UVEngines
    List has task parameters defined in UVTask object
    This function extracts time, freq, Antenna1, Antenna2
    """
    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    if not isinstance(sources, np.ndarray):
        raise TypeError("sources must be a numpy array")

    freq = input_uv.freq_array[0, :]  # units.Hz

    telescope = Telescope(input_uv.telescope_name,
                          EarthLocation.from_geocentric(*input_uv.telescope_location, unit='m'),
                          beam_list)

    if len(beam_list) > 1 and beam_dict is not None:
        raise ValueError('beam_dict must be supplied if beam_list has more than one element.')

    times = input_uv.time_array

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

    baselines = []
    print('Generating Baselines')
    for count, antnum1 in enumerate(input_uv.ant_1_array):
        antnum2 = input_uv.ant_2_array[count]
        index1 = np.where(input_uv.antenna_numbers == antnum1)[0][0]
        index2 = np.where(input_uv.antenna_numbers == antnum2)[0][0]
        baselines.append(Baseline(antennas[index1], antennas[index2]))

    baselines = np.array(baselines)

    blts_index = np.arange(input_uv.Nblts)
    frequency_index = np.arange(input_uv.Nfreqs)
    source_index = np.arange(len(sources))
    print('Making Meshgrid')
    blts_ind, freq_ind, source_ind = np.meshgrid(blts_index, frequency_index, source_index)
    print('Raveling')
    blts_ind = blts_ind.ravel()
    freq_ind = freq_ind.ravel()
    source_ind = source_ind.ravel()

    uvtask_list = []
    print('Making Tasks')
    print('Number of tasks:', len(blts_ind))

    if progsteps or progbar:
        count = 0
        tot = len(blts_ind)
        if progbar:
            pbar = progressbar.ProgressBar(maxval=tot).start()
        else:
            pbar = utils.progsteps(maxval=tot)

    for (bl, freqi, t, source, blti, fi) in izip(baselines[blts_ind],
                                                 freq[freq_ind], times[blts_ind],
                                                 sources[source_ind], blts_ind,
                                                 freq_ind):

        task = UVTask(source, t, freqi, bl, telescope)
        task.uvdata_index = (blti, 0, fi)    # 0 = spectral window index
        uvtask_list.append(task)

        if progbar or progsteps:
            count += 1
            pbar.update(count)

    if progbar:
        pbar.finish()
    return uvtask_list


@profile
def initialize_uvdata(uvtask_list):
    # writing 4 pol polarization files now

    # Assume all tasks have the same telescope.
    #   Enforce this generally!

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
    freq_0 = uvtask_list[0].freq
    for task in uvtask_list:
        if not task.source == source_0:
            continue
        task_freqs.append(task.freq)

        if task.freq == freq_0:
            task_bls.append(task.baseline)
            task_times.append(task.time)
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
    antpos_ecef = uvutils.ECEF_from_ENU(antpos.T, *uv_obj.telescope_location_lat_lon_alt).T - uv_obj.telescope_location
    uv_obj.antenna_positions = antpos_ecef
    uv_obj.ant_1_array = np.array(ant_1_array)
    uv_obj.ant_2_array = np.array(ant_2_array)
    uv_obj.time_array = np.array(task_times)
    uv_obj.uvw_array = np.array(task_uvw)
    uv_obj.baseline_array = uv_obj.antnums_to_baseline(ant_1_array, ant_2_array)
    uv_obj.Nbls = np.unique(uv_obj.baseline_array).size
    if uv_obj.Nfreqs == 1:
        uv_obj.channel_width = 1.  # Hz
    else:
        uv_obj.channel_width = np.diff(freqs)[0]

    if uv_obj.Ntimes == 1:
        uv_obj.integration_time = 1.  # Second
    else:
        uv_obj.integration_time = np.diff(np.unique(task_times))[0]
    uv_obj.set_lsts_from_time_array()
    uv_obj.zenith_ra = uv_obj.lst_array
    uv_obj.zenith_dec = np.repeat(uv_obj.telescope_location_lat_lon_alt[0], uv_obj.Nblts)  # Latitude
    uv_obj.object_name = 'zenith'
    uv_obj.set_drift()
    uv_obj.vis_units = 'Jy'
    uv_obj.polarization_array = np.array([-5, -6, -7, -8])
    uv_obj.spw_array = np.array([0])
    uv_obj.freq_array = np.array([freqs])

    uv_obj.Nspws = uv_obj.spw_array.size
    uv_obj.Npols = uv_obj.polarization_array.size

    uv_obj.data_array = np.zeros((uv_obj.Nblts, uv_obj.Nspws, uv_obj.Nfreqs, uv_obj.Npols), dtype=np.complex)
    uv_obj.flag_array = np.zeros((uv_obj.Nblts, uv_obj.Nspws, uv_obj.Nfreqs, uv_obj.Npols), dtype=bool)
    uv_obj.nsample_array = np.ones_like(uv_obj.data_array, dtype=float)
    uv_obj.history = 'UVSim'

    uv_obj.check()

    return uv_obj


@profile
def serial_gather(uvtask_list, uv_out):
    """
        Initialize uvdata object, loop over uvtask list, acquire visibilities,
        and add to uvdata object.
    """
    for task in uvtask_list:
        blt_ind, spw_ind, freq_ind = task.uvdata_index
        uv_out.data_array[blt_ind, spw_ind, freq_ind, :] += task.visibility_vector

    return uv_out


@profile
def create_mock_catalog(time, arrangement='zenith', array_location=None, Nsrcs=None, zen_ang=None, save=False):
    """
        Create mock catalog with test sources at zenith.

        arrangment = Choose test point source pattern (default = 1 source at zenith)
        Keywords:
            Nsrcs = Number of sources to put at zenith (ignored for other source arrangements)
            array_location = EarthLocation object.
            zen_ang = For off-zenith and triangle arrangements, how far from zenith to place sources. (deg)
            save = Save mock catalog as npz file.
        Accepted arrangements:
            'triangle' = Three point sources forming a triangle around the zenith
            'cross'    = An asymmetric cross
            'horizon'  = A single source on the horizon   ## TODO
            'zenith'   = Some number of sources placed at the zenith.
            'off-zenith' = A single source off zenith
            'long-line' = Horizon to horizon line of point sources
            'hera_text' = Spell out HERA around the zenith

    """

    if array_location is None:
        array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                       height=1073.)
    freq = (150e6 * units.Hz)

    if arrangement not in ['off-zenith', 'zenith', 'cross', 'triangle', 'long-line', 'hera_text']:
        raise KeyError("Invalid mock catalog arrangement" + str(arrangement))

    if arrangement == 'off-zenith':
        if zen_ang is None:
            zen_ang = 5.0  # Degrees
        Nsrcs = 1
        alts = [90. - zen_ang]
        azs = [90.]   # 0 = North pole, 90. = East pole
        fluxes = [1.0]

    if arrangement == 'triangle':
        Nsrcs = 3
        if zen_ang is None:
            zen_ang = 3.0
        alts = [90. - zen_ang, 90. - zen_ang, 90. - zen_ang]
        azs = [0., 120., 240.]
        fluxes = [1.0, 1.0, 1.0]

    if arrangement == 'cross':
        Nsrcs = 4
        alts = [88., 90., 86., 82.]
        azs = [270., 0., 90., 135.]
        fluxes = [5., 4., 1.0, 2.0]

    if arrangement == 'zenith':
        if Nsrcs is None:
            Nsrcs = 1
        alts = np.ones(Nsrcs) * 90.
        azs = np.zeros(Nsrcs, dtype=float)
        fluxes = np.ones(Nsrcs) * 1 / Nsrcs
        # Divide total Stokes I intensity among all sources
        # Test file has Stokes I = 1 Jy
    if arrangement == 'long-line':
        if 'Nsrcs' is None:
            Nsrcs = 10
        fluxes = np.ones(Nsrcs, dtype=float)
        zas = np.linspace(-85, 85, Nsrcs)
        alts = 90. - zas
        azs = np.zeros(Nsrcs, dtype=float)
        inds = np.where(alts > 90.0)
        azs[inds] = 180.
        alts[inds] = 90. + zas[inds]

    if arrangement == 'hera_text':

        azs = np.array([-254.055, -248.199, -236.310, -225.000, -206.565,
                        -153.435, -123.690, -111.801, -105.945, -261.870,
                        -258.690, -251.565, -135.000, -116.565, -101.310,
                        -98.130, 90.000, 90.000, 90.000, 90.000, 90.000,
                        -90.000, -90.000, -90.000, -90.000, -90.000,
                        -90.000, 81.870, 78.690, 71.565, -45.000, -71.565,
                        -78.690, -81.870, 74.055, 68.199, 56.310, 45.000,
                        26.565, -26.565, -45.000, -56.310, -71.565])

        zas = np.array([7.280, 5.385, 3.606, 2.828, 2.236, 2.236, 3.606,
                        5.385, 7.280, 7.071, 5.099, 3.162, 1.414, 2.236,
                        5.099, 7.071, 7.000, 6.000, 5.000, 3.000, 2.000,
                        1.000, 2.000, 3.000, 5.000, 6.000, 7.000, 7.071,
                        5.099, 3.162, 1.414, 3.162, 5.099, 7.071, 7.280,
                        5.385, 3.606, 2.828, 2.236, 2.236, 2.828, 3.606, 6.325])

        alts = 90. - zas
        Nsrcs = zas.size
        fluxes = np.ones_like(azs)

    catalog = []

    source_coord = SkyCoord(alt=Angle(alts, unit=units.deg), az=Angle(azs, unit=units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    for si in range(Nsrcs):
        catalog.append(Source('src' + str(si), ra[si], dec[si], freq, [fluxes[si], 0, 0, 0]))
    if rank == 0 and save:
        np.savez('mock_catalog', ra=ra.rad, dec=dec.rad)

    catalog = np.array(catalog)
    return catalog


def run_serial_uvsim(input_uv, beam_list, catalog=None, Nsrcs=3):
    """Run uvsim."""

    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    time = Time(input_uv.time_array[0], scale='utc', format='jd')
    if catalog is None:
        catalog = create_mock_catalog(time, arrangement='zenith')

    uvtask_list = uvdata_to_task_list(input_uv, catalog, beam_list)
    uvdata_out = initialize_uvdata(uvtask_list)

    for task in uvtask_list:
        engine = UVEngine(task)
        task.visibility_vector = engine.make_visibility()

    uvdata_out = serial_gather(uvtask_list, uvdata_out)

    return uvdata_out


@profile
def run_uvsim(input_uv, beam_list, catalog=None, Nsrcs=None, mock_arrangement='zenith'):
    """Run uvsim."""
    if not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    # The Head node will initialize our simulation
    # Read input file and make uvtask list
    uvtask_list = []
    if rank == 0:
        print('Nblts:', input_uv.Nblts)
        print('Nfreqs:', input_uv.Nfreqs)

        time = Time(input_uv.time_array[0], scale='utc', format='jd')
        if catalog is None:
            array_loc = EarthLocation.from_geocentric(*input_uv.telescope_location, unit='m')
            if Nsrcs is not None:
                print("Nsrcs:", Nsrcs)
                catalog = create_mock_catalog(time, arrangement=mock_arrangement, array_location=array_loc, Nsrcs=Nsrcs)
            else:
                catalog = create_mock_catalog(time, arrangement=mock_arrangement, array_location=array_loc)

        uvtask_list = uvdata_to_task_list(input_uv, catalog, beam_list)
        uv_container = initialize_uvdata(uvtask_list)
        # To split into PUs make a list of lists length NPUs
        print("Splitting Task List")
        uvtask_list = np.array_split(uvtask_list, Npus)
        uvtask_list = [list(tl) for tl in uvtask_list]

        print("Sending Tasks To Processing Units")
    # Scatter the task list among all available PUs
    local_task_list = comm.scatter(uvtask_list, root=0)
    if rank == 0:
        print("Tasks Received. Begin Calculations.")
    summed_task_dict = {}

    if rank == 0:
        if progsteps or progbar:
            count = 0
            tot = len(local_task_list)
            if progbar:
                pbar = progressbar.ProgressBar(maxval=tot).start()
            else:
                pbar = utils.progsteps(maxval=tot)

        for count, task in enumerate(local_task_list):
            engine = UVEngine(task)
            if task.uvdata_index not in summed_task_dict.keys():
                summed_task_dict[task.uvdata_index] = task
            if summed_task_dict[task.uvdata_index].visibility_vector is None:
                summed_task_dict[task.uvdata_index].visibility_vector = engine.make_visibility()
            else:
                summed_task_dict[task.uvdata_index].visibility_vector += engine.make_visibility()

            if progbar or progsteps:
                pbar.update(count)

        if progbar or progsteps:
            pbar.finish()
    else:
        for task in local_task_list:
            engine = UVEngine(task)
            if task.uvdata_index not in summed_task_dict.keys():
                summed_task_dict[task.uvdata_index] = task
            if summed_task_dict[task.uvdata_index].visibility_vector is None:
                summed_task_dict[task.uvdata_index].visibility_vector = engine.make_visibility()
            else:
                summed_task_dict[task.uvdata_index].visibility_vector += engine.make_visibility()

    if rank == 0:
        print("Calculations Complete.")

    # All the sources in this summed list are foobar-ed
    # Source are summed over but only have 1 name
    # Some source may be correct
    summed_local_task_list = summed_task_dict.values()
    # gather all the finished local tasks into a list of list of len NPUs
    # gather is a blocking communication, have to wait for all PUs
    full_tasklist = comm.gather(summed_local_task_list, root=0)

    # Concatenate the list of lists into a flat list of tasks
    if rank == 0:
        uvtask_list = sum(full_tasklist, [])
        uvdata_out = serial_gather(uvtask_list, uv_container)

        return uvdata_out
