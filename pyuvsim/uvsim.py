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
#

class Source(object):
    def __init__(self):
        self.freq = None   # Hz, float
        self.stokes = None  # [I, Q, U, V] in RA Dec
        self.ra = None
        self.dec = None
        self.epoch = None

        self.coherency_radec = np.array([[self.stokes[0] + self.stokes[1],
                                          self.stokes[2] - 1j * self.stokes[3]],
                                         [self.stokes[2] + 1j * self.stokes[3],
                                          self.stokes[0] - self.stokes[1]]])

    def coherency_calc(self, time, array_location):
        # 2x2 matrix giving electric field correlation in Jy.
        # Specified as coherency in ra/dec basis, but must be rotated into local az/za.
        coherency_local = # calculate coherency from coherency_radec in local az/za

        return coherency_local

    def az_za_calc(self, time, array_location):
        # calculate az_za direction of source at current time and array location
        # 2-element tuple with (az, za)
        az_za = calc_az_za(time, array_location, self.ra, self.dec, self.epoch)
        return az_za

    def pos_lmn(self, time, array_location):
        # calculate direction cosines of source at current time and array location
        az_za = az_za_calc(time, array_location)

        pos_l = cos(az_za[0]) * sin(az_za[1])
        pos_m = sin(az_za[0]) * sin(az_za[1])
        pos_n = cos(az_za[1])
        return (pos_l, pos_m, pos_n)

    # debate: will an object ever be _re_-used?
    # answer, that is not our initial intent. Remake objects for each t,f etc
    # new plan from 2/1: Source objects describe sources in celestial coords
    # but have methods to convert those coords to the local az/za frame


class Array(object):
    def __init__(self, array_location, beam_list):
        # array location in ECEF
        self.array_location = array_location

        # list of UVBeam objects, length of number of unique beams
        self.beam_list = beam_list


class Antenna(object):
    def __init__(self, enu_position, beam_id):
        # ENU position in meters relative to the array_location
        self.pos_enu
        # index of beam for this antenna from array.beam_list
        self.beam_id

    def get_beam_jones(array, source_lmn, frequency):
        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in Az/ZA
        # return array.beam_list[self.beam_id].get_direction_jones(source_lmn, frequency)
        return np.array([[1, 0], [0, 1]])

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


class UVTask(object):
    # holds all the information necessary to calculate a single src, t, f, bl, array
    # need the array because we need an array location for mapping to locat az/za
    def __init__(self, source, time, freq, baseline, array):
        self.time = time
        self.freq = freq
        self.source = source.update(self)

    def check(self):
        # make sure all input objects are syncronized
        self.source.time = self.time
        self.source.freq = self.freq


class Jones(object):
    # Holds a single jones matrix and coordinate system, and methods for rotating.
    def __init__(self):
        self.jones = None    # (2,2) matrix
        self.coord_sys = None  # [ radec, altaz, EN ]

    def rotate2altaz(self):
        if self.coord_sys == 'radec':
        if self.coord_sys == 'EN':

    def rotate2radec(self):
        if self.coord_sys == 'altaz':
        if self.coord_sys == 'EN':


class UVEngine(object):
    # inputs x,y,z,flux,baseline(u,v,w), time, freq
    # x,y,z in same coordinate system as uvws

    def __init__(self, task_array):   # task_array  = list of tuples (source,time,freq,uvw)
        self.rank
        self.tasks = [UVTask(t) for t in task_array]
        # construct self based on MPI input

    def calculate_beams(self):
        # calculate beam pierce for every source and antenna
        # implicitly recalculating for every baseline
        for task in self.tasks:
            source = task.source
            source_lmn = source.pos_lmn
            beam1_jones = antenna1.get_beam_jones(task.array, source_lmn, task.frequency)
            if antenna1.beam_id == antenna2.beam_id:
                return beam1_jones, beam1_jones
            else:
                beam2_jones = antenna2.get_beam_jones(task.array, source_lmn, task.frequency)
                return beam1_jones, beam2_jones

    def calculate_sky_model(self):
        # convert list of fluxes and spectral indices (or whatever)
        for task in self.tasks:
            task.source.calc()  # update anything about the source depending on the location, time, freq
            # calculate local xyz coords
            task.source.s(array_location)
            # calculate electric field coherency (electric field arriving at ant)

    # Debate --- Do we allow for baseline-defined beams, or stick with just antenna beams?
    #   This would necessitate the Mueller matrix formalism.
    #   As long as we stay modular, it should be simple to redefine things.

    def apply_beam(self, jonesL, jonesR, jones_coord_sys):
        # Supply jones matrices and their coordinate. Knows current coords of visibilities, applies rotations.
        # for every antenna calculate the apparent jones flux
        # Beam --> Takes coherency matrix alt/az to ENU
        self.apparent_coherency = []
        for task in self.tasks:
            baseline = task.baseline
            source = task.source
            # coherency is a 2x2 matrix
            # (Ex^2 conj(Ex)Ey, conj(Ey)Ex Ey^2)
            # where x and y vectors along the local za/az coordinates.
            beam1_jones, beam2_jones = calculate_beams()
            this_apparent_coherency = np.dot(beam1_jones,
                                             source.coherency_calc(self.task.time,
                                                                   self.task.array.array_location))
            this_apparent_coherency = np.dot(this_apparent_coherency,
                                             (beam2_jones.conj().T))

            self.apparent_coherency.append(this_apparent_coherency)

    def make_visibility(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.fringe = np.exp(-2j * np.pi * np.dot(self.task.baseline.uvw,
                                                  self.task.source.pos_lmn(self.task.time,
                                                                           self.task.array.array_location)))
        self.vij = self.apparent_coherency * self.fringe

# objects that are defined on input
# class Array(object):
#    array_location (maybe an astropy observer class)
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
