import numpy as np
import astropy.units
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

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

    def coherency_calc(self, time, array_location):
        """
        Calculate the local coherency in az/za basis for this source at a time & location.

        The coherency is a 2x2 matrix giving electric field correlation in Jy.
        It's specified on the object as a coherency in the ra/dec basis,
        but must be rotated into local az/za.

        Args:
            time: astropy Time object
            array_location: astropy EarthLocation object

        Returns:
            local coherency in az/za basis
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(array_location, EarthLocation):
            raise ValueError('array_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=array_location))

        # First need to calculate the sin & cos of the parallactic angle
        # See Meeus's astronomical algorithms eq 14.1
        # also see Astroplan.observer.parallactic_angle method
        time.location = array_location
        lst = time.sidereal_time('mean')
        hour_angle = (lst - self.ra).rad
        sinX = np.sin(hour_angle)
        cosX = np.tan(array_location.lat) * np.cos(self.dec) - np.sin(self.dec) * np.cos(hour_angle)

        rotation_matrix = np.array([[cosX, sinX], [-sinX, cosX]])

        coherency_local = np.einsum('ab,bc,cd->ad', rotation_matrix.T,
                                    self.coherency_radec, rotation_matrix)

        return coherency_local

    def az_za_calc(self, time, array_location):
        """
        calculate the azimuth & zenith angle for this source at a time & location

        Args:
            time: astropy Time object
            array_location: astropy EarthLocation object

        Returns:
            (azimuth, zenith_angle) in radians
        """
        if not isinstance(time, Time):
            raise ValueError('time must be an astropy Time object. '
                             'value was: {t}'.format(t=time))

        if not isinstance(array_location, EarthLocation):
            raise ValueError('array_location must be an astropy EarthLocation object. '
                             'value was: {al}'.format(al=array_location))

        source_coord = SkyCoord(self.ra, self.dec, frame='icrs', equinox=self.epoch)

        source_altaz = source_coord.transform_to(AltAz(obstime=time, location=array_location))

        az_za = (source_altaz.az.rad, source_altaz.zen.rad)
        return az_za

    def pos_lmn(self, time, array_location):
        """
        calculate the direction cosines of this source at a time & location

        Args:
            time: astropy Time object
            array_location: astropy EarthLocation object

        Returns:
            (l, m, n) direction cosine values
        """
        # calculate direction cosines of source at current time and array location
        az_za = self.az_za_calc(time, array_location)

        pos_l = np.cos(az_za[0]) * np.sin(az_za[1])
        pos_m = np.sin(az_za[0]) * np.sin(az_za[1])
        pos_n = np.cos(az_za[1])
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
        self.pos_enu = enu_position
        # index of beam for this antenna from array.beam_list
        self.beam_id = beam_id

    def get_beam_jones(self, array, source_lmn, frequency):
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
        self.source = source
        self.baseline = baseline
        self.array = array


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
        beam1_jones = baseline.antenna1.get_beam_jones(self.task.array, source.pos_lmn, self.task.freq)
        beam2_jones = baseline.antenna2.get_beam_jones(self.task.array, source.pos_lmn, self.task.freq)
        this_apparent_coherency = np.dot(beam1_jones,
                                         source.coherency_calc(self.task.time,
                                                               self.task.array.array_location))
        this_apparent_coherency = np.dot(this_apparent_coherency,
                                         (beam2_jones.conj().T))

        self.apparent_coherency = this_apparent_coherency

    def make_visibility(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.apply_beam()

        fringe = np.exp(-2j * np.pi * np.dot(self.task.baseline.uvw,
                                             self.task.source.pos_lmn(self.task.time,
                                                                      self.task.array.array_location)))
        vij = self.apparent_coherency * fringe

        # need to reshape to be [xx, yy, xy, yx]
        vis_vector = [vij[0, 0], vij[1, 1], vij[0, 1], vij[1, 0]]

        return vis_vector

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
