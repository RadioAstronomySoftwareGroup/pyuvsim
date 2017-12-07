#execution pseudo code

#if rank==0
    # load stuff
    # broadcast
    # gather
    # sum sources
    # save file
#else:
    # rankN
    # load array details
    # create object
    # recv
    # put details into object
    # object.calculate()
    # format results
    # send, return to recv

#TODO:
# loading array details, observation settings
# what does the MPI traffic look like?
#    who does the formatting to input mpi parameters into uvengine
#    ANS: each MPI process takes a list of inputs, creates a list of tasks, and
#      hands them to the UVEngine
# unitttests

# FUTURE:
#   Baseline-dependent beams via Mueller matrix formalism

# read source catalog, generate Source objects,
#    set Source.calc to whatever
# 30 Nov 2017 - started source object, made _Blank_ test file
#    next steps: add basic parameter test for Source, flux calculation, ENU direction, coherence
#       -- HW Adam add some unittest boilerplate, next time testing Source attributes exist
#    after that, Baseline object, ECEF to ENU, empty beam hook
#    GOAL: 1 source on 1 baseline. -  first important unittest
# 7 Dec 2017 - Debated Jones vs Mueller and coordinate systems. Began defining Jones class and apply_beam method.
#    next steps: add in jones coordinate transformations and make unit tests for it.
#     refs: Smirnov Series 2011 papers. Nunhokee 2017
#    Future tie-in: UVBeam to get jones matrix option. Open github issue. Request support for coordinate transformations.
#    Zac --> Make us some diagrams of coord systems for different parts of the transformation.

class Source(object):
    def __init__(self):
        self.freq = None   # Hz, float
        self.coherence_radec = None  # Jy, ndarray, shape=(2,2) dtype=complex128
        self.coherence_local = None  # Jy, ndarray, shape=(2,2) dtype=complex128; In local ENU coordinates.
        self.time = None  # 
        self.array_location = None
        self.s = None
        self.ra = None
        self.dec = None
        self.polarization_angle = None  #pol angle in ra-dec
        self.flux_calc_done = False
        self.update_done = False
        self.epoch = None
    def update(self,task):
        #update with the current simulation time, freq etc
        self.freq = task.freq
        self.time = task.time
        self.array_location = task.array_location
    def flux_calc():
        self.coherence  #2x2 matrix giving electric field correlation in Jy.
                        # Specified as coherence in ra/dec basis, but must be rotated into local coords.
    def xyz_calc():
        #calculate xyz direction of source at current time and array location
        self.s = calc_ENU_direction(self.time,self.array_location,(self.ra,self.dec,self.epoch))
    def calc():
        self.flux_calc()
        self.xyz_calc()
    #debate: will an object ever be _re_-used?
    # answer, that is not our initial intent. Remake objects for each t,f etc


class Baseline(object):
    def __init__(self):
        self.antenna1 = None
        self.antenna2 = None


class UVTask(object):
    #holds all the information necessary to calculate a single src, t,f,bl
    def __init__(self,source,time,freq,uvw):
        self.time=time
        self.freq=freq
        self.source = source.update(self)
    def check(self):
        #make sure all input objects are syncronized
        self.source.time = self.time
        self.source.freq = self.freq

class Jones(object):
    ## Holds a single jones matrix and coordinate system, and methods for rotating.
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
 #inputs x,y,z,flux,baseline(u,v,w), time, freq
 #x,y,z in same coordinate system as uvws

    def __init__(self,array_details):
        self.rank
        #construct self based on MPI input
    def calculate_beams(self):
        #calculate beam pierce for every source and antenna
        # implicitly recalculating for every baseline
        for task in self.tasks:
            source = task.source
            task.baseline.antenna1.beam.calc_beam_jones(source)
            task.baseline.antenna2.beam.calc_beam_jones(source)
    def calculate_sky_model(self):
        #convert list of fluxes and spectral indices (or whatever)
        for task in self.tasks:

            task.source.calc()#update anything about the source depending on the location, time, freq
            #calculate local xyz coords
            task.source.s(array_location)
            #calculate electric field coherence (electric field arriving at ant)

    ## Debate --- Do we allow for baseline-defined beams, or stick with just antenna beams?
       ## This would necessitate the Mueller matrix formalism.
       ## As long as we stay modular, it should be simple to redefine things.

    def apply_beam(self,jonesL,jonesR,jones_coord_sys):
        # Supply jones matrices and their coordinate. Knows current coords of visibilities, applies rotations.
        #for every antenna calculate the apparent jones flux
        ## Beam --> Takes coherence matrix alt/az to ENU
        self.fluxes = []
        for task in self.tasks:
            baseline = task.baseline
            source = task.source
            #coherence is a 2x2 matrix
            # (Ex^2 ExEy, EyEx Ey^2)
            # where x and y are in the local alt/az coordinate system.
            x = np.dot(baseline.antenna1.beam.jones, source.coherence)
            x = np.dot(x,(baseline.antenna2.beam.jones.conj().T))
            
            self.fluxes.append(np.dot(x,(baseline.antenna2.beam.jones.conj().T)))
    def sim(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.fringe = np.exp(-2j*np.pi * np.dot(self.baseline_uvw,self.source_pos_lmn))
        self.vij = self.apparent_jones * 
        self.vij *= np.exp(-2j*np.pi)
        self.vij *= np.dot(self.baseline_uvw,self.source_pos_lmn)


#objects that are defined on input
#class SkyModel(object):
#       """ """
##can be derived from sources, a healpix map, or an image
##(internally represented as lists of sources/components)
##hands per freq, time, direction inputs to a broadcaster
#class Spectrum(object):
#
#class SkySelection(object):
#
#class Compute(object):
#
#class Observation(object):
##defines user settable parameters
##(start time, stop time, pointing, integration time)
#class Instrument(object):
##limits of the instrument
##available ranges of parameters
##(frequencies, pointing ranges, integration times)
##coordinate transforms
#    def Location(object):
#        #location of observatory, timekeeping?
#class Geometry(object):
#    #takes source positions and observation positions and calculates xyz
#    #could also calculate baselines
#class Beam(object)
