#execution psuedo code

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
# unitttests


class UVEngine(object):
 #inputs x,y,z,flux,baseline(u,v,w), time, freq
 #x,y,z in same coordinate system as uvws

    def __init__(self,array_details):
        self.rank
        #construct self based on MPI input
    def calculate_beams(self):
        for task in self.tasks:
            source = task.source
            task.baseline.antenna1.beam.calc_beam_jones(source)
            task.baseline.antenna2.beam.calc_beam_jones(source)
    def calculate_sky_model(self):
        #convert list of fluxes and spectral indices (or whatever)
        for task in self.tasks:
            #calculate local xyz coords
            task.source.local_position(array_location)
            #calculate electric field coherence (electric field arriving at ant)

    def flux(self):
        #for every antenna calculate the apparent jones flux
        self.fluxes = []
        for task in self.tasks:
            baseline = task.baseline
            source = task.source
            #coherence is a 2x2 matrix
            # (Ex^2 ExEy, EyEx Ey^2)
            # where x and y are in the local ENU coordinate system
            #
            x = np.dot(baseline.antenna1.beam.jones * source.coherence)
            self.fluxes.append(np.dot(x,(baseline.antenna2.beam.jones.conj().T))
    def sim(self):
        # dimensions?
        # polarization, ravel index (freq,time,source)
        self.fringe = np.exp(-2j*np.pi * np.dot(self.baseline_uvw,self.source_pos_lmn))
        self.vij = self.apparent_jones *
        self.vij *= np.exp(-2j*np.pi
        self.vij *= np.dot(self.baseline_uvw,self.source_pos_lmn))


#objects that are defined on input
class SkyModel(object):
#can be derived from sources, a healpix map, or an image
#(internally represented as lists of sources/components)
#hands per freq, time, direction inputs to a broadcaster
    class Spectrum(object)
    class SkySelection(object)

class Compute(object)

class Observation(object)
#defines user settable parameters
#(start time, stop time, pointing, integration time)
class Instrument(object)
#limits of the instrument
#available ranges of parameters
#(frequencies, pointing ranges, integration times)
#coordinate transforms
    def Location(object):
        #location of observatory, timekeeping?
class Geometry(object):
    #takes source positions and observation positions and calculates xyz
    #could also calculate baselines
class Beam(object)
