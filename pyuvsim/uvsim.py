class UVEngine(object)
 #inputs x,y,z,flux,baseline(u,v,w), time, freq
 #x,y,z in same coordinate system as uvws



#objects that are defined on input
class SkyModel(object)
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
