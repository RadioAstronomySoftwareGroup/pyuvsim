
Comparison with other simulators
================================

pyuvsim is designed to calculate visibilities from point sources and a given primary beam model. It does this through direct calculation of the full-polarization Radio Interferometer Measurement Equation (RIME), fully parallelized across a set of MPI processes. As such, it has relatively poor performance compared with several other simulators. Deciding which tool to use will ultimately depend on what factors are important to the user. We recommend pyuvsim as a reference tool to verify the results of other simulators, and for producing limited datasets for comparison with observations (particularly on long baselines). Users will need access to a lot of memory and a large number of processor cores to take advantage of pyuvsim's accuracy and precision.

On this page, we compare pyuvsim's features and performance to several other visibility simulators.

Precision Radio Interferometry Simulator (PRISim)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A python package for generating full-sky simulations.

Advantages (over pyuvsim):

1. A bit faster.
2. Well-tested through years of use.
3. Support for a range of sky models including point sources, diffuse emission, and spectral cubes
4. Support for tracking as well as transit telescopes 
5. Parallelized over baselines, frequencies, or sky model/catalog
6. Has been tested using imaging with CASA

Disadvantages:

1. MPI overhead is large if time axis in simulation is long
2. Supports only one polarization at a time
3. Does not support non-identical antenna patterns

Fast Holographic Deconvolution (FHD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FHD is a data calibration and modeling framework written in IDL by researchers at the University of Washington. It calculates sky models

Advantages:

1. Fast.
2. Well-tested through years of use.
3. Support for diffuse foreground models.

Disadvantages:

1. Written in the proprietary IDL language.
2. Designed for the Murchison Widefield Array (MWA), and has limited support for other instruments.
3. The discrete convolution in Fourier Space can introduce aliasing and ringing artifacts at a level relevant to 21cm cosmology.


Common Astronomy Software Applications (CASA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CASA is a set of tools published by the National Radio Astronomical Observatory (NRAO), and has tools available for visibility simulation.

Advantages:

1. Faster, and less
2. Established in the field already and very well-documented.

Disadvantages:

1. Limited support for user-defined primary beam models.
2. Lacks full-sky coverage.
3. Unparallelized.
