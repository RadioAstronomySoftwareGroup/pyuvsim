
Comparison with other simulators
================================

pyuvsim is designed to calculate visibilities from point sources and a given primary beam model. It does this through direct calculation of the full-polarization Radio Interferometer Measurement Equation (RIME), parallelized across a set of MPI processes. Our design approach prioritizes precision, design clarity, flexible usefulness over speed and efficiency. As such, it has relatively slow speed and large memory usage compared with several other simulators. Deciding which tool to use will ultimately depend on what factors are important to the user. We recommend pyuvsim as a reference tool to verify the results of other simulators, and for producing limited datasets for comparison with observations (particularly on long baselines). Users will need access to a lot of memory and a large number of processor cores to take advantage of pyuvsim's accuracy and precision, ideally on a dedicated high-performance computing cluster.

On this page, we compare pyuvsim's features and performance to some other visibility simulators.

Precision Radio Interferometry Simulator (PRISim)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A python package for generating full-sky simulations.

Advantages:

1. Faster due to vectorized calculations.
2. Well-tested through years of use.
3. Support for a range of sky models including point sources, diffuse emission, and spectral cubes.
4. Support for tracking as well as transit telescopes.
5. Parallelized over baselines, frequencies, or sky model/catalog.
6. Has been tested using imaging with CASA.

Disadvantages:

1. MPI overhead is large if time axis in simulation is long.
2. Limited to accurately simulating the total intensity (stokes I) radio interferometer measurement equation, which is accurate if the polarized sky emission is negligible and the off-diagonal terms in the beam pattern (e.g. XY or YX for linear polarization) are negligible.
3. Does not support non-identical antenna patterns.

Fast Holographic Deconvolution (FHD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FHD is a data calibration and modeling framework written in IDL by researchers
by researchers at the University of Washington, originally developed for analyzing
EoR data from the Murchison Widefield Array (MWA) in Western Australia.
It calculates sky models by convolving a high-resolution psf with the exact DFT
of point sources.

Advantages:

1. Optimized for smooth-spectrum sources.
2. Well-tested through years of use.
3. Support for diffuse foreground models.
4. Support for discrete pointings.
5. Supports non-identical antenna patterns (can be memory intensive if there are many different patterns).
6. Uses IDL's native parallelization efficiently.


Disadvantages:

1. Written in the proprietary IDL language.
2. The discrete convolution in Fourier Space can introduce aliasing and ringing artifacts at a level relevant to 21cm cosmology.
3. Memory usage scales with the size of the uv-plane. It is very fast for compact 21 cm cosmology arrays, but for higher resolution arrays with sparse uv coverage the memory scaling can be poor.
4. Primary beam motion on the sky is limited to the "snapshot" size, not the time step size.
5. Much slower for simulating complex spectral structure.

.. image:: fhd_uvsim_compare.png
    :width: 600
    :alt: Delay-spectrum comparison of FHD and pyuvsim point source simulations.

The figure above shows a comparison of a `pyuvsim` simulation of several point sources with a corresponding FHD simulation. The visibilities have been Fourier-transformed along the frequency axis (a "delay transform") with a Blackman-Harris window function, squared, and then averaged in time. This forms a so-called *delay spectrum*, which is proportional to the power spectrum modes along the line of sight. For smooth-spectrum foregrounds, the power at high delays should be very low. The FHD simulation (shown in orange) introduces excess power at high delays, which can be a real problem for 21cm cosmological research.

Common Astronomy Software Applications (CASA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CASA is a set of tools published by the National Radio Astronomical Observatory (NRAO), and has tools available for visibility simulation. Fully parallelized by partitioning Measurement Set (MS) files into smaller tasks and distributed with MPI across any number of available processing units. Non-trivial parallelization also available via OpenMP for shared memory computations on a single node.

Advantages:

1. Uses compiled C code internally, which is faster.
2. Established in the field already and very well-documented.
3. OpenMP utilizes shared memory on a single node if the calculation an be decomposed efficiently. MPI for all other parallel processing needs.
4. Support for a source component lists and FITS image source models.

Disadvantages:

1. Limited support for user-defined primary beam models.
2. Internal UVW rotation is known to be incorrect, affecting coherence far from the phase center (CASA helpdesk ticket 2291, listed as closed but apparently not fixed).
3. In its default (and fastest) mode of operation, point sources are gridded to pixel locations so an FFT can be performed. This pixel-scale imprecision can introduce point source subtraction errors that are significant to 21cm cosmology experiments [CTROTT2012]_.
4. Full direction-dependent Jones matrices can only be simulated if the beam times sky model calculation is carried out in separate software [JAGANNATHAN17]_.
5. Does not support non-identical antenna beam patterns.


.. [CTROTT2012]
   Trott, Cathryn M., Randall B. Wayth, and Steven J. Tingay. "The impact of point-source subtraction residuals on 21 cm epoch of reionization estimation." The Astrophysical Journal 757.1 (2012): 101.

.. [JAGANNATHAN17]
   Jagannathan, P., et al. "Direction-dependent Corrections in Polarimetric Radio Imaging. I. Characterizing the Effects of the Primary Beam on Full-Stokes Imaging." The Astronomical Journal 154.2 (2017): 56.
