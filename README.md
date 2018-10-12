# pyuvsim
PYUVsim is a comprehensive simulation package for radio interferometers in python.
[![Build Status](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim.svg?branch=master)](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim)
[![Coverage Status](https://coveralls.io/repos/github/RadioAstronomySoftwareGroup/pyuvsim/badge.svg?branch=master)](https://coveralls.io/github/RadioAstronomySoftwareGroup/pyuvsim?branch=master)


A number of analysis tools are available to simulate the output of a radio interferometer (CASA, OSCAR, FHD, PRISim, et al), however each makes numerical approximations to enable speed ups.  The PYUVsim goal is to provide a simulated instrument output which emphasizes accuracy and extensibility.

 # Motivation and Approach
The two primary pyuvsim goals are interferometer simulation accuracy at the level of precision necessary for 21cm cosmology science. Key elements of this approach include:
 1. High level of test coverage including accuracy (design goal is 97%).
 2. Include analytic tests in unittests.
 3. Comparison with external simulations.
 4. Design for scalability across many cpus.

 # Installation
 * Bleeding edge: `git clone https://github.com/RadioAstronomySoftwareGroup/pyuvsim`
 * pip: `pip install pyuvsim`
 * PYUVSim is mainly intended to run on clusters running the linux operating system

 ## Dependencies
  * `numpy`, `astropy`, `scipy`, `mpi4py`, `pyyaml`, `six`, `pyuvdata`
  * optionally `line_profiler` if you want to do profiling (support for profiling is built in)

 # Inputs
 Adam L.
 The simulator requires specification of telescope, sky model, and observation details. These variables are set via an input text file formatted with yaml.  The observation details can optionally defined instead by inputting a data file.
 * Telescope definition file
 * Observation definition file
 * Sky model definition file

 # Quick start guide
  Adam L.
 How to run a basic simulation.
1. Install from github or pip.
2. Use included template simulation files
3. run with mpi
4. run with profiling (?)


 # How to contribute
 Contributions to this package to add new features or address any of the
 issues in the [issue log](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/issues) are very welcome.
 Please submit improvements as pull requests against the repo after verifying that
 the existing tests pass and any new code is well covered by unit tests.

 Bug reports or feature requests are also very welcome, please add them to the
 issue log after verifying that the issue does not already exist.
 Comments on existing issues are also welcome.

 # Versioning Approach
We use a `generation.major.minor` format.
 * generation - certified comparisons with other simulations
 * major - substantial package changes released on a sub-yearly cycle
 * minor - fixes that shouldn't change any outputs
