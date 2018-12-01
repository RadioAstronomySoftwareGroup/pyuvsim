# pyuvsim

[![Build Status](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim.svg?branch=master)](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim)
[![Coverage Status](https://coveralls.io/repos/github/RadioAstronomySoftwareGroup/pyuvsim/badge.svg?branch=master)](https://coveralls.io/github/RadioAstronomySoftwareGroup/pyuvsim?branch=master)

PYUVsim is a comprehensive simulation package for radio interferometers in python.

A number of analysis tools are available to simulate the output of a radio interferometer (CASA, OSCAR, FHD, PRISim, et al), however each makes numerical approximations to enable speed ups.  The PYUVsim goal is to provide a simulated instrument output which emphasizes accuracy and extensibility.

# Motivation and Approach
The two primary pyuvsim goals are interferometer simulation accuracy at the level of precision necessary for 21cm cosmology science. Key elements of this approach include:

1. High level of test coverage including accuracy (design goal is 97%).
2. Include analytic tests in unittests.
3. Comparison with external simulations.
4. Design for scalability across many cpus.

# Installation
* Bleeding edge: `git clone https://github.com/RadioAstronomySoftwareGroup/pyuvsim`
<!---* pip: `pip install pyuvsim`-->
* PYUVSim is mainly intended to run on clusters running the linux operating system

## Dependencies
* `numpy`, `astropy`, `scipy`, `mpi4py`, `pyyaml`, `six`, `pyuvdata`
* optionally `line_profiler` if you want to do profiling (support for profiling is built in)

# Inputs

A simulation requires sets of times, frequencies, source positions and brightnesses, antenna positions, and direction-dependent primary beam responses. PYUVsim specifies times, frequencies, and array configuration via a UVData object (from the pyuvdata package), source positions and brightnesses via Source objects, and primary beams either through UVBeam or AnalyticBeam objects.

* All sources are treated as point sources, with flux specified in Stokes parameters and position in right ascension / declination in the International Celestial Reference Frame (equivalently, in J2000 epoch).
* Primary beams are specified as full electric field components, and are interpolated in angle and frequency. This allows for an exact Jones matrix to be constructed for each desired source position.
* Multiple beam models may be used throughout the array, allowing for more complex instrument responses to be modeled.

These input objects may be made from a data file or from a set of `yaml` configuration files. See [Running a simulation](https://pyuvsim.readthedocs.io/en/latest/usage.html).

# Quick start guide
Example `obsparam` configuration files may be found in the `reference_simulations` directory.

1. Install from github or pip.
2. Run off of a parameter file with 20 MPI ranks:
```
mpirun -n 20 python run_param_pyuvsim.py reference_simulations/obsparam_1.1.yaml
```

# Documentation
Documentation on how to run simulations and developer API documentation is hosted on [ReadTheDocs](https://pyuvsim.readthedocs.io).

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

* Generation - Release combining multiple new physical effects and or major computational improvements.
Testing: Backed by unittests, internal model validation, and significant external comparison.
* Major - Adds new physical effect or major computational improvement. Small number of improvements with each release.
Testing: Backed by unittests, internal model validation and limited external comparison.
* Minor - Bug fixes and small improvements not expected to change physical model.
Testing: Backed by unittests

## Some helpful definitions
* __Physical effects__: things like polarization effects, noise, ionospheric modeling, or nonterrestrial observing positions.
* __Major computational improvement__:  Support for new catalog types (e.g, diffuse maps), new analysis tools, changes to parallelization scheme
* __Small improvements__: Better documentation or example code, outer framework redesign.
