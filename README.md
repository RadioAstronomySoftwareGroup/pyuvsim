# pyuvsim

[![Build Status](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim.svg?branch=master)](https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim)
[![Coverage Status](https://coveralls.io/repos/github/RadioAstronomySoftwareGroup/pyuvsim/badge.svg?branch=master)](https://coveralls.io/github/RadioAstronomySoftwareGroup/pyuvsim?branch=master)

pyuvsim is a comprehensive simulation package for radio interferometers in python.

A number of analysis tools are available to simulate the output of a radio
interferometer (CASA, OSKAR, FHD, PRISim, et al), however each makes numerical
approximations to enable speed ups. pyuvsim's goal is to provide a simulated
instrument output which emphasizes accuracy and extensibility, and can represent the most
general simulator design.

A comparison to other simulators may be found at [ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/comparison.html).

# Motivation and Approach
pyuvsim's two primary goals are interferometer simulation accuracy at the level of precision necessary for 21cm cosmology science, and maximum flexibility in use cases. Key elements of this approach include:

1. High level of test coverage including accuracy (design goal is 97%).
2. Include analytic tests in unittests.
3. Comparison with external simulations.
4. Design for scalability across many cpus.
5. Fully-polarized instrument response, floating-point source position accuracy, full-sky field of view, and exact antenna positions.
6. Support for varied beam models across the array.
7. Defining a clear, user-friendly standard for simulation design.

# Installation
* For simple installation, the latest stable release is available via pip (`pip install pyuvsim`)
* To install the development version: Clone the repository using `git clone https://github.com/RadioAstronomySoftwareGroup/pyuvsim`, navigate into the pyuvsim directory and run
```
pip install -r requirements.txt
pip install pyuvdata
pip install .
```
* `pyuvsim` is intended to run on clusters running the linux operating system.


**Note**
The `mpi4py` module is installed as a wrapper around an existing installation of MPI. The easiest way to install it is to use anaconda, which will install a compatible version of MPI and configure mpi4py to use it:
```
conda install -c conda-forge mpi4py
```


## Dependencies
* `numpy`, `astropy`, `scipy`, `mpi4py`, `pyyaml`, `six`, `pyuvdata`
* optionally `line_profiler` if you want to do profiling (support for profiling is built in)
* `h5py` is required to support the new `uvh5` data format.

# Inputs

A simulation requires sets of times, frequencies, source positions and brightnesses, antenna positions, and direction-dependent primary beam responses. pyuvsim specifies times, frequencies, and array configuration via a UVData object (from the pyuvdata package), source positions and brightnesses via Source objects, and primary beams either through UVBeam or AnalyticBeam objects.

* All sources are treated as point sources, with flux specified in Stokes parameters and position in right ascension / declination in the International Celestial Reference Frame (equivalently, in J2000 epoch).
* Primary beams are specified as full electric field components, and are interpolated in angle and frequency. This allows for an exact Jones matrix to be constructed for each desired source position.
* Multiple beam models may be used throughout the array, allowing for more complex instrument responses to be modeled.

These input objects may be made from a data file or from a set of `yaml` configuration files. See [Running a simulation](https://pyuvsim.readthedocs.io/en/latest/usage.html).

# Outputs

Data from a simulation run are written out to a file in any format accessible with `pyuvdata`. This includes:

* uvfits
* MIRIAD
* uvh5

When read into a UVData object, the `history` string will contain information on the pyuvsim and pyuvdata versions used for that run (including the latest git hash, if available), and details on the catalog used.

# Quick start guide
Example `obsparam` configuration files may be found in the `reference_simulations` directory.

1. Install from github or pip.
2. Run off of a parameter file with 20 MPI ranks:
```
mpirun -n 20 python run_param_pyuvsim.py reference_simulations/obsparam_1.1.yaml
```

# Documentation
Documentation on how to run simulations and developer API documentation is hosted on [ReadTheDocs](https://pyuvsim.readthedocs.io).

# Testing

`pyuvsim` uses the `nose` package for unit testing. If you've cloned the source into a directory `pyuvsim/`, you may verify it as follows:

1. Install `nose` from anaconda or pip.
2. Run the nosetests from `pyuvsim/`
```
nosetests pyuvsim/tests/
```

You will need to have all dependencies installed.

# Where to find Support

Please feel free to submit new issues to the [issue log](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/issues) to request new features, document new bugs, or ask questions.

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
