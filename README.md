# pyuvsim

![](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/workflows/Tests/badge.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/RadioAstronomySoftwareGroup/pyuvsim/badge.svg?branch=master)](https://coveralls.io/github/RadioAstronomySoftwareGroup/pyuvsim?branch=master)
[![codecov](https://codecov.io/gh/RadioAstronomySoftwareGroup/pyuvsim/branch/master/graph/badge.svg)](https://codecov.io/gh/RadioAstronomySoftwareGroup/pyuvsim)

pyuvsim is a comprehensive simulation package for radio interferometers in python.

A number of analysis tools are available to simulate the output of a radio
interferometer (CASA, OSKAR, FHD, PRISim, et al), however each makes numerical
approximations to enable speed ups. pyuvsim's goal is to provide a simulated
instrument output which emphasizes accuracy and extensibility, and can represent the most
general simulator design.

A comparison to other simulators may be found at [ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/comparison.html).

## pyuvsim, the Interferometer Simulator of Record
pyuvsim's  primary goal is to be an interferometer simulator accurate at the level of
precision necessary for 21cm cosmology science,

1. High level of test coverage including accuracy (design goal is 97%).
2. Testing against analytic calculations, monitored by continuous integration (see memo #XXX)
3. Comparison with external simulations with standardized reference simulations
## Usability and extensibility
A secondary goal is a community simulation environment which provides well documented and flexible code to support a diversity of use cases.
Key elements of this approach include:
1. Design for scalability across many cpus.
2. Defining a clear, user-friendly standard for simulation design.
3. Documentation of analytic validation and reference simulations

# Physical Instrumental Effects
Each addition of new physics is validated against analytic calculations and included in a new reference simulation. Physics that have been included or are on the roadmap.
1. Fully-polarized instrument response (complete)
1. Polarized sources (analytic testing ~90% )
1. Floating-point source position accuracy (complete)
1. Full-sky field of view (complete)
1. Exact antenna positions.  (complete)
1. Varied beam models across the array (complete, tested against analytic)
1. Diffuse emission (complete, tested against analytic, paper in prep)
1. Arbitrary spectrum (complete)
1. Non-terrestrial observatories (Lunar observatory complete)
1. Time domain sources (TODO)
1. Ionospheric scintillation (TODO)



## Installation
Simple installation via pip is available for users, developers should follow
the directions under [Developer Installation](#developer-installation) below.

A user-installation is achieved simply with `pip install pyuvsim`, or to get the
bleeding-edge: `pip install https://github.com/RadioAstronomySoftwareGroup/pyuvsim`.

By default, `mpi` capabilities are not enabled -- many of the utilities provided
in `pyuvsim` do not require it. To use the simulator within `pyuvsim`, you
should install `pyuvsim` with  `pip install pyuvsim[sim]`. Note that the
`pyuvsim` simulator is intended to run on clusters running the linux operating
system, but we do test against Mac OSX as well.

There are a few more optional dependencies for `pyuvsim` which enable some features,
such as `astropy_healpix` to use healpix based sky catalogs or healpix beams,
`python-casacore` for writing out measurement sets and `lunarsky` for simulating telescopes
on the moon. If you would like these tools as well as the full simulator, install
`pyuvsim` with `pip install pyuvsim[all]` (or use the `[healpix]`, `[casa]` or `[moon]`
options to only get the dependencies for each of those functionalities).

If you wish to manage dependencies manually read on.

### Dependencies
If you are using `conda` to manage your environment, you may wish to install the
following packages before installing `pyuvsim`:

Required:

* astropy>=4.0
* numpy>=1.15
* psutil
* pyradiosky>=0.1.0
* pyuvdata>=2.1.5
* pyyaml>=5.1
* scipy>1.0.1
* setuptools_scm

Optional:

* astropy-healpix (for using healpix based sky catalogs or beams)
* mpi4py>=3.0.0 (for actually running simulations)
* lunarsky (for simulating telescopes on the moon)
* python-casacore >= 3.1.0 (for writing CASA measurement sets)

### Developer Installation
If you are developing `pyuvsim`, you will need to download and install the
repository using `git clone https://github.com/RadioAstronomySoftwareGroup/pyuvsim.git`.

If you use conda, you may wish to use a fresh environment, in which case you can
use the included `environment.yaml` file to make a conda environment with all
the extra dependencies required for testing/development as well as the
standard ones using `conda env create -f environment.yml`. If you do not wish to
make a fresh dedicated environment, you should verify that the environment you
are using contains the packages listed in that file.
Then do a developer install of pyuvsim using `pip install -e .` (or
`pip install --no-deps -e .` if you do not want pip to install any missing
requirements).

If you do not use conda, after downloading the repository, install using
`pip install -e .[dev]` to install all the extra dependencies required for
testing/development as well as the standard ones.

Finally, install the pre-commit hook using `pre-commit install` to help prevent
committing code that does not meet our style guidelines.


## Inputs

A simulation requires sets of times, frequencies, source positions and brightnesses, antenna positions, and direction-dependent primary beam responses. pyuvsim specifies times, frequencies, and array configuration via a UVData object (from the pyuvdata package), source positions and brightnesses via Source objects, and primary beams either through UVBeam or AnalyticBeam objects.

* All sources are treated as point sources, with flux specified in Stokes parameters and position in right ascension / declination in the International Celestial Reference Frame (equivalently, in J2000 epoch).
* Primary beams are specified as full electric field components, and are interpolated in angle and frequency. This allows for an exact Jones matrix to be constructed for each desired source position.
* Multiple beam models may be used throughout the array, allowing for more complex instrument responses to be modeled.

These input objects may be made from a data file or from a set of `yaml` configuration files. See [Running a simulation](https://pyuvsim.readthedocs.io/en/latest/usage.html).

## Outputs

Data from a simulation run are written out to a file in any format accessible with `pyuvdata`. This includes:

* uvfits
* MIRIAD
* uvh5

When read into a UVData object, the `history` string will contain information on the pyuvsim and pyuvdata versions used for that run (including the latest git hash, if available), and details on the catalog used.

## Quick start guide
Example `obsparam` configuration files may be found in the `reference_simulations` directory.

1. Install from github or pip.
2. Run off of a parameter file with 20 MPI ranks:
```
mpirun -n 20 run_param_pyuvsim.py reference_simulations/obsparam_1.1.yaml
```

## Documentation
Documentation on how to run simulations and developer API documentation is hosted on [ReadTheDocs](https://pyuvsim.readthedocs.io).

## Testing

`pyuvsim` uses the `pytest` package for unit testing. If you've cloned the source into a directory `pyuvsim/`, you may verify it as follows:

1. Install `pytest` from anaconda or pip.
2. Run the pytest from `pyuvsim/`
```
pytest
```
You can alternatively run ```python -m pytest pyuvsim``` or ```python setup.py test```.
You will need to have all dependencies installed.

Some tests are run in parallel using the mpi4py module. Those tests have a decorator ``pytest.mark.parallel(n)``` where ```n``` is an integer giving the number
of parallel processes to run the test on. To temporarily disable parallel tests, run pytest with the option ```--nompi```.

## Where to find Support

Please feel free to submit new issues to the [issue log](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/issues) to request new features, document new bugs, or ask questions.

## How to contribute
Contributions to this package to add new features or address any of the
issues in the [issue log](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/issues) are very welcome, as are bug reports or feature requests.
Please see our [guide on contributing](.github/CONTRIBUTING.md)

## Versioning Approach
We use a `generation.major.minor` format.

* Generation - Release combining multiple new physical effects and or major computational improvements.
Testing: Backed by unittests, internal model validation, and significant external comparison.
* Major - Adds new physical effect or major computational improvement. Small number of improvements with each release.
Testing: Backed by unittests, internal model validation and limited external comparison.
* Minor - Bug fixes and small improvements not expected to change physical model.
Testing: Backed by unittests

### Some helpful definitions
* __Physical effects__: things like polarization effects, noise, ionospheric modeling, or nonterrestrial observing positions.
* __Major computational improvement__:  Support for new catalog types (e.g, diffuse maps), new analysis tools, changes to parallelization scheme
* __Small improvements__: Better documentation or example code, outer framework redesign.


# Maintainers
pyuvsim is maintained by the RASG Managers, which currently include:

 - Adam Beardsley (Arizona State University)
 - Bryna Hazelton (University of Washington)
 - Daniel Jacobs (Arizona State University)
 - Paul La Plante (University of California, Berkeley)
 - Jonathan Pober (Brown University)

Please use the channels discussed in the [guide on contributing](.github/CONTRIBUTING.md)
for code-related discussions. You can contact us privately if needed at
[rasgmanagers@gmail.com](mailto:rasgmanagers@gmail.com).
