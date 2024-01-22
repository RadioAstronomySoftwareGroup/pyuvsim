# Changelog

## [Unreleased]

### Changed
- Add tracking of the beam `freq_interp_kind` to the BeamList object since
it is moving from a UVBeam attribute to a parameter to the `UVBeam.interp` method
in future pyuvdata versions. Also ensure compatibility with current and future
pyuvdata versions.

## [1.2.6] - 2022-07-17

### Changed
- Require pyradiosky >= 0.2.0
- Increase the minimum compatible version for lunarsky to 0.2.1 to fix astropy
deprecation warnings.

### Fixed
- Bugs where shared memory resources were not properly freed.
- Bugs when running with pyradiosky>=0.2, mostly related to the new `skycoord` parameter
on pyradiosky's `SkyModel` object and calls to functions removed from pyradiosky.
- Better handling of lunarsky imports, required for the newest version of pyradiosky

## [1.2.5] - 2022-06-01

### Changed
- Initial ordering of blt axis in `initialize_uvdata_from_params`
(unchanged output by default).
- Use future array shapes on UVBeam objects.
- Update dependency requirements to: pyuvdata >= 2.2.10,<2.3, pyradiosky>=0.1.0,<0.2,
numpy>=1.19, scipy>=1.3
- Update optional dependency requirements to: astropy-healpix>=0.6, lunarsky>=0.1.2,
python-casacore>=3.3.1

### Fixed
- A major bug introduced between v1.2.1 and v1.2.2 that caused errors in simulations
when the source list was large enough that it needed to be split among processing units.

### Added
- `reorder_kw` option to `initialize_uvdata_from_params` function.
- `check_kw` option to `initialize_uvdata_from_params` function.
- `select:` parameter to obsparam file definition for `telescope:`.

## [1.2.4] - 2022-06-01

### Added
- A check that the beam basis vectors are aligned with the azimuth and zenith angle in
each pixel.
- A new parameter, `beam_interp_check`, to `run_uvsim` and `run_uvdata_uvsim` to allow
the check that the beam covers the interpolation location to be turned off. The default
behavior is to turn off the check if the beam covers the full sky horizon to horizon.

### Changed
- Require pyuvdata >= 2.2.8

## [1.2.3] - 2022-05-10

### Added
- Added a `return_beams` parameter to the `initialize_uvdata_from_params` function to
have it return the beam list and dict. Defaults to True currently, but will default to
False starting in version 1.4.
- Added a `return_catname` parameter to the `initialize_catalog_from_params` function to
have it return the catalog name. Defaults to True currently, but will default to False
starting in version 1.4.
- An option to the `run_param_pyuvsim.py` script and the `run_uvsim` and
`run_uvdata_uvsim` functions to allow users to keep the output from nonzero ranks for
debugging purposes.

### Changed
- Fix auto visibilities to be real on file write-out by default if the pyuvdata
version is 2.2.7 or greater.
- Updated the astropy requirement to >= 5.0.4
- Dropped support for python 3.7

## [1.2.2] - 2022-02-22

### Added
- A `return_beams` parameter to `simsetup.initialize_uvdata_from_params` to allow beam
information to be returned. Currently defaults to True, will default to False starting
in version 1.4.
- A `return_catname` parameter to `simsetup.initialize_catalog_from_params` to allow
the catalog name to be returned. Currently defaults to True, will default to False
starting in version 1.4.
- A `filetype` parameter to `simsetup.initialize_catalog_from_params` to allow users to
specify the catalog filetype.
- Ensure that the version numbers for pyuvsim, pyradiosky and pyuvdata are written to
the history of the output UVData files.
- New `check_consistency` method for `BeamList` objects.
- Ability to simulate UVData objects where Nblts != Nbls * Ntimes. Currently only supported with direct `run_uvdata_uvsim` calls.

## Changed
- Updated the astropy requirement to >= 5.0.4
- Dropped support for python 3.7

### Fixed
- Input `uvdata.blt_order` forced to be (time, baseline) before a simulation is run. Attempts to reorder output uvdata object if input ordering was different.

## [1.2.1] - 2021-10-13

### Added
- Support for writing out measurement set files.
- Support for unit tests parallelized with MPI.
- Require that future changes not drastically increase runtime for current capabilities.

### Changed
- Require pyuvdata >= 2.1.5
- Use future array shapes on UVData objects.
- Use the replacements for the deprecated `SkyModel.source_cuts` method if the pyradiosky
version is new enough to have them.
- Use remote memory access to collect finished visibility data, without serialization.

### Fixed
- `x_orientation` is now checked on UVBeam objects and flowed into the output UVData files.
- `SkyModel.name` is now coerced to an array.
- Bug MPI-enabled tests which caused failure on tests that wouldn't pass in serial mode.
- Fix bugs in reading catalogs using pyradiosky for `skyh5` files and files with unusual extension names.
- Corrects the distribution of random points for the random mock catalog.
- pixel interpolation was defaulting to az_za_simple for all beams, breaking healpix-coord UVBeams.

## [1.2.0] - 2020-7-20

### Added
- Diffuse models in mock catalogs, via the analytic_diffuse module.
- quantity_shared_bcast function to allow objects derived from astropy.units.Quantity to use shared memory broadcasting.
- SkyModelData class to replace recarray conversion in pyradiosky.
- quiet keyword for run_uvsim, to suppress stdout printing.
- Support for Moon-based observing -- keyword "world: moon" in telescope config.
- Option to pass along interpolating spline order to UVBeam from teleconfig file.
- Scripts for running / verifying reference simulations.
- Benchmarking tools.

### Changed
- Switch to a new definition for the Counter class, without threading.
- Cleaned up unit tests
- Use `at_frequencies` method to enable support for all pyradiosky spectral types.
- Only do coherency calculation when the time changes
- Only do beam eval when time, freq, or beam type changes.
- The definition of the Airy beam now uses the exact value of c, not 3e8.

### Fixed
- Keep a UVBeam with more than two frequencies for tests, so the (default) cubic interpolation works.
- Use pytest hooks to ensure the profiler tests run last.
- Switched to using a remote-memory-access based counter class. Appears to have fixed bug in the Counter test.
- Ensure that source positions, coherency matrices, and Jones matrices are updated at the right times.
- Error early if the task list is too long for gather.

### Deprecated
- Support for pyradiosky versions <0.1.0.

## [1.1.2] - 2020-2-14

### Added
- BeamList class for handling the set of beams, and their string representations.
- Support for individual shape parameters for analytic beams.

### Fixed
- Removed warning catch for astropy GLEAM read in tests -- new versions don't raise the warning.

## [1.1.1] - 2019-11-14

### Added
- A function for checking the memory usage on each Node

### Changed
- Replaced init_uvdata_out function with complete_uvdata
- init_uvdata_out function is more modular.
- unit tests check error messages raised.
- Polarization selection is now allowed in setup. Will break if incorrect polarization is used in pyuvsim.
- Added functions to read healpix maps and added support for frequency axis

### Fixed
- Fixed a serious bug with improper conversion of Angle to rad/deg in simsetup.py
- skymodel splitting was not working correctly to avoid running out of memory.
- Flush stderr in excepthook, so stack trace is printed when an exception is raised in MPI.
- No longer calling UVBeam.interp with freq_interp_kind, for updated pyuvdata.
- Circular-installation bug
- Bug in how integration_time array is set in setup
- Bug in how filepath is interpreted in setup

## [1.1.0] - 2019-6-14

### Added
- A parallelized counter that can track progress across all ranks.
- shared memory broadcast, which will create a shared memory window for a given object on each node
- Function to generate UVData object from keywords, optionally saving config files.
- Optionally, can skip having telescope_config file if there is no beam list.

### Changed
- UVBeam frequency interpolation is cubic spline by default, not linear.
- Tasks are split over time, frequency, and baseline only (in that order).
- Tasks are split over the source axis if the estimated memory footprint exceeds available resources.
- The source class is replaced with SkyModel, which supports vectorized coordinate transformation and beam interpolation.

### Fixed
- MPI.Init -> MPI.Init_thread(), for threaded applications.
- Progress steps now update in real time, accurately reflecting job progress.
- Analytic visibility calculation tests now also check analytic beams.
- Analytic beams are redefined so that off-diagonal Jones matrix terms are zero.
- parameter dicts are not modified by functions using them.

### Deprecated
- Coarse horizon cuts are no longer performed. Should be restored in a future version.

## [1.0.0] - 2019-5-10

### Added
- More detailed comparison to other simulators in the documentation.
- Option to convert analytic beams from efield to pstokes
- Enabled Gaussian beams to be defined with power law widths.
- Enabled Gaussian beams to be defined from antenna diameter and have chromaticity
- Checking the beam kernel width makes sense for Airy beams

### Changed
- Pending deprecation of gaussian beams defined from sigma parameter

### Fixed
- pip installation instructions


## [0.2.3] - 2019-2-17

## [0.2.2] - 2019-2-8

### Added
- Support for only simulating one baseline per redundant group

### Fixed
- Using the correct Bessel function in the definition of the Airy beam


## [0.2.1] - 2018-12-10

### Added
- Support for miriad and uvh5 output
- stdout printing explicitly disabled for non-root ranks
- Flux cut options with source catalogs
- Coarse horizon cut for reduced computational load.
- Selection keys in obsparam files (valid options for UVData.select)
- New batch scripts for profiling jobs on a SLURM-managed cluster, and plotting the results.

### Changed
- Rename write_uvfits to write_uvdata
- Fixed memory overload bug due to local task indices creation step.
- Moved parameter simulation functions from wrapper to uvsim
- MPI3 enabled in travis.
- Line profiling may be done by running `pyuvsim.profiling.set_profiler`.
- Which functions to profile is set by listing them by name, instead of by function decorators.
- Improved test coverage

### Deprecated
- Simulation directly off of uvfits files.


## [0.2.0] - 2018-10-26

### Changed
- Tasks are generated in place on each rank.
- UVEngines are reused by the worker loops.
- Alt/Az quantities are only recalculated for each source when the time changes.
- Use `reuse_splines` option on UVBeams to reduce interpolation time.
- Beams are scattered as strings, not objects.


## [0.1.0] - 2018-10-24
