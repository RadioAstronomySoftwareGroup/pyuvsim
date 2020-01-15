# Changelog

## [Unreleased]

### Added
- Support for individual shape parameters for analytic beams.

### Changed

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
