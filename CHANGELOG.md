# Changelog

## [Unreleased]

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
