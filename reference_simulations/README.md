# Reference Simulations

``pyuvsim`` reference simulations are intended as a reference point for comparison
with existing simulators, and as a benchmark for regression testing. As the capability
of `pyuvsim` has improved, our reference simulations have grown in complexity.
Each generation replaces the previous set as the validation standard for new `pyuvsim`
versions. For explicit documentation of the reference simulations and discussion and
sample analysis of simulation output, see
[first_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/first_generation) and
[second_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/second_generation).
Currently the first generation configuration files are stored in
[data](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/src/pyuvsim/data).
The second generation configuration files still exist in
[second_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/second_generation)
but are out of date and will not run without additional changes. Note that the README is also
very outdated.

For regression testing, only the lightweight first generation simulations are integrated into
the testing pipeline. First generation reference simulations can be run, downloaded, and compared
locally through `pytest` in addition to the scripts offered in `reference_simulations`. See
[Benchmarking](https://pyuvsim.readthedocs.io/en/latest/developers.html#benchmarking)
for documentation of running the reference simulations using `pytest`.

The current approach to downloading large files uses [astropy](https://github.com/astropy/astropy)
functionality to download files from all current download points.
We have configured downloading to use the astropy cache location for the package pyuvsim. To run
the reference simulations, all necessary data files should be installed into the cache directory
using `download_data_files`.

The current approach to storing historical reference simulations is the
[Brown Digital Repository (BDR)](https://repository.library.brown.edu/studio/collections/bdr:wte2qah8/).
We aim to consistently upload first generation reference simulations to the BDR across
code changes for regression testing.

For convenience, command line tools are provided to help run the reference simulations and
compare their results to the latest set posted to the BDR:

 - `download_data_files`
   Downloads necessary beam and catalog files. The script requires
   `pyradiosky`, `astropy`, and `astroquery`. Sample script usage:
   ```
   # see argparse help menu for sample usage
   download_data_files --help
   # download all data files by default if no argument given
   # files downloaded to the cache
   download_data_files
   # download specific data files by keyword
   download_data_files keyword1 keyword2
   ```
   Currently downloadable data files by keyword:
    - gleam: GLEAM extragalactic source catalog from Vizier, using `astroquery`, and save it to a
    VOTable file. GLEAM is used for several reference simulations.
    - mwa: mwa uvbeam file from [here](https://github.com/MWATelescope/mwa_pb).
    - healpix: gsm 2016 nside 128 healpix map saved as skyh5
    [here](https://repository.library.brown.edu/studio/item/bdr:eafzyycj/).
 - `download_ref_sims`
   Run this to download the latest reference simulation data from the BDR using `astropy` and
   `requests`. This cli tool is only used in pytest to perform comparisons for output consistency,
   and is not necessary to run the reference simulations.
   Sample script usage:
   ```
   # see argparse help menu for sample usage
   download_ref_sims --help
   # download all latest reference by default if no argument given
   # files downloaded to the cache
   download_ref_sims
   # download specific reference simulation files by keyword
   download_ref_sims keyword1 keyword2
   ```
   Currently downloadable reference simulation files by keyword:
    - 1.1_baseline_number
    - 1.2_time_axis
    - 1.3_frequency_axis
    - 1.4_source_axis
    - 1.5_uvbeam
    - 1.6_healpix
    - 1.7_multi_beam
    - 1.8_lunar

The first generation reference simulations run in only a couple minutes on a single core. The
simulations run from the top directory with:
```
# Running with 20 MPI processing units, but 1 is sufficient for first generation
# This should soon become an entrypoint
# 1.1_baseline_number used as example, works for all 1st gen sims with data files downloaded
mpirun -n 20 python scripts/run_param_pyuvsim.py reference_simulations/first_generation/obsparam_ref_1.1_baseline_number.yaml
```
The second generation reference simulations can take a few more hours, and as such should ideally
be run on some form of computing cluster. Running the second generation reference simulations
interactively works well but is not preferable. The following scripts are designed for a cluster
managed by Simple Linux Utility Resource Manager (SLURM). These scripts are nonfunctional
as they were last updated in 2020 and are out of date. Some of the second generation reference
simulations also require some manual file finding and placing to run.

 - run_ref_sims.sh
        To run all of the reference simulations, run:
        ```
            ./run_ref_sims.sh obsparam_*
        ```
        This launches SLURM jobs, described by the script `jobscript.sh`.

 - jobscript.sh
        This sets up an output directory, named by the current date, into which output files,
        slurm output files, profiling data, and records of which job IDs correspond with each
        obsparam file are saved.

More information about how to run the simulations may be found at
[ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/usage.html).
