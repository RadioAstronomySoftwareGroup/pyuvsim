# Reference Simulations

``pyuvsim`` reference simulations are intended as a reference point for comparison
with existing simulators, and as a benchmark for regression testing. As the capability
of `pyuvsim` has improved, our reference simulations have grown in complexity.
Each generation is detailed in its own folder and replaces the previous set as the
validation standard for new `pyuvsim` versions. For explicit documentation of the reference
simulations and discussion and sample analysis of simulation output, see
[first_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/first_generation) and
[second_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/second_generation).

For regression testing, only the lightweight first generation simulations are integrated into
the testing pipeline. First generation reference simulations can be ran, downloaded, and compared
locally through `pytest` in addition to the scripts offered in `reference_simulations`. See
[Benchmarking](https://pyuvsim.readthedocs.io/en/latest/developers.html#benchmarking)
for documentation of running the reference simulations using `pytest`.

The current approach to downloading large files is [pooch](https://github.com/fatiando/pooch).
`pooch` is used to download files to the cache from all current download points. We have configured
downloading to use `pooch.os_cache(pyuvsim)` which should place files in a pyuvsim folder in the
default cache directory of your os. To run the reference simulatiions, all necessary data files
should be installed into the cache directory.

The current approach to storing historical reference simulations is the [Brown Digital Repository (BDR)](https://repository.library.brown.edu/studio/collections/bdr:wte2qah8/).
We aim to consistently upload first generation reference simulations to the BDR across
code changes for regression testing.

For convenience, several scripts are provided to help run the reference simulations and
compare their results to the latest set posted to the BDR:

 - `download_data_files.py`
   An extensible script to download necessary beam and catalog files. The script requires
   `pyradiosky`, `pooch`, and `astroquery`. Sample script usage:
   ```
   # see argpare help menu for sample usage
   python3 download_data_files.py --help
   # download all data files by default if no argument given
   # files downloaded to the cache
   python3 download_data_files.py
   # download specific data files by keyword
   python3 download_data_files.py keyword1 keyword2
   ```
   Currently downloadable data files by keyword:
    - gleam: GLEAM extragalactic source catalog from Vizier, using `astroquery`, and save it to a VOTable file. GLEAM, treated as flat-spectrum, is used for several reference simulations.
    - mwa: mwa uvbeam file from [here](https://github.com/MWATelescope/mwa_pb).
    - healpix: "healpix": gsm 2016 nside 128 healpix map saved as skyh5 [here](https://repository.library.brown.edu/studio/item/bdr:eafzyycj/).
 
 - `download_ref_sim.py`
   Run this first to download the reference sim data from the BDR using `pooch` (requires
   `pooch`). The new data can be compared to these files to check for output consistency.
   The files will by default be downloaded to the current working directory in a folder named
   `latest_ref_sims`. The folder can also be placed in `pooch.os_cache(pyuvsim)` by passing
   the necessary flag. Sample script usage:
   ```
   # see argpare help menu for sample usage
   python3 download_ref_sim.py --help
   # download all latest reference by default if no argument given
   # files downloaded to the cache
   python3 download_ref_sim.py
   # download specific reference simulation files by keyword
   python3 download_ref_sim.py keyword1 keyword2
   ```
   Currently downloadable reference simulation files by keyword:
    - 1.1_uniform
    - 1.1_gauss
    - 1.1_mwa
    - 1.2_uniform
    - 1.2_gauss
    - 1.3_uniform
    - 1.3_gauss
 - `compare_with_last.py`
   Given the paths to the latest output files (uvh5), this will compare the data in those files
   to the corresponding files in `latest_ref_data`. 

At this time, the following scripts are designed for a cluster managed by Simple Linux Utility Resource
Manager (SLURM): (TODO: fix and adapt these scripts as they are broken and outdated)

 - run_ref_sims.sh
        To run all of the reference simulations, run:
        ```
            ./run_ref_sims.sh obsparam_*
        ```
        This launches SLURM jobs, described by the script `jobscript.sh`.

 - jobscript.sh
        This sets up an output directory, named by the current date, into which output files,
        slurm output files, profiling data, and records of which job IDs correspond with each obsparam file.

More information about how to run the simulations may be found at
[ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/usage.html).
