# Reference Simulations

``pyuvsim`` reference simulations are intended as a reference point for comparison
with existing simulators. As the capability of `pyuvsim` has improved, our
reference simulations have grown in complexity. Each generation is detailed in its own
folder and replaces the previous set as the validation standard for new `pyuvsim`
versions.

For convenience, several scripts are provided to help run the reference simulations and
compare their results to a standard set (hosted on a publicly-accessible Google drive).
At this time, these scripts are designed for a cluster managed by Simple Linux Utility Resource
Manager (SLURM).


 - run_ref_sims.sh
        To run all of the reference simulations, run:
        ```
            ./run_ref_sims.sh obsparam_*
        ```
        This launches SLURM jobs, described by the script `jobscript.sh`.

 - jobscript.sh
        This sets up an output directory, named by the current date, into which output files,
        slurm output files, profiling data, and records of which job IDs correspond with each obsparam file.

 - download_latest_data.py
        Run this first to download the reference sim data from the google drive. The new data will be
        compared to these files. They will end up in a directory `latest_ref_data`.

 - compare_with_last.py
        Given the paths to the latest output files (uvh5), this will compare the data in those files
        to the corresponding files in `latest_ref_data`.

 - get_gleam.py
        A script to download the GLEAM extragalactic source catalog from Vizier, using `astroquery`,
        and save it to a VOTable file. GLEAM, treated as flat-spectrum, is used for several reference simulations.


More information about how to run the simulations may be found at
[ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/usage.html).
