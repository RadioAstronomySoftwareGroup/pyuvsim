# Reference Simulations

``pyuvsim`` reference simulations are intended as a reference point for comparison
with existing simulators. As the capability of `pyuvsim` has improved, our
reference simulations have grown in complexity. Each set is detailed in its own
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
        This setups up an output directory, named by the current date, into which output files,
   slurm output files, profiling data, and records of which job IDs correspond with each obsparam file.
 - compare_with_last.py
        Given the paths to the latest output files (uvh5), this will compare the data in those files
   to the corresponding files in `latest_ref_data`.
 - download_latest_data.py
        Run this first to download the reference sim data from the google drive. The new data should be
   compared with these files. They will end up in a directory `latest_ref_data`.


More information about how to run the simulations may be found at
[ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/usage.html).
