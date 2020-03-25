# Benchmarking tools

This folder contains helper scripts to run large simulation jobs and record memory usage,
runtime, and line-by-line profiling information. Any changes to `pyuvsim` that may affect
overall performance should run one of these jobs to ensure that performance is not affected.

Currently, the default benchmarking simulation simulates a full-frequency HEALPix sky model,
an array of 500 baselines, 256 frequencies, 21 times, and analytic (Gaussian) beams.

## Usage
These can only be used effectively if `pyuvsim` is installed in "developer" mode. That is,
the repository is cloned locally and installed in conda using (from the outermost directory of
the repo):
```
    conda develop .
```
or in pip using
```
    pip install -e
```

The file `setting.yaml` contains the parameters describing the simulation. This includes the scales
of the various simulation axes (baselines, times, frequencies) as well as the resources. The current
configuration uses a HEALPix sky model, so the Nside parameter of that map must also be given.

The script `run_benchmarking.py` makes the configuration files and SLURM job script, and can submit the SLURM
job or update the BENCHMARKS.log log file with the latest run. To make the configuration files, run the following
from the benchmarking directory:
```
    python run_benchmarking.py settings.yaml --init_config
```
By default, this will make a directory with the current date and put the configuration and output files into it.

To submit the job (or initialize and submit, if the config hasn't been made yet):
```
    python run_benchmarking.py settings.yaml --submit
```
This will start the SLURM job.

When finished, you can update the log file BENCHMARKS.log by running the following:
```
    python run_benchmarking.py settings.yaml --update_log
```

Note that BENCHMARKS.log is tracked by git, and is intended as a record of performance over time. Updates to this
can be included with PRs for performance-affecting changes.

## Other Scripts

The `analyze_runtimes.py` script attempts to estimate how the runtime scales with different axis combinations.
Ideally, this will give a formula for the approximate runtime given a set of parameters. This is still experimental.
