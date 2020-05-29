For Developers
==============

The Radio Astronomy Software Group welcomes contributions from anyone. Full contribution guidelines and code of conduct may be found on the `RASG Github Pages <https://radioastronomysoftwaregroup.github.io>`_.

This page discusses more specific points to help developers contribute to ``pyuvsim``.

Computing Resources
------------------

Instrument simulation is a computationally-demanding task, and ``pyuvsim`` is designed to run a brute-force calculation that requires much computing power. There is a continual effort to improve in performance, but for the foreseeable future the best results can be obtained using large, distributed computing clusters.

Much of the development and testing of ``pyuvsim`` was done on the Oscar computing cluster at Brown University's Center for Computation and Visualization (CCV_). This cluster is managed with the Simple Linux Utility Resource Manager (SLURM). Though ``pyuvsim`` can be run on most Linux or Mac-OS systems, its operation with large numbers of times/baselines/frequencies/sources or especially complicated beam objects has not been tested well. Testing within other cluster management environments is encouraged.

.. _CCV: https://docs.ccv.brown.edu/oscar/


Running and Comparing Reference Simulations
-------------------------------------------

The *reference simulations* comprise a set of standard simulations with well-defined, simple parameters that most instrument simulation tools should be capable of running. As such, they represent a standard point of comparison between ``pyuvsim`` and other simulators. They are also used to track that ``pyuvsim`` produces consistent results over time.

For some types of pull requests, we ask that the reference simulations be re-run with the new code to confirm that any changes to results are detected. The older set of files are stored in a Google drive, and are downloadable using scripts in the ``reference_simulations`` directory. If a change to the reference simulations is expected, it should be documented in the CHANGELOG and the new simulation file given to the project managers to store on the Google drive.

More information on checking the reference simulations can be found in the README in the ``reference_simulations`` directory in the repository.

We will periodically issue a new set of reference simulations as ``pyuvsim`` gains more capabilities.

For more details, see `reference_simulations/README.md <https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/master/reference_simulations>`_.

Benchmarking
------------

The ``benchmarking`` directory contains tools to test the runtime and memory usage of large simulations. There is no requirement to check benchmarks for pull requests, but it's a good idea to make sure changes don't drastically alter the runtime. The file BENCHMARKS.log keeps a record of performance over time.

The README file in the ``benchmarking`` directory gives more details on how to do benchmarking.

Note that the benchmarking scripts are designed only for SLURM systems.

For more details, see `benchmarking/README.md <https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/master/benchmarking>`_.
